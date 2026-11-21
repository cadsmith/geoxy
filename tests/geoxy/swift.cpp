/***************************************************************************************************
 *  @file      tests/geoxy/swift.cpp
 *  @brief     Tests for the swift kernels: a line, a ray or an endless line met with an ellipse, a
 *             polyline or a Bézier curve, each crossing the double nearest the exact one, and every
 *             answer the thorough way's to the bit or none.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "draws.hpp"
#include "expansion.hpp"
#include "real.hpp"
#include "swift/lane.hpp"
#include "swift/met.hpp"
#include "thorough.hpp"

#include "geoxy/angle.hpp"
#include "geoxy/at.hpp"
#include "geoxy/bezier/cubic.hpp"
#include "geoxy/bezier/quadratic.hpp"
#include "geoxy/curve.hpp"
#include "geoxy/derivative.hpp"
#include "geoxy/domain.hpp"
#include "geoxy/ellipse.hpp"
#include "geoxy/intersect.hpp"
#include "geoxy/line.hpp"
#include "geoxy/matrix.hpp"
#include "geoxy/point.hpp"
#include "geoxy/polyline.hpp"
#include "geoxy/ray.hpp"
#include "geoxy/tolerance.hpp"
#include "geoxy/vector.hpp"
#include "geoxy/xline.hpp"

#include <doctest/doctest.h>

#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <optional>
#include <type_traits>
#include <variant>
#include <vector>

namespace Geoxy
{
  namespace
  {
    /**
     *  A case of the curves audit's corpus: two curves, the tolerance they were met within, and
     *  where they cross as its exact rational reference finds them, each rounded once.
     */
    struct SwiftOracle
    {
      int                Id     = 0;   /**< The case's number in the corpus.           */
      Curve              First  = {};  /**< The straight curve.                        */
      Curve              Second = {};  /**< The other.                                 */
      double             Within = 0.0; /**< The tolerance's length.                    */
      std::vector<Point> Points = {};  /**< Where they cross, each the double nearest. */
    };

    /**
     *  A double from its bits: a number of a case, written so where it lies near a named constant
     *  that a literal would be taken for.
     *
     *  @param  bits The bits.
     *  @return      The double.
     */
    constexpr auto SwiftBits(std::uint64_t bits) -> double
    {
      return std::bit_cast<double>(bits);
    }

    /**
     *  Whether two meetings are the same to the bit, hit for hit and stretch for stretch.
     *
     *  @param  one   One.
     *  @param  other The other.
     *  @return       Whether they are.
     */
    auto SwiftSame(const Intersections &one, const Intersections &other) -> bool
    {
      const auto hit = [](const Hit &lhs, const Hit &rhs) -> bool
      {
        return lhs.At     == rhs.At
            && lhs.First  == rhs.First
            && lhs.Second == rhs.Second
            && lhs.Kind   == rhs.Kind;
      };
      return std::ranges::equal(one.Hits, other.Hits, hit)
          && one.Overlaps.size() == other.Overlaps.size()
          && one.Spent == other.Spent;
    }

    /**
     *  Two curves met by the swift kernels and by the thorough way.
     */
    struct SwiftMet
    {
      bool Sure   = false;         /**< Whether a kernel was sure. */
      bool Agrees = true;          /**< Whether it found what the thorough
                                        way finds, to the bit.             */
      Intersections Swiftly  = {}; /**< What the kernel found.       */
      Intersections Thorough = {}; /**< What the thorough way found. */
    };

    /**
     *  Two curves met by the swift kernels, within the length Intersect would use, and by the
     *  thorough way where the kernels are sure.
     *
     *  @param  curves    The two curves.
     *  @param  tolerance The tolerance asked for.
     *  @return           What each found.
     */
    auto SwiftBoth(const Curves &curves, const Tolerance &tolerance) -> SwiftMet
    {
      SwiftMet                    met;
      const std::optional<double> length = Reachable(curves, tolerance);
      if (!length)
      {
        return met;
      }
      met.Sure = Swift::Met(curves, {.Length = *length, .Radians = tolerance.Radians}, met.Swiftly);
      if (met.Sure)
      {
        Thorough(curves, tolerance, met.Thorough);
        met.Agrees = SwiftSame(met.Swiftly, met.Thorough);
      }
      return met;
    }

    /**
     *  A straight curve and a Bézier curve, the pair the exact check below takes.
     */
    struct Pairing
    {
      const Curve       *Straight = nullptr; /**< The straight curve.                */
      std::vector<Point> Poles    = {};      /**< The Bézier curve's control points. */
    };

    /**
     *  A Bézier curve's control points.
     *
     *  @param  curve The curve.
     *  @return       Its control points; none for any other curve.
     */
    auto SwiftPoles(const Curve &curve) -> std::vector<Point>
    {
      if (const auto *quadratic = std::get_if<Bezier::Quadratic>(&curve))
      {
        return {quadratic->Start, quadratic->Control, quadratic->End};
      }
      if (const auto *cubic = std::get_if<Bezier::Cubic>(&curve))
      {
        return {cubic->Start, cubic->Out, cubic->In, cubic->End};
      }
      return {};
    }

    /**
     *  A polynomial in Bernstein form at a number, by de Casteljau's algorithm, exactly.
     *
     *  @param  values Its coefficients.
     *  @param  along  Where.
     *  @return        The value.
     */
    auto SwiftValued(std::vector<Expansion::Number> values, const Expansion::Number &along)
        -> Expansion::Number
    {
      const Expansion::Number rest = Expansion::Difference(Expansion::Of(1.0), along);
      for (std::size_t round = 1; round < values.size(); ++round)
      {
        for (std::size_t index = 0; index + round < values.size(); ++index)
        {
          values[index] = Expansion::Sum(Expansion::Product(values[index], rest),
                                         Expansion::Product(values[index + 1], along));
        }
      }
      return values.front();
    }

    /**
     *  The double nearest a quotient of two numbers held exactly.
     *
     *  @param  top   The numerator.
     *  @param  below The denominator, not zero.
     *  @return       The double.
     */
    auto Nearest(const Expansion::Number &top, const Expansion::Number &below) -> double
    {
      return Expansion::Rounded({.Top = top, .Bottom = {}, .Root = {}, .Below = below}, 1);
    }

    /**
     *  Whether a crossing of a straight curve and a Bézier curve is exact, checked in exact
     *  arithmetic: the line's equation on the curve changes sign between the places halfway to the
     *  neighbours of the hit's place on the curve, so the true root lies between them and the place
     *  is the double nearest it; and at both, the curve's point and the line's place round to the
     *  hit's, as they then do at the root between, the curve too short a way across a step of its
     *  parameter to turn back.
     *
     *  @param  pairing The two curves.
     *  @param  hit     The hit.
     *  @param  first   Whether the straight curve is the first.
     *  @return         Whether it is exact.
     */
    auto SwiftProved(const Pairing &pairing, const Hit &hit, bool first) -> bool
    {
      using Expansion::Difference;
      using Expansion::Of;
      using Expansion::Product;
      using Expansion::Sum;
      const Curve                     &line  = *pairing.Straight;
      const double                     place = first ? hit.Second.Along : hit.First.Along;
      const double                     along = first ? hit.First.Along : hit.Second.Along;
      Point                            start = {};
      std::array<Expansion::Number, 2> way   = {};
      if (const auto *segment = std::get_if<Line>(&line))
      {
        start = segment->Start;
        way   = {Difference(Of(segment->End.X), Of(start.X)),
                 Difference(Of(segment->End.Y), Of(start.Y))};
      }
      else
      {
        const Ray ray =
            std::holds_alternative<Ray>(line)
                ? std::get<Ray>(line)
                : Ray{.Start = std::get<Xline>(line).Through, .Way = std::get<Xline>(line).Way};
        start = ray.Start;
        way   = {Of(ray.Way.X), Of(ray.Way.Y)};
      }
      std::vector<Expansion::Number> abscissas;
      std::vector<Expansion::Number> ordinates;
      std::vector<Expansion::Number> equation;
      for (const Point pole : pairing.Poles)
      {
        abscissas.push_back(Of(pole.X));
        ordinates.push_back(Of(pole.Y));
        equation.push_back(Difference(Product(way[0], Difference(Of(pole.Y), Of(start.Y))),
                                      Product(way[1], Difference(Of(pole.X), Of(start.X)))));
      }
      const auto half = [place](double toward) -> Expansion::Number
      {
        return Expansion::Scaled(Sum(Of(place), Of(std::nextafter(place, toward))), 0.5);
      };
      std::array<Expansion::Number, 2> ends  = {half(0.0), half(1.0)};
      const int                        below = Expansion::Sign(SwiftValued(equation, ends[0]));
      if (below == 0 || below != -Expansion::Sign(SwiftValued(equation, ends[1])))
      {
        return false;
      }
      const Expansion::Number square = Sum(Product(way[0], way[0]), Product(way[1], way[1]));
      const auto              round  = [&](const Expansion::Number &end) -> std::array<double, 3>
      {
        const Expansion::Number x     = SwiftValued(abscissas, end);
        const Expansion::Number y     = SwiftValued(ordinates, end);
        const Expansion::Number reach = Sum(Product(Difference(x, Of(start.X)), way[0]),
                                            Product(Difference(y, Of(start.Y)), way[1]));
        return {Nearest(x, Of(1.0)), Nearest(y, Of(1.0)), Nearest(reach, square)};
      };
      /**
       *  Halved about the root until the point and the line's place round alike at both ends.
       */
      for (int step = 0; step < 24; ++step)
      {
        const std::array<double, 3> low  = round(ends[0]);
        const std::array<double, 3> high = round(ends[1]);
        if (low == high)
        {
          return low == std::array<double, 3>{hit.At.X, hit.At.Y, along};
        }
        const Expansion::Number middle = Expansion::Scaled(Sum(ends[0], ends[1]), 0.5);
        ends.at(Expansion::Sign(SwiftValued(equation, middle)) == below ? 0 : 1) = middle;
      }
      return false;
    }

    /**
     *  Whether the kernels' answer for a straight curve and a Bézier curve, where it differs from
     *  the thorough way's, is the better: as many hits, of the same kinds, each proved exact.
     *
     *  @param  curves The two curves.
     *  @param  met    What each found.
     *  @return        Whether it is.
     */
    auto SwiftBetter(const Curves &curves, const SwiftMet &met) -> bool
    {
      const bool           first   = SwiftPoles(curves[0]).empty();
      const Pairing        pairing = {.Straight = &curves.at(first ? 0 : 1).get(),
                                      .Poles    = SwiftPoles(curves.at(first ? 1 : 0))};
      const Intersections &swift   = met.Swiftly;
      return !pairing.Poles.empty()
          && swift.Hits.size() == met.Thorough.Hits.size()
          && std::ranges::equal(swift.Hits, met.Thorough.Hits,
                                [](const Hit &lhs, const Hit &rhs) -> bool
                                { return lhs.Kind == rhs.Kind; })
          && std::ranges::all_of(swift.Hits, [&pairing, first](const Hit &hit) -> bool
                                 { return SwiftProved(pairing, hit, first); });
    }

    /**
     *  Whether the kernels are sure of two curves, within the length Intersect would use.
     *
     *  @param  first     One curve.
     *  @param  second    The other.
     *  @param  tolerance The tolerance.
     *  @return           Whether they are.
     */
    auto SwiftCertain(const Curve &first, const Curve &second, const Tolerance &tolerance) -> bool
    {
      Intersections meeting;
      return Swift::Met({first, second}, Usable({first, second}, tolerance), meeting);
    }

    /**
     *  A double from -1 to 1.
     *
     *  @param  draws Where it comes from.
     *  @return       The double.
     */
    auto SwiftUnit(Draws &draws) -> double
    {
      return draws.Near(0.0, 1.0);
    }

    /**
     *  A choice among some.
     *
     *  @param  draws Where it comes from.
     *  @param  count How many.
     *  @return       Which, from 0.
     */
    auto SwiftPick(Draws &draws, std::uint64_t count) -> std::uint64_t
    {
      return draws.Next() % count;
    }

    /**
     *  A straight curve along a way through a point: a line about it, a ray from behind it, or an
     *  endless line through it.
     *
     *  @param  draws Where the choices come from.
     *  @param  spot  The point.
     *  @param  way   The way, of some length.
     *  @return       The curve.
     */
    auto AlongWay(Draws &draws, Point spot, Vector way) -> Curve
    {
      const Vector unit   = way / Length(way);
      const double before = 0.2 + std::abs(SwiftUnit(draws));
      const double after  = 0.2 + std::abs(SwiftUnit(draws));
      switch (SwiftPick(draws, 4))
      {
      case 0:
        return Ray{.Start = spot - (unit * before),
                   .Way   = unit * std::ldexp(1.0, 2 - static_cast<int>(SwiftPick(draws, 5)))};
      case 1:
        return Xline{.Through = spot - (unit * before), .Way = unit};
      default:
        return Line{.Start = spot - (unit * before), .End = spot + (unit * after)};
      }
    }

    /**
     *  A straight curve between two random points in the unit square's double.
     *
     *  @param  draws Where the numbers come from.
     *  @return       The curve.
     */
    auto SwiftStraight(Draws &draws) -> Curve
    {
      const Point  spot = {.X = 1.5 * SwiftUnit(draws), .Y = 1.5 * SwiftUnit(draws)};
      const Vector way  = {.X = SwiftUnit(draws), .Y = SwiftUnit(draws)};
      return AlongWay(draws, spot, Length(way) > 0.0 ? way : Vector{.X = 1.0, .Y = 0.0});
    }

    /**
     *  An ellipse at random: skewed, mirrored, flat, a circle at times, whole or an arc.
     *
     *  @param  draws Where the numbers come from.
     *  @return       The ellipse.
     */
    auto SwiftOval(Draws &draws) -> Ellipse
    {
      const double wide = 0.2 + (0.8 * std::abs(SwiftUnit(draws)));
      const double high = SwiftPick(draws, 8) == 0
                            ? wide * std::ldexp(1.0, -7 - static_cast<int>(SwiftPick(draws, 12)))
                            : 0.05 + (0.95 * std::abs(SwiftUnit(draws)));
      const double turn = 3.0 * SwiftUnit(draws);
      const double skew = SwiftPick(draws, 3) == 0 ? 0.6 * SwiftUnit(draws) : 0.0;
      Ellipse      oval = {
          .Map   = {.X      = {.X = wide * std::cos(turn), .Y = wide * std::sin(turn)},
                    .Y      = {.X = -high * std::sin(turn + skew), .Y = high * std::cos(turn + skew)},
                    .Origin = {.X = 0.5 * SwiftUnit(draws), .Y = 0.5 * SwiftUnit(draws)}},
          .Start = {},
          .Sweep = Turn};
      if (SwiftPick(draws, 4) == 0)
      {
        std::swap(oval.Map.X, oval.Map.Y);
      }
      if (SwiftPick(draws, 8) == 0)
      {
        oval.Map.Y = Perpendicular(oval.Map.X);
      }
      if (SwiftPick(draws, 2) == 0)
      {
        oval.Start = Degrees(std::round(180.0 * SwiftUnit(draws)));
        oval.Sweep = SwiftPick(draws, 2) == 0 ? Degrees(std::round(350.0 * SwiftUnit(draws)))
                                              : Radians(6.0 * SwiftUnit(draws));
      }
      return oval;
    }

    /**
     *  A polyline at random: two to six vertices, a bulge on two segments in three, closed at
     *  times.
     *
     *  @param  draws Where the numbers come from.
     *  @return       The polyline.
     */
    auto Bends(Draws &draws) -> Polyline
    {
      Polyline          polyline = {.Vertices = {}, .Closed = SwiftPick(draws, 3) == 0};
      const std::size_t count    = 2 + SwiftPick(draws, 5);
      for (std::size_t index = 0; index < count; ++index)
      {
        const double size  = std::pow(10.0, (1.8 * std::abs(SwiftUnit(draws))) - 3.0 + 0.6);
        const double bulge = SwiftPick(draws, 3) == 0 ? 0.0 : std::copysign(size, SwiftUnit(draws));
        polyline.Vertices.push_back(
            {.Position = {.X = SwiftUnit(draws), .Y = SwiftUnit(draws)}, .Bulge = bulge});
      }
      return polyline;
    }

    /**
     *  A quadratic or cubic Bézier curve at random.
     *
     *  @param  draws Where the numbers come from.
     *  @return       The curve.
     */
    auto SwiftBent(Draws &draws) -> Curve
    {
      const auto point = [&draws]() -> Point
      {
        return {.X = SwiftUnit(draws), .Y = SwiftUnit(draws)};
      };
      if (SwiftPick(draws, 3) == 0)
      {
        return Bezier::Quadratic{.Start = point(), .Control = point(), .End = point()};
      }
      return Bezier::Cubic{.Start = point(), .Out = point(), .In = point(), .End = point()};
    }

    /**
     *  A curve moved and scaled: by a power of two and a little more, then far off.
     *
     *  @param  curve The curve.
     *  @param  scale The scale.
     *  @param  by    How far it moves.
     *  @return       The curve, moved.
     */
    auto SwiftMoved(const Curve &curve, double scale, Vector by) -> Curve
    {
      const auto shifted = [scale, by](Point point) -> Point
      {
        return {.X = (point.X * scale) + by.X, .Y = (point.Y * scale) + by.Y};
      };
      return std::visit(
          [&shifted, scale](const auto &shape) -> Curve
          {
            using Shape = std::decay_t<decltype(shape)>;
            if constexpr (std::is_same_v<Shape, Line>)
            {
              return Line{.Start = shifted(shape.Start), .End = shifted(shape.End)};
            }
            else if constexpr (std::is_same_v<Shape, Ray>)
            {
              return Ray{.Start = shifted(shape.Start), .Way = shape.Way * scale};
            }
            else if constexpr (std::is_same_v<Shape, Xline>)
            {
              return Xline{.Through = shifted(shape.Through), .Way = shape.Way * scale};
            }
            else if constexpr (std::is_same_v<Shape, Ellipse>)
            {
              Ellipse moved    = shape;
              moved.Map.X      = shape.Map.X * scale;
              moved.Map.Y      = shape.Map.Y * scale;
              moved.Map.Origin = shifted(shape.Map.Origin);
              return moved;
            }
            else if constexpr (std::is_same_v<Shape, Polyline>)
            {
              Polyline moved = shape;
              for (Vertex &vertex : moved.Vertices)
              {
                vertex.Position = shifted(vertex.Position);
              }
              return moved;
            }
            else if constexpr (std::is_same_v<Shape, Bezier::Quadratic>)
            {
              return Bezier::Quadratic{.Start   = shifted(shape.Start),
                                       .Control = shifted(shape.Control),
                                       .End     = shifted(shape.End)};
            }
            else if constexpr (std::is_same_v<Shape, Bezier::Cubic>)
            {
              return Bezier::Cubic{.Start = shifted(shape.Start),
                                   .Out   = shifted(shape.Out),
                                   .In    = shifted(shape.In),
                                   .End   = shifted(shape.End)};
            }
            else
            {
              return shape;
            }
          },
          curve);
    }

    /**
     *  A straight curve near a curve at an edge of what the kernels decide: along its tangent and a
     *  few tolerances off, through it at a slight angle, through one of its ends, or through it
     *  plainly.
     *
     *  @param  draws Where the numbers come from.
     *  @param  curve The curve.
     *  @param  near  The tolerance's length.
     *  @return       The straight curve.
     */
    auto SwiftEdged(Draws &draws, const Curve &curve, double near) -> Curve
    {
      const Range  domain = Domain(curve);
      const auto   kind   = SwiftPick(draws, 4);
      const double end    = SwiftPick(draws, 2) == 0 ? domain.Start : domain.End;
      const double along =
          kind == 2
              ? end
              : domain.Start + ((domain.End - domain.Start) * (0.5 + (0.45 * SwiftUnit(draws))));
      const Point spot = At(curve, along);
      Vector      way  = Derivative(curve, along);
      if (!(Length(way) > 0.0))
      {
        way = {.X = 1.0, .Y = 0.0};
      }
      if (kind == 0)
      {
        const double off = near * 16.0 * SwiftUnit(draws);
        return AlongWay(draws, spot + (Perpendicular(way) * (off / Length(way))), way);
      }
      const double angle = kind == 1
                             ? std::ldexp(SwiftUnit(draws), -static_cast<int>(SwiftPick(draws, 30)))
                             : 3.0 * SwiftUnit(draws);
      const Vector turned = {.X = (way.X * std::cos(angle)) - (way.Y * std::sin(angle)),
                             .Y = (way.X * std::sin(angle)) + (way.Y * std::cos(angle))};
      return AlongWay(draws, spot, turned);
    }

    /**
     *  A pair of a straight curve and another, and the tolerance asked for.
     */
    struct SwiftPosed
    {
      Curve            First     = {}; /**< The first curve.         */
      Curve            Second    = {}; /**< The second.              */
      Geoxy::Tolerance Tolerance = {}; /**< The tolerance asked for. */
    };

    /**
     *  A straight curve and a curve the maker draws, in one of the ways the kernels must answer as
     *  the thorough way does: at random or at an edge, far off or scaled, under a coarse angle at
     *  times, either way round.
     *
     *  @tparam MakeT The maker's kind.
     *  @param  draws Where the numbers come from.
     *  @param  make  What draws the other curve.
     *  @return       The pair.
     */
    template <typename MakeT>
    auto SwiftPosing(Draws &draws, const MakeT &make) -> SwiftPosed
    {
      const Curve  other  = make(draws);
      const bool   coarse = SwiftPick(draws, 4) == 0;
      const double near =
          coarse ? std::pow(10.0, -3.0 - (9.0 * std::abs(SwiftUnit(draws)))) : 2e-10;
      const bool   wide    = SwiftPick(draws, 6) == 0;
      const double radians = wide ? std::pow(10.0, -6.0 * std::abs(SwiftUnit(draws))) : 1e-10;
      const Curve  line =
          SwiftPick(draws, 2) == 0 ? SwiftStraight(draws) : SwiftEdged(draws, other, near);
      const auto   kind  = SwiftPick(draws, 6);
      const double scale = std::ldexp(kind == 0 ? 1.37 : 1.0,
                                      kind == 0 ? static_cast<int>(SwiftPick(draws, 80)) - 40 : 0);
      const Vector by    = {.X = kind == 1 ? std::ldexp(SwiftUnit(draws), 26) : 0.0,
                            .Y = kind == 1 ? 1e8 : 0.0};
      const bool   swapped = SwiftPick(draws, 2) == 0;
      const Curve  first   = SwiftMoved(line, scale, by);
      const Curve  second  = SwiftMoved(other, scale, by);
      return {.First     = swapped ? second : first,
              .Second    = swapped ? first : second,
              .Tolerance = {.Length = near * scale, .Radians = radians}};
    }

    /**
     *  Pairs of a straight curve and a curve the maker draws, as SwiftPosing draws them: how many
     *  the kernels were sure of, how many they got wrong, and how many they answered better, where
     *  the thorough way's crossing is not the exact one.
     *
     *  @tparam MakeT  What draws the curve.
     *  @param  make   The maker.
     *  @param  rounds How many pairs.
     *  @return        How many the kernels were sure of, how many disagreed, and how many were
     *                 proved better.
     */
    template <typename MakeT>
    auto SwiftDrawn(const MakeT &make, int rounds) -> std::array<int, 3>
    {
      Draws              draws{0x5EED5EED};
      std::array<int, 3> count = {};
      for (int round = 0; round < rounds; ++round)
      {
        const SwiftPosed posed   = SwiftPosing(draws, make);
        const Curves     curves  = {posed.First, posed.Second};
        const SwiftMet   met     = SwiftBoth(curves, posed.Tolerance);
        const bool       agrees  = met.Agrees || SwiftBetter(curves, met);
        count[0]                += static_cast<int>(met.Sure);
        count[1]                += static_cast<int>(!agrees);
        count[2]                += static_cast<int>(!met.Agrees);
      }
      return count;
    }
  } // namespace

  TEST_CASE("swift: the curves corpus, every crossing the double nearest the exact one")
  {
    const std::vector<SwiftOracle> corpus = {
        SwiftOracle{
            .Id    = 88,
            .First = Line{.Start = {.X = -0x1.c47108a239f60p-1, .Y = -0x1.e97519afd1afcp-2},
                          .End   = {.X = 0x1.e1bbf82148004p-2, .Y = -0x1.eae7d13496ea4p-2}},
            .Second =
                Ellipse{.Map = {.X      = {.X = -0x1.4dcfb76307485p-1, .Y = 0x1.6b4048054226bp-2},
                                .Y      = {.X = -0x1.697a5e0c860a3p-4, .Y = -0x1.4c2e9702cb63ep-3},
                                .Origin = {.X = -0x1.eac50862db768p-4, .Y = -0x1.880d7b92732fcp-2}},
                        .Start = Degrees(-75.0),
                        .Sweep = Degrees(-280.0)},
            .Within = 0x1.b7cdfd9d7bdbbp-33,
            .Points = {Point{.X = -0x1.4bbda8f6c6531p-2, .Y = -0x1.ea0e54f50ce06p-2},
                       Point{.X = 0x1.6ebd879a67898p-2, .Y = -0x1.eac912bcff603p-2}}},
        SwiftOracle{
            .Id     = 94,
            .First  = Line{.Start = {.X = 0x1.7356444406d44p-2, .Y = 0x1.b79e841b60430p-2},
                           .End   = {.X = -0x1.2ec0405e80b1cp-1, .Y = -0x1.bb9fccaf16004p-1}},
            .Second = Ellipse{.Map = {.X = {.X = 0x1.5f26810b2589ap-1, .Y = 0x1.1da9a40486268p-1},
                                      .Y = {.X = -0x1.07431936746cdp-1, .Y = 0x1.439d5391852cbp-1},
                                      .Origin = {.X = 0x1.23b3bf0d81578p-4,
                                                 .Y = 0x1.20d57d14185d0p-3}},
                              .Start = Degrees(81.0),
                              .Sweep = Degrees(322.0)},
            .Within = 0x1.b7cdfd9d7bdbbp-33,
            .Points = {Point{.X = -0x1.95007168ce6d1p-2, .Y = -0x1.3373a1641d13cp-1}}},
        SwiftOracle{
            .Id     = 155,
            .First  = Line{.Start = {.X = -0x1.56a7a8157e6d0p-1, .Y = -0x1.da969d5309808p-2},
                           .End   = {.X = 0x1.471108432135cp-1, .Y = -0x1.967bdf6fc1f88p-2}},
            .Second = Ellipse{.Map = {.X = {.X = -0x1.50a8401122af8p-5, .Y = 0x1.a62605dfbd242p-1},
                                      .Y = {.X = 0x1.f5f973dbc4397p-1, .Y = 0x1.90513443792fdp-5},
                                      .Origin = {.X = 0x1.12e84e43072bap-2,
                                                 .Y = -0x1.3e0cf56b5f664p+0}},
                              .Start = Degrees(0.0),
                              .Sweep = Degrees(360.0)},
            .Within = 0x1.b7cdfd9d7bdbbp-33,
            .Points = {Point{.X = 0x1.cce02b9598894p-3, .Y = -0x1.ac06dd6c0feb1p-2},
                       Point{.X = 0x1.d1a68c81c4e01p-3, .Y = -0x1.abe7c9ee0355dp-2}}},
        SwiftOracle{
            .Id     = 178,
            .First  = Line{.Start = {.X = -0x1.7d78403ec0183p+26, .Y = -0x1.01bc38b77a86dp+24},
                           .End   = {.X = -0x1.7d783ff3a006cp+26, .Y = -0x1.01bc371fcd8a0p+24}},
            .Second = Ellipse{.Map = {.X = {.X = -0x1.9a529d54df769p-2, .Y = 0x1.d6d6862990144p-3},
                                      .Y = {.X = -0x1.851633f7a46eep-3, .Y = -0x1.53142473e95afp-2},
                                      .Origin = {.X = -0x1.7d784026ae9e0p+26,
                                                 .Y = -0x1.01bc38298fb3bp+24}},
                              .Start = Degrees(-24.0),
                              .Sweep = Degrees(360.0)},
            .Within = 0x1.0000000000000p-21,
            .Points = {Point{.X = -0x1.7d784033cca1ep+26, .Y = -0x1.01bc387c0d393p+24},
                       Point{.X = -0x1.7d784016cae7ep+26, .Y = -0x1.01bc37dea4ae0p+24}}},
        SwiftOracle{
            .Id    = 184,
            .First = Line{.Start = {.X = -0x1.47309e18e7256p-31, .Y = 0x1.a5a7c15eea974p-31},
                          .End   = {.X = 0x1.e024110f7e6c0p-33, .Y = -0x1.88e48fc407fe0p-32}},
            .Second =
                Ellipse{
                    .Map   = {.X      = {.X = 0x1.f36b7f8c74595p-32, .Y = 0x1.26e2ebd4ffc2bp-32},
                              .Y      = {.X = -0x1.36c745d81b7d1p-32, .Y = 0x1.072ad4cd48a4ap-31},
                              .Origin = {.X = 0x1.bd88fed38e4c0p-35, .Y = -0x1.29573710591a0p-34}},
                    .Start = Degrees(-22.0),
                    .Sweep = Degrees(164.0)},
            .Within = 0x1.b7cdfd9d7bdbbp-63,
            .Points = {Point{.X = -0x1.4bbdfa6b77c17p-32, .Y = 0x1.8d661ca80e412p-32}}},
        SwiftOracle{
            .Id    = 196,
            .First = Line{.Start = {.X = -0x1.d51216402865bp+22, .Y = 0x1.3a26e93fc4203p+22},
                          .End   = {.X = 0x1.0073cad4ebc26p+24, .Y = 0x1.75e2025a0a8a8p+24}},
            .Second =
                Ellipse{.Map = {.X = {.X = -0x1.53f20457c6c8dp+24, .Y = 0x1.00ccd1583bc1cp+21},
                                .Y = {.X = -0x1.539268d679045p+20, .Y = -0x1.c184420153bb7p+23},
                                .Origin = {.X = 0x1.b88c3238fd674p+24, .Y = 0x1.9b99b86716390p+24}},
                        .Start = Degrees(-105.0),
                        .Sweep = Degrees(360.0)},
            .Within = 0x1.8fc504816f007p-7,
            .Points = {Point{.X = 0x1.4080efaf6ab68p+23, .Y = 0x1.29db6ef78073fp+24}}},
        SwiftOracle{.Id    = 209,
                    .First = Line{.Start = {.X = 0x1.0d12d4cd89988p-3, .Y = -0x1.8d446ff378b2ep-1},
                                  .End   = {.X = 0x1.7b5f60527e850p-2, .Y = 0x1.afe7653d51170p-3}},
                    .Second =
                        Bezier::Quadratic{
                            .Start   = {.X = 0x1.19860db21bbecp-2, .Y = 0x1.81b59092a4bc4p-2},
                            .Control = {.X = 0x1.b62a5f7195b92p-1, .Y = -0x1.af32b3a24ad60p-4},
                            .End     = {.X = 0x1.c6aa56ff286a8p-3, .Y = 0x1.e79a62bffbb60p-5}},
                    .Within = 0x1.b7cdfd9d7bdbbp-33,
                    .Points = {Point{.X = 0x1.4fad8c7d70e5ap-2, .Y = 0x1.1cec458c71c73p-5}}},
        SwiftOracle{.Id    = 210,
                    .First = Line{.Start = {.X = 0x1.6ca437bf60010p-3, .Y = 0x1.6f4b720b2eda0p-2},
                                  .End = {.X = -0x1.3f8661ea4a092p-1, .Y = -0x1.c71f043101a2ep-1}},
                    .Second =
                        Bezier::Quadratic{
                            .Start   = {.X = -0x1.88e06698c984ep-1, .Y = -0x1.d443def0fc8cap-1},
                            .Control = {.X = 0x1.c3547a1e7fa20p-4, .Y = -0x1.8f23063f3bf78p-2},
                            .End     = {.X = 0x1.237a6b99aa99ap-1, .Y = -0x1.5a36ff3bbc7d0p-1}},
                    .Within = 0x1.b7cdfd9d7bdbbp-33,
                    .Points = {Point{.X = -0x1.2379255407360p-1, .Y = -0x1.9b7d8530a1567p-1}}},
        SwiftOracle{.Id    = 279,
                    .First = Line{.Start = {.X = 0x1.fb84b8e146660p-3, .Y = -0x1.8d6073dff54acp-1},
                                  .End = {.X = -0x1.c9b5df3c7f442p-1, .Y = -0x1.a327bf48528b4p-2}},
                    .Second =
                        Bezier::Quadratic{
                            .Start   = {.X = 0x1.97a7a9bc26df6p-5, .Y = -0x1.6c6d6524a4d17p-1},
                            .Control = {.X = -0x1.5753600db552fp-2, .Y = -0x1.2da095fd62de4p-1},
                            .End     = {.X = -0x1.709f29db409f9p-1, .Y = -0x1.dc83decf48deap-2}},
                    .Within = 0x1.b7cdfd9d7bdbbp-33,
                    .Points = {Point{.X = -0x1.573c07a699bf9p-2, .Y = -0x1.2d7c2021c3bfcp-1}}},
        SwiftOracle{.Id    = 301,
                    .First = Line{.Start = {.X = 0x1.7d78401ebe35fp+26, .Y = 0x1.f700b2b0d6aeap+25},
                                  .End = {.X = 0x1.7d783fff1fef3p+26, .Y = 0x1.f700b31128be6p+25}},
                    .Second =
                        Bezier::Quadratic{
                            .Start   = {.X = 0x1.7d783fff1dbc0p+26, .Y = 0x1.f700b2f319a33p+25},
                            .Control = {.X = 0x1.7d78403f4b3eap+26, .Y = 0x1.f700b2ea6401fp+25},
                            .End     = {.X = 0x1.7d784029334d2p+26, .Y = 0x1.f700b2efee271p+25}},
                    .Within = 0x1.0000000000000p-21,
                    .Points = {Point{.X = 0x1.7d784009722f3p+26, .Y = 0x1.f700b2f1b7813p+25}}},
        SwiftOracle{.Id = 309,
                    .First =
                        Line{.Start = {.X = 0x1.fde01bcd44910p-34, .Y = -0x1.bb5bef55ba598p-33},
                             .End   = {.X = 0x1.5b741091356dcp-32, .Y = 0x1.fe104a1779baep-31}},
                    .Second =
                        Bezier::Quadratic{
                            .Start   = {.X = -0x1.ca6a541dfa432p-31, .Y = 0x1.ab8b99863dff8p-31},
                            .Control = {.X = 0x1.f5ce122f034bcp-31, .Y = 0x1.73b11a4f3cfaap-31},
                            .End     = {.X = -0x1.2c19335a9a8c0p-34, .Y = 0x1.4e5661b84dc0cp-32}},
                    .Within = 0x1.b7cdfd9d7bdbbp-63,
                    .Points = {Point{.X = 0x1.fee6e91af41f2p-33, .Y = 0x1.f4c6b85075c54p-32},
                               Point{.X = 0x1.19a7810e79968p-32, .Y = 0x1.4458bd46c7b2dp-31}}},
        SwiftOracle{
            .Id    = 337,
            .First = Line{.Start = {.X = 0x1.c7a4771b39060p-1, .Y = 0x1.df330c2ba0a5ap-1},
                          .End   = {.X = -0x1.2400118ba7358p-1, .Y = -0x1.a1a4f0b6db5f0p-2}},
            .Second =
                Bezier::Cubic{.Start = {.X = 0x1.7312303b884bcp-1, .Y = -0x1.1dffd936e1490p-4},
                              .Out   = {.X = 0x1.dd8a324ae5d40p-4, .Y = 0x1.b6080a826a498p-2},
                              .In    = {.X = -0x1.288786a6f7314p-1, .Y = -0x1.3073627e7c940p-4},
                              .End   = {.X = -0x1.aa7859a20fbc8p-1, .Y = -0x1.e069c9768e0f2p-1}},
            .Within = 0x1.b7cdfd9d7bdbbp-33,
            .Points = {Point{.X = -0x1.5e34486bda790p-5, .Y = 0x1.3e00306542a6fp-4}}},
        SwiftOracle{
            .Id    = 338,
            .First = Line{.Start = {.X = 0x1.df9b78a682da2p-1, .Y = 0x1.c6ceb9fe489f4p-2},
                          .End   = {.X = 0x1.f335d6e613d38p-3, .Y = -0x1.2e3c10a3ce2f6p-1}},
            .Second =
                Bezier::Cubic{.Start = {.X = -0x1.a34be62f6a2c8p-2, .Y = 0x1.adb12bb50961ap-1},
                              .Out   = {.X = 0x1.f9bf1617ac340p-1, .Y = 0x1.23ecf09b5095cp-2},
                              .In    = {.X = 0x1.b19b0c5536340p-1, .Y = -0x1.af475b3ffe9c8p-3},
                              .End   = {.X = 0x1.600967ed0ecf4p-1, .Y = 0x1.d2df4f76897c0p-2}},
            .Within = 0x1.b7cdfd9d7bdbbp-33,
            .Points = {Point{.X = 0x1.88cab0905c3f0p-1, .Y = 0x1.873c5f8d0c097p-3},
                       Point{.X = 0x1.7fca05de976ddp-1, .Y = 0x1.517b1d19c2c3cp-3}}},
        SwiftOracle{
            .Id     = 403,
            .First  = Line{.Start = {.X = -0x1.c836f59c42340p-6, .Y = -0x1.c06c301cc77c8p-3},
                           .End   = {.X = -0x1.7d2565df26f84p-2, .Y = -0x1.b13ec7635189cp-2}},
            .Second = Bezier::Cubic{.Start = {.X = 0x1.88d951ec7c37bp-2, .Y = 0x1.9c557c69aaa0bp-6},
                                    .Out = {.X = 0x1.acffa5edfe1d6p-5, .Y = -0x1.5e43853c0536cp-3},
                                    .In  = {.X = -0x1.1d514a869a9dcp-2, .Y = -0x1.78823eb62a41ep-2},
                                    .End = {.X = -0x1.37fd35f049764p-1,
                                            .Y = SwiftBits(0xbfe212e0e40ee24a)}},
            .Within = 0x1.b7cdfd9d7bdbbp-33,
            .Points = {Point{.X = -0x1.cf1a8da7539ffp-4, .Y = -0x1.13f05916f906ap-2}}},
        SwiftOracle{
            .Id    = 430,
            .First = Line{.Start = {.X = -0x1.7d78400d3dd35p+26, .Y = 0x1.bc76a0b1b0e2dp+20},
                          .End   = {.X = -0x1.7d7840009eb38p+26, .Y = 0x1.bc76b6c07c13ap+20}},
            .Second =
                Bezier::Cubic{.Start = {.X = -0x1.7d783fd230991p+26, .Y = 0x1.bc76a5c9e3c01p+20},
                              .Out   = {.X = -0x1.7d784000940d9p+26, .Y = 0x1.bc76b8725e2dep+20},
                              .In    = {.X = -0x1.7d78402e4c220p+26, .Y = 0x1.bc76a7ea51fecp+20},
                              .End   = {.X = -0x1.7d78401282205p+26, .Y = 0x1.bc76b4c1ce906p+20}},
            .Within = 0x1.0000000000000p-21,
            .Points = {Point{.X = -0x1.7d784004d1006p+26, .Y = 0x1.bc76af6b02ad6p+20}}},
        SwiftOracle{
            .Id    = 434,
            .First = Line{.Start = {.X = 0x1.553cf58f79698p-32, .Y = -0x1.64db481279e1cp-31},
                          .End   = {.X = -0x1.e9a9108218f38p-31, .Y = 0x1.01ed4d76ade74p-31}},
            .Second =
                Bezier::Cubic{.Start = {.X = -0x1.e7f989ea894d8p-33, .Y = -0x1.7852ac555075ap-31},
                              .Out   = {.X = 0x1.5c95f76217bfcp-32, .Y = 0x1.c3492f371eda0p-34},
                              .In    = {.X = 0x1.4d6083894a6c0p-31, .Y = 0x1.a28a318ea6360p-33},
                              .End   = {.X = -0x1.9611cac490b90p-34, .Y = -0x1.6083d0cb0d4d4p-32}},
            .Within = 0x1.b7cdfd9d7bdbbp-63,
            .Points = {Point{.X = 0x1.27be80126f5fbp-37, .Y = -0x1.9497c85320838p-32},
                       Point{.X = -0x1.1b6847d5e7b37p-34, .Y = -0x1.4a0491d947156p-32}}},
        SwiftOracle{.Id     = 594,
                    .First  = Line{.Start = {.X = -0x1.7ed78be26220ep-1, .Y = 0x1.949938358e6eap-1},
                                   .End   = {.X = 0x1.ba4c91f2b2390p-4, .Y = -0x1.4eac73bed49cap-1}},
                    .Second = Polyline{.Vertices = {{.Position = {.X = -0x1.807b379693696p-1,
                                                                  .Y = 0x1.554a1c323bc40p-3},
                                                     .Bulge    = -0x1.73a116048e25dp+0},
                                                    {.Position = {.X = -0x1.4967750fa2cbap-1,
                                                                  .Y = -0x1.aa62d7ba276ecp-1},
                                                     .Bulge    = 0x0.0p+0}},
                                       .Closed   = false},
                    .Within = 0x1.b7cdfd9d7bdbbp-33,
                    .Points = {Point{.X = -0x1.a14caf9785149p-2, .Y = 0x1.bab5b8c091de2p-3},
                               Point{.X = 0x1.39c46bfc55346p-8, .Y = -0x1.eb0aec67b0710p-2}}},
        SwiftOracle{
            .Id    = 596,
            .First = Line{.Start = {.X = 0x1.b1c08f1faf5c8p-1, .Y = -0x1.8a0f3fa889560p-3},
                          .End   = {.X = 0x1.b7bb6a4abf4ccp-2, .Y = -0x1.3b61c5612fec4p-2}},
            .Second =
                Polyline{.Vertices =
                             {{.Position = {.X = 0x1.a2c1788ebebc4p-1, .Y = 0x1.b5bdd3cec1426p-1},
                               .Bulge    = 0x0.0p+0},
                              {.Position = {.X = -0x1.f2990cf2f6a78p-1, .Y = -0x1.26451727f6026p-1},
                               .Bulge    = 0x0.0p+0},
                              {.Position = {.X = -0x1.e23843b704ca6p-1, .Y = -0x1.1049eca54e3b0p-2},
                               .Bulge    = 0x0.0p+0},
                              {.Position = {.X = 0x1.d12491803dc08p-3, .Y = -0x1.56479a643cbcep-1},
                               .Bulge    = 0x0.0p+0},
                              {.Position = {.X = 0x1.ebf6e173042a0p-5, .Y = -0x1.3d8612f868d70p-4},
                               .Bulge    = -0x1.87a3f14d7c140p-5},
                              {.Position = {.X = 0x1.b4bf360269fbcp-1, .Y = -0x1.f37141f0271a6p-1},
                               .Bulge    = 0x0.0p+0}},
                         .Closed = true},
            .Within = 0x1.b7cdfd9d7bdbbp-33,
            .Points = {Point{.X = 0x1.ad133bc890c4ep-1, .Y = -0x1.8f3c4f249d18ap-3}}},
        SwiftOracle{
            .Id    = 663,
            .First = Line{.Start = {.X = -0x1.3bd6fe59cbd8cp-1, .Y = 0x1.7089456c452d2p-1},
                          .End   = {.X = -0x1.f88075c2195aep-1, .Y = 0x1.fab6fb06edda4p-1}},
            .Second =
                Polyline{.Vertices =
                             {{.Position = {.X = -0x1.4b0abea8da9a4p-1, .Y = 0x1.0de03f64e413ep+0},
                               .Bulge    = 0x0.0p+0},
                              {.Position = {.X = -0x1.e32482bde4e08p-1, .Y = 0x1.032bdf8563ea3p+0},
                               .Bulge    = 0x1.7ee7d2ce80489p-3},
                              {.Position = {.X = -0x1.241b41be479aep-1, .Y = 0x1.7ab7f83e28c9bp-1},
                               .Bulge    = 0x0.0p+0}},
                         .Closed = false},
            .Within = 0x1.b7cdfd9d7bdbbp-33,
            .Points = {Point{.X = -0x1.90ad4feacc8f7p-1, .Y = 0x1.aeac0eca0bddap-1},
                       Point{.X = -0x1.9134925b186a1p-1, .Y = 0x1.af0f1fa7c35a6p-1}}},
        SwiftOracle{
            .Id    = 685,
            .First = Line{.Start = {.X = -0x1.7d784035a5a84p+26, .Y = -0x1.127767db6a4a4p+26},
                          .End   = {.X = -0x1.7d783fc8abce1p+26, .Y = -0x1.127767e067dd5p+26}},
            .Second =
                Polyline{
                    .Vertices =
                        {{.Position = {.X = -0x1.7d7840349de72p+26, .Y = -0x1.1277678da895fp+26},
                          .Bulge    = 0x1.83fd0e0ba9d00p-3},
                         {.Position = {.X = -0x1.7d783fec4977ep+26, .Y = -0x1.127767f2506f0p+26},
                          .Bulge    = -0x1.38ba136cd2204p+0},
                         {.Position = {.X = -0x1.7d78400d7eee0p+26, .Y = -0x1.1277679e3da8fp+26},
                          .Bulge    = 0x1.f2a31ae21febcp-1}},
                    .Closed = false},
            .Within = 0x1.0000000000000p-21,
            .Points = {Point{.X = -0x1.7d78403097d22p+26, .Y = -0x1.127767dba58a8p+26},
                       Point{.X = -0x1.7d7840065d7e3p+26, .Y = -0x1.127767dd9498dp+26}}},
        SwiftOracle{
            .Id    = 692,
            .First = Line{.Start = {.X = 0x1.e29cc7327b578p-31, .Y = 0x1.c0438ecb9a086p-31},
                          .End   = {.X = -0x1.2aa51459c6d62p-31, .Y = 0x1.ebfe75cc5de18p-31}},
            .Second =
                Polyline{
                    .Vertices =
                        {{.Position = {.X = -0x1.74003be86e4e4p-31, .Y = 0x1.ee8de19667be6p-31},
                          .Bulge    = -0x1.07e1f79029f40p+0},
                         {.Position = {.X = 0x1.f848df5315a40p-32, .Y = -0x1.53e93bb0f7964p-31},
                          .Bulge    = 0x1.a2391a1ae3c28p-1},
                         {.Position = {.X = 0x1.5cd03aff3af06p-31, .Y = 0x1.8bae3efe01920p-31},
                          .Bulge    = 0x0.0p+0},
                         {.Position = {.X = 0x1.b13a75b3aaa28p-31, .Y = 0x1.31ef10c478b04p-31},
                          .Bulge    = 0x0.0p+0}},
                    .Closed = false},
            .Within = 0x1.b7cdfd9d7bdbbp-63,
            .Points = {Point{.X = 0x1.3f482ffd5ec22p-31, .Y = 0x1.c967f9b418e41p-31}}},
        SwiftOracle{
            .Id    = 704,
            .First = Line{.Start = {.X = -0x1.a8f5b307faf60p+24, .Y = -0x1.2e89906750f23p+25},
                          .End   = {.X = -0x1.daad694beac35p+20, .Y = 0x1.87252626787bap+25}},
            .Second =
                Polyline{
                    .Vertices =
                        {{.Position = {.X = -0x1.328993fb70ef9p+22, .Y = -0x1.25249b6ef1e2fp+24},
                          .Bulge    = 0x0.0p+0},
                         {.Position = {.X = -0x1.22340229d48fdp+24, .Y = 0x1.7a0d6c84a8969p+25},
                          .Bulge    = 0x0.0p+0},
                         {.Position = {.X = -0x1.1b233ff6fe6edp+22, .Y = -0x1.45bbdd1322c53p+24},
                          .Bulge    = 0x0.0p+0}},
                    .Closed = false},
            .Within = 0x1.8fc504816f007p-7,
            .Points = {Point{.X = -0x1.71ee98c315b0ap+23, .Y = 0x1.da78ab019d015p+23},
                       Point{.X = -0x1.717e55187f4e3p+23, .Y = 0x1.dc02af217ba69p+23}}},
        SwiftOracle{
            .Id     = 844,
            .First  = Ray{.Start = {.X = -0x1.b8c9b04257f26p-1, .Y = -0x1.e57cf9f9b6820p-4},
                          .Way   = {.X = 0x1.e7f9807d4a287p-1, .Y = 0x1.8923d8253dad0p-4}},
            .Second = Ellipse{.Map = {.X = {.X = 0x1.9837e9562051ep-4, .Y = 0x1.c679c307ba2d1p-1},
                                      .Y = {.X = -0x1.491137f2adbd8p-2, .Y = 0x1.279304960b385p-5},
                                      .Origin = {.X = 0x1.4e681f71d9708p-4,
                                                 .Y = 0x1.079010e3e55bep-2}},
                              .Start = Degrees(168.0),
                              .Sweep = Degrees(-325.0)},
            .Within = 0x1.b7cdfd9d7bdbbp-33,
            .Points = {Point{.X = -0x1.032310a60daf9p-2, .Y = -0x1.d58124276f662p-5},
                       Point{.X = 0x1.79b7fc49bc219p-2, .Y = 0x1.5ca91c0c1b6f0p-8}}},
        SwiftOracle{
            .Id     = 846,
            .First  = Ray{.Start = {.X = 0x1.82f2ac4bd5ed0p-1, .Y = -0x1.51c1f3189160fp-1},
                          .Way   = {.X = -0x1.c237dd1157050p-1, .Y = 0x1.8b9aa0f16d6abp-1}},
            .Second = Ellipse{.Map = {.X = {.X = 0x1.bdb2ab8e1dc5dp-3, .Y = 0x1.4a67f66268bcap-3},
                                      .Y = {.X = -0x1.c5a8427da4938p-4, .Y = 0x1.31fa68c57499cp-3},
                                      .Origin = {.X = -0x1.1176a783a27b8p-4,
                                                 .Y = -0x1.12c52fc57494cp-3}},
                              .Start = Degrees(50.0),
                              .Sweep = Degrees(-205.0)},
            .Within = 0x1.b7cdfd9d7bdbbp-33,
            .Points = {Point{.X = 0x1.4a1469af293fcp-3, .Y = -0x1.1908c6d62912ap-3}}},
        SwiftOracle{
            .Id     = 913,
            .First  = Ray{.Start = {.X = 0x1.a404c1c6775e0p-2, .Y = 0x1.60ce00c4deff0p+0},
                          .Way   = {.X = -0x1.e69b6cf1e8d98p-3, .Y = -0x1.0850cd60170b5p+0}},
            .Second = Ellipse{.Map = {.X = {.X = -0x1.ea0d73c9a3425p-2, .Y = 0x1.c528e65b1ad0fp-4},
                                      .Y = {.X = -0x1.e4553c41d8adbp-3, .Y = -0x1.05e1be2176400p+0},
                                      .Origin = {.X = 0x1.643ea196e3bc3p-1,
                                                 .Y = 0x1.b7590208f0890p-2}},
                              .Start = Degrees(0.0),
                              .Sweep = Degrees(360.0)},
            .Within = 0x1.b7cdfd9d7bdbbp-33,
            .Points = {Point{.X = 0x1.bcdf9ec848741p-3, .Y = 0x1.14519dcfdba73p-1},
                       Point{.X = 0x1.bacfd75a06174p-3, .Y = 0x1.121441f38b5a9p-1}}},
        SwiftOracle{
            .Id     = 936,
            .First  = Ray{.Start = {.X = -0x1.7d783ff20b035p+26, .Y = -0x1.67fea8bd88b7cp+25},
                          .Way   = {.X = -0x1.a85c753c1868ep-2, .Y = -0x1.8bd7ae9fb2028p-4}},
            .Second = Ellipse{.Map = {.X = {.X = -0x1.855c11fdc8adap-3, .Y = 0x1.e28fea23ffdd3p-1},
                                      .Y = {.X = -0x1.f14ac7d0af686p-1, .Y = -0x1.913ea0a78b6bap-3},
                                      .Origin = {.X = -0x1.7d784014e45bap+26,
                                                 .Y = -0x1.67fea8ae27d24p+25}},
                              .Start = Degrees(-49.0),
                              .Sweep = Degrees(360.0)},
            .Within = 0x1.0000000000000p-21,
            .Points = {Point{.X = -0x1.7d78404d27b4fp+26, .Y = -0x1.67fea8e8074fcp+25}}},
        SwiftOracle{.Id    = 942,
                    .First = Ray{.Start = {.X = 0x1.c1dba09a837a0p-34, .Y = -0x1.d1ea026c59610p-33},
                                 .Way   = {.X = 0x1.bff036fd33a78p-33, .Y = 0x1.e6f4ae9166f30p-34}},
                    .Second =
                        Ellipse{
                            .Map = {.X = {.X = 0x1.8e7cafb671553p-33, .Y = 0x1.4fac8936eeacap-31},
                                    .Y = {.X = -0x1.916bafc31e06dp-31, .Y = 0x1.dc895eee22ceep-33},
                                    .Origin = {.X = 0x1.089dd4fab9480p-37,
                                               .Y = 0x1.870dabe80326cp-33}},
                            .Start = Degrees(131.0),
                            .Sweep = Degrees(360.0)},
                    .Within = 0x1.b7cdfd9d7bdbbp-63,
                    .Points = {Point{.X = 0x1.a14563af4d44bp-31, .Y = 0x1.3f0f4b53f4118p-33}}},
        SwiftOracle{.Id    = 954,
                    .First = Ray{.Start = {.X = 0x1.7bd8af7cc1468p+25, .Y = -0x1.42b8192857573p+24},
                                 .Way = {.X = -0x1.90b2422436169p+25, .Y = 0x1.26fd3e9e3b3a0p+24}},
                    .Second =
                        Ellipse{
                            .Map = {.X = {.X = -0x1.46750f631fa5dp+24, .Y = 0x1.87dc83c24ea77p+23},
                                    .Y = {.X = -0x1.f8982c0f84689p+23, .Y = -0x1.a45fdc268cb8bp+24},
                                    .Origin = {.X = 0x1.b05090a0f9ea7p+24,
                                               .Y = -0x1.e566dc8754acbp+22}},
                            .Start = Degrees(76.0),
                            .Sweep = Degrees(267.0)},
                    .Within = 0x1.8fc504816f007p-7,
                    .Points = {Point{.X = 0x1.7a1fac9747c1fp+25, .Y = -0x1.41736df25c710p+24}}},
        SwiftOracle{
            .Id    = 969,
            .First = Ray{.Start = {.X = -0x1.24338b6e10a8ap-1, .Y = -0x1.1e0fc015e2992p+0},
                         .Way   = {.X = 0x1.a95b6cd8b69dcp-3, .Y = 0x1.98c285ccf8e07p+0}},
            .Second =
                Bezier::Quadratic{.Start   = {.X = 0x1.ffecf6da6147ap-1, .Y = 0x1.fee0dd594d528p-1},
                                  .Control = {.X = 0x1.4a48e84afecf4p-2, .Y = 0x1.d47f1deb727c0p-6},
                                  .End = {.X = -0x1.193a23c965a8ap-1, .Y = 0x1.436729614c2d0p-2}},
            .Within = 0x1.b7cdfd9d7bdbbp-33,
            .Points = {Point{.X = -0x1.8f23c4e3f9ddfp-2, .Y = 0x1.1806ca57dc3fep-2}}},
        SwiftOracle{
            .Id    = 972,
            .First = Ray{.Start = {.X = -0x1.8d6072f9f1622p-1, .Y = -0x1.bd393e75a6e80p-4},
                         .Way   = {.X = 0x1.356e14293ed7fp-1, .Y = 0x1.08c0dce7d40f4p-2}},
            .Second =
                Bezier::Quadratic{
                    .Start   = {.X = SwiftBits(0xbfe20cd5624bd8fc), .Y = 0x1.8d26772cd0584p-2},
                    .Control = {.X = 0x1.6997164530e68p-3, .Y = 0x1.8669c87488c70p-4},
                    .End     = {.X = -0x1.00190816b8f18p-3, .Y = -0x1.69ca1cf11aadcp-2}},
            .Within = 0x1.b7cdfd9d7bdbbp-33,
            .Points = {Point{.X = -0x1.4ceca9ed462f2p-3, .Y = 0x1.3af6e95db701bp-3}}},
        SwiftOracle{.Id    = 1037,
                    .First = Ray{.Start = {.X = -0x1.7b4ec227c3e00p-5, .Y = -0x1.d94c7f1255148p-2},
                                 .Way   = {.X = -0x1.bf346ad98deecp-3, .Y = 0x1.4137114109130p-3}},
                    .Second =
                        Bezier::Quadratic{
                            .Start   = {.X = 0x1.3fe346fb80357p-2, .Y = -0x1.700d20b1ac2a1p-1},
                            .Control = {.X = -0x1.4b28d9710e498p-3, .Y = -0x1.84b4fc14d5816p-2},
                            .End     = {.X = -0x1.45414cdb57639p-1, .Y = -0x1.447ebb3cd7f9ep-5}},
                    .Within = 0x1.b7cdfd9d7bdbbp-33,
                    .Points = {Point{.X = -0x1.4ae416161ea1cp-3, .Y = -0x1.8485043d2782dp-2}}},
        SwiftOracle{.Id = 1060,
                    .First =
                        Ray{.Start = {.X = -0x1.7d783fc7cdbe0p+26, .Y = -0x1.c0a404192dc86p+24},
                            .Way   = {.X = -0x1.252cf455b6db6p+0, .Y = -0x1.77c1fc417998ap+0}},
                    .Second =
                        Bezier::Quadratic{
                            .Start   = {.X = -0x1.7d783ff52a739p+26, .Y = -0x1.c0a405c7b540cp+24},
                            .Control = {.X = -0x1.7d783ff94fb75p+26, .Y = -0x1.c0a406375214dp+24},
                            .End     = {.X = -0x1.7d784046a7845p+26, .Y = -0x1.c0a404fcd7bc2p+24}},
                    .Within = 0x1.0000000000000p-21,
                    .Points = {Point{.X = -0x1.7d784016af170p+26, .Y = -0x1.c0a405ad938b9p+24}}},
        SwiftOracle{.Id    = 1066,
                    .First = Ray{.Start = {.X = 0x1.497e066d936d4p-30, .Y = -0x1.de8237a7b3470p-33},
                                 .Way = {.X = -0x1.bdc9b3f1757acp-30, .Y = -0x1.0d3df63afafe4p-32}},
                    .Second =
                        Bezier::Quadratic{
                            .Start   = {.X = -0x1.a88417d86a450p-33, .Y = -0x1.e42cd379cb7eap-31},
                            .Control = {.X = -0x1.50bfbcbe2bb8ap-31, .Y = 0x1.546a308e51506p-31},
                            .End     = {.X = 0x1.7a8158005872cp-32, .Y = 0x1.5c881fe9a7ecep-31}},
                    .Within = 0x1.b7cdfd9d7bdbbp-63,
                    .Points = {Point{.X = -0x1.407214e030062p-32, .Y = -0x1.e6a4694e25567p-32}}},
        SwiftOracle{.Id    = 1080,
                    .First = Ray{.Start = {.X = 0x1.2a05a427ba2adp+22, .Y = 0x1.6e80d1c2325e0p+25},
                                 .Way = {.X = -0x1.9d925d5b959cfp+23, .Y = -0x1.674d5cab59ed4p+25}},
                    .Second =
                        Bezier::Quadratic{
                            .Start   = {.X = -0x1.6b94d8c5b6a20p+23, .Y = 0x1.3616f6efa6b6cp+24},
                            .Control = {.X = -0x1.5d516cde3adf7p+25, .Y = 0x1.9f052c979eb22p+24},
                            .End     = {.X = 0x1.b3db9b927553ep+25, .Y = -0x1.a59a7e0af17c0p+25}},
                    .Within = 0x1.8fc504816f007p-7,
                    .Points = {Point{.X = -0x1.030765f087901p+23, .Y = 0x1.8038c7496ab3ep+20}}},
        SwiftOracle{
            .Id    = 1092,
            .First = Ray{.Start = {.X = 0x1.02b0004123be8p-1, .Y = -0x1.2714b899d79ccp-2},
                         .Way   = {.X = -0x1.f971a8178412bp-1, .Y = 0x1.161f8f771fa4fp-1}},
            .Second =
                Bezier::Cubic{.Start = {.X = -0x1.15241398a7c30p-4, .Y = -0x1.902b791efa180p-1},
                              .Out   = {.X = 0x1.dedc72a6fb2b4p-1, .Y = 0x1.fd818345ff542p-1},
                              .In    = {.X = 0x1.b533c64d90c36p-1, .Y = -0x1.f93729c4be4ccp-2},
                              .End   = {.X = -0x1.4c4098642d920p-1, .Y = -0x1.605547ac8e058p-1}},
            .Within = 0x1.b7cdfd9d7bdbbp-33,
            .Points = {Point{.X = 0x1.95528126c5e58p-2, .Y = -0x1.d2d8bdd015d4bp-3},
                       Point{.X = 0x1.40650e0dc2acap-2, .Y = -0x1.75620a67bd438p-3}}},
        SwiftOracle{
            .Id    = 1093,
            .First = Ray{.Start = {.X = SwiftBits(0xbfe276451a1fee28), .Y = 0x1.f340ee811f13cp-1},
                         .Way   = {.X = 0x1.8bf7746e54051p-1, .Y = -0x1.5a2b3d6af79c4p+0}},
            .Second =
                Bezier::Cubic{.Start = {.X = -0x1.507713426ea88p-1, .Y = -0x1.f23911e5f01e0p-5},
                              .Out   = {.X = -0x1.4bdf484f386e8p-3, .Y = 0x1.28baf5178db7ep-1},
                              .In    = {.X = -0x1.2304afce3a650p-3, .Y = -0x1.c3cc941c242e6p-1},
                              .End   = {.X = 0x1.463c6348d8284p-1, .Y = -0x1.ee99b81c995a8p-2}},
            .Within = 0x1.b7cdfd9d7bdbbp-33,
            .Points = {Point{.X = 0x1.2a4445bff4760p-2, .Y = -0x1.15fcd80224b69p-1}}},
        SwiftOracle{
            .Id    = 1161,
            .First = Ray{.Start = {.X = 0x1.1b5aa78acfb40p-5, .Y = -0x1.edc21a31d734cp-1},
                         .Way   = {.X = -0x1.04ff96f79d4dcp-1, .Y = 0x1.0de3fa39d3a3ep-1}},
            .Second =
                Bezier::Cubic{.Start = {.X = 0x1.9775f2b36d625p-3, .Y = -0x1.22b2dee1142c2p+0},
                              .Out   = {.X = -0x1.2cc6ba64286ffp-4, .Y = -0x1.b4a6a1670875dp-1},
                              .In    = {.X = -0x1.625cc218ddfa3p-2, .Y = -0x1.2405a4226b509p-1},
                              .End   = {.X = -0x1.3ce32092e274ep-1, .Y = -0x1.27058be8a1d10p-2}},
            .Within = 0x1.b7cdfd9d7bdbbp-33,
            .Points = {Point{.X = -0x1.adada678710efp-3, .Y = -0x1.6c5daa8a5ab7cp-1}}},
        SwiftOracle{
            .Id    = 1183,
            .First = Ray{.Start = {.X = -0x1.7d7840514e457p+26, .Y = -0x1.1ee0bae125941p+23},
                         .Way   = {.X = 0x1.59bfb10f4799cp+0, .Y = 0x1.78db4298d4275p+0}},
            .Second =
                Bezier::Cubic{.Start = {.X = -0x1.7d7840342ef2cp+26, .Y = -0x1.1ee0b9be728e1p+23},
                              .Out   = {.X = -0x1.7d78403956585p+26, .Y = -0x1.1ee0b9e262424p+23},
                              .In    = {.X = -0x1.7d783fcc2790dp+26, .Y = -0x1.1ee0b8445f884p+23},
                              .End   = {.X = -0x1.7d7840206e712p+26, .Y = -0x1.1ee0b889921c7p+23}},
            .Within = 0x1.0000000000000p-21,
            .Points = {Point{.X = -0x1.7d78402de9e45p+26, .Y = -0x1.1ee0b9ac891aap+23},
                       Point{.X = -0x1.7d78400bc370dp+26, .Y = -0x1.1ee0b882c0e44p+23}}},
        SwiftOracle{
            .Id    = 1190,
            .First = Ray{.Start = {.X = -0x1.b7b3930d21b69p-31, .Y = -0x1.d2878e4b956b8p-32},
                         .Way   = {.X = 0x1.4097f745ab636p-30, .Y = 0x1.d10fdb5823b29p-31}},
            .Second =
                Bezier::Cubic{.Start = {.X = 0x1.fd42966af636cp-32, .Y = -0x1.d907134486e58p-32},
                              .Out   = {.X = 0x1.2f5ba023b8fbep-31, .Y = -0x1.d8babb0cdd520p-35},
                              .In    = {.X = -0x1.5fc4857dd63a8p-33, .Y = 0x1.f8ee90e5161b0p-32},
                              .End   = {.X = -0x1.787ec27901648p-31, .Y = 0x1.15a26cc60aefcp-32}},
            .Within = 0x1.b7cdfd9d7bdbbp-63,
            .Points = {Point{.X = 0x1.243e735af7a06p-35, .Y = 0x1.8b9e86b8552e9p-33}}},
        SwiftOracle{
            .Id    = 1202,
            .First = Ray{.Start = {.X = -0x1.2aaada8739dcdp+25, .Y = 0x1.7b12f8385c896p+25},
                         .Way   = {.X = 0x1.36ddaa74e441ap+24, .Y = -0x1.82aec816c8173p+25}},
            .Second =
                Bezier::Cubic{.Start = {.X = -0x1.7caee5753aa7ep+24, .Y = -0x1.b83e817c7177dp+25},
                              .Out   = {.X = -0x1.48b2bb602b5edp+25, .Y = 0x1.6e9588942671ep+25},
                              .In    = {.X = -0x1.367b502086557p+23, .Y = -0x1.15667d3d2c6a2p+25},
                              .End   = {.X = -0x1.13811134fb4a5p+22, .Y = -0x1.44d31754c5ea3p+24}},
            .Within = 0x1.8fc504816f007p-7,
            .Points = {Point{.X = -0x1.51404da9aa9e2p+23, .Y = -0x1.3462a755c42b0p+24}}},
        SwiftOracle{.Id     = 1351,
                    .First  = Ray{.Start = {.X = -0x1.89ce180ff2d61p-1, .Y = -0x1.562ca752177c0p-3},
                                  .Way   = {.X = 0x1.2302506dbfe3ap-2, .Y = 0x1.457e76cc082c0p-3}},
                    .Second = Polyline{.Vertices = {{.Position = {.X = -0x1.95820f14422d8p-3,
                                                                  .Y = 0x1.53c9f96ba65e4p-1},
                                                     .Bulge    = 0x0.0p+0},
                                                    {.Position = {.X = -0x1.48256acee7350p-1,
                                                                  .Y = -0x1.c55e65bbd876ap-1},
                                                     .Bulge    = 0x0.0p+0}},
                                       .Closed   = false},
                    .Within = 0x1.b7cdfd9d7bdbbp-33,
                    .Points = {Point{.X = -0x1.7cf0f2a092328p-2, .Y = 0x1.c2bdbc41c382fp-5}}},
        SwiftOracle{
            .Id    = 1353,
            .First = Ray{.Start = {.X = -0x1.3ec988411e370p-3, .Y = -0x1.5e8be462d7984p+0},
                         .Way   = {.X = 0x1.792b1eccbdd78p-2, .Y = 0x1.1d63fe3c651aap+0}},
            .Second =
                Polyline{.Vertices =
                             {{.Position = {.X = 0x1.56680b809b044p-1, .Y = -0x1.7c480372d6f6ap-1},
                               .Bulge    = -0x1.fb97639b1d79cp-2},
                              {.Position = {.X = 0x1.c1616b4bfb648p-2, .Y = 0x1.4e2719eb7f6a2p-1},
                               .Bulge    = 0x1.1ffe8fba0053ap+0},
                              {.Position = {.X = -0x1.5ff8f817f8d4cp-1, .Y = 0x1.22f438e3a717ap-1},
                               .Bulge    = -0x1.3f62729758c16p-1},
                              {.Position = {.X = 0x1.b41402def7978p-1, .Y = -0x1.8d9ef21f8cf60p-1},
                               .Bulge    = 0x1.34c9f4620b528p-2},
                              {.Position = {.X = 0x1.2f66b23cf5fccp-1, .Y = -0x1.709c8eba623c0p-1},
                               .Bulge    = 0x1.551e6c5fe5798p+0},
                              {.Position = {.X = 0x1.9423ed27edc32p-1, .Y = -0x1.a75a16007fec8p-1},
                               .Bulge    = 0x0.0p+0}},
                         .Closed = true},
            .Within = 0x1.b7cdfd9d7bdbbp-33,
            .Points = {Point{.X = 0x1.d78695f3ed325p-3, .Y = -0x1.9c5bfb4bfbb74p-3},
                       Point{.X = 0x1.c087d4da3ccdbp-2, .Y = 0x1.b5cc02d9ecf65p-2}}},
        SwiftOracle{
            .Id    = 1421,
            .First = Ray{.Start = {.X = -0x1.3f0d7ed1ef04ep-1, .Y = -0x1.d8a25b5cce540p-4},
                         .Way   = {.X = 0x1.87b3dc1e3844cp-1, .Y = -0x1.79ffe9b4ead70p-3}},
            .Second =
                Polyline{.Vertices =
                             {{.Position = {.X = -0x1.5035bf91429e6p-1, .Y = -0x1.c3fe40b4f38d6p-3},
                               .Bulge    = 0x0.0p+0},
                              {.Position = {.X = -0x1.7f24b38b448a9p-2, .Y = -0x1.ea6611c1884f8p-4},
                               .Bulge    = 0x1.842414ac31c12p-3},
                              {.Position = {.X = 0x1.7a45a8eb8df22p-3, .Y = -0x1.040f1ffca6994p-2},
                               .Bulge    = 0x0.0p+0}},
                         .Closed = false},
            .Within = 0x1.b7cdfd9d7bdbbp-33,
            .Points = {Point{.X = -0x1.de77c81fce5b6p-2, .Y = -0x1.3957f47397846p-3},
                       Point{.X = -0x1.be5f538c14230p-4, .Y = -0x1.ea5d217b33c9cp-3},
                       Point{.X = -0x1.b81e260403b51p-4, .Y = -0x1.eb1e46b9eff85p-3}}},
        SwiftOracle{
            .Id    = 1443,
            .First = Ray{.Start = {.X = -0x1.7d78404514fd4p+26, .Y = 0x1.ee9da92d822d0p+25},
                         .Way   = {.X = 0x1.348477ed803c0p+0, .Y = -0x1.7e82ec04aa4b7p-1}},
            .Second =
                Polyline{
                    .Vertices =
                        {{.Position = {.X = -0x1.7d7840024693ep+26, .Y = 0x1.ee9da8f725fe6p+25},
                          .Bulge    = 0x1.55cb83ed7e972p+0},
                         {.Position = {.X = -0x1.7d783fd9c3df0p+26, .Y = 0x1.ee9da898a1524p+25},
                          .Bulge    = -0x1.1dfe227edeef0p+0},
                         {.Position = {.X = -0x1.7d783fd478bedp+26, .Y = 0x1.ee9da8c1ef576p+25},
                          .Bulge    = 0x0.0p+0},
                         {.Position = {.X = -0x1.7d783fc837f9bp+26, .Y = 0x1.ee9da93c4e439p+25},
                          .Bulge    = 0x0.0p+0}},
                    .Closed = false},
            .Within = 0x1.0000000000000p-21,
            .Points = {Point{.X = -0x1.7d78400cde55fp+26, .Y = 0x1.ee9da8e7d01bcp+25},
                       Point{.X = -0x1.7d783fe29f7cep+26, .Y = 0x1.ee9da8b36f752p+25}}},
        SwiftOracle{
            .Id    = 1450,
            .First = Ray{.Start = {.X = -0x1.11c4a2e5efa30p-32, .Y = 0x1.e203219af5850p-31},
                         .Way   = {.X = 0x1.f01f260de3812p-32, .Y = -0x1.af34bdb989f7ep-31}},
            .Second =
                Polyline{
                    .Vertices =
                        {{.Position = {.X = -0x1.08493cb9b3aa8p-32, .Y = 0x1.8c866d3fa0798p-31},
                          .Bulge    = -0x1.2687457dc2500p-2},
                         {.Position = {.X = 0x1.23e8cf0162c72p-31, .Y = 0x1.3f89618996860p-34},
                          .Bulge    = 0x1.59897cb24107cp+0},
                         {.Position = {.X = -0x1.31de2d4704afcp-31, .Y = -0x1.76a4b8fbf7b3cp-32},
                          .Bulge    = -0x1.ae76dbffd7b60p-4},
                         {.Position = {.X = 0x1.8648416cad934p-31, .Y = 0x1.a8318c01dfee0p-31},
                          .Bulge    = 0x1.80ff754712c00p-1}},
                    .Closed = true},
            .Within = 0x1.b7cdfd9d7bdbbp-63,
            .Points = {Point{.X = -0x1.cf2dc8c5ae3d9p-33, .Y = 0x1.bd5a3b35639d7p-31},
                       Point{.X = -0x1.490261a1f57d6p-33, .Y = 0x1.830ba4e0ec05dp-31},
                       Point{.X = -0x1.ebf201a88953fp-34, .Y = 0x1.5ef5add0755edp-31},
                       Point{.X = 0x1.2c8137f0b0db1p-34, .Y = 0x1.6589f3b4e27ebp-32}}},
        SwiftOracle{
            .Id    = 1462,
            .First = Ray{.Start = {.X = -0x1.cca019f9d4719p+20, .Y = -0x1.f893e42f7ad6cp+24},
                         .Way   = {.X = -0x1.17d92d5b2b0fap+21, .Y = 0x1.3dd92e66df721p+25}},
            .Second =
                Polyline{
                    .Vertices =
                        {{.Position = {.X = -0x1.40a60d27c9614p+23, .Y = 0x1.dcf208cd59146p+24},
                          .Bulge    = 0x0.0p+0},
                         {.Position = {.X = 0x1.5fa0fe691ff47p+24, .Y = -0x1.e0aacce71b677p+24},
                          .Bulge    = 0x0.0p+0},
                         {.Position = {.X = -0x1.cd39052f4a801p+25, .Y = 0x1.66b0a2ae777f2p+25},
                          .Bulge    = 0x0.0p+0},
                         {.Position = {.X = 0x1.26a2a32d10f25p+25, .Y = -0x1.6b8b085d49746p+25},
                          .Bulge    = 0x0.0p+0},
                         {.Position = {.X = 0x1.616f7e0076cedp+25, .Y = -0x1.d4ec190ab976cp+23},
                          .Bulge    = -0x1.5184b14707558p+0}},
                    .Closed = false},
            .Within = 0x1.8fc504816f007p-7,
            .Points = {Point{.X = -0x1.9154ec9b1eea0p+21, .Y = -0x1.d0614761479a6p+22},
                       Point{.X = -0x1.977b4a4155f04p+21, .Y = -0x1.98801d23902d5p+22},
                       Point{.X = -0x1.279eceb1bc98dp+22, .Y = 0x1.3b4c2d42625b6p+24}}},
        SwiftOracle{
            .Id    = 1560,
            .First = Xline{.Through = {.X = 0x1.d4e61f20ac2dcp-3, .Y = 0x1.e5981876dc92ep-2},
                           .Way     = {.X = 0x1.cad762bee1424p-1, .Y = 0x1.c659c0636d3bep-2}},
            .Second =
                Ellipse{.Map   = {.X      = {.X = -0x1.4ca93a31c88f4p-3, .Y = 0x1.0c1596d7f3604p-3},
                                  .Y      = {.X = -0x1.6e9dba39a0339p-2, .Y = -0x1.c6ed688ac1cdap-2},
                                  .Origin = {.X = 0x1.1ccfbcc0acb90p-2, .Y = 0x1.5a6c3e04f9fa6p-2}},
                        .Start = Degrees(-96.0),
                        .Sweep = Degrees(360.0)},
            .Within = 0x1.b7cdfd9d7bdbbp-33,
            .Points = {Point{.X = 0x1.bd9ef67c7e7edp-5, .Y = 0x1.8d186d45d98fap-2},
                       Point{.X = 0x1.57acfcf72533cp-1, .Y = SwiftBits(0x3fe62ea16a006b00)}}},
        SwiftOracle{
            .Id    = 1561,
            .First = Xline{.Through = {.X = 0x1.44d947c8a4698p-2, .Y = -0x1.e691c8ab2cb40p-4},
                           .Way     = {.X = 0x1.eca090d14ff78p-2, .Y = 0x1.c0dbad78ef166p-1}},
            .Second =
                Ellipse{.Map   = {.X      = {.X = -0x1.2f1b58d35c101p-1, .Y = 0x1.280b938437f1bp-4},
                                  .Y      = {.X = -0x1.1d9ebb2cf71e6p-4, .Y = -0x1.246ed77282189p-1},
                                  .Origin = {.X = 0x1.b882708495e7ap-2, .Y = 0x1.9e466b4b93358p-3}},
                        .Start = Degrees(-143.0),
                        .Sweep = Degrees(-320.0)},
            .Within = 0x1.b7cdfd9d7bdbbp-33,
            .Points = {Point{.X = 0x1.9ebae81176683p-3, .Y = -0x1.4fbbba413d90ap-2}}},
        SwiftOracle{
            .Id    = 1629,
            .First = Xline{.Through = {.X = 0x1.66253ed2f6112p-2, .Y = 0x1.7265ee1dfa2ecp-3},
                           .Way     = {.X = 0x1.fbcd0af30e575p-4, .Y = 0x1.fc0cd3b0a620ap-1}},
            .Second =
                Ellipse{.Map   = {.X      = {.X = -0x1.e1d5f6d6f8997p-3, .Y = 0x1.ddaf6489294eep-6},
                                  .Y      = {.X = 0x1.fc406b14d5777p-5, .Y = 0x1.00557dc1e3109p-1},
                                  .Origin = {.X = 0x1.43c6376c0d591p-1, .Y = 0x1.0fb52f1aec287p-1}},
                        .Start = Degrees(0.0),
                        .Sweep = Degrees(360.0)},
            .Within = 0x1.b7cdfd9d7bdbbp-33,
            .Points = {Point{.X = 0x1.96a1736c9e669p-2, .Y = 0x1.1ea2aa3f3577bp-1},
                       Point{.X = 0x1.96e73a6d1cbddp-2, .Y = 0x1.1fb9e9502a93bp-1}}},
        SwiftOracle{
            .Id     = 1651,
            .First  = Xline{.Through = {.X = -0x1.7d78401fe4307p+26, .Y = -0x1.7cc87487254b1p+25},
                            .Way     = {.X = -0x1.dbd357d74c742p-1, .Y = 0x1.7a11a84abc68bp-2}},
            .Second = Ellipse{.Map = {.X = {.X = 0x1.0f9bf5ba9d883p-1, .Y = 0x1.2812bc97136e6p-2},
                                      .Y = {.X = -0x1.4c6ca778941fap-4, .Y = 0x1.30f4f2675105ap-3},
                                      .Origin = {.X = -0x1.7d783febf9810p+26,
                                                 .Y = -0x1.7cc8747539a65p+25}},
                              .Start = Degrees(-179.0),
                              .Sweep = Degrees(360.0)},
            .Within = 0x1.0000000000000p-21,
            .Points = {Point{.X = -0x1.7d784003221d7p+26, .Y = -0x1.7cc8749dfedb5p+25},
                       Point{.X = -0x1.7d78400e3d6a9p+26, .Y = -0x1.7cc874952bb38p+25}}},
        SwiftOracle{
            .Id    = 1658,
            .First = Xline{.Through = {.X = 0x1.5d07a3a968864p-33, .Y = -0x1.a104c1f0cdd30p-32},
                           .Way     = {.X = 0x1.dd52cbee65f4fp-31, .Y = 0x1.7270808d5429dp-32}},
            .Second =
                Ellipse{
                    .Map   = {.X      = {.X = 0x1.29b029ecab412p-31, .Y = 0x1.1aeec01ec6e42p-32},
                              .Y      = {.X = -0x1.948c43465f8e7p-32, .Y = 0x1.a9a55a63fb951p-31},
                              .Origin = {.X = 0x1.99afbcf2df07cp-32, .Y = -0x1.e98f22007f638p-32}},
                    .Start = Degrees(-159.0),
                    .Sweep = Degrees(-54.0)},
            .Within = 0x1.b7cdfd9d7bdbbp-63,
            .Points = {Point{.X = -0x1.0220b257df683p-32, .Y = -0x1.24734266aa462p-31}}},
        SwiftOracle{
            .Id    = 1670,
            .First = Xline{.Through = {.X = 0x1.912fd2865ec07p+18, .Y = -0x1.0b98a081a52b7p+23},
                           .Way     = {.X = -0x1.ec2815b67f07bp+24, .Y = 0x1.8b039f89232f6p+25}},
            .Second =
                Ellipse{.Map = {.X      = {.X = 0x1.d7d32c5d6fa1ap+23, .Y = 0x1.1b74dd7d67d8bp+24},
                                .Y      = {.X = -0x1.87b35dbbb1cecp+23, .Y = 0x1.4600371f0dfabp+23},
                                .Origin = {.X = 0x1.8a93a583b55d2p+22, .Y = 0x1.b67fd009a5156p+22}},
                        .Start = Degrees(111.0),
                        .Sweep = Degrees(141.0)},
            .Within = 0x1.8fc504816f007p-7,
            .Points = {Point{.X = 0x1.aef2a6bba6942p+21, .Y = -0x1.a46a25476fb80p+23},
                       Point{.X = -0x1.585aa18669419p+23, .Y = 0x1.314c8350786f9p+23}}},
        SwiftOracle{.Id = 1683,
                    .First =
                        Xline{.Through = {.X = 0x1.6eeb7dd7d4960p-5, .Y = 0x1.f7d472a59a07ep-2},
                              .Way     = {.X = 0x1.95f99fbeed099p-3, .Y = 0x1.f5d6f1f914abfp-1}},
                    .Second =
                        Bezier::Quadratic{
                            .Start   = {.X = -0x1.a7bcabce141e4p-2, .Y = -0x1.d1247a34cca80p-3},
                            .Control = {.X = 0x1.f8fb3e771ff26p-1, .Y = -0x1.5ced93985f29ap-1},
                            .End     = {.X = -0x1.3e3159fea3954p-1, .Y = -0x1.8dc963f7cf630p-1}},
                    .Within = 0x1.b7cdfd9d7bdbbp-33,
                    .Points = {Point{.X = -0x1.a286e51183ac7p-3, .Y = -0x1.7ad4bba7d2673p-1},
                               Point{.X = -0x1.f22705c7f051dp-4, .Y = -0x1.52bc006542415p-2}}},
        SwiftOracle{.Id = 1684,
                    .First =
                        Xline{.Through = {.X = 0x1.428e16269c7cep-2, .Y = -0x1.2251e889f8138p-4},
                              .Way     = {.X = 0x1.aef91c2f21570p-1, .Y = 0x1.146aada988d10p-1}},
                    .Second =
                        Bezier::Quadratic{
                            .Start   = {.X = -0x1.075ecf902e990p-3, .Y = -0x1.517ecdbe0a4b0p-1},
                            .Control = {.X = 0x1.7aff86daf2c00p-6, .Y = -0x1.3235bc0648720p-3},
                            .End     = {.X = 0x1.f4e2fccd4c6a6p-1, .Y = 0x1.3d53cd26f447ap-1}},
                    .Within = 0x1.b7cdfd9d7bdbbp-33,
                    .Points = {Point{.X = 0x1.28371e65ddf72p-3, .Y = -0x1.70eef5a12d41ap-3}}},
        SwiftOracle{
            .Id    = 1753,
            .First = Xline{.Through = {.X = -0x1.281d45f36512cp-2, .Y = 0x1.6888a485c6b68p-3},
                           .Way     = {.X = -0x1.a2fb12ffb504dp-1, .Y = 0x1.2646fc87abda2p-1}},
            .Second =
                Bezier::Quadratic{
                    .Start   = {.X = -0x1.82f4a675eb863p-3, .Y = 0x1.af0e8bcf37ec5p-4},
                    .Control = {.X = -0x1.4feed40bf6abbp-1, .Y = SwiftBits(0x3fdbbe1a416d1133)},
                    .End     = {.X = -0x1.1f6c1127ed61cp+0, .Y = 0x1.8666a09e450e1p-1}},
            .Within = 0x1.b7cdfd9d7bdbbp-33,
            .Points = {Point{.X = -0x1.4fdcbd0150a21p-1, .Y = SwiftBits(0x3fdbc150b177e668)}}},
        SwiftOracle{.Id = 1775,
                    .First =
                        Xline{.Through = {.X = 0x1.7d783ffe8b3e6p+26, .Y = -0x1.d04d78e27080cp+22},
                              .Way     = {.X = -0x1.96dc86d8cc35dp-1, .Y = 0x1.36d162b4d8101p-1}},
                    .Second =
                        Bezier::Quadratic{
                            .Start   = {.X = 0x1.7d783fd05305ep+26, .Y = -0x1.d04d7a244d776p+22},
                            .Control = {.X = 0x1.7d784011b3c0fp+26, .Y = -0x1.d04d7814beb95p+22},
                            .End     = {.X = 0x1.7d784021a201fp+26, .Y = -0x1.d04d75dddc0a0p+22}},
                    .Within = 0x1.0000000000000p-21,
                    .Points = {Point{.X = 0x1.7d783ff93153ap+26, .Y = -0x1.d04d78a107f14p+22}}},
        SwiftOracle{.Id = 1797,
                    .First =
                        Xline{.Through = {.X = 0x1.f12007e1e58bap+23, .Y = 0x1.2d6332a7030e7p+22},
                              .Way     = {.X = 0x1.94304689ea8a8p+24, .Y = 0x1.a3394bc64b46cp+25}},
                    .Second =
                        Bezier::Quadratic{
                            .Start   = {.X = -0x1.34f537ec6463cp+25, .Y = 0x1.4b068cd3adbc4p+25},
                            .Control = {.X = -0x1.7471cefbcfc2dp+25, .Y = 0x1.5d5459acd2600p+24},
                            .End     = {.X = 0x1.c07f0aebae966p+25, .Y = 0x1.c0106ce62d6b9p+25}},
                    .Within = 0x1.8fc504816f007p-7,
                    .Points = {Point{.X = 0x1.2aaac3f604416p+25, .Y = 0x1.8f6b3565a7762p+25}}},
        SwiftOracle{
            .Id    = 1807,
            .First = Xline{.Through = {.X = 0x1.41bc966f023d8p-3, .Y = -0x1.8383566724468p-3},
                           .Way     = {.X = 0x1.2d669b1fe8fffp-1, .Y = 0x1.9de2cbac8c135p-1}},
            .Second =
                Bezier::Cubic{.Start = {.X = -0x1.13463353350e4p-1, .Y = -0x1.c3f429be2944cp-1},
                              .Out   = {.X = 0x1.8ea9bbf373d78p-1, .Y = 0x1.24dc300dc1b04p-2},
                              .In    = {.X = -0x1.11b0904abefd0p-2, .Y = -0x1.dbe51ad60ade0p-2},
                              .End   = {.X = 0x1.51257585eca16p-1, .Y = 0x1.bd5bcb1c62a82p-1}},
            .Within = 0x1.b7cdfd9d7bdbbp-33,
            .Points = {Point{.X = 0x1.250f9f8db95dcp-4, .Y = -0x1.3a0dc22fcd861p-2},
                       Point{.X = 0x1.6e299adccdcf2p-3, .Y = -0x1.4681d270377bbp-3}}},
        SwiftOracle{
            .Id    = 1808,
            .First = Xline{.Through = {.X = -0x1.046947fcc9048p-2, .Y = -0x1.0b8e13ca652fep-2},
                           .Way     = {.X = -0x1.0a8ec11e05a9cp-5, .Y = 0x1.ffba980d77f6bp-1}},
            .Second =
                Bezier::Cubic{.Start = {.X = -0x1.f3f21df33c69cp-1, .Y = 0x1.a55b8168e8fecp-1},
                              .Out   = {.X = -0x1.da5fcdc99ffd0p-3, .Y = -0x1.6d8a17c1123c8p-2},
                              .In    = {.X = -0x1.d93602ba671fep-1, .Y = 0x1.d8a2d493aaa98p-1},
                              .End   = {.X = 0x1.24937d9bcf076p-1, .Y = -0x1.05aaff06de7d0p-4}},
            .Within = 0x1.b7cdfd9d7bdbbp-33,
            .Points = {Point{.X = -0x1.186ceaa729006p-2, .Y = 0x1.5b35010673c5ep-2}}},
        SwiftOracle{
            .Id    = 1877,
            .First = Xline{.Through = {.X = -0x1.1b87ed41ebee0p-5, .Y = 0x1.557e638ea6850p-2},
                           .Way     = {.X = -0x1.e15d95d8f1e7ep-2, .Y = 0x1.c3e70505b2628p-1}},
            .Second =
                Bezier::Cubic{.Start = {.X = -0x1.0ca8a7c2e49c1p-5, .Y = 0x1.50ff7070cac6fp-2},
                              .Out   = {.X = -0x1.27b5eccb4fbebp-3, .Y = 0x1.14303fe6714c9p-1},
                              .In    = {.X = -0x1.064111bc9c966p-2, .Y = 0x1.7fd837ecb0809p-1},
                              .End   = {.X = -0x1.78c766fd3ab89p-2, .Y = 0x1.eb77a04b22ff5p-1}},
            .Within = 0x1.b7cdfd9d7bdbbp-33,
            .Points = {Point{.X = -0x1.9a2c2517190c0p-3, .Y = 0x1.4a0217ff9daa9p-1}}},
        SwiftOracle{
            .Id    = 1899,
            .First = Xline{.Through = {.X = 0x1.7d7840119241dp+26, .Y = -0x1.094b1053ad08ap+26},
                           .Way     = {.X = 0x1.eb094694b7833p-1, .Y = 0x1.220444aead522p-2}},
            .Second =
                Bezier::Cubic{.Start = {.X = 0x1.7d78403bde83ep+26, .Y = -0x1.094b102f6b02ap+26},
                              .Out   = {.X = 0x1.7d783fce99214p+26, .Y = -0x1.094b10619af78p+26},
                              .In    = {.X = 0x1.7d783fcdee747p+26, .Y = -0x1.094b107262f95p+26},
                              .End   = {.X = 0x1.7d7840439777bp+26, .Y = -0x1.094b105210a1bp+26}},
            .Within = 0x1.0000000000000p-21,
            .Points = {Point{.X = 0x1.7d783fea9e4dfp+26, .Y = -0x1.094b105f2ddc1p+26}}},
        SwiftOracle{
            .Id    = 1906,
            .First = Xline{.Through = {.X = -0x1.472932f732234p-32, .Y = -0x1.d93aac50f9348p-33},
                           .Way     = {.X = -0x1.36c3d69fa3fe3p-31, .Y = 0x1.96e6dfe6071b2p-31}},
            .Second =
                Bezier::Cubic{.Start = {.X = 0x1.0ea93d8a5e076p-31, .Y = 0x1.5227cd51a9e2cp-31},
                              .Out   = {.X = -0x1.f856d2ea9eb10p-34, .Y = -0x1.d5661f9159f20p-32},
                              .In    = {.X = 0x1.824a3d58c2adep-31, .Y = -0x1.e8a5bb45c1372p-31},
                              .End   = {.X = -0x1.1dc714563c3ccp-31, .Y = -0x1.a95882e044920p-35}},
            .Within = 0x1.b7cdfd9d7bdbbp-63,
            .Points = {Point{.X = -0x1.5ea882d0dd638p-32, .Y = -0x1.9bb24f43321acp-33}}},
        SwiftOracle{
            .Id    = 1918,
            .First = Xline{.Through = {.X = -0x1.818c096c514a3p+24, .Y = -0x1.dd6f527d02a66p+22},
                           .Way     = {.X = -0x1.a2964f0017408p+25, .Y = 0x1.96d1c57995464p+24}},
            .Second =
                Bezier::Cubic{.Start = {.X = 0x1.2e90098ee3130p+24, .Y = -0x1.8c6512cfa5460p+25},
                              .Out   = {.X = 0x1.40c83c01121b2p+25, .Y = -0x1.5f3335241fb88p+24},
                              .In    = {.X = -0x1.76a611e84e62dp+25, .Y = -0x1.8d69dd1734e7cp+23},
                              .End   = {.X = 0x1.bf194ae5992d0p+25, .Y = 0x1.325b9b0f89700p+24}},
            .Within = 0x1.8fc504816f007p-7,
            .Points = {Point{.X = 0x1.e37921d8ca44dp+23, .Y = -0x1.a82ee0b62c6c5p+24}}},
        SwiftOracle{
            .Id    = 2067,
            .First = Xline{.Through = {.X = -0x1.2d38d6619cc9ep-2, .Y = -0x1.9fc1584940b5ap-2},
                           .Way     = {.X = -0x1.cac2ec0434fc6p-1, .Y = 0x1.c6ac61109fc0ap-2}},
            .Second =
                Polyline{
                    .Vertices =
                        {{.Position = {.X = -0x1.21832a45971fep-1, .Y = -0x1.48be030cdf800p-10},
                          .Bulge    = 0x0.0p+0},
                         {.Position = {.X = 0x1.4417d502d5484p-2, .Y = -0x1.b3e893bbaafb8p-1},
                          .Bulge    = 0x0.0p+0},
                         {.Position = {.X = -0x1.4a36455175490p-4, .Y = -0x1.8ebf0853ba900p-3},
                          .Bulge    = 0x1.1b6e7fe897640p+0}},
                    .Closed = true},
            .Within = 0x1.b7cdfd9d7bdbbp-33,
            .Points = {Point{.X = 0x1.8bb7da2626480p-3, .Y = -0x1.4b893ba1ad3a0p-1},
                       Point{.X = 0x1.7f104ccd1a839p-7, .Y = -0x1.1d7a647ce526ep-1}}},
        SwiftOracle{
            .Id    = 2068,
            .First = Xline{.Through = {.X = 0x1.7337f4f7f8ac8p-2, .Y = -0x1.c78598fcff230p-3},
                           .Way     = {.X = -0x1.c9e6734bf16dep-3, .Y = 0x1.f309ddac276b6p-1}},
            .Second =
                Polyline{.Vertices =
                             {{.Position = {.X = 0x1.e2ad58fea2370p-4, .Y = 0x1.04e6d3a31f500p-3},
                               .Bulge    = 0x0.0p+0},
                              {.Position = {.X = 0x1.f3d2c0bc77558p-2, .Y = 0x1.842e6af86b0bap-1},
                               .Bulge    = 0x1.07cafc3445db4p+0},
                              {.Position = {.X = -0x1.d6296305d10c8p-3, .Y = -0x1.733d4fe0d4ba8p-1},
                               .Bulge    = 0x0.0p+0},
                              {.Position = {.X = 0x1.8fa95f63775c0p-5, .Y = 0x1.35ef9db1a8c68p-3},
                               .Bulge    = 0x0.0p+0}},
                         .Closed = false},
            .Within = 0x1.b7cdfd9d7bdbbp-33,
            .Points = {Point{.X = 0x1.e3794dc30f740p-3, .Y = 0x1.50b2683ce9a76p-2},
                       Point{.X = 0x1.dbd484d05355ap-4, .Y = 0x1.b3f7f1d1bfae7p-1}}},
        SwiftOracle{
            .Id    = 2137,
            .First = Xline{.Through = {.X = -0x1.a3b6030d6dbf4p-2, .Y = 0x1.6c81d7414a642p-2},
                           .Way     = {.X = 0x1.dff2c8d2be18bp-1, .Y = 0x1.649d6d58135dep-2}},
            .Second =
                Polyline{.Vertices =
                             {{.Position = {.X = -0x1.ee84499c6dfe1p-2, .Y = 0x1.ba0577bdb1260p-2},
                               .Bulge    = 0x0.0p+0},
                              {.Position = {.X = -0x1.02eb0b186f118p-2, .Y = 0x1.3f94f6c85ac82p-1},
                               .Bulge    = 0x1.353e37a3f0b61p-2},
                              {.Position = {.X = 0x1.ee0fe8c1b396cp-1, .Y = 0x1.134279f4da5b9p+0},
                               .Bulge    = 0x0.0p+0}},
                         .Closed = false},
            .Within = 0x1.b7cdfd9d7bdbbp-33,
            .Points = {Point{.X = 0x1.b258566bfd11fp-2, .Y = 0x1.54e68259d0791p-1},
                       Point{.X = 0x1.b498e280395a0p-2, .Y = 0x1.55519b58716e0p-1}}},
        SwiftOracle{
            .Id    = 2159,
            .First = Xline{.Through = {.X = 0x1.7d78402221d4dp+26, .Y = 0x1.84d8a1ac06a91p+22},
                           .Way     = {.X = 0x1.fa54c369d7b62p-5, .Y = 0x1.ff0565df83217p-1}},
            .Second =
                Polyline{.Vertices =
                             {{.Position = {.X = 0x1.7d78403e72e77p+26, .Y = 0x1.84d89ff748a85p+22},
                               .Bulge    = 0x0.0p+0},
                              {.Position = {.X = 0x1.7d783feb66790p+26, .Y = 0x1.84d8a5e34cbc1p+22},
                               .Bulge    = 0x1.71c2fa086e090p-1},
                              {.Position = {.X = 0x1.7d78404494c74p+26, .Y = 0x1.84d8a33e8a936p+22},
                               .Bulge    = -0x1.8f4ce654d82e0p-4},
                              {.Position = {.X = 0x1.7d78400b41c38p+26, .Y = 0x1.84d8a72f4b943p+22},
                               .Bulge    = 0x0.0p+0}},
                         .Closed = false},
            .Within = 0x1.0000000000000p-21,
            .Points = {Point{.X = 0x1.7d7840226c055p+26, .Y = 0x1.84d8a1f6e714fp+22},
                       Point{.X = 0x1.7d784022b3ff1p+26, .Y = 0x1.84d8a23f8b82cp+22},
                       Point{.X = 0x1.7d7840256df7fp+26, .Y = 0x1.84d8a4fffb6f0p+22}}},
        SwiftOracle{
            .Id    = 2167,
            .First = Xline{.Through = {.X = 0x1.cc3e63aacd6b0p-32, .Y = -0x1.a14794fc7e0c0p-37},
                           .Way     = {.X = -0x1.edf924bba9a1bp-31, .Y = 0x1.0d541904099acp-32}},
            .Second =
                Polyline{
                    .Vertices =
                        {{.Position = {.X = -0x1.f877c41a440eap-31, .Y = 0x1.2bf5785694c46p-31},
                          .Bulge    = -0x1.ba4d04d182128p-3},
                         {.Position = {.X = 0x1.ece259f06c8f0p-31, .Y = 0x1.ddda25a32a7b6p-31},
                          .Bulge    = 0x1.3bcb697b31db0p+0},
                         {.Position = {.X = -0x1.f57888199a7e4p-32, .Y = 0x1.9810764422b90p-34},
                          .Bulge    = 0x1.337324b87ca84p-1},
                         {.Position = {.X = 0x1.ef53a77e81bfap-31, .Y = 0x1.68e17760fbf64p-32},
                          .Bulge    = -0x1.140c0852b9ff4p-2},
                         {.Position = {.X = 0x1.2e5a1257b2410p-31, .Y = -0x1.0cd3a194762a8p-32},
                          .Bulge    = -0x1.6dbd13f8d6ce8p+0},
                         {.Position = {.X = 0x1.91189b2dec70cp-31, .Y = 0x1.8346ef2c5319ep-31},
                          .Bulge    = 0x1.cd65dfd2b9668p-1}},
                    .Closed = false},
            .Within = 0x1.b7cdfd9d7bdbbp-63,
            .Points = {Point{.X = 0x1.93bb87a0a0d31p-31, .Y = -0x1.aeca4db135b8bp-34},
                       Point{.X =  0x1.5020f46f1ad1ap-31, .Y = -0x1.1b5a1e97993acp-34},
                       Point{.X = -0x1.9beed8a91ea8bp-36, .Y =  0x1.ddca8a86c438ap-34},
                       Point{.X = -0x1.393038c918c97p-31, .Y = 0x1.1b304060ddb2bp-32}}},
        SwiftOracle{
            .Id    = 2179,
            .First = Xline{.Through = {.X = -0x1.a74f9065a9863p+24, .Y = -0x1.dc4cdec4fa65bp+23},
                           .Way     = {.X = 0x1.ac9f25fb9a5c0p+25, .Y = 0x1.6a9fec08aa8eep+24}},
            .Second =
                Polyline{
                    .Vertices =
                        {{.Position = {.X = -0x1.3bca7e2df0bd7p+25, .Y = 0x1.cc87e182368e6p+25},
                          .Bulge    = -0x1.fa2a8c54b8b05p-1},
                         {.Position = {.X = 0x1.e7e5abae45bb1p+22, .Y = -0x1.846a39ee3d547p+23},
                          .Bulge    = 0x1.2a06aa390d22ap+0},
                         {.Position = {.X = -0x1.e41ac2281e57ap+20, .Y = -0x1.2650abfc16cdep+25},
                          .Bulge    = 0x0.0p+0},
                         {.Position = {.X = -0x1.4f3989cae010ap+24, .Y = 0x1.478d6f6f89b63p+25},
                          .Bulge    = -0x1.97fcb587424fep-1},
                         {.Position = {.X = -0x1.484c50ae61fd8p+24, .Y = 0x1.1c96a77fe73a1p+25},
                          .Bulge    = 0x0.0p+0}},
                    .Closed = false},
            .Within = 0x1.8fc504816f007p-7,
            .Points = {Point{.X = -0x1.221d3dbf2acfcp+23, .Y = -0x1.e1c80613a5cb2p+22},
                       Point{.X = 0x1.646f4a617bea0p+24, .Y = 0x1.6ec49d9032bd9p+22}}},
    };
    int sure = 0;
    for (const SwiftOracle &oracle : corpus)
    {
      for (const bool swapped : {false, true})
      {
        const Curve        &first     = swapped ? oracle.Second : oracle.First;
        const Curve        &second    = swapped ? oracle.First : oracle.Second;
        const Tolerance     tolerance = {.Length = oracle.Within, .Radians = 1e-10};
        const Intersections meeting   = Intersect(first, second, tolerance);
        std::vector<Point>  found;
        for (const Hit &hit : meeting.Hits)
        {
          CHECK(hit.Kind == Contact::Cross);
          found.push_back(hit.At);
        }
        const auto order = [](Point lhs, Point rhs) -> bool
        {
          return lhs.X < rhs.X || (Exactly(lhs.X, rhs.X) && lhs.Y < rhs.Y);
        };
        std::vector<Point> expected = oracle.Points;
        std::ranges::sort(found, order);
        std::ranges::sort(expected, order);
        INFO("case " << oracle.Id);
        CHECK(found == expected);
        sure += SwiftCertain(first, second, tolerance) ? 1 : 0;
      }
    }
    MESSAGE("sure of " << sure << " of " << 2 * corpus.size());
    CHECK(sure >= static_cast<int>(corpus.size()));
  }

  TEST_CASE("swift: a margin or a value that is not a size is near")
  {
    const Range  unit     = {.Start = 0.0, .End = 1.0};
    const double nan      = std::numeric_limits<double>::quiet_NaN();
    const double infinity = std::numeric_limits<double>::infinity();
    CHECK(Swift::Within(unit, 0.5, nan)      == Swift::Side::Near);
    CHECK(Swift::Within(unit, 0.5, infinity) == Swift::Side::Near);
    CHECK(Swift::Within(unit, 0.5, -1.0)     == Swift::Side::Near);
    CHECK(Swift::Within(unit, nan, 0.1)      == Swift::Side::Near);
    CHECK(Swift::Within(unit, 0.5, 0.1)      == Swift::Side::In);
    CHECK(Swift::Within(unit, 1.5, 0.1)      == Swift::Side::Out);
    CHECK(!Swift::Held({1.0, 1e300}));
    CHECK(Swift::Held({1.0, -2.0}));
  }

  TEST_CASE("swift: numbers past the range, ellipses flat, tiny or swept oddly: not sure")
  {
    const Tolerance tolerance = {.Length = 2e-10, .Radians = 1e-10};
    const Curve     line      = Line{.Start = {.X = -2, .Y = 0.25}, .End = {.X = 2, .Y = 0.5}};
    const Ellipse   oval      = {.Map = Scaling(1, 0.5), .Start = {}, .Sweep = Turn};
    CHECK(SwiftCertain(line, oval, tolerance));
    CHECK(!SwiftCertain(Line{.Start = {.X = -2, .Y = 1e300}, .End = {.X = 2, .Y = 0.5}}, oval,
                        tolerance));
    CHECK(!SwiftCertain(Line{.Start = {.X = -2, .Y = 0.25}, .End = {.X = 1e300, .Y = 0.5}}, oval,
                        tolerance));
    CHECK(!SwiftCertain(Ray{.Start = {.X = -2, .Y = 0.25}, .Way = {.X = 1e-300, .Y = 0}}, oval,
                        tolerance));
    CHECK(!SwiftCertain(Xline{.Through = {.X = -2, .Y = 0.25}, .Way = {.X = 1e-300, .Y = 0}}, oval,
                        tolerance));
    Ellipse far    = oval;
    far.Map.Origin = {.X = 1e300, .Y = 0};
    CHECK(!SwiftCertain(line, far, tolerance));
    CHECK(!SwiftCertain(line, Ellipse{.Map = Scaling(0x1p-150), .Start = {}, .Sweep = Turn},
                        tolerance));
    CHECK(
        !SwiftCertain(line, Ellipse{.Map = Scaling(1, 0.5), .Start = {}, .Sweep = {}}, tolerance));
    CHECK(!SwiftCertain(line, Ellipse{.Map = Scaling(1, 0.5), .Start = {}, .Sweep = Degrees(720)},
                        tolerance));
    CHECK(!SwiftCertain(
        line,
        Ellipse{.Map   = {.X = {.X = 1, .Y = 0}, .Y = {.X = 1, .Y = 1e-13}, .Origin = {}},
                .Start = {},
                .Sweep = Turn},
        tolerance));
    CHECK(
        !SwiftCertain(Ray{.Start = {.X = -0x1p-98, .Y = 0x1p-100}, .Way = {.X = 0x1p-199, .Y = 0}},
                      Ellipse{.Map = Scaling(0x1p-99), .Start = {}, .Sweep = Turn}, tolerance));
  }

  TEST_CASE("swift: a crossing of an ellipse with a coordinate exactly zero: the exact way's")
  {
    const Tolerance tolerance = {.Length = 2e-10, .Radians = 1e-10};
    const Curve     level     = Line{.Start = {.X = -2, .Y = 0}, .End = {.X = 2, .Y = 0}};
    const Curve     upright   = Line{.Start = {.X = 0, .Y = -2}, .End = {.X = 0, .Y = 2}};
    const Ellipse   round     = {.Map = Scaling(1, 0.5), .Start = Degrees(90), .Sweep = Turn};
    CHECK(!SwiftCertain(level, round, tolerance));
    CHECK(!SwiftCertain(upright,
                        Ellipse{.Map   = {.X      = {.X = 1, .Y = 0},
                                          .Y      = {.X = 0, .Y = 0.5},
                                          .Origin = {.X = 0.25, .Y = 0.25}},
                                .Start = {},
                                .Sweep = Turn},
                        tolerance));
  }

  TEST_CASE("swift: polylines out of range, of no segment, or bulged oddly: not sure")
  {
    const Tolerance tolerance = {.Length = 2e-10, .Radians = 1e-10};
    const Curve     line      = Line{.Start = {.X = -2, .Y = 0.25}, .End = {.X = 2, .Y = 0.5}};
    const Polyline  lone      = {.Vertices = {{.Position = {.X = 0, .Y = 1}, .Bulge = 0}},
                                 .Closed   = false};
    CHECK(!SwiftCertain(line, lone, tolerance));
    const Polyline huge = {.Vertices = {{.Position = {.X = 0, .Y = 1}, .Bulge = 0},
                                        {.Position = {.X = 0, .Y = 1e300}, .Bulge = 0}},
                           .Closed   = false};
    CHECK(!SwiftCertain(line, huge, tolerance));
    const Polyline odd = {.Vertices = {{.Position = {.X = 0, .Y = 1},
                                        .Bulge    = std::numeric_limits<double>::quiet_NaN()},
                                       {.Position = {.X = 0, .Y = -1}, .Bulge = 0}},
                          .Closed   = false};
    CHECK(!SwiftCertain(line, odd, tolerance));
  }

  TEST_CASE("swift: a crossing of a polyline with a coordinate exactly zero: the exact way's")
  {
    const Tolerance tolerance = {.Length = 2e-10, .Radians = 1e-10};
    const Curve     level     = Line{.Start = {.X = -2, .Y = 0}, .End = {.X = 2, .Y = 0}};
    const Curve     upright   = Line{.Start = {.X = 0, .Y = -2}, .End = {.X = 0, .Y = 2}};
    const Polyline  post      = {.Vertices = {{.Position = {.X = 0.5, .Y = -1}, .Bulge = 0},
                                              {.Position = {.X = 0.5, .Y = 1}, .Bulge = 0}},
                                 .Closed   = false};
    CHECK(!SwiftCertain(level, post, tolerance));
    const Polyline arch = {.Vertices = {{.Position = {.X = -1, .Y = -1}, .Bulge = 0.5},
                                        {.Position = {.X = -1, .Y = 1}, .Bulge = 0}},
                           .Closed   = false};
    CHECK(!SwiftCertain(level, arch, tolerance));
    const Polyline aside = {.Vertices = {{.Position = {.X = -1, .Y = -0.5}, .Bulge = 0.5},
                                         {.Position = {.X = 1, .Y = -0.5}, .Bulge = 0}},
                            .Closed   = false};
    CHECK(!SwiftCertain(upright, aside, tolerance));
  }

  TEST_CASE("swift: a short segment crossing a line all but square may lie along: not sure")
  {
    /**
     *  The thorough way takes two runs to lie along each other where the stretch of one that the
     *  other's projection covers is longer than the tolerance and its ends within it of the other's
     *  line: so it may for a run crossing another all but square, both short.
     */
    const Tolerance tolerance = {.Length = 2e-10, .Radians = 1e-10};
    const Curve     level     = Line{.Start = {.X = -1, .Y = 0.125}, .End = {.X = 1, .Y = 0.125}};
    const Polyline  stub      = {
        .Vertices = {{.Position = {.X = 0.25, .Y = 0.125 - 3e-9}, .Bulge = 0},
                     {.Position = {.X = 0.25 + 2e-10, .Y = 0.125 + 3e-9}, .Bulge = 0}},
        .Closed   = false};
    CHECK(!SwiftCertain(level, stub, tolerance));
    const Curve    post = Line{.Start = {.X = 0.25, .Y = 0.125 - 5e-9},
                               .End   = {.X = 0.25 + 2e-10, .Y = 0.125 + 5e-9}};
    const Polyline flat = {.Vertices = {{.Position = {.X = -1, .Y = 0.125}, .Bulge = 0},
                                        {.Position = {.X = 1, .Y = 0.125}, .Bulge = 0}},
                           .Closed   = false};
    CHECK(!SwiftCertain(post, flat, tolerance));
  }

  TEST_CASE("swift: a crossing near a line's end or a polyline's vertex: not sure")
  {
    const Tolerance tolerance = {.Length = 2e-10, .Radians = 1e-10};

    /**
     *  A line ending a few tolerances past a crossing, a straight segment's or a bulge's.
     */
    const Curve stops = Line{.Start = {.X = -1, .Y = 0.125}, .End = {.X = 0.5 + 3e-9, .Y = 0.125}};
    const Polyline post = {.Vertices = {{.Position = {.X = 0.5, .Y = -1}, .Bulge = 0},
                                        {.Position = {.X = 0.5, .Y = 1}, .Bulge = 0}},
                           .Closed   = false};
    CHECK(!SwiftCertain(stops, post, tolerance));
    const Polyline bow = {.Vertices = {{.Position = {.X = 0, .Y = -0.875}, .Bulge = -0.5},
                                       {.Position = {.X = 0, .Y = 1.125}, .Bulge = 0}},
                          .Closed   = false};
    const Curve    shy =
        Line{.Start = {.X = -1, .Y = 0.125}, .End = {.X = At(bow, 0.5).X + 3e-9, .Y = 0.125}};
    CHECK(!SwiftCertain(shy, bow, tolerance));

    /**
     *  Through a bulge's vertex, near no segment's crossing, or a line's end above a segment.
     */
    const Polyline cap = {.Vertices = {{.Position = {.X = 0, .Y = 0}, .Bulge = 1},
                                       {.Position = {.X = 2, .Y = 0}, .Bulge = 0}},
                          .Closed   = false};
    CHECK(!SwiftCertain(Line{.Start = {.X = -0.5, .Y = -1}, .End = {.X = 0.5, .Y = 1}}, cap,
                        tolerance));
    const Curve    ground = Line{.Start = {.X = -1, .Y = 0.125}, .End = {.X = 1, .Y = 0.125}};
    const Polyline above  = {.Vertices = {{.Position = {.X = 0, .Y = 0.125 + 1e-9}, .Bulge = 0},
                                          {.Position = {.X = 1, .Y = 1}, .Bulge = 0}},
                             .Closed   = false};
    CHECK(!SwiftCertain(ground, above, tolerance));
    const Curve    from = Line{.Start = {.X = 0, .Y = 0.125}, .End = {.X = 1, .Y = 0.125}};
    const Polyline over = {.Vertices = {{.Position = {.X = -1, .Y = 0.125 + 1e-9}, .Bulge = 0},
                                        {.Position = {.X = 1, .Y = 0.125 + 2.5e-9}, .Bulge = 0}},
                           .Closed   = false};
    CHECK(!SwiftCertain(from, over, tolerance));
    const Curve    twice   = Polyline{.Vertices = {{.Position = {.X = 0.5, .Y = -1}, .Bulge = 0},
                                                   {.Position = {.X = 0.5, .Y = -1}, .Bulge = 0},
                                                   {.Position = {.X = 0.5, .Y = 1}, .Bulge = 0}},
                                      .Closed   = false};
    const SwiftMet stutter = SwiftBoth({ground, twice}, tolerance);
    CHECK(stutter.Sure);
    CHECK(stutter.Agrees);
  }

  TEST_CASE("swift: Bézier curves out of range, exactly zero, all but square, or closed")
  {
    const Tolerance tolerance = {.Length = 2e-10, .Radians = 1e-10};
    const Curve     line      = Line{.Start = {.X = -2, .Y = 0.25}, .End = {.X = 2, .Y = 0.5}};
    const Curve     bent      = Bezier::Quadratic{
        .Start = {.X = -1, .Y = 1}, .Control = {.X = 0, .Y = 1e300}, .End = {.X = 1, .Y = 1}};
    CHECK(!SwiftCertain(line, bent, tolerance));
    const Curve level = Line{.Start = {.X = -2, .Y = 0}, .End = {.X = 2, .Y = 0}};
    const Curve cup   = Bezier::Quadratic{
        .Start = {.X = -1, .Y = 1}, .Control = {.X = 0, .Y = -1}, .End = {.X = 1, .Y = 1}};
    CHECK(!SwiftCertain(level, cup, tolerance));

    /**
     *  A short curve all but square to a line, as a short segment may lie along it.
     */
    const Curve raised = Line{.Start = {.X = -1, .Y = 0.125}, .End = {.X = 1, .Y = 0.125}};
    const Curve tiny   = Bezier::Quadratic{.Start   = {.X = 0.25, .Y = 0.125 - 3e-9},
                                           .Control = {.X = 0.25 + 1e-10, .Y = 0.125},
                                           .End     = {.X = 0.25 + 2e-10, .Y = 0.125 + 3e-9}};
    CHECK(!SwiftCertain(raised, tiny, tolerance));

    /**
     *  A curve that closes on itself has no chord to lie along a line by.
     */
    const Curve    loop   = Bezier::Cubic{.Start = {.X = 0, .Y = 0},
                                          .Out   = {.X = 2, .Y = 2},
                                          .In    = {.X = -2, .Y = 2},
                                          .End   = {.X = 0, .Y = 0}};
    const Curve    across = Line{.Start = {.X = -2, .Y = 0.5}, .End = {.X = 2, .Y = 0.75}};
    const SwiftMet met    = SwiftBoth({across, loop}, tolerance);
    CHECK(met.Sure);
    CHECK(met.Agrees);
    CHECK(met.Swiftly.Hits.size() == 2);
  }

  TEST_CASE("swift: lines and ellipses as met thoroughly, to the bit")
  {
    const int rounds = Rounds(40000);
    const auto [sure, differ, better] =
        SwiftDrawn([](Draws &draws) -> Curve { return SwiftOval(draws); }, rounds);
    MESSAGE("sure of " << sure << " of " << rounds);
    CHECK(differ == 0);
    CHECK(better == 0);
    CHECK(sure > rounds / 3);
  }

  TEST_CASE("swift: lines and polylines as met thoroughly, to the bit")
  {
    const int rounds = Rounds(40000);
    const auto [sure, differ, better] =
        SwiftDrawn([](Draws &draws) -> Curve { return Bends(draws); }, rounds);
    MESSAGE("sure of " << sure << " of " << rounds);
    CHECK(differ == 0);
    CHECK(better == 0);
    CHECK(sure > rounds / 3);
  }

  TEST_CASE("swift: lines and Bézier curves as met thoroughly, to the bit")
  {
    const int rounds                  = Rounds(40000);
    const auto [sure, differ, better] = SwiftDrawn(SwiftBent, rounds);
    MESSAGE("sure of "
            << sure
            << " of "
            << rounds
            << ", proved better than the thorough way's in "
            << better);
    CHECK(differ == 0);
    CHECK(sure > rounds / 3);
  }

  TEST_CASE(
      "swift: a crossing a step of a double before an arc's start is left to the thorough way")
  {
    /**
     *  A line across a polyline's arc at its very start, the numbers near 5e8 and the arc a
     *  hundredth long, met with a tolerance of 1e-18, far under the numbers' grain as Intersect
     *  never asks: the crossing found lies a step of a double before the start, which is not on the
     *  arc by more than its place may be rounded, and nothing is told.
     */
    const Curve bow = Polyline{
        .Vertices = {{.Position = {.X = 0x1.3c8c25969ca7ep+23, .Y = 0x1.b8e97c0ad68b9p+28},
                      .Bulge    = -0x1.0ebf138507b2cp-1},
                     {.Position = {.X = 0x1.3c8c259447eb4p+23, .Y = 0x1.b8e97c0ad1814p+28}}}};
    const Curve     line  = Line{.Start = {.X = 0x1.3c8c25964db43p+23, .Y = 0x1.b8e97c0ae9b44p+28},
                                 .End   = {.X = 0x1.3c8c2596eb9b7p+23, .Y = 0x1.b8e97c0ac362ep+28}};
    const Tolerance under = {.Length = 0x1.4dd3500b2ba37p-60, .Radians = 0x1.12e0be826d695p-30};
    Intersections   met;
    CHECK_FALSE(Swift::Met({bow, line}, under, met));
    CHECK(met.Hits.empty());
  }
} // namespace Geoxy
