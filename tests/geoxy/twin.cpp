/***************************************************************************************************
 *  @file      tests/geoxy/twin.cpp
 *  @brief     Tests for the twin kernel: two quadratic or cubic Bézier curves met, each crossing
 *             the double nearest the exact one, and every answer the thorough way's to the bit or a
 *             better one proved exact.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "swift/twin.hpp"

#include "draws.hpp"
#include "expansion.hpp"
#include "real.hpp"
#include "swift/met.hpp"
#include "thorough.hpp"

#include "geoxy/at.hpp"
#include "geoxy/bezier/cubic.hpp"
#include "geoxy/bezier/quadratic.hpp"
#include "geoxy/curve.hpp"
#include "geoxy/derivative.hpp"
#include "geoxy/intersect.hpp"
#include "geoxy/point.hpp"
#include "geoxy/tolerance.hpp"
#include "geoxy/vector.hpp"

#include <doctest/doctest.h>

#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <utility>
#include <variant>
#include <vector>

namespace Geoxy
{
  namespace
  {
    /**
     *  A case of the curves audit's Bézier corpus: two curves, the tolerance they were met within,
     *  and where they cross as its exact rational reference finds them, each rounded once.
     */
    struct TwinOracle
    {
      int                Id     = 0;   /**< The case's number in the corpus.           */
      Curve              First  = {};  /**< The first curve.                           */
      Curve              Second = {};  /**< The second.                                */
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
    constexpr auto TwinBits(std::uint64_t bits) -> double
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
    auto TwinSame(const Intersections &one, const Intersections &other) -> bool
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
     *  Two curves met by the twin kernel and by the thorough way.
     */
    struct TwinMet
    {
      bool Sure   = false;         /**< Whether the kernel was sure. */
      bool Agrees = true;          /**< Whether it found what the thorough
                                        way finds, to the bit.              */
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
    auto TwinBoth(const Curves &curves, const Tolerance &tolerance) -> TwinMet
    {
      TwinMet                     met;
      const std::optional<double> length = Reachable(curves, tolerance);
      if (!length)
      {
        return met;
      }
      met.Sure = Swift::Met(curves, {.Length = *length, .Radians = tolerance.Radians}, met.Swiftly);
      if (met.Sure)
      {
        Thorough(curves, tolerance, met.Thorough);
        met.Agrees = TwinSame(met.Swiftly, met.Thorough);
      }
      return met;
    }

    /**
     *  Whether the kernels are sure of two curves, within the length Intersect would use.
     *
     *  @param  first     One curve.
     *  @param  second    The other.
     *  @param  tolerance The tolerance.
     *  @return           Whether they are.
     */
    auto TwinCertain(const Curve &first, const Curve &second, const Tolerance &tolerance) -> bool
    {
      Intersections meeting;
      return Swift::Met({first, second}, Usable({first, second}, tolerance), meeting);
    }

    /**
     *  A Bézier curve's control points.
     *
     *  @param  curve The curve.
     *  @return       Its control points; none for any other curve.
     */
    auto TwinPoles(const Curve &curve) -> std::vector<Point>
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
    auto TwinValued(std::vector<Expansion::Number> values, const Expansion::Number &along)
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
     *  A Bézier curve held exactly: its control points' x and y, its hodograph's, and a bound on
     *  its second derivative.
     */
    struct ExactBezier
    {
      std::array<std::vector<Expansion::Number>, 2> Poles = {};  /**< x and y.               */
      std::array<std::vector<Expansion::Number>, 2> Ways  = {};  /**< The hodograph's.       */
      double                                        Bend  = 0.0; /**< The most |B''| may be. */
    };

    /**
     *  A Bézier curve held exactly.
     *
     *  @param  poles Its control points.
     *  @return       It.
     */
    auto Exacted(const std::vector<Point> &poles) -> ExactBezier
    {
      ExactBezier exact;
      const auto  order = static_cast<double>(poles.size() - 1);
      for (std::size_t index = 0; index < poles.size(); ++index)
      {
        exact.Poles[0].push_back(Expansion::Of(poles[index].X));
        exact.Poles[1].push_back(Expansion::Of(poles[index].Y));
        if (index > 0)
        {
          for (std::size_t axis = 0; axis < 2; ++axis)
          {
            exact.Ways.at(axis).push_back(Expansion::Scaled(
                Expansion::Difference(exact.Poles.at(axis)[index], exact.Poles.at(axis)[index - 1]),
                order));
          }
        }
        if (index > 1)
        {
          const Vector second =
              (poles[index] - poles[index - 1]) - (poles[index - 1] - poles[index - 2]);
          exact.Bend = std::max(exact.Bend, 2.0 * order * (order - 1.0) * Length(second));
        }
      }
      return exact;
    }

    /**
     *  A place on each curve, exactly.
     */
    using Exacts = std::array<Expansion::Number, 2>;

    /**
     *  A box of places, exactly: the least place on each curve and the most.
     */
    struct Cell
    {
      Exacts Low  = {}; /**< The least. */
      Exacts High = {}; /**< The most.  */
    };

    /**
     *  The gap between two curves, A(s) - B(t), exactly.
     *
     *  @param  curves The curves.
     *  @param  place  The places s and t.
     *  @return        The gap's x and y.
     */
    auto Separation(const std::array<ExactBezier, 2> &curves, const Exacts &place) -> Exacts
    {
      return {Expansion::Difference(TwinValued(curves[0].Poles[0], place[0]),
                                    TwinValued(curves[1].Poles[0], place[1])),
              Expansion::Difference(TwinValued(curves[0].Poles[1], place[0]),
                                    TwinValued(curves[1].Poles[1], place[1]))};
    }

    /**
     *  The middle of two numbers, exactly.
     *
     *  @param  low  One.
     *  @param  high The other.
     *  @return      The middle.
     */
    auto Middle(const Expansion::Number &low, const Expansion::Number &high) -> Expansion::Number
    {
      return Expansion::Scaled(Expansion::Sum(low, high), 0.5);
    }

    /**
     *  Whether a box holds a root of the gap, by the theorem of Poincaré and Miranda on Y F, Y the
     *  inverse of the Jacobian at the box's middle, in doubles: along each edge where s is fixed,
     *  the first coordinate of Y F keeps one sign, opposite on the two, as its value at the edge's
     *  middle, worked out exactly, outweighs what its slope there and its bend may move it along
     *  the half edge; and likewise the second coordinate on the edges where t is fixed.
     *
     *  @param  curves The curves.
     *  @param  cell   The box.
     *  @return        Whether it is shown to hold one.
     */
    auto Holds(const std::array<ExactBezier, 2> &curves, const Cell &cell) -> bool
    {
      const Exacts middle = {Middle(cell.Low[0], cell.High[0]), Middle(cell.Low[1], cell.High[1])};
      const Exacts half   = {
          Expansion::Scaled(Expansion::Difference(cell.High[0], cell.Low[0]), 0.5),
          Expansion::Scaled(Expansion::Difference(cell.High[1], cell.Low[1]), 0.5)};
      std::array<std::array<double, 2>, 2> slopes = {};
      for (std::size_t which = 0; which < 2; ++which)
      {
        for (std::size_t axis = 0; axis < 2; ++axis)
        {
          slopes.at(which).at(axis) =
              Expansion::Approximate(TwinValued(curves.at(which).Ways.at(axis), middle.at(which)));
        }
      }
      const auto [ax, ay]                                = slopes[0];
      const auto [bx, by]                                = slopes[1];
      const double                               det     = (bx * ay) - (ax * by);
      const std::array<std::array<double, 2>, 2> inverse = {std::array{-by / det, bx / det},
                                                            std::array{-ay / det, ax / det}};
      for (std::size_t row = 0; row < 2; ++row)
      {
        const std::size_t  other = 1 - row;
        const auto         y     = inverse.at(row);
        const ExactBezier &moved = curves.at(other);
        std::array<int, 2> signs = {};
        for (std::size_t end = 0; end < 2; ++end)
        {
          Exacts place                = middle;
          place.at(row)               = end == 0 ? cell.Low.at(row) : cell.High.at(row);
          const Exacts            gap = Separation(curves, place);
          const Expansion::Number value =
              Expansion::Sum(Expansion::Scaled(gap[0], y[0]), Expansion::Scaled(gap[1], y[1]));
          const Expansion::Number slope =
              Expansion::Sum(Expansion::Scaled(TwinValued(moved.Ways[0], middle.at(other)), y[0]),
                             Expansion::Scaled(TwinValued(moved.Ways[1], middle.at(other)), y[1]));
          const double            bend  = (std::abs(y[0]) + std::abs(y[1])) * moved.Bend;
          const Expansion::Number reach = Expansion::Sum(
              Expansion::Product(Expansion::Scaled(slope, Expansion::Sign(slope)), half.at(other)),
              Expansion::Scaled(Expansion::Product(half.at(other), half.at(other)), bend));
          signs.at(end) = Expansion::Sign(value);
          if (Expansion::Sign(
                  Expansion::Difference(Expansion::Scaled(value, Expansion::Sign(value)), reach))
              <= 0)
          {
            return false;
          }
        }
        if (signs[0] != -signs[1])
        {
          return false;
        }
      }
      return true;
    }

    /**
     *  The double nearest a number held exactly.
     *
     *  @param  number The number.
     *  @return        The double.
     */
    auto Nearest(const Expansion::Number &number) -> double
    {
      return Expansion::Rounded(
          {.Top = number, .Bottom = {}, .Root = {}, .Below = Expansion::Of(1.0)}, 1);
    }

    /**
     *  Whether a hit of two Bézier curves is exact, checked in exact arithmetic: the box of the
     *  places halfway to the neighbours of the hit's places holds a root, so each place is the
     *  double nearest the root's; and, halved about the root until the first curve's point at both
     *  ends of the box's stretch rounds alike, that point is the hit's, the curve too short a way
     *  across the box to turn back.
     *
     *  @param  curves The two curves.
     *  @param  hit    The hit.
     *  @return        Whether it is exact.
     */
    auto TwinProved(const Curves &curves, const Hit &hit) -> bool
    {
      const std::array<ExactBezier, 2> exact = {Exacted(TwinPoles(curves[0])),
                                                Exacted(TwinPoles(curves[1]))};
      const auto halfway                     = [](double place, double toward) -> Expansion::Number
      {
        return Middle(Expansion::Of(place), Expansion::Of(std::nextafter(place, toward)));
      };
      const double along  = hit.First.Along;
      const double beside = hit.Second.Along;
      Cell         cell   = {.Low  = {halfway(along, -1.0), halfway(beside, -1.0)},
                             .High = {halfway(along, 2.0), halfway(beside, 2.0)}};
      if (!Holds(exact, cell))
      {
        return false;
      }
      for (int step = 0; step < 48; ++step)
      {
        const auto point = [&exact](const Expansion::Number &place) -> Point
        {
          return {.X = Nearest(TwinValued(exact[0].Poles[0], place)),
                  .Y = Nearest(TwinValued(exact[0].Poles[1], place))};
        };
        const Point low  = point(cell.Low[0]);
        const Point high = point(cell.High[0]);
        if (low == high)
        {
          return low == hit.At;
        }
        const Exacts middle = {Middle(cell.Low[0], cell.High[0]),
                               Middle(cell.Low[1], cell.High[1])};
        bool         found  = false;
        for (std::size_t quarter = 0; quarter < 4 && !found; ++quarter)
        {
          Cell part                                          = cell;
          ((quarter & 1U) == 0 ? part.High[0] : part.Low[0]) = middle[0];
          ((quarter & 2U) == 0 ? part.High[1] : part.Low[1]) = middle[1];
          if (Holds(exact, part))
          {
            cell  = part;
            found = true;
          }
        }
        if (!found)
        {
          return false;
        }
      }
      return false;
    }

    /**
     *  Whether the kernel's answer for two Bézier curves, where it differs from the thorough way's,
     *  is the better: as many hits, each a crossing, each proved exact.
     *
     *  @param  curves The two curves.
     *  @param  met    What each found.
     *  @return        Whether it is.
     */
    auto TwinBetter(const Curves &curves, const TwinMet &met) -> bool
    {
      return met.Swiftly.Hits.size() == met.Thorough.Hits.size()
          && std::ranges::all_of(met.Thorough.Hits,
                                 [](const Hit &hit) -> bool { return hit.Kind == Contact::Cross; })
          && std::ranges::all_of(met.Swiftly.Hits, [&curves](const Hit &hit) -> bool
                                 { return TwinProved(curves, hit); });
    }

    /**
     *  A double from -1 to 1.
     *
     *  @param  draws Where it comes from.
     *  @return       The double.
     */
    auto TwinUnit(Draws &draws) -> double
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
    auto TwinPick(Draws &draws, std::uint64_t count) -> std::uint64_t
    {
      return draws.Next() % count;
    }

    /**
     *  A quadratic or cubic Bézier curve at random in the square of side two about the origin.
     *
     *  @param  draws Where the numbers come from.
     *  @return       The curve.
     */
    auto TwinBent(Draws &draws) -> Curve
    {
      const auto point = [&draws]() -> Point
      {
        return {.X = TwinUnit(draws), .Y = TwinUnit(draws)};
      };
      if (TwinPick(draws, 3) == 0)
      {
        return Bezier::Quadratic{.Start = point(), .Control = point(), .End = point()};
      }
      return Bezier::Cubic{.Start = point(), .Out = point(), .In = point(), .End = point()};
    }

    /**
     *  A quadratic or cubic Bézier curve through a point along a way at its middle, bent to one
     *  side by an amount.
     *
     *  @param  draws Where the choices come from.
     *  @param  spot  The point and the way.
     *  @param  bend  How far its ends lie to the side, signed.
     *  @return       The curve.
     */
    auto TwinThrough(Draws &draws, const std::pair<Point, Vector> &spot, double bend) -> Curve
    {
      const auto [point, way] = spot;
      const Vector along      = way * ((0.3 + (0.5 * std::abs(TwinUnit(draws)))) / Length(way));
      const Vector side       = Perpendicular(along) * (bend / Length(along));
      if (TwinPick(draws, 2) == 0)
      {
        return Bezier::Quadratic{
            .Start = point - along + side, .Control = point - side, .End = point + along + side};
      }
      return Bezier::Cubic{.Start = point - along + side,
                           .Out   = point - (along / 3.0) - (side / 3.0),
                           .In    = point + (along / 3.0) - (side / 3.0),
                           .End   = point + along + side};
    }

    /**
     *  A curve moved and scaled: by a power of two and a little more, then far off.
     *
     *  @param  curve The curve.
     *  @param  scale The scale.
     *  @param  by    How far it moves.
     *  @return       The curve, moved.
     */
    auto TwinMoved(const Curve &curve, double scale, Vector by) -> Curve
    {
      const auto shifted = [scale, by](Point point) -> Point
      {
        return {.X = (point.X * scale) + by.X, .Y = (point.Y * scale) + by.Y};
      };
      if (const auto *quadratic = std::get_if<Bezier::Quadratic>(&curve))
      {
        return Bezier::Quadratic{.Start   = shifted(quadratic->Start),
                                 .Control = shifted(quadratic->Control),
                                 .End     = shifted(quadratic->End)};
      }
      const auto &cubic = std::get<Bezier::Cubic>(curve);
      return Bezier::Cubic{.Start = shifted(cubic.Start),
                           .Out   = shifted(cubic.Out),
                           .In    = shifted(cubic.In),
                           .End   = shifted(cubic.End)};
    }

    /**
     *  A Bézier curve near another at an edge of what the kernel decides: through it at a slight
     *  angle, along its tangent a few tolerances off, through one of its ends, or a copy of it
     *  moved a little.
     *
     *  @param  draws Where the numbers come from.
     *  @param  curve The other.
     *  @param  near  The tolerance's length.
     *  @return       The curve.
     */
    auto TwinEdged(Draws &draws, const Curve &curve, double near) -> Curve
    {
      const auto   kind  = TwinPick(draws, 4);
      const double end   = TwinPick(draws, 2) == 0 ? 0.0 : 1.0;
      const double along = kind == 2 ? end : 0.5 + (0.4 * TwinUnit(draws));
      const Point  spot  = At(curve, along);
      Vector       way   = Derivative(curve, along);
      if (!(Length(way) > 0.0))
      {
        way = {.X = 1.0, .Y = 0.0};
      }
      if (kind == 3)
      {
        return TwinMoved(curve, 1.0,
                         {.X = near * 16.0 * TwinUnit(draws), .Y = near * 16.0 * TwinUnit(draws)});
      }
      const double angle  = kind == 0
                              ? std::ldexp(TwinUnit(draws), -static_cast<int>(TwinPick(draws, 30)))
                              : 3.0 * TwinUnit(draws);
      const Vector turned = {.X = (way.X * std::cos(angle)) - (way.Y * std::sin(angle)),
                             .Y = (way.X * std::sin(angle)) + (way.Y * std::cos(angle))};
      const double off    = kind == 1 ? near * 16.0 * TwinUnit(draws) : 0.0;
      const Point  moved  = spot + (Perpendicular(way) * (off / Length(way)));
      const double bend =
          kind == 1 ? std::copysign(0.3 * std::abs(TwinUnit(draws)), off) : 0.3 * TwinUnit(draws);
      return TwinThrough(draws, {moved, kind == 1 ? way : turned}, bend);
    }

    /**
     *  A pair of Bézier curves and the tolerance asked for.
     */
    struct TwinPosed
    {
      Curve            First     = {}; /**< The first curve.         */
      Curve            Second    = {}; /**< The second.              */
      Geoxy::Tolerance Tolerance = {}; /**< The tolerance asked for. */
    };

    /**
     *  Two Bézier curves in one of the ways the kernel must answer as the thorough way does: at
     *  random or at an edge, far off or scaled, under a coarse angle at times, either way round.
     *
     *  @param  draws Where the numbers come from.
     *  @return       The pair.
     */
    auto TwinPosing(Draws &draws) -> TwinPosed
    {
      const Curve  one    = TwinBent(draws);
      const bool   coarse = TwinPick(draws, 4) == 0;
      const double near = coarse ? std::pow(10.0, -3.0 - (9.0 * std::abs(TwinUnit(draws)))) : 2e-10;
      const bool   wide = TwinPick(draws, 6) == 0;
      const double radians = wide ? std::pow(10.0, -6.0 * std::abs(TwinUnit(draws))) : 1e-10;
      const Curve  other = TwinPick(draws, 2) == 0 ? TwinBent(draws) : TwinEdged(draws, one, near);
      const auto   kind  = TwinPick(draws, 6);
      const double scale = std::ldexp(kind == 0 ? 1.37 : 1.0,
                                      kind == 0 ? static_cast<int>(TwinPick(draws, 80)) - 40 : 0);
      const Vector by    = {.X = kind == 1 ? std::ldexp(TwinUnit(draws), 26) : 0.0,
                            .Y = kind == 1 ? 1e8 : 0.0};
      const bool   swapped = TwinPick(draws, 2) == 0;
      const Curve  first   = TwinMoved(one, scale, by);
      const Curve  second  = TwinMoved(other, scale, by);
      return {.First     = swapped ? second : first,
              .Second    = swapped ? first : second,
              .Tolerance = {.Length = near * scale, .Radians = radians}};
    }
  } // namespace

  TEST_CASE("twin: the Bézier corpus, every crossing the double nearest the exact one")
  {
    const std::vector<TwinOracle> corpus = {
        TwinOracle{
            .Id = 186,
            .First =
                Bezier::Quadratic{.Start = {.X = 0x1.0c20ff46dc0e0p-2, .Y = -0x1.0d7a52966bc7ap-1},
                                  .Control = {.X = 0x1.2c5be8dd3c748p-1, .Y = 0x1.109395dd2d67ep-1},
                                  .End = {.X = -0x1.ee1009a596120p-2, .Y = 0x1.3fe7e0fb7bdd0p-4}},
            .Second =
                Bezier::Quadratic{
                    .Start   = {.X = 0x1.ae0ac4c8a60e8p-1, .Y = -0x1.60601ef8a886ap-1},
                    .Control = {.X = -0x1.18c6e4b9402f0p-4, .Y = 0x1.36c1e7ffec6d2p-1},
                    .End     = {.X = -0x1.295d9b291744cp-2, .Y = 0x1.30e2166fa3350p-3}},
            .Within = 0x1.b7cdfd9d7bdbbp-33,
            .Points = {Point{.X = -0x1.27ee674002e84p-2, .Y = 0x1.36be1bf8603f5p-3},
                       Point{.X = 0x1.1a6cbb415ceaep-5, .Y = 0x1.b74afe6fad9d1p-3},
                       Point{.X = 0x1.558b29da64105p-2, .Y = -0x1.67a224f75036bp-5}}},
        TwinOracle{.Id = 67,
                   .First =
                       Bezier::Quadratic{
                           .Start   = {.X = 0x1.e4091b40e5bc0p-1, .Y = -0x1.bfe67c729c728p-2},
                           .Control = {.X = 0x1.9eb4d66d31800p-9, .Y = -0x1.9346f986bfa9ep-1},
                           .End     = {.X = 0x1.a019be5329aa4p-2, .Y = 0x1.da57930d3464cp-1}},
                   .Second =
                       Bezier::Quadratic{
                           .Start   = {.X = 0x1.bb57891e3e782p-1, .Y = -0x1.f1163dbed9792p-1},
                           .Control = {.X = -0x1.810aa4e569f5ap-1, .Y = 0x1.f9d7df07c7302p-1},
                           .End     = {.X = 0x1.ba8bd7d6868dep-1, .Y = 0x1.ead45bada66d0p-2}},
                   .Within = 0x1.b7cdfd9d7bdbbp-33,
                   .Points = {Point{.X = 0x1.57af0bad67855p-2, .Y = 0x1.2b19734c198a7p-1},
                              Point{.X = 0x1.5d30834094276p-2, .Y = -0x1.19e6d46535544p-2},
                              Point{.X = 0x1.c42742d53a0fbp-2, .Y = -0x1.b05a39b709057p-2}}},
        TwinOracle{
            .Id = 1059,
            .First =
                Bezier::Quadratic{
                    .Start   = {.X = 0x1.a45a2a5b80bf6p-1, .Y = 0x1.4887b7bf724dcp-2},
                    .Control = {.X = -0x1.55312e0e06a30p-1, .Y = -0x1.4b9d84b6ccf4cp-1},
                    .End     = {.X = 0x1.b10c5b4d7a080p-7, .Y = 0x1.5e07766308f2cp-1}},
            .Second =
                Bezier::Cubic{.Start = {.X = -0x1.e9c755db70d48p-1, .Y = 0x1.1882c6dbe50a4p-2},
                              .Out   = {.X = 0x1.8f7c00eb49280p-4, .Y = 0x1.7f4b2649353f8p-2},
                              .In    = {.X = 0x1.9d2caeb6bf468p-1, .Y = -0x1.4df8c67a99386p-1},
                              .End   = {.X = -0x1.3c420b7f3603ep-1, .Y = 0x1.00fb80faacf7ep-1}},
            .Within = 0x1.b7cdfd9d7bdbbp-33,
            .Points = {Point{.X = -0x1.83a811093f038p-3, .Y = 0x1.5f5989c97e508p-3},
                       Point{.X = 0x1.3553701063749p-3, .Y = -0x1.6672d67bd54e9p-5},
                       Point{.X = 0x1.64f59d542c6fcp-3, .Y = -0x1.1f621f94d15c2p-5},
                       Point{.X = -0x1.86136c3dd8456p-3, .Y = 0x1.52c5fba9bc949p-3}}},
        TwinOracle{
            .Id = 1227,
            .First =
                Bezier::Quadratic{
                    .Start   = {.X = 0x1.6833ffdb6a938p-2, .Y = -0x1.1cc04ca5a9222p-1},
                    .Control = {.X = -0x1.4c8472a5c4aa0p-3, .Y = 0x1.a618e59398272p-1},
                    .End     = {.X = 0x1.7de7f1d05a5e4p-2, .Y = -0x1.e8c55aee24f00p-7}},
            .Second =
                Bezier::Cubic{.Start = {.X = 0x1.0a2d9ff312d06p-1, .Y = -0x1.d26fb1d9f6218p-2},
                              .Out   = {.X = -0x1.c2e974fad9c50p-2, .Y = -0x1.746647c3cea5ep-1},
                              .In    = {.X = 0x1.9222ef42db5e0p-5, .Y = -0x1.f8b79da57fee8p-1},
                              .End   = {.X = 0x1.6639796bb3078p-2, .Y = 0x1.11727a2256facp-2}},
            .Within = 0x1.b7cdfd9d7bdbbp-33,
            .Points = {Point{.X = 0x1.36fb2e5a01e3cp-2, .Y = 0x1.7460d416d28d8p-4},
                       Point{.X = 0x1.ba61ffc72b3e0p-3, .Y = -0x1.6a22cbc094417p-3},
                       Point{.X = 0x1.57cd0968c56efp-2, .Y = -0x1.06a6f368ca278p-1}}},
        TwinOracle{
            .Id    = 2302,
            .First = Bezier::Cubic{.Start = {.X = 0x1.7feebd53d71d8p-2, .Y = 0x1.b89ac84865160p-5},
                                   .Out = {.X = -0x1.06c17dfc43630p-2, .Y = -0x1.7084fc4e5d758p-2},
                                   .In  = {.X = -0x1.b316176ff3a48p-2, .Y = 0x1.3b59b34313f98p-3},
                                   .End = {.X = 0x1.d1dee243ebc54p-1, .Y = 0x1.3a4397958592ap-1}},
            .Second =
                Bezier::Cubic{.Start = {.X = -0x1.3903c83d41cecp-1, .Y = 0x1.d97ad3851bc56p-1},
                              .Out   = {.X = -0x1.7d3086408ed1cp-2, .Y = 0x1.9c0419da57100p-4},
                              .In    = {.X = 0x1.965868ae84258p-1, .Y = -0x1.b7102fd803f96p-1},
                              .End   = {.X = -0x1.2cfb006b0da94p-2, .Y = 0x1.80c9ca1e43650p-1}},
            .Within = 0x1.b7cdfd9d7bdbbp-33,
            .Points = {Point{.X = 0x1.cba7dc3032d82p-5, .Y = 0x1.80c70589cb051p-3},
                       Point{.X = 0x1.48534a101044ep-3, .Y = -0x1.eeadbb8ae84dap-5},
                       Point{.X = 0x1.0eb045274cd99p-4, .Y = -0x1.6f4b49c361dfap-4},
                       Point{.X = -0x1.44b92ec86a65dp-4, .Y = 0x1.4a0a2d3392c4ap-5}}},
        TwinOracle{
            .Id = 2335,
            .First =
                Bezier::Cubic{.Start = {.X = -0x1.8d3547b9f8ec0p-3, .Y = -0x1.dd3ed627f4ebcp-2},
                              .Out   = {.X = -0x1.3e15a38fc4b24p-1, .Y = -0x1.4f4114a3075bap-1},
                              .In    = {.X = 0x1.771a8cd8b7f0cp-1, .Y = -0x1.ad701b824bfb8p-2},
                              .End   = {.X = -0x1.f919146337e2cp-1, .Y = -0x1.970bde5884fecp-2}},
            .Second =
                Bezier::Cubic{.Start = {.X = -0x1.5e24f6ed389c8p-3, .Y = -0x1.8e9f3585b9c1cp-1},
                              .Out   = {.X = -0x1.bf631bc2ac698p-1, .Y = -0x1.9061b8fc01c00p-3},
                              .In    = {.X = 0x1.c8b0964df5934p-1, .Y = 0x1.5fd8a00aa94c0p-3},
                              .End   = {.X = -0x1.53e10f87d4740p-2, .Y = -0x1.4dbb8235ee318p-1}},
            .Within = 0x1.b7cdfd9d7bdbbp-33,
            .Points = {Point{.X = -0x1.46baf6f69a556p-3, .Y = -0x1.0ee8199431252p-1},
                       Point{.X = -0x1.813b2c79cb323p-4, .Y = -0x1.e90374c837e4ep-2},
                       Point{.X = -0x1.53b8b640c77e1p-2, .Y = -0x1.b4417ba3da7a7p-2}}},
        TwinOracle{
            .Id    = 3078,
            .First = Bezier::Cubic{.Start = {.X = -0x1.b3889061a9e0ep-1, .Y = 0x1.937a66453a5e4p-1},
                                   .Out = {.X = -0x1.6506ae8f7f600p-5, .Y = -0x1.6886e8b5fc058p-1},
                                   .In  = {.X = -0x1.1441eecc1e150p+1, .Y = 0x1.960356bacb5a2p-1},
                                   .End = {.X = 0x1.1e46fa78abb32p-1, .Y = -0x1.ae31bac27f85ep-3}},
            .Second =
                Bezier::Cubic{.Start = {.X = -0x1.d1a56e016e940p-1, .Y = -0x1.7407362398b00p-6},
                              .Out   = {.X = -0x1.46a218075e400p-7, .Y = 0x1.b5f46c0f188dcp-1},
                              .In    = {.X = 0x1.095030fd7ac00p-8, .Y = 0x1.e6bddd3e64192p-1},
                              .End   = {.X = -0x1.e7ab75ecabd56p-1, .Y = -0x1.9b1c249b5a384p-1}},
            .Within = 0x1.b7cdfd9d7bdbbp-33,
            .Points = {Point{.X = -0x1.cf76271c206f9p-2, .Y = 0x1.d019f5e7ab2cbp-4},
                       Point{.X = -0x1.4eb12eebaeea1p-1, .Y = 0x1.d4131b70c824fp-3},
                       Point{.X = -0x1.712dbac17cc33p-1, .Y = 0x1.4c2bd0d3b5c41p-3},
                       Point{.X = -0x1.9aa15b0b5becap-1, .Y = 0x1.5215abf5cdce9p-4}}},
        TwinOracle{
            .Id    = 3065,
            .First = Bezier::Cubic{.Start = {.X = 0x1.05649329a9b3ap-1, .Y = 0x1.5c8240978b4e0p-2},
                                   .Out = {.X = -0x1.6b867c6a3f02bp+0, .Y = -0x1.1a0096c4ecdb7p+0},
                                   .In  = {.X = 0x1.ef9fe8dc54fd2p-1, .Y = 0x1.8bfe6a9b89c50p-4},
                                   .End = {.X = -0x1.14dcac50c6d1bp+0, .Y = -0x1.d57f0de4d89c0p-2}},
            .Second =
                Bezier::Cubic{.Start = {.X = -0x1.435557d307fb0p-2, .Y = -0x1.e8e3d4eb8b79cp-1},
                              .Out   = {.X = 0x1.902ec4b673770p-3, .Y = 0x1.378f3bdc8d6c2p-1},
                              .In    = {.X = -0x1.fa50c67ef99e0p-1, .Y = 0x1.74d26b94a175ep-1},
                              .End   = {.X = -0x1.d699ad092d0f8p-3, .Y = 0x1.fe243878a4b6cp-1}},
            .Within = 0x1.b7cdfd9d7bdbbp-33,
            .Points = {Point{.X = -0x1.7704beb3dfc83p-3, .Y = -0x1.d01e8c94b4ac1p-3},
                       Point{.X = -0x1.7e52188a743aap-3, .Y = -0x1.4d402c26a6330p-2},
                       Point{.X = -0x1.814fdf457b9e0p-3, .Y = -0x1.63129cf591682p-2}}},
        TwinOracle{
            .Id = 2029,
            .First =
                Bezier::Quadratic{
                    .Start   = {.X = -0x1.846930b3bb084p-2, .Y = -0x1.b14703f620b80p-5},
                    .Control = {.X = 0x1.e49e709d6a0b8p-3, .Y = -0x1.e3cf3e4984630p-1},
                    .End     = {.X = 0x1.2098d09a8a300p-3, .Y = 0x1.5c196dec6d7b8p-1}},
            .Second =
                Bezier::Cubic{.Start = {.X = -0x1.e9334987e37e0p-1, .Y = -0x1.d6012b417a75ap-1},
                              .Out   = {.X = 0x1.87cffec940ba1p+0, .Y = 0x1.96768967e7688p+0},
                              .In    = {.X = -0x1.7152c855f9712p-1, .Y = -0x1.000bf6ccc9010p+0},
                              .End   = {.X = 0x1.176e8fbdc3c58p-1, .Y = 0x1.9d577faf70c1ap-1}},
            .Within = 0x1.b7cdfd9d7bdbbp-33,
            .Points = {Point{.X = 0x1.3907285950ff1p-3, .Y = 0x1.957e29e7bf41ap-3},
                       Point{.X = 0x1.33575d9f3c6e5p-3, .Y = 0x1.d90cc97bbbb86p-4},
                       Point{.X = 0x1.36ad6461d3d26p-3, .Y = 0x1.400bf72866fe3p-3},
                       Point{.X = -0x1.0157f883cfc88p-2, .Y = -0x1.c2720c42d5fdfp-3}}},
        TwinOracle{
            .Id = 2771,
            .First =
                Bezier::Cubic{.Start = {.X = 0x1.ffffeda19e8cfp+19, .Y = 0x1.aa353e20bd19dp+19},
                              .Out   = {.X = 0x1.0000117901af4p+20, .Y = 0x1.aa351e0833bd7p+19},
                              .In    = {.X = 0x1.000004428ad4cp+20, .Y = 0x1.aa355105f424dp+19},
                              .End   = {.X = 0x1.0000098c59ef4p+20, .Y = 0x1.aa352dae01b95p+19}},
            .Second =
                Bezier::Cubic{.Start = {.X = 0x1.00000adbea351p+20, .Y = 0x1.aa353e010f094p+19},
                              .Out   = {.X = 0x1.00000107752e3p+20, .Y = 0x1.aa3535baed0aap+19},
                              .In    = {.X = 0x1.00000e1ea8e25p+20, .Y = 0x1.aa351fbf6363cp+19},
                              .End   = {.X = 0x1.fffff5b3099f2p+19, .Y = 0x1.aa354d4e7598fp+19}},
            .Within = 0x1.0000000000000p-27,
            .Points = {Point{.X = 0x1.000005759f467p+20, .Y = 0x1.aa35330e527f8p+19},
                       Point{.X = 0x1.00000694f67cfp+20, .Y = 0x1.aa3533a215f0ep+19},
                       Point{.X = 0x1.0000081fcb1cdp+20, .Y = 0x1.aa353a62bf5a2p+19},
                       Point{.X = 0x1.0000083e6def4p+20, .Y = 0x1.aa353aa1cc169p+19}}},
        TwinOracle{
            .Id = 2807,
            .First =
                Bezier::Cubic{.Start = {.X = -0x1.7d783fcfd4c5cp+26, .Y = -0x1.039c1aac197dap+26},
                              .Out   = {.X = -0x1.7d783fda04b15p+26, .Y = -0x1.039c1aedd987bp+26},
                              .In    = {.X = -0x1.7d784006eed8bp+26, .Y = -0x1.039c1ab5bfda0p+26},
                              .End   = {.X = -0x1.7d783ff3bad40p+26, .Y = -0x1.039c1ad03e096p+26}},
            .Second =
                Bezier::Cubic{.Start = {.X = -0x1.7d783fca506bdp+26, .Y = -0x1.039c1ac8ed902p+26},
                              .Out   = {.X = -0x1.7d78403f9e1fap+26, .Y = -0x1.039c1abe9f96dp+26},
                              .In    = {.X = -0x1.7d7840262098ap+26, .Y = -0x1.039c1b12b7b06p+26},
                              .End   = {.X = -0x1.7d78403166d19p+26, .Y = -0x1.039c1ade9a9b2p+26}},
            .Within = 0x1.0000000000000p-21,
            .Points = {Point{.X = -0x1.7d783ff7e359bp+26, .Y = -0x1.039c1aca32397p+26},
                       Point{.X = -0x1.7d783ff527e72p+26, .Y = -0x1.039c1ac9bd2f4p+26},
                       Point{.X = -0x1.7d783fdaa031bp+26, .Y = -0x1.039c1ac8109bbp+26}}},
        TwinOracle{
            .Id = 1780,
            .First =
                Bezier::Quadratic{
                    .Start   = {.X = 0x1.7d783fe4e169cp+26, .Y = -0x1.705f8e000ff99p+26},
                    .Control = {.X = 0x1.7d784008a2b3bp+26, .Y = -0x1.705f8ddd2ed9ep+26},
                    .End     = {.X = 0x1.7d783fd83e832p+26, .Y = -0x1.705f8e0c94525p+26}},
            .Second =
                Bezier::Cubic{.Start = {.X = 0x1.7d784000e5a02p+26, .Y = -0x1.705f8dcfe2e56p+26},
                              .Out   = {.X = 0x1.7d784033a2dc5p+26, .Y = -0x1.705f8df697444p+26},
                              .In    = {.X = 0x1.7d783fe17e34bp+26, .Y = -0x1.705f8dda09a34p+26},
                              .End   = {.X = 0x1.7d783fd75fe68p+26, .Y = -0x1.705f8e1bd3aeep+26}},
            .Within = 0x1.0000000000000p-21,
            .Points = {Point{.X = 0x1.7d783fdebff6ap+26, .Y = -0x1.705f8e0634ee0p+26},
                       Point{.X = 0x1.7d783ff18d06fp+26, .Y = -0x1.705f8df3b6e1fp+26},
                       Point{.X = 0x1.7d783ff1b2034p+26, .Y = -0x1.705f8df3a0d55p+26}}},
        TwinOracle{.Id = 746,
                   .First =
                       Bezier::Quadratic{
                           .Start   = {.X = -0x1.7d783fc7f655ep+26, .Y = 0x1.3475db6f8338ap+25},
                           .Control = {.X = -0x1.7d78402fca1b9p+26, .Y = 0x1.3475dbf350927p+25},
                           .End     = {.X = -0x1.7d784002acc18p+26, .Y = 0x1.3475dc07b8f36p+25}},
                   .Second =
                       Bezier::Quadratic{
                           .Start   = {.X = -0x1.7d783ff612d7ap+26, .Y = 0x1.3475dbfb72a00p+25},
                           .Control = {.X = -0x1.7d7840176781ep+26, .Y = 0x1.3475db4e9e24dp+25},
                           .End     = {.X = -0x1.7d78400968fc8p+26, .Y = 0x1.3475dc06c78d2p+25}},
                   .Within = 0x1.0000000000000p-21,
                   .Points = {Point{.X = -0x1.7d784009b10a3p+26, .Y = 0x1.3475dc030c36cp+25},
                              Point{.X = -0x1.7d78400c49a64p+26, .Y = 0x1.3475dbdc3e9e6p+25},
                              Point{.X = -0x1.7d784001a11d1p+26, .Y = 0x1.3475dbc419f5dp+25}}},
        TwinOracle{
            .Id = 2867,
            .First =
                Bezier::Cubic{.Start = {.X = -0x1.6d59fc680d86cp+29, .Y = -0x1.00a5e3aff9ef8p+27},
                              .Out   = {.X = 0x1.828341be61d26p+29, .Y = 0x1.691232433781cp+28},
                              .In    = {.X = 0x1.33b2f6401a5f2p+29, .Y = -0x1.213cdd0ca5af8p+27},
                              .End   = {.X = -0x1.fc76605c94f7ep+29, .Y = 0x1.4db6887924620p+29}},
            .Second =
                Bezier::Cubic{.Start = {.X = 0x1.e58d5919d4442p+29, .Y = -0x1.e45b956e40ce0p+29},
                              .Out   = {.X = 0x1.5dd7a71b4ca58p+28, .Y = 0x1.ed44259512fa0p+29},
                              .In    = {.X = -0x1.f9ee76a62130ap+29, .Y = 0x1.99bdf2cf3b9a0p+25},
                              .End   = {.X = -0x1.815f1a7e66a3cp+29, .Y = 0x1.2ff9567be6aecp+29}},
            .Within = 0x1.b7cdfd9d7bdbbp-3,
            .Points = {Point{.X = -0x1.8aa79b5bcf6fbp+29, .Y = 0x1.170777b288eedp+29},
                       Point{.X = -0x1.49bfc473cd5c4p+28, .Y = 0x1.6a85d4ec7e51cp+28},
                       Point{.X =  0x1.6fb242085a44dp+27, .Y = 0x1.79cf6d30fbb10p+27},
                       Point{.X = 0x1.168800e8aaa11p+28, .Y = 0x1.fa41dca4b4a77p+26}}},
        TwinOracle{
            .Id = 2862,
            .First =
                Bezier::Cubic{.Start = {.X = 0x1.af2753825ce00p-31, .Y = 0x1.c8286bbf163a0p-34},
                              .Out   = {.X = -0x1.f9d20c9df1c90p-31, .Y = 0x1.b6a7b6ce4a448p-32},
                              .In    = {.X = -0x1.8eff2efbcf684p-32, .Y = -0x1.5e16712b24466p-31},
                              .End   = {.X = -0x1.272ed3eb1b140p-33, .Y = 0x1.aaef08e5893d0p-31}},
            .Second =
                Bezier::Cubic{.Start = {.X = 0x1.080a4e4e4bb32p-31, .Y = 0x1.3794d9b87f618p-31},
                              .Out   = {.X = 0x1.6188f8b40f582p-31, .Y = 0x1.275b3422a74c6p-31},
                              .In    = {.X = -0x1.c9b37999885b2p-31, .Y = -0x1.bbfd2433162d0p-31},
                              .End   = {.X = -0x1.0bd34a226f838p-32, .Y = 0x1.d961aab65ebd8p-32}},
            .Within = 0x1.b7cdfd9d7bdbbp-63,
            .Points = {Point{.X = -0x1.a394e28d268c5p-32, .Y = 0x1.3d9a8564d86d7p-35},
                       Point{.X = -0x1.a3ea00665aed0p-32, .Y = 0x1.21de940787413p-35},
                       Point{.X = 0x1.2e9ba2884da20p-33, .Y = 0x1.5d88243cb58ffp-33}}},
        TwinOracle{.Id = 872,
                   .First =
                       Bezier::Quadratic{
                           .Start   = {.X = 0x1.0e2c0b264cf68p-29, .Y = 0x1.419006d221410p-30},
                           .Control = {.X = 0x1.845e7263b2f5ap-29, .Y = -0x1.976fa67de9c4bp-29},
                           .End     = {.X = 0x1.fc8fb36ab9c33p-35, .Y = 0x1.19c2b945db461p-33}},
                   .Second =
                       Bezier::Quadratic{
                           .Start   = {.X = 0x1.4c2f6645d59f4p-31, .Y = -0x1.9c37190f48440p-38},
                           .Control = {.X = -0x1.16965671b40c7p-29, .Y = -0x1.51af606cedbaap-32},
                           .End     = {.X = 0x1.e394f6973dcc9p-29, .Y = -0x1.bab6328f1b6dbp-31}},
                   .Within = 0x1.b4d1a36891c79p-61,
                   .Points = {Point{.X = 0x1.2803201624e9ep-29, .Y = -0x1.71c6c8949c419p-31},
                              Point{.X = 0x1.5f907aec85404p-31, .Y = -0x1.085df5904bfc5p-31},
                              Point{.X = 0x1.f13807b8ccedcp-33, .Y = -0x1.f40855da3dd90p-35}}},
        TwinOracle{
            .Id = 1882,
            .First =
                Bezier::Quadratic{
                    .Start   = {.X = 0x1.79951c9149f10p+19, .Y = -0x1.4bfad9ba26f8ap+25},
                    .Control = {.X = -0x1.482a710b7ffc3p+25, .Y = 0x1.1f0ea96e8a60cp+25},
                    .End     = {.X = 0x1.44b1757d330e5p+25, .Y = 0x1.1a49b87e1e738p+24}},
            .Second =
                Bezier::Cubic{.Start = {.X = 0x1.3e397af59a563p+25, .Y = 0x1.c137050c170d8p+25},
                              .Out   = {.X = 0x1.09117b099d881p+25, .Y = 0x1.012eb4387a11ap+22},
                              .In    = {.X = -0x1.6df026c0af0fap+24, .Y = 0x1.b664929e3abd3p+24},
                              .End   = {.X = -0x1.736b039e3c54dp+21, .Y = 0x1.8d38705595811p+23}},
            .Within = 0x1.8fc504816f007p-7,
            .Points = {Point{.X = -0x1.a9bc30188bedep+22, .Y = 0x1.f824aa1c721bep+23},
                       Point{.X = -0x1.259477bb12535p+19, .Y = 0x1.30491ef7236e0p+24},
                       Point{.X = 0x1.7199f19ed5221p+23, .Y = 0x1.5146436633973p+24}}},
        TwinOracle{
            .Id    = 2635,
            .First = Bezier::Cubic{.Start = {.X = 0x1.b32f518600d3cp-2, .Y = -0x1.d1d9bc55b350cp-2},
                                   .Out = {.X = -0x1.507da7efb64aep-1, .Y = -0x1.ea3ffc7d7b408p-1},
                                   .In  = {.X = 0x1.3341acf7d8d3cp-2, .Y = -0x1.cea9f79d3568cp-1},
                                   .End = {.X = 0x1.7c392232f5f92p-1, .Y = 0x1.a5e86b4d31b84p-1}},
            .Second =
                Bezier::Cubic{
                    .Start = {.X = -0x1.de7d25f99bf69p-3, .Y = -0x1.eb6a8f60fe75bp-1},
                    .Out   = {.X = -0x1.52dd6a66a289fp-7, .Y = -0x1.5fe1399eadcfep-1},
                    .In    = {.X = 0x1.b375a8945d8f1p-3, .Y = -0x1.a86990ec669e0p-2},
                    .End   = {.X = TwinBits(0x3fdbdb6abdb5d983), .Y = -0x1.2194ef9e3bcd4p-3}},
            .Within = 0x1.b7cdfd9d7bdbbp-33,
            .Points = {Point{.X = 0x1.a88f58400e45dp-4, .Y = -0x1.1862a0e31c32ep-1},
                       Point{.X = 0x1.9df1e9e1bea75p-4, .Y = -0x1.1a023a30e60d1p-1},
                       Point{.X = -0x1.67a96b93f0856p-7, .Y = -0x1.6029b72d6a818p-1}}},
        TwinOracle{
            .Id = 1592,
            .First =
                Bezier::Quadratic{
                    .Start   = {.X = 0x1.3f13caecfa780p-3, .Y = 0x1.cacff3fae4404p-2},
                    .Control = {.X = -0x1.855bd6d470dc0p-6, .Y = -0x1.2223c75b5a180p-4},
                    .End     = {.X = 0x1.5409e945e6cd4p-2, .Y = -0x1.2ba92fa88fb78p-2}},
            .Second =
                Bezier::Cubic{.Start = {.X = 0x1.bce7597032a0ep-4, .Y = 0x1.4f330c4c20667p-1},
                              .Out = {.X = 0x1.9ab35bb31d4e9p-4, .Y = TwinBits(0x3fd4600dcce964e6)},
                              .In  = {.X = 0x1.76d4708bf5b7fp-4, .Y = -0x1.261977caf9f50p-6},
                              .End = {.X = 0x1.514a97fabbdd0p-4, .Y = -0x1.6ac1444990fc3p-2}},
            .Within = 0x1.b7cdfd9d7bdbbp-33,
            .Points = {Point{.X = 0x1.88592ac504f29p-4, .Y = 0x1.33a0a91118fe7p-3},
                       Point{.X = 0x1.887e9058f0072p-4, .Y = 0x1.366ec6ca87d0bp-3}}},
        TwinOracle{
            .Id = 554,
            .First =
                Bezier::Quadratic{.Start   = {.X = 0x1.7b6389a447060p-3, .Y = 0x1.65a4f7fdfe180p-6},
                                  .Control = {.X = 0x1.5238d015b0c84p-1, .Y = 0x1.943b0708ebbd8p-3},
                                  .End = {.X = -0x1.75cddf12d0502p-1, .Y = 0x1.81e209c232aa8p-3}},
            .Second =
                Bezier::Quadratic{.Start   = {.X = 0x1.4fda40dee0bc8p-2, .Y = 0x1.dd904e78911f0p-6},
                                  .Control = {.X = 0x1.7514a6f8d584ep-7, .Y = 0x1.2bd21b3f67170p-2},
                                  .End = {.X = -0x1.5bd4da903f1b4p-2, .Y = 0x1.3c1288b10ed44p-4}},
            .Within = 0x1.b7cdfd9d7bdbbp-33,
            .Points = {Point{.X = 0x1.c6f6b25cee459p-10, .Y = 0x1.6267446a2801ep-3},
                       Point{.X = 0x1.f23403ba12929p-9, .Y = 0x1.621962d887508p-3},
                       Point{.X = 0x1.1ec47320478d1p-2, .Y = 0x1.0f11215d59e4ep-4}}},
        TwinOracle{
            .Id    = 2577,
            .First = Bezier::Cubic{.Start = {.X = 0x1.d92d3392586c8p-1, .Y = -0x1.d9f6ecdbc5d6ep-1},
                                   .Out   = {.X = -0x1.96a958b709878p-2, .Y = 0x1.1fc42836dfd1cp-2},
                                   .In    = {.X = 0x1.bc6fa99d89be0p-1, .Y = -0x1.363e33aca6b00p-6},
                                   .End = {.X = -0x1.ed12528679780p-4, .Y = 0x1.3f8492802e680p-7}},
            .Second =
                Bezier::Cubic{.Start = {.X = 0x1.3b566cedb4e14p-1, .Y = -0x1.2ec4e592a450ap-3},
                              .Out   = {.X = 0x1.9c6601b04d372p-2, .Y = 0x1.388300c7ecb3ap-4},
                              .In    = {.X = 0x1.7c41fa16a9a98p-3, .Y = 0x1.483fba7a44facp-4},
                              .End   = {.X = -0x1.2111a09bfb24fp-5, .Y = -0x1.1729cf071fe5ep-3}},
            .Within = 0x1.b7cdfd9d7bdbbp-33,
            .Points = {Point{.X = 0x1.b594b8e3a5f16p-3, .Y = 0x1.e2920bbc7b59fp-7},
                       Point{.X = 0x1.2c0da0965e4fdp-2, .Y = 0x1.7b447f2c27371p-6},
                       Point{.X = 0x1.2c785b31332dap-2, .Y = 0x1.7b264a2839a9cp-6}}},
        TwinOracle{
            .Id = 1935,
            .First =
                Bezier::Quadratic{
                    .Start   = {.X = 0x1.68e957923235cp-1, .Y = 0x1.7bdda4a0a36d0p-1},
                    .Control = {.X = -0x1.9f54f2ff1960cp-2, .Y = -0x1.ad7b00fb5aac0p-1},
                    .End     = {.X = -0x1.0d17464a250f4p-1, .Y = -0x1.414f7a97e01b4p-2}},
            .Second =
                Bezier::Cubic{.Start = {.X = 0x1.89837ac91db66p-2, .Y = -0x1.7e8606b0e5956p-4},
                              .Out   = {.X = -0x1.994d836551b05p-6, .Y = -0x1.5f900979d5b8bp-4},
                              .In    = {.X = -0x1.d9b9b43612ec5p-3, .Y = -0x1.2598c67fc89b8p-2},
                              .End   = {.X = -0x1.e0a76ae12fb45p-3, .Y = -0x1.646010ff66968p-1}},
            .Within = 0x1.b7cdfd9d7bdbbp-33,
            .Points = {Point{.X = -0x1.3e8e7f62d3753p-4, .Y = -0x1.e8ef203c81712p-3},
                       Point{.X = -0x1.3b9645be05152p-4, .Y = -0x1.e77d0aa9068d5p-3}}},
        TwinOracle{.Id = 1,
                   .First =
                       Bezier::Quadratic{
                           .Start   = {.X = 0x1.324b4ca059680p-1, .Y = -0x1.648b51cd7ccf0p-2},
                           .Control = {.X = -0x1.09bc11674782cp-1, .Y = 0x1.4a4c15f846d90p-1},
                           .End     = {.X = -0x1.8a2d0bda95184p-1, .Y = -0x1.855630602aba0p-5}},
                   .Second =
                       Bezier::Quadratic{
                           .Start   = {.X = 0x1.da460bffb0ff0p-4, .Y = 0x1.c8052e79c24d0p-3},
                           .Control = {.X = -0x1.6b098738d891cp-1, .Y = 0x1.0b9288a93f328p-2},
                           .End     = {.X = 0x1.44f1009f00000p-10, .Y = 0x1.678bb223d516cp-1}},
                   .Within = 0x1.b7cdfd9d7bdbbp-33,
                   .Points = {}},
        TwinOracle{
            .Id    = 2041,
            .First = Bezier::Cubic{.Start = {.X = 0x1.f6e5311300778p-3, .Y = -0x1.a0c7d654498e0p-2},
                                   .Out = {.X = -0x1.ee62ce294ce7cp-2, .Y = -0x1.2bed22a5b1646p-1},
                                   .In  = {.X = -0x1.f97dbfd132c00p-9, .Y = 0x1.ea78146ce0034p-1},
                                   .End = {.X = -0x1.33b4542850470p-1, .Y = 0x1.8eaf4e97e6d70p-3}},
            .Second =
                Bezier::Cubic{.Start = {.X = 0x1.15ec5ef88de88p-2, .Y = -0x1.cecf2fd39bbc8p-3},
                              .Out   = {.X = 0x1.9edaaf38e8f60p-2, .Y = -0x1.0540908f16254p-1},
                              .In    = {.X = 0x1.3a95374e2a91cp-1, .Y = -0x1.450df8b04ffc8p-3},
                              .End   = {.X = -0x1.9e387c4a1b8c0p-6, .Y = 0x1.3eddf8a8aa036p-1}},
            .Within = 0x1.b7cdfd9d7bdbbp-33,
            .Points = {}},
        TwinOracle{
            .Id    = 6176,
            .First = Bezier::Cubic{.Start = {.X = -0x1.d895bde28560cp-1, .Y = 0x1.2bb8905e2f52bp+0},
                                   .Out   = {.X = 0x1.66a4c7d78aceap-1, .Y = -0x1.142adc23a81d0p+0},
                                   .In    = {.X = -0x1.55178bd399e57p+0, .Y = 0x1.abd5651ab046cp+0},
                                   .End = {.X = 0x1.ba99ae3948b61p-2, .Y = -0x1.5339e0ec1734ep-1}},
            .Second =
                Bezier::Cubic{.Start = {.X = 0x1.1ad7aa9cfa7e8p-2, .Y = 0x1.5bee0980a6c02p-1},
                              .Out   = {.X = 0x1.05657e12b99b8p-2, .Y = 0x1.86ec8b39c751cp-1},
                              .In    = {.X = -0x1.d3cbdb6916ac8p-1, .Y = -0x1.0ba41421644f0p-2},
                              .End   = {.X = -0x1.506acddd8cbd8p-1, .Y = -0x1.8ab81cf5187e0p-1}},
            .Within = 0x1.b7cdfd9d7bdbbp-33,
            .Points = {Point{.X = -0x1.029535ddf326bp-2, .Y = 0x1.cd9309022831ep-3},
                       Point{.X = -0x1.00693fa0190e4p-2, .Y = 0x1.d27e391608924p-3},
                       Point{.X = -0x1.fb8009479c3a4p-3, .Y = 0x1.d881f4d49cc40p-3}}},
    };
    int sure = 0;
    for (const TwinOracle &oracle : corpus)
    {
      for (const bool swapped : {false, true})
      {
        const Curve        &first     = swapped ? oracle.Second : oracle.First;
        const Curve        &second    = swapped ? oracle.First : oracle.Second;
        const Tolerance     tolerance = {.Length = oracle.Within, .Radians = 0x1.b7cdfd9d7bdbbp-34};
        const Intersections meeting   = Intersect(first, second, tolerance);
        std::vector<Point>  found;
        for (const Hit &hit : meeting.Hits)
        {
          CHECK(hit.Kind == Contact::Cross);
          CHECK(TwinProved({first, second}, hit));
          found.push_back(hit.At);
        }
        const auto order = [](Point lhs, Point rhs) -> bool
        {
          return lhs.X < rhs.X || (Exactly(lhs.X, rhs.X) && lhs.Y < rhs.Y);
        };
        std::vector<Point> expected = oracle.Points;
        std::ranges::sort(found, order);
        std::ranges::sort(expected, order);
        INFO("case " << oracle.Id << (swapped ? ", swapped" : ""));
        CHECK(found == expected);
        sure += TwinCertain(first, second, tolerance) ? 1 : 0;
      }
    }
    MESSAGE("sure of " << sure << " of " << 2 * corpus.size());
    CHECK(sure == static_cast<int>(2 * corpus.size()));
  }

  TEST_CASE("twin: Bézier curves as met thoroughly, to the bit, or proved better")
  {
    const int          rounds = Rounds(20000);
    Draws              draws{0x7717D5EEDU};
    std::array<int, 3> count = {};
    for (int round = 0; round < rounds; ++round)
    {
      const TwinPosed posed   = TwinPosing(draws);
      const Curves    curves  = {posed.First, posed.Second};
      const TwinMet   met     = TwinBoth(curves, posed.Tolerance);
      count[0]               += static_cast<int>(met.Sure);
      count[1]               += static_cast<int>(!met.Agrees);
      if (!met.Agrees)
      {
        INFO("round " << round);
        CHECK(TwinBetter(curves, met));
      }
      count[2] += static_cast<int>(met.Sure && !met.Swiftly.Hits.empty());
    }
    MESSAGE("sure of "
            << count[0]
            << " of "
            << rounds
            << ", crossing in "
            << count[2]
            << ", proved better than the thorough way's in "
            << count[1]);
    CHECK(count[0] > rounds / 3);
    CHECK(count[2] > rounds / 8);
  }

  TEST_CASE("twin: numbers past the range, curves of no size or both straight: not sure")
  {
    const Tolerance tolerance = {.Length = 2e-10, .Radians = 1e-10};
    const Curve     arch      = Bezier::Quadratic{
        .Start = {.X = -1, .Y = 0}, .Control = {.X = 0, .Y = 2}, .End = {.X = 1, .Y = 0}};
    const Curve dip = Bezier::Cubic{.Start = {.X = -1, .Y = 1},
                                    .Out   = {.X = -0.5, .Y = -1},
                                    .In    = {.X = 0.5, .Y = -1},
                                    .End   = {.X = 1, .Y = 1}};
    CHECK(TwinCertain(arch, dip, tolerance));
    CHECK(!TwinCertain(arch,
                       Bezier::Cubic{.Start = {.X = -1, .Y = 1},
                                     .Out   = {.X = -0.5, .Y = 1e300},
                                     .In    = {.X = 0.5, .Y = -1},
                                     .End   = {.X = 1, .Y = 1}},
                       tolerance));
    const Curve point = Bezier::Quadratic{
        .Start = {.X = 0, .Y = 1}, .Control = {.X = 0, .Y = 1}, .End = {.X = 0, .Y = 1}};
    CHECK(!TwinCertain(arch, point, tolerance));
    const Curve level = Bezier::Quadratic{
        .Start = {.X = -1, .Y = 0.5}, .Control = {.X = 0, .Y = 0.5}, .End = {.X = 1, .Y = 0.5}};
    const Curve upright = Bezier::Cubic{.Start = {.X = 0.25, .Y = -1},
                                        .Out   = {.X = 0.25, .Y = 0},
                                        .In    = {.X = 0.25, .Y = 0.5},
                                        .End   = {.X = 0.25, .Y = 1}};
    CHECK(!TwinCertain(level, upright, tolerance));
    CHECK(TwinCertain(level, dip, tolerance));
  }

  TEST_CASE("twin: curves that touch, lie along each other, cross at an end or a cusp: not sure")
  {
    const Tolerance tolerance = {.Length = 2e-10, .Radians = 1e-10};
    const Curve     arch      = Bezier::Quadratic{
        .Start = {.X = -1, .Y = 0}, .Control = {.X = 0, .Y = 2}, .End = {.X = 1, .Y = 0}};
    const Curve cap = Bezier::Quadratic{
        .Start = {.X = -1, .Y = 2}, .Control = {.X = 0, .Y = 0}, .End = {.X = 1, .Y = 2}};
    CHECK(!TwinCertain(arch, cap, tolerance));
    CHECK(!TwinCertain(arch, arch, tolerance));
    const Curve shifted = Bezier::Quadratic{.Start   = {.X = -1, .Y = 1e-10},
                                            .Control = {.X = 0, .Y = 2 + 1e-10},
                                            .End     = {.X = 1, .Y = 1e-10}};
    CHECK(!TwinCertain(arch, shifted, tolerance));
    const Curve raised = Bezier::Cubic{.Start = {.X = -1, .Y = 0},
                                       .Out   = {.X = -1.0 / 3.0, .Y = 4.0 / 3.0},
                                       .In    = {.X = 1.0 / 3.0, .Y = 4.0 / 3.0},
                                       .End   = {.X = 1, .Y = 0}};
    CHECK(!TwinCertain(arch, raised, tolerance));
    /**
     *  Two low arches a little farther apart than the kernel keeps what does not cross: pieces
     *  short enough to tell them apart all along are more than it halves.
     */
    const Curve low = Bezier::Quadratic{
        .Start = {.X = -1, .Y = 0}, .Control = {.X = 0, .Y = 0.2}, .End = {.X = 1, .Y = 0}};
    const Curve above = Bezier::Quadratic{.Start   = {.X = -1, .Y = 2e-7},
                                          .Control = {.X = 0, .Y = 0.2 + 2e-7},
                                          .End     = {.X = 1, .Y = 2e-7}};
    CHECK(Intersect(low, above, tolerance).Hits.empty());
    CHECK(!TwinCertain(low, above, tolerance));
    const Curve post = Bezier::Quadratic{
        .Start = {.X = -1, .Y = -1}, .Control = {.X = -1, .Y = 0}, .End = {.X = -1, .Y = 1}};
    CHECK(!TwinCertain(arch, post, tolerance));
    const Curve next = Bezier::Cubic{.Start = {.X = 1, .Y = 0},
                                     .Out   = {.X = 1.5, .Y = -1},
                                     .In    = {.X = 2, .Y = 1},
                                     .End   = {.X = 3, .Y = 0}};
    CHECK(!TwinCertain(arch, next, tolerance));
    const Curve cusp    = Bezier::Cubic{.Start = {.X = -1, .Y = 0},
                                        .Out   = {.X = 1, .Y = 1},
                                        .In    = {.X = -1, .Y = 1},
                                        .End   = {.X = 1, .Y = 0}};
    const Curve through = Bezier::Quadratic{
        .Start = {.X = 0, .Y = 0}, .Control = {.X = 0, .Y = 0.75}, .End = {.X = 0, .Y = 1.5}};
    CHECK(!TwinCertain(cusp, through, tolerance));
  }

  TEST_CASE("twin: nine crossings, more than the kernel is sure of")
  {
    const Tolerance tolerance = {.Length = 2e-10, .Radians = 1e-10};
    const Curve     one       = Bezier::Cubic{.Start = {.X = -0.7, .Y = -0.8},
                                              .Out   = {.X = -0.7 / 3.0, .Y = 4},
                                              .In    = {.X = 0.7 / 3.0, .Y = -4},
                                              .End   = {.X = 0.7, .Y = 0.8}};
    const Curve     other     = Bezier::Cubic{.Start = {.X = -0.8, .Y = -0.7},
                                              .Out   = {.X = 4, .Y = -0.7 / 3.0},
                                              .In    = {.X = -4, .Y = 0.7 / 3.0},
                                              .End   = {.X = 0.8, .Y = 0.7}};
    CHECK(Intersect(one, other, tolerance).Hits.size() == 9);
    CHECK(!TwinCertain(one, other, tolerance));
  }

  TEST_CASE("twin: a shallow crossing, the second curve hooked back near it")
  {
    /**
     *  An arch crossed at a shallow angle by a cubic whose end curls back over the crossing: the
     *  stretch of the arch within the tolerance of the cubic runs some way, and a piece of the
     *  cubic the thorough way weighs against it for lying along may reach the end, where the curl
     *  comes within its reach.
     */
    const Tolerance tolerance = {.Length = 2e-10, .Radians = 1e-10};
    const Curve     arch      = Bezier::Quadratic{
        .Start = {.X = -1, .Y = 0}, .Control = {.X = 0, .Y = 0.002}, .End = {.X = 1, .Y = 0}};
    const Curve near =
        Bezier::Cubic{.Start = {.X = -0x1.999999999999ap-1, .Y = 0x1.a36e2eb1c432cp-13},
                      .Out   = {.X = 1, .Y = 0x1.0624dd2f1a9fcp-9},
                      .In    = {.X = 1, .Y = 0x1.a1cac083126e9p-5},
                      .End   = {.X = -0x1.2b998b4b75686p-1, .Y = 0x1.5907b3b0660b5p-11}};
    CHECK(Intersect(arch, near, tolerance).Hits.size() == 1);
    CHECK(!TwinCertain(arch, near, tolerance));
    const Curve clear =
        Bezier::Cubic{.Start = {.X = -0x1.999999999999ap-1, .Y = -0x1.6f0068db8bac8p-10},
                      .Out   = {.X = 1, .Y = 0x1.0624dd2f1a9fcp-8},
                      .In    = {.X = 1, .Y = 0x1.a1cac083126e9p-5},
                      .End   = {.X = -0x1.7bb57f6268633p-2, .Y = 0x1.c475efdb16173p-11}};
    const TwinMet met = TwinBoth({arch, clear}, tolerance);
    CHECK(met.Sure);
    CHECK(met.Agrees);
    CHECK(met.Swiftly.Hits.size() == 1);
  }

  TEST_CASE("twin: curves touching where no halving falls: not sure")
  {
    const Tolerance tolerance = {.Length = 2e-10, .Radians = 1e-10};
    const Curve     arch      = Bezier::Quadratic{
        .Start = {.X = -1, .Y = 0}, .Control = {.X = 0, .Y = 2}, .End = {.X = 1, .Y = 0}};
    const double third = 1.0 / 3.0;
    const Point  spot  = At(arch, third);
    const Vector way   = Derivative(arch, third) / Length(Derivative(arch, third));
    const Vector side  = Perpendicular(way) * 0.3;
    const Curve  under = Bezier::Quadratic{.Start   = spot - (way * 0.5) - side,
                                           .Control = spot + side,
                                           .End     = spot + (way * 0.5) - side};
    CHECK(!TwinCertain(arch, under, tolerance));
  }
} // namespace Geoxy
