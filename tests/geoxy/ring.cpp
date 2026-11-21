/***************************************************************************************************
 *  @file      tests/geoxy/ring.cpp
 *  @brief     Tests for the ring kernel: an arc or an ellipse met with a Bézier curve or another
 *             conic, each crossing the double nearest the exact one, and every answer the thorough
 *             way's to the bit or none.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "swift/ring.hpp"

#include "draws.hpp"
#include "swift/met.hpp"
#include "swift/section.hpp"
#include "thorough.hpp"
#include "waves.hpp"

#include "geoxy/angle.hpp"
#include "geoxy/arc.hpp"
#include "geoxy/at.hpp"
#include "geoxy/bezier/cubic.hpp"
#include "geoxy/bezier/quadratic.hpp"
#include "geoxy/box.hpp"
#include "geoxy/curve.hpp"
#include "geoxy/derivative.hpp"
#include "geoxy/domain.hpp"
#include "geoxy/ellipse.hpp"
#include "geoxy/intersect.hpp"
#include "geoxy/point.hpp"
#include "geoxy/spline.hpp"
#include "geoxy/tolerance.hpp"
#include "geoxy/vector.hpp"

#include <doctest/doctest.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <type_traits>
#include <utility>
#include <variant>
#include <vector>

namespace Geoxy
{
  namespace
  {
    /**
     *  Whether two meetings are the same to the bit, hit for hit and stretch for stretch.
     *
     *  @param  one   One.
     *  @param  other The other.
     *  @return       Whether they are.
     */
    auto RingSame(const Intersections &one, const Intersections &other) -> bool
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
    struct RingMet
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
    auto RingBoth(const Curves &curves, const Tolerance &tolerance) -> RingMet
    {
      RingMet                     met;
      const std::optional<double> length = Reachable(curves, tolerance);
      if (!length)
      {
        return met;
      }
      met.Sure = Swift::Met(curves, {.Length = *length, .Radians = tolerance.Radians}, met.Swiftly);
      if (met.Sure)
      {
        Thorough(curves, tolerance, met.Thorough);
        met.Agrees = RingSame(met.Swiftly, met.Thorough);
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
    auto RingCertain(const Curve &first, const Curve &second, const Tolerance &tolerance) -> bool
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
    auto RingUnit(Draws &draws) -> double
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
    auto RingPick(Draws &draws, std::uint64_t count) -> std::uint64_t
    {
      return draws.Next() % count;
    }

    /**
     *  An ellipse at random: skewed, mirrored, flat, a circle at times, whole or an arc.
     *
     *  @param  draws Where the numbers come from.
     *  @return       The ellipse.
     */
    auto RingOval(Draws &draws) -> Curve
    {
      const double wide = 0.2 + (0.8 * std::abs(RingUnit(draws)));
      const double high = RingPick(draws, 8) == 0
                            ? wide * std::ldexp(1.0, -3 - static_cast<int>(RingPick(draws, 8)))
                            : 0.05 + (0.95 * std::abs(RingUnit(draws)));
      const double turn = 3.0 * RingUnit(draws);
      const double skew = RingPick(draws, 3) == 0 ? 0.6 * RingUnit(draws) : 0.0;
      Ellipse      oval = {
          .Map   = {.X      = {.X = wide * std::cos(turn), .Y = wide * std::sin(turn)},
                    .Y      = {.X = -high * std::sin(turn + skew), .Y = high * std::cos(turn + skew)},
                    .Origin = {.X = 0.5 * RingUnit(draws), .Y = 0.5 * RingUnit(draws)}},
          .Start = {},
          .Sweep = Turn};
      if (RingPick(draws, 4) == 0)
      {
        std::swap(oval.Map.X, oval.Map.Y);
      }
      if (RingPick(draws, 8) == 0)
      {
        oval.Map.Y = Perpendicular(oval.Map.X);
      }
      if (RingPick(draws, 2) == 0)
      {
        oval.Start = Degrees(std::round(180.0 * RingUnit(draws)));
        oval.Sweep = RingPick(draws, 2) == 0 ? Degrees(std::round(350.0 * RingUnit(draws)))
                                             : Radians(6.0 * RingUnit(draws));
      }
      return oval;
    }

    /**
     *  An arc at random: whole at times, its angles in degrees or radians.
     *
     *  @param  draws Where the numbers come from.
     *  @return       The arc.
     */
    auto Round(Draws &draws) -> Curve
    {
      const Point  center = {.X = 0.5 * RingUnit(draws), .Y = 0.5 * RingUnit(draws)};
      const double radius = 0.1 + (0.9 * std::abs(RingUnit(draws)));
      if (RingPick(draws, 4) == 0)
      {
        return Arc{.Center = center, .Radius = radius, .Start = {}, .Sweep = Turn};
      }
      return RingPick(draws, 2) == 0 ? Arc{.Center = center,
                                           .Radius = radius,
                                           .Start  = Degrees(std::round(180.0 * RingUnit(draws))),
                                           .Sweep  = Degrees(std::round(350.0 * RingUnit(draws)))}
                                     : Arc{.Center = center,
                                           .Radius = radius,
                                           .Start  = Radians(3.0 * RingUnit(draws)),
                                           .Sweep  = Radians(6.0 * RingUnit(draws))};
    }

    /**
     *  An arc or an ellipse at random.
     *
     *  @param  draws Where the numbers come from.
     *  @return       The curve.
     */
    auto Conic(Draws &draws) -> Curve
    {
      return RingPick(draws, 3) == 0 ? Round(draws) : RingOval(draws);
    }

    /**
     *  A quadratic or cubic Bézier curve at random in the square of side two about the origin.
     *
     *  @param  draws Where the numbers come from.
     *  @return       The curve.
     */
    auto RingBent(Draws &draws) -> Curve
    {
      const auto point = [&draws]() -> Point
      {
        return {.X = RingUnit(draws), .Y = RingUnit(draws)};
      };
      if (RingPick(draws, 3) == 0)
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
    auto RingThrough(Draws &draws, const std::pair<Point, Vector> &spot, double bend) -> Curve
    {
      const auto [point, way] = spot;
      const Vector along      = way * ((0.3 + (0.5 * std::abs(RingUnit(draws)))) / Length(way));
      const Vector side       = Perpendicular(along) * (bend / Length(along));
      if (RingPick(draws, 2) == 0)
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
    auto RingMoved(const Curve &curve, double scale, Vector by) -> Curve
    {
      const auto shifted = [scale, by](Point point) -> Point
      {
        return {.X = (point.X * scale) + by.X, .Y = (point.Y * scale) + by.Y};
      };
      return std::visit(
          [&shifted, scale](const auto &shape) -> Curve
          {
            using Shape = std::decay_t<decltype(shape)>;
            if constexpr (std::is_same_v<Shape, Arc>)
            {
              return Arc{.Center = shifted(shape.Center),
                         .Radius = shape.Radius * scale,
                         .Start  = shape.Start,
                         .Sweep  = shape.Sweep};
            }
            else if constexpr (std::is_same_v<Shape, Ellipse>)
            {
              Ellipse moved    = shape;
              moved.Map.X      = shape.Map.X * scale;
              moved.Map.Y      = shape.Map.Y * scale;
              moved.Map.Origin = shifted(shape.Map.Origin);
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
            else if constexpr (std::is_same_v<Shape, Spline>)
            {
              Spline moved = shape;
              for (Point &pole : moved.Poles)
              {
                pole = shifted(pole);
              }
              return moved;
            }
            else
            {
              return shape;
            }
          },
          curve);
    }

    /**
     *  A Bézier curve near a curve at an edge of what the kernel decides: through it at a slight
     *  angle, along its tangent a few tolerances off, through one of its ends, or through it
     *  plainly.
     *
     *  @param  draws Where the numbers come from.
     *  @param  curve The curve.
     *  @param  near  The tolerance's length.
     *  @return       The Bézier curve.
     */
    auto RingEdged(Draws &draws, const Curve &curve, double near) -> Curve
    {
      const auto   kind  = RingPick(draws, 4);
      const double end   = RingPick(draws, 2) == 0 ? 0.0 : 1.0;
      const double along = kind == 2 ? end : 0.5 + (0.45 * RingUnit(draws));
      const Point  spot  = At(curve, along);
      Vector       way   = Derivative(curve, along);
      if (!(Length(way) > 0.0))
      {
        way = {.X = 1.0, .Y = 0.0};
      }
      const double angle  = kind == 0
                              ? std::ldexp(RingUnit(draws), -static_cast<int>(RingPick(draws, 30)))
                              : 3.0 * RingUnit(draws);
      const Vector turned = {.X = (way.X * std::cos(angle)) - (way.Y * std::sin(angle)),
                             .Y = (way.X * std::sin(angle)) + (way.Y * std::cos(angle))};
      const double off    = kind == 1 ? near * 16.0 * RingUnit(draws) : 0.0;
      const Point  moved  = spot + (Perpendicular(way) * (off / Length(way)));
      const double bend   = 0.3 * RingUnit(draws);
      return RingThrough(draws, {moved, kind == 1 ? way : turned}, bend);
    }

    /**
     *  An arc or an ellipse near a curve at an edge of what the kernel decides: through one of its
     *  points at an angle, slight or not, or along it a few tolerances off, or through its end.
     *
     *  @param  draws Where the numbers come from.
     *  @param  curve The curve.
     *  @param  near  The tolerance's length.
     *  @return       The conic.
     */
    auto Ringing(Draws &draws, const Curve &curve, double near) -> Curve
    {
      const auto   kind   = RingPick(draws, 4);
      const double end    = RingPick(draws, 2) == 0 ? 0.0 : 1.0;
      const Range  domain = Domain(curve);
      const double along =
          domain.Start
          + ((domain.End - domain.Start) * (kind == 2 ? end : 0.5 + (0.45 * RingUnit(draws))));
      const Point spot = At(curve, along);
      Vector      way  = Derivative(curve, along);
      if (!(Length(way) > 0.0))
      {
        way = {.X = 1.0, .Y = 0.0};
      }
      const double angle  = kind == 0
                              ? std::ldexp(RingUnit(draws), -static_cast<int>(RingPick(draws, 30)))
                              : 3.0 * RingUnit(draws);
      const Vector turned = {.X = (way.X * std::cos(angle)) - (way.Y * std::sin(angle)),
                             .Y = (way.X * std::sin(angle)) + (way.Y * std::cos(angle))};
      const Vector unit   = (kind == 1 ? way : turned) / Length(kind == 1 ? way : turned);
      const double off    = kind == 1 ? near * 16.0 * RingUnit(draws) : 0.0;
      const double radius = 0.2 + std::abs(RingUnit(draws));
      const Vector inward = Perpendicular(unit) * (RingPick(draws, 2) == 0 ? radius : -radius);
      const Point  center = spot + (Perpendicular(unit) * off) + inward;
      const Vector toward = spot - center;
      const double start  = std::atan2(toward.Y, toward.X) - (0.5 + std::abs(RingUnit(draws)));
      if (RingPick(draws, 2) == 0)
      {
        return Arc{.Center = center,
                   .Radius = Length(toward),
                   .Start  = Radians(start),
                   .Sweep  = Radians(2.0 + RingUnit(draws))};
      }
      const double ratio = 0.5 + (0.5 * std::abs(RingUnit(draws)));
      const Vector major = toward;
      const Vector minor = Perpendicular(toward) * ratio;
      return Ellipse{.Map   = {.X = major, .Y = minor, .Origin = center},
                     .Start = Radians(-0.5 - std::abs(RingUnit(draws))),
                     .Sweep = Radians(2.0 + RingUnit(draws))};
    }

    /**
     *  A pair of curves and the tolerance asked for.
     */
    struct RingPosed
    {
      Curve            First     = {};    /**< The first curve.         */
      Curve            Second    = {};    /**< The second.              */
      Geoxy::Tolerance Tolerance = {};    /**< The tolerance asked for. */
      bool             Afield    = false; /**< Whether a spline moved 1e8 off: there the thorough
                                               way's polish misses the double nearest a crossing
                                               the kernel proves, as its own test shows. */
    };

    /**
     *  What the other curve of a pair is: a Bézier curve met with a conic, a conic met with a
     *  conic, or a conic met with a spline.
     */
    enum class Partner : std::uint8_t
    {
      Bent,   /**< A Bézier curve with a conic.                         */
      Conic,  /**< A conic with a conic.                                */
      Spline, /**< A conic with a spline.                               */
      Rated,  /**< A conic with a rational one.                         */
      Broad,  /**< A conic with a spline of degree four to eight.       */
      Heavy,  /**< A conic with a rational one of degree four to eight. */
    };

    /**
     *  A spline as RingPosing draws one: of degree two or three, or four to eight for a broad pair,
     *  its weights from a quarter to three for a rational one.
     *
     *  @param  draws Where the numbers come from.
     *  @param  kind  What the pair is.
     *  @return       The spline.
     */
    auto RingSplined(Draws &draws, Partner kind) -> Curve
    {
      const bool broad  = kind == Partner::Broad || kind == Partner::Heavy;
      Spline     spline = broad ? Tests::Splined(draws, 4, 4) : Tests::Splined(draws, 2, 1);
      if (kind == Partner::Rated || kind == Partner::Heavy)
      {
        spline.Weights.resize(spline.Poles.size());
        for (double &weight : spline.Weights)
        {
          weight =
              std::ldexp(1.0 + (0.5 * RingUnit(draws)), static_cast<int>(RingPick(draws, 3)) - 1);
        }
      }
      return spline;
    }

    /**
     *  A pair in one of the ways the kernel must answer as the thorough way does: at random or at
     *  an edge, far off or scaled, under a coarse angle at times, either way round.
     *
     *  @param  draws Where the numbers come from.
     *  @param  kind  What the pair is.
     *  @return       The pair.
     */
    auto RingPosing(Draws &draws, Partner kind) -> RingPosed
    {
      const bool   strung = kind != Partner::Bent && kind != Partner::Conic;
      const Curve  one    = strung ? RingSplined(draws, kind) : Conic(draws);
      const bool   coarse = RingPick(draws, 4) == 0;
      const double near = coarse ? std::pow(10.0, -3.0 - (9.0 * std::abs(RingUnit(draws)))) : 2e-10;
      const bool   wide = RingPick(draws, 6) == 0;
      const double radians = wide ? std::pow(10.0, -6.0 * std::abs(RingUnit(draws))) : 1e-10;
      const bool   random  = RingPick(draws, 2) == 0;
      const auto   drawn   = [&draws, &one, kind, random, near]() -> Curve
      {
        if (kind == Partner::Bent)
        {
          return random ? RingBent(draws) : RingEdged(draws, one, near);
        }
        return random ? Conic(draws) : Ringing(draws, one, near);
      };
      const Curve  other  = drawn();
      const auto   moving = RingPick(draws, 6);
      const double scale = std::ldexp(moving == 0 ? 1.37 : 1.0,
                                      moving == 0 ? static_cast<int>(RingPick(draws, 80)) - 40 : 0);
      const Vector by    = {.X = moving == 1 ? std::ldexp(RingUnit(draws), 26) : 0.0,
                            .Y = moving == 1 ? 1e8 : 0.0};
      const bool   swapped = RingPick(draws, 2) == 0;
      const Curve  first   = RingMoved(one, scale, by);
      const Curve  second  = RingMoved(other, scale, by);
      return {.First     = swapped ? second : first,
              .Second    = swapped ? first : second,
              .Tolerance = {.Length = near * scale, .Radians = radians},
              .Afield    = strung && moving == 1};
    }

    /**
     *  Pairs as RingPosing draws them: how many the kernel was sure of, how many it got wrong, and
     *  how many of those it was sure of cross.
     *
     *  @param  kind   What the pairs are.
     *  @param  rounds How many pairs.
     *  @return        How many it was sure of, how many disagreed, how many crossing.
     */
    auto RingDrawn(Partner kind, int rounds) -> std::array<int, 3>
    {
      constexpr std::array<std::uint64_t, 6> Seeds = {0xC0B1CU, 0xE11195EU, 0x5B11E5U,
                                                      0x2A7EDU, 0xB10ADU,   0x4EA7EU};
      Draws                                  draws{Seeds.at(static_cast<std::size_t>(kind))};
      std::array<int, 3>                     count = {};
      for (int round = 0; round < rounds; ++round)
      {
        const RingPosed posed   = RingPosing(draws, kind);
        const Curves    curves  = {posed.First, posed.Second};
        const RingMet   met     = RingBoth(curves, posed.Tolerance);
        count[0]               += static_cast<int>(met.Sure);
        count[1]               += static_cast<int>(!met.Agrees && !posed.Afield);
        count[2]               += static_cast<int>(met.Sure && !met.Swiftly.Hits.empty());
        if (!met.Agrees && !posed.Afield)
        {
          INFO("round " << round);
          CHECK(met.Agrees);
        }
      }
      return count;
    }
  } // namespace

  TEST_CASE("ring: arcs and ellipses with Bézier curves as met thoroughly, to the bit")
  {
    const int rounds                    = Rounds(20000);
    const auto [sure, differ, crossing] = RingDrawn(Partner::Bent, rounds);
    MESSAGE("sure of " << sure << " of " << rounds << ", crossing in " << crossing);
    CHECK(differ == 0);
    CHECK(sure > rounds / 3);
  }

  TEST_CASE("ring: arcs and ellipses with polynomial splines as met thoroughly, to the bit")
  {
    const int rounds                    = Rounds(20000);
    const auto [sure, differ, crossing] = RingDrawn(Partner::Spline, rounds);
    MESSAGE("sure of " << sure << " of " << rounds << ", crossing in " << crossing);
    CHECK(differ == 0);
    CHECK(sure > rounds / 3);
  }

  TEST_CASE("ring: arcs and ellipses with rational splines as met thoroughly, to the bit")
  {
    const int rounds                    = Rounds(20000);
    const auto [sure, differ, crossing] = RingDrawn(Partner::Rated, rounds);
    MESSAGE("sure of " << sure << " of " << rounds << ", crossing in " << crossing);
    CHECK(differ == 0);
    CHECK(sure > rounds / 3);
  }

  TEST_CASE("ring: arcs and ellipses with splines of degree four to eight as met thoroughly")
  {
    const int rounds                    = Rounds(20000);
    const auto [sure, differ, crossing] = RingDrawn(Partner::Broad, rounds);
    MESSAGE("sure of " << sure << " of " << rounds << ", crossing in " << crossing);
    CHECK(differ == 0);
    CHECK(sure > rounds / 3);
  }

  TEST_CASE("ring: arcs and ellipses with rational splines of degree four to eight, to the bit")
  {
    const int rounds                    = Rounds(20000);
    const auto [sure, differ, crossing] = RingDrawn(Partner::Heavy, rounds);
    MESSAGE("sure of " << sure << " of " << rounds << ", crossing in " << crossing);
    CHECK(differ == 0);
    CHECK(sure > rounds / 3);
  }

  TEST_CASE("ring: a spline far from nothing met with the double nearest each crossing")
  {
    /**
     *  An arc and a cubic spline 1e8 off, met once: each span's control points taken about its
     *  first pole keep the digits the kernel needs; the crossing is mpmath's in 80 digits, which
     *  the thorough way's polish misses by hundreds of steps of a double.
     */
    const Curve arc = Arc{.Center = {.X = -0x1.a3e97f616e1c8p+23, .Y = 0x1.7d783ffdc586ap+26},
                          .Radius = 0x1.c60430d84408ap-3,
                          .Start  = Radians(-0x1.632a6c43fc862p-2),
                          .Sweep  = Radians(0x1.2a7273bf0d59cp+1)};
    const Curve spline =
        Spline{.Degree  = 3,
               .Knots   = {0, 0, 0, 0, 1, 2, 3, 4, 4, 4, 4},
               .Poles   = {{.X = -0x1.a3e97f631cb34p+23, .Y = 0x1.7d784009ddfcbp+26},
                           {.X = -0x1.a3e97f1b178e8p+23, .Y = 0x1.7d78400da6023p+26},
                           {.X = -0x1.a3e97ee2f26a2p+23, .Y = 0x1.7d783ffad1052p+26},
                           {.X = -0x1.a3e97e7c20a89p+23, .Y = 0x1.7d783ffd0ca74p+26},
                           {.X = -0x1.a3e97e1764542p+23, .Y = 0x1.7d78401950541p+26},
                           {.X = -0x1.a3e97dc5ebde4p+23, .Y = 0x1.7d784001e1a7ap+26},
                           {.X = -0x1.a3e97d666c27ep+23, .Y = 0x1.7d78401927ea3p+26}},
               .Weights = {}};
    for (const bool swapped : {false, true})
    {
      Intersections met;
      const Curves  curves = {swapped ? spline : arc, swapped ? arc : spline};
      REQUIRE(Swift::Met(curves,
                         Usable(curves, {.Length = 0x1.7663573df2dbfp-40, .Radians = 1e-10}), met));
      REQUIRE(met.Hits.size() == 1);
      const Hit &hit = met.Hits.front();
      CHECK(hit.At == Point{.X = -0x1.a3e97f10facc1p+23, .Y = 0x1.7d784007c7b4bp+26});
      CHECK((swapped ? hit.First.Along : hit.Second.Along) == 0x1.0516fd37a9135p-1);
    }
  }

  TEST_CASE("ring: arcs and ellipses with each other as met thoroughly, to the bit")
  {
    const int rounds                    = Rounds(20000);
    const auto [sure, differ, crossing] = RingDrawn(Partner::Conic, rounds);
    MESSAGE("sure of " << sure << " of " << rounds << ", crossing in " << crossing);
    CHECK(differ == 0);
    CHECK(sure > rounds / 3);
  }
} // namespace Geoxy

namespace Geoxy
{
  TEST_CASE("ring: a circle and an arch cross plainly; two arcs are the chord form's")
  {
    const Tolerance tolerance = {.Length = 2e-10, .Radians = 1e-10};
    const Curve     arch      = Bezier::Quadratic{
        .Start = {.X = -1, .Y = 0}, .Control = {.X = 0, .Y = 2}, .End = {.X = 1, .Y = 0}};
    const Curve circle =
        Arc{.Center = {.X = 0, .Y = 0.5}, .Radius = 0.7, .Start = {}, .Sweep = Turn};
    CHECK(RingCertain(circle, arch, tolerance));
    CHECK(RingCertain(arch, circle, tolerance));
    CHECK(!RingCertain(circle, circle, tolerance));
    const RingMet met = RingBoth({circle, arch}, tolerance);
    CHECK(met.Agrees);
    CHECK(met.Swiftly.Hits.size() == 2);
  }

  TEST_CASE("ring: a crossing with a coordinate exactly zero: the thorough way's")
  {
    const Tolerance tolerance = {.Length = 2e-10, .Radians = 1e-10};
    const Curve     circle    = Arc{.Center = {}, .Radius = 1, .Start = {}, .Sweep = Turn};
    const Curve     upright   = Bezier::Quadratic{
        .Start = {.X = 0, .Y = 0.25}, .Control = {.X = 0, .Y = 1}, .End = {.X = 0, .Y = 1.5}};
    CHECK(!RingCertain(circle, upright, tolerance));
    const Intersections met = Intersect(circle, upright, tolerance);
    REQUIRE(met.Hits.size() == 1);
    CHECK(met.Hits[0].At    == Point{.X = 0, .Y = 1});
  }

  TEST_CASE("ring: conics the kernel does not take, and Bézier curves out of its range")
  {
    const Tolerance tolerance = {.Length = 2e-10, .Radians = 1e-10};
    const Curve     arch      = Bezier::Quadratic{
        .Start = {.X = -1, .Y = 0}, .Control = {.X = 0, .Y = 2}, .End = {.X = 1, .Y = 0}};
    const Curve inside =
        Arc{.Center = {.X = 0, .Y = 0.5}, .Radius = -0.7, .Start = {}, .Sweep = Turn};
    CHECK(!RingCertain(inside, arch, tolerance));
    const Curve coiled = Ellipse{
        .Map   = {.X = {.X = 0.7, .Y = 0}, .Y = {.X = 0, .Y = 0.5}, .Origin = {.X = 0, .Y = 0.5}},
        .Start = {},
        .Sweep = Degrees(720)};
    CHECK(!RingCertain(coiled, arch, tolerance));
    const Curve flat = Ellipse{
        .Map = {.X = {.X = 0.7, .Y = 0}, .Y = {.X = 0.7, .Y = 1e-14}, .Origin = {.X = 0, .Y = 0.5}},
        .Start = {},
        .Sweep = Turn};
    CHECK(!RingCertain(arch, flat, tolerance));
    const Curve oval = Ellipse{
        .Map   = {.X = {.X = 0.7, .Y = 0}, .Y = {.X = 0, .Y = 0.3}, .Origin = {.X = 0, .Y = 0.5}},
        .Start = {},
        .Sweep = Turn};
    CHECK(RingCertain(oval, arch, tolerance));
    CHECK(!RingCertain(oval, coiled, tolerance));
    const Curve weighted = Spline{.Degree = 2,
                                  .Knots  = {0, 0, 0, 1, 1, 1},
                                  .Poles  = {{.X = -1, .Y = 0}, {.X = 0, .Y = 2}, {.X = 1, .Y = 0}},
                                  .Weights = {1, 2, 1}};
    const RingMet rated  = RingBoth({oval, weighted}, tolerance);
    CHECK(rated.Sure);
    CHECK(rated.Agrees);
    const Curve plain = Spline{.Degree  = 2,
                               .Knots   = {0, 0, 0, 1, 1, 1},
                               .Poles   = {{.X = -1, .Y = 0}, {.X = 0, .Y = 2}, {.X = 1, .Y = 0}},
                               .Weights = {}};
    CHECK(RingCertain(oval, plain, tolerance));
    const double far    = 0x1p45;
    const Curve  coarse = Spline{.Degree  = 2,
                                 .Knots   = {far, far, far, far + 1, far + 1, far + 1},
                                 .Poles   = {{.X = -1, .Y = 0}, {.X = 0, .Y = 2}, {.X = 1, .Y = 0}},
                                 .Weights = {}};
    CHECK(!RingCertain(oval, coarse, tolerance));
    const Curve wild = Bezier::Quadratic{
        .Start = {.X = -1, .Y = 0}, .Control = {.X = 0, .Y = 1e300}, .End = {.X = 1, .Y = 0}};
    CHECK(!RingCertain(oval, wild, tolerance));
  }

  namespace
  {
    /**
     *  The cubic Bézier curve through four points at four places, by solving for its control
     *  points.
     *
     *  @param  points The points.
     *  @param  places Their places, from 0 to 1, all different.
     *  @return        The curve.
     */
    auto Interpolated(const std::array<Point, 4> &points, const std::array<double, 4> &places)
        -> Curve
    {
      std::array<std::array<double, 6>, 4> rows = {};
      for (std::size_t row = 0; row < 4; ++row)
      {
        const double along = places.at(row);
        const double rest  = 1.0 - along;
        rows.at(row)       = {rest * rest * rest,         3.0 * rest * rest * along,
                              3.0 * rest * along * along, along * along * along,
                              points.at(row).X,           points.at(row).Y};
      }
      for (std::size_t pivot = 0; pivot < 4; ++pivot)
      {
        for (std::size_t row = 0; row < 4; ++row)
        {
          if (row != pivot)
          {
            const double times = rows.at(row).at(pivot) / rows.at(pivot).at(pivot);
            for (std::size_t column = 0; column < 6; ++column)
            {
              rows.at(row).at(column) -= times * rows.at(pivot).at(column);
            }
          }
        }
      }
      std::array<Point, 4> poles = {};
      for (std::size_t row = 0; row < 4; ++row)
      {
        poles.at(row) = {.X = rows.at(row)[4] / rows.at(row).at(row),
                         .Y = rows.at(row)[5] / rows.at(row).at(row)};
      }
      return Bezier::Cubic{.Start = poles[0], .Out = poles[1], .In = poles[2], .End = poles[3]};
    }

    /**
     *  The unit circle's point at an angle in degrees.
     *
     *  @param  degrees The angle.
     *  @return         The point.
     */
    auto On(double degrees) -> Point
    {
      const Vector way = Direction(Degrees(degrees));
      return {.X = way.X, .Y = way.Y};
    }
  } // namespace

  TEST_CASE("ring: rational splines the kernel does not cut")
  {
    const Tolerance tolerance = {.Length = 2e-10, .Radians = 1e-10};
    const Curve     oval      = Ellipse{
        .Map   = {.X = {.X = 0.7, .Y = 0}, .Y = {.X = 0, .Y = 0.3}, .Origin = {.X = 0, .Y = 0.5}},
        .Start = {},
        .Sweep = Turn};
    const std::vector<Point> poles = {{.X = -1, .Y = 0}, {.X = 0, .Y = 2}, {.X = 1, .Y = 0}};
    const Spline             faint = {
        .Degree  = 2,
        .Knots   = {0, 0, 0, 1e-310, 1, 1, 1},
        .Poles   = {{.X = -1, .Y = 0}, {.X = -0.5, .Y = 2}, {.X = 0.5, .Y = 2}, {.X = 1, .Y = 0}},
        .Weights = {1, 2, 1, 1}};
    const Spline light = {
        .Degree = 2, .Knots = {0, 0, 0, 1, 1, 1}, .Poles = poles, .Weights = {1, 1e-250, 1}};
    const Spline slight = {
        .Degree = 2, .Knots = {0, 0, 0, 1, 1, 1}, .Poles = poles, .Weights = {1, 1e-20, 1}};
    const Spline heavy = {
        .Degree = 2, .Knots = {0, 0, 0, 1, 1, 1}, .Poles = poles, .Weights = {1, 1e250, 1}};
    const Spline plain = {.Degree = 2, .Knots = {0, 0, 0, 1, 1, 1}, .Poles = poles, .Weights = {}};
    REQUIRE(Valid(faint));
    const Box around = {.Min = {.X = -2, .Y = -2}, .Max = {.X = 2, .Y = 2}};
    const Box beside = {.Min = {.X = 3, .Y = -2}, .Max = {.X = 4, .Y = 2}};
    for (const Spline &spline : {faint, light, slight, heavy})
    {
      CHECK(!Swift::Weighed(spline, around));
      CHECK(!RingCertain(oval, spline, tolerance));
    }
    const auto weighed = Swift::Weighed(plain, around);
    const auto cut     = Swift::Sectioned(plain);
    REQUIRE(weighed);
    REQUIRE(cut);
    CHECK(weighed->size() == cut->size());
    const auto missed = Swift::Weighed(plain, beside);
    REQUIRE(missed);
    CHECK(missed->empty());
    CHECK(!Swift::Sectioned(slight));
    /**
     *  A polynomial span of degree four whose knots' fraction falls below the normal doubles.
     */
    const Spline wisp = {.Degree  = 4,
                         .Knots   = {0, 0, 0, 0, 0, 1e-310, 1, 1, 1, 1, 1},
                         .Poles   = {{.X = -1, .Y = 0},
                                     {.X = -0.75, .Y = 2},
                                     {.X = -0.25, .Y = 2},
                                     {.X =  0.25, .Y = 2},
                                     {.X =  0.75, .Y = 2},
                                     {.X = 1, .Y = 0}},
                         .Weights = {}};
    REQUIRE(Valid(wisp));
    CHECK(!Swift::Weighed(wisp, around));
  }

  TEST_CASE("ring: a rational span whose control points pass the kernels' range: not sure")
  {
    /**
     *  The middle control point w (P - B) is 2^161 times 2^41 a coordinate, past 2^200, though
     *  every pole and weight is within it and the slack below every weight.
     */
    const double side = 0x1p41;
    const Curve  circle =
        Arc{.Center = {.X = 0, .Y = side}, .Radius = side, .Start = {}, .Sweep = Turn};
    const Spline heavy = {
        .Degree  = 2,
        .Knots   = {0, 0, 0, 1, 1, 1},
        .Poles   = {{.X = 0, .Y = 0}, {.X = side, .Y = side}, {.X = 0, .Y = 2 * side}},
        .Weights = {0x1p160, 0x1p161, 0x1p160}};
    const Box  around  = {.Min = {.X = -side, .Y = -side}, .Max = {.X = 2 * side, .Y = 3 * side}};
    const auto weighed = Swift::Weighed(heavy, around);
    REQUIRE(weighed);
    REQUIRE(weighed->size() == 1);
    CHECK(!(std::abs((*weighed)[0].Poles[1].X) <= 0x1p200));
    const Curve               shape = heavy;
    Intersections             into;
    const std::optional<bool> sure =
        Swift::Ringed({circle, shape}, {.Length = 0x1p-10, .Radians = 1e-10}, into);
    REQUIRE(sure);
    CHECK(!*sure);
    CHECK(into.Hits.empty());
  }

  TEST_CASE("ring: crossings the thorough way may take for one: not sure")
  {
    const Tolerance tolerance = {.Length = 1e-9, .Radians = 1e-10};
    const Curve     circle    = Arc{.Center = {}, .Radius = 1, .Start = {}, .Sweep = Turn};
    /**
     *  Three crossings, the middle one at the place halfway between the outer two, and one more
     *  near the start.
     */
    const Curve middle =
        Interpolated({On(100), On(60), On(20), Point{.X = 1.6, .Y = -1.2}}, {0.2, 0.4, 0.6, 1.0});
    CHECK(Intersect(circle, middle, tolerance).Hits.size() == 4);
    CHECK(!RingCertain(circle, middle, tolerance));
    /**
     *  Three crossings, the last at the place half a curve past halfway between the first two, and
     *  one more near the start.
     */
    const Curve round =
        Interpolated({On(150), On(80), On(-40), Point{.X = 1.6, .Y = -1.5}}, {0.2, 0.4, 0.8, 1.0});
    CHECK(Intersect(circle, round, tolerance).Hits.size() == 4);
    CHECK(!RingCertain(circle, round, tolerance));
    /**
     *  A loop whose double point lies on a circle: two crossings a rounding apart, away from the
     *  axes, where a coordinate all but zero would settle nothing.
     */
    const Curve loop = Bezier::Cubic{.Start = {.X = -0.75, .Y = 0.125},
                                     .Out   = {.X = 2.25, .Y = 2.125},
                                     .In    = {.X = -1.75, .Y = 2.125},
                                     .End   = {.X = 1.25, .Y = 0.125}};
    double      low  = 0.1;
    double      high = 0.3;
    for (int step = 0; step < 60; ++step)
    {
      const double along                      = low + ((high - low) / 2.0);
      (At(loop, along).X < 0.25 ? low : high) = along;
    }
    const Point twice = At(loop, low);
    const Curve ring  = Arc{
        .Center = {.X = twice.X, .Y = twice.Y - 0.25}, .Radius = 0.25, .Start = {}, .Sweep = Turn};
    CHECK(Intersect(ring, loop, tolerance).Hits.size() == 4);
    CHECK(!RingCertain(ring, loop, tolerance));
  }
} // namespace Geoxy
