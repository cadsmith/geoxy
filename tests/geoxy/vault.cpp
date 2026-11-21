/***************************************************************************************************
 *  @file      tests/geoxy/vault.cpp
 *  @brief     Tests for the parabola kernel: a quadratic Bézier curve met with a rational spline or
 *             one past a cubic on the curve's parabola, each crossing the double nearest the exact
 *             one, and every answer the thorough way's to the bit or none.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "swift/vault.hpp"

#include "draws.hpp"
#include "swift/met.hpp"
#include "thorough.hpp"
#include "waves.hpp"

#include "geoxy/at.hpp"
#include "geoxy/bezier/cubic.hpp"
#include "geoxy/bezier/quadratic.hpp"
#include "geoxy/curve.hpp"
#include "geoxy/derivative.hpp"
#include "geoxy/domain.hpp"
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
#include <utility>
#include <variant>
#include <vector>

namespace Geoxy
{
  namespace
  {
    /**
     *  Whether two meetings are the same to the bit, hit for hit.
     *
     *  @param  one   One.
     *  @param  other The other.
     *  @return       Whether they are.
     */
    auto VaultSame(const Intersections &one, const Intersections &other) -> bool
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
     *  Two curves met by the swift kernels, within the length Intersect would use, and by the
     *  thorough way where they are sure: whether they were, and whether the two agree to the bit.
     *
     *  @param  curves    The two curves.
     *  @param  tolerance The tolerance asked for.
     *  @return           Whether the kernels were sure, whether they agreed, and whether they found
     *                    a crossing.
     */
    auto VaultBoth(const Curves &curves, const Tolerance &tolerance) -> std::array<bool, 3>
    {
      const std::optional<double> length = Reachable(curves, tolerance);
      Intersections               swiftly;
      if (!length
          || !Swift::Met(curves, {.Length = *length, .Radians = tolerance.Radians}, swiftly))
      {
        return {false, true, false};
      }
      Intersections thorough;
      Thorough(curves, tolerance, thorough);
      return {true, VaultSame(swiftly, thorough), !swiftly.Hits.empty()};
    }

    /**
     *  What the parabola kernel answers for two curves, within the length Intersect would use.
     *
     *  @param  first     One curve.
     *  @param  second    The other.
     *  @param  tolerance The tolerance.
     *  @return           Whether it was sure, nothing where it does not take them, and the hits.
     */
    auto Vaulting(const Curve &first, const Curve &second, const Tolerance &tolerance)
        -> std::pair<std::optional<bool>, Intersections>
    {
      Intersections             into;
      const Curves              curves = {first, second};
      const std::optional<bool> sure   = Swift::Vaulted(curves, Usable(curves, tolerance), into);
      return {sure, into};
    }

    /**
     *  A double from -1 to 1.
     *
     *  @param  draws Where it comes from.
     *  @return       The double.
     */
    auto VaultUnit(Draws &draws) -> double
    {
      return draws.Near(0.0, 1.0);
    }

    /**
     *  A whole number below a count.
     *
     *  @param  draws Where it comes from.
     *  @param  count The count.
     *  @return       The number.
     */
    auto VaultPick(Draws &draws, std::uint64_t count) -> std::uint64_t
    {
      return draws.Next() % count;
    }

    /**
     *  A spline the kernel takes: rational of degree two to eight, its weights from an eighth to
     *  eight times each other, or polynomial of degree four to eight.
     *
     *  @param  draws Where the numbers come from.
     *  @return       The spline.
     */
    auto VaultSplined(Draws &draws) -> Spline
    {
      const bool rated  = VaultPick(draws, 3) != 0;
      Spline     spline = rated ? Tests::Splined(draws, 2, 6) : Tests::Splined(draws, 4, 4);
      spline.Weights.clear();
      if (rated)
      {
        for (std::size_t index = 0; index < spline.Poles.size(); ++index)
        {
          spline.Weights.push_back(std::ldexp(1.0 + (0.5 * VaultUnit(draws)),
                                              static_cast<int>(VaultPick(draws, 7)) - 3));
        }
      }
      return spline;
    }

    /**
     *  A quadratic Bézier curve near a spline: at random, through one of its points at an angle,
     *  slight or not, along its tangent a few tolerances off, or through one of its ends.
     *
     *  @param  draws  Where the numbers come from.
     *  @param  spline The spline.
     *  @param  near   The tolerance's length.
     *  @return        The curve.
     */
    auto VaultArched(Draws &draws, const Spline &spline, double near) -> Curve
    {
      const auto kind = VaultPick(draws, 5);
      if (kind == 4)
      {
        return Bezier::Quadratic{.Start   = {.X = VaultUnit(draws), .Y = VaultUnit(draws)},
                                 .Control = {.X = VaultUnit(draws), .Y = VaultUnit(draws)},
                                 .End     = {.X = VaultUnit(draws), .Y = VaultUnit(draws)}};
      }
      const Range  domain = Domain(spline);
      const double end    = VaultPick(draws, 2) == 0 ? 0.0 : 1.0;
      const double share  = kind == 2 ? end : 0.5 + (0.45 * VaultUnit(draws));
      const double along  = domain.Start + ((domain.End - domain.Start) * share);
      const Point  spot   = At(spline, along);
      Vector       way    = Derivative(spline, along);
      if (!(Length(way) > 0.0))
      {
        way = {.X = 1.0, .Y = 0.0};
      }
      const double angle = kind == 0
                             ? std::ldexp(VaultUnit(draws), -static_cast<int>(VaultPick(draws, 30)))
                             : 3.0 * VaultUnit(draws);
      const Vector turned = {.X = (way.X * std::cos(angle)) - (way.Y * std::sin(angle)),
                             .Y = (way.X * std::sin(angle)) + (way.Y * std::cos(angle))};
      const Vector toward = kind == 1 ? way : turned;
      const double off    = kind == 1 ? near * 16.0 * VaultUnit(draws) : 0.0;
      const Point  moved  = spot + (Perpendicular(way) * (off / Length(way)));
      const Vector reach  = toward * ((0.3 + (0.5 * std::abs(VaultUnit(draws)))) / Length(toward));
      const Vector side   = Perpendicular(reach) * (0.3 * VaultUnit(draws) / Length(reach));
      const Vector ending = kind == 2 ? reach + side : Vector{};
      return Bezier::Quadratic{.Start   = moved - reach + side - ending,
                               .Control = moved - side - ending,
                               .End     = moved + reach + side - ending};
    }

    /**
     *  A curve scaled by a power of two and moved.
     *
     *  @param  curve A quadratic Bézier curve or a spline.
     *  @param  scale The scale.
     *  @param  by    How far it moves.
     *  @return       The curve, moved.
     */
    auto VaultMoved(const Curve &curve, double scale, Vector by) -> Curve
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
      Spline moved = std::get<Spline>(curve);
      for (Point &pole : moved.Poles)
      {
        pole = shifted(pole);
      }
      return moved;
    }

    /**
     *  A pair of curves and the tolerance asked for.
     */
    struct VaultPosed
    {
      Curve            First     = {}; /**< The first curve.         */
      Curve            Second    = {}; /**< The second.              */
      Geoxy::Tolerance Tolerance = {}; /**< The tolerance asked for. */
    };

    /**
     *  A pair drawn at random or at the kernel's edges, scaled and moved now and then, under a
     *  coarse tolerance at times, either way round.
     *
     *  @param  draws Where the numbers come from.
     *  @return       The pair.
     */
    auto VaultPosing(Draws &draws) -> VaultPosed
    {
      const Spline spline = VaultSplined(draws);
      const bool   coarse = VaultPick(draws, 4) == 0;
      const double near =
          coarse ? std::pow(10.0, -3.0 - (9.0 * std::abs(VaultUnit(draws)))) : 2e-10;
      const bool   wide    = VaultPick(draws, 6) == 0;
      const double radians = wide ? std::pow(10.0, -6.0 * std::abs(VaultUnit(draws))) : 1e-10;
      const Curve  arch    = VaultArched(draws, spline, near);
      const auto   moving  = VaultPick(draws, 6);
      const double scale   = std::ldexp(
          moving == 0 ? 1.37 : 1.0, moving == 0 ? static_cast<int>(VaultPick(draws, 80)) - 40 : 0);
      const Vector by      = {.X = moving == 1 ? std::ldexp(VaultUnit(draws), 19) : 0.0,
                              .Y = moving == 1 ? 1e6 : 0.0};
      const Curve  one     = VaultMoved(spline, scale, by);
      const Curve  other   = VaultMoved(arch, scale, by);
      const bool   swapped = VaultPick(draws, 2) == 0;
      return {.First     = swapped ? other : one,
              .Second    = swapped ? one : other,
              .Tolerance = {.Length = near * scale, .Radians = radians}};
    }

    /**
     *  Pairs as VaultPosing draws them: how many the kernels were sure of, how many they got wrong,
     *  and how many of those they were sure of cross.
     *
     *  @param  rounds How many pairs.
     *  @return        How many it was sure of, how many disagreed, how many crossing.
     */
    auto VaultDrawn(int rounds) -> std::array<int, 3>
    {
      Draws              draws{0xA4C4EDU};
      std::array<int, 3> count = {};
      for (int round = 0; round < rounds; ++round)
      {
        const VaultPosed posed = VaultPosing(draws);
        const auto [sure, agrees, crossed] =
            VaultBoth({posed.First, posed.Second}, posed.Tolerance);
        count[0] += static_cast<int>(sure);
        count[1] += static_cast<int>(!agrees);
        count[2] += static_cast<int>(crossed);
        if (!agrees)
        {
          INFO("round " << round);
          CHECK(agrees);
        }
      }
      return count;
    }

    /**
     *  A Bézier curve raised to degree four, the same curve: each round C_i = (i / (n + 1)) P_i-1 +
     *  (1 - i / (n + 1)) P_i.
     *
     *  @param  poles Its control points, of degree two or three.
     *  @return       Its five control points at degree four.
     */
    auto Raised(const std::vector<Point> &poles) -> std::vector<Point>
    {
      std::vector<Point> raised = poles;
      for (std::size_t degree = poles.size() - 1; degree < 4; ++degree)
      {
        std::vector<Point> next   = {raised.front()};
        const auto         higher = static_cast<double>(degree + 1);
        for (std::size_t index = 1; index <= degree; ++index)
        {
          const double share = static_cast<double>(index) / higher;
          next.push_back({.X = (share * raised[index - 1].X) + ((1.0 - share) * raised[index].X),
                          .Y = (share * raised[index - 1].Y) + ((1.0 - share) * raised[index].Y)});
        }
        next.push_back(raised.back());
        raised = next;
      }
      return raised;
    }

    /**
     *  y = 1 - x^2 + x (x^2 - 0.25) over x from a start to an end as a polynomial spline of degree
     *  four, one span from 0 to 1: the cubic's control ordinates from its values and slopes at the
     *  ends, raised once. It crosses y = 1 - x^2 at x = -0.5, 0 and 0.5.
     *
     *  @param  start Where x starts.
     *  @param  end   Where it ends.
     *  @return       The spline.
     */
    auto Graphed(double start, double end) -> Curve
    {
      const auto graph = [](double x) -> double
      {
        return 1.0 - (x * x) + (x * ((x * x) - 0.25));
      };
      const auto slope = [](double x) -> double
      {
        return (-2.0 * x) + (3.0 * x * x) - 0.25;
      };
      const double third = (end - start) / 3.0;
      return Spline{.Degree = 4,
                    .Knots  = {0, 0, 0, 0, 0, 1, 1, 1, 1, 1},
                    .Poles =
                        Raised({{.X = start, .Y = graph(start)},
                                {.X = start + third, .Y = graph(start) + (third * slope(start))},
                                {.X = end - third  , .Y = graph(end) - (third * slope(end))    },
                                {.X = end, .Y = graph(end)}}),
                    .Weights = {}};
    }
  } // namespace

  TEST_CASE("vault: quadratic Bézier curves with rational splines and splines past a cubic")
  {
    const int rounds                    = Rounds(20000);
    const auto [sure, differ, crossing] = VaultDrawn(rounds);
    MESSAGE("sure of " << sure << " of " << rounds << ", crossing in " << crossing);
    CHECK(differ == 0);
    CHECK(sure > rounds / 4);
  }

  TEST_CASE("vault: pairs the kernel does not take")
  {
    const Tolerance tolerance = {.Length = 2e-10, .Radians = 1e-10};
    const Curve     arch      = Bezier::Quadratic{
        .Start = {.X = -1, .Y = 0}, .Control = {.X = 0, .Y = 2}, .End = {.X = 1, .Y = 0}};
    const std::vector<Point> poles = {{.X = -1, .Y = 1}, {.X = 0, .Y = -1}, {.X = 1, .Y = 1}};
    const Curve              plain =
        Spline{.Degree = 2, .Knots = {0, 0, 0, 1, 1, 1}, .Poles = poles, .Weights = {}};
    const Curve rated =
        Spline{.Degree = 2, .Knots = {0, 0, 0, 1, 1, 1}, .Poles = poles, .Weights = {1, 2, 1}};
    const Curve cubic = Bezier::Cubic{.Start = {.X = -1, .Y = 1},
                                      .Out   = {.X = 0, .Y = -1},
                                      .In    = {.X = 0.5, .Y = -1},
                                      .End   = {.X = 1, .Y = 1}};
    CHECK(!Vaulting(arch, plain, tolerance).first);
    CHECK(!Vaulting(cubic, rated, tolerance).first);
    CHECK(!Vaulting(arch, arch, tolerance).first);
    CHECK(!Vaulting(rated, cubic, tolerance).first);
    const auto [taken, met] = Vaulting(rated, arch, tolerance);
    REQUIRE(taken);
    CHECK(*taken);
    CHECK(met.Hits.size() == 2);
    CHECK(VaultBoth({rated, arch}, tolerance) == std::array{true, true, true});
    /**
     *  A Bézier curve whose control points lie on a line has no parabola; one past the kernels'
     *  range, none it takes; a spline of degree nine is cut by no kernel.
     */
    const Curve flat = Bezier::Quadratic{
        .Start = {.X = -1, .Y = 0}, .Control = {.X = 0, .Y = 0}, .End = {.X = 1, .Y = 0}};
    const Curve wild = Bezier::Quadratic{
        .Start = {.X = -1, .Y = 0}, .Control = {.X = 0, .Y = 1e300}, .End = {.X = 1, .Y = 0}};
    std::vector<double> knots(10, 0.0);
    knots.resize(20, 1.0);
    std::vector<Point> nine;
    nine.reserve(10);
    for (int index = 0; index < 10; ++index)
    {
      nine.push_back({.X = (index / 4.5) - 1.0, .Y = index % 2 == 0 ? 1.0 : -1.0});
    }
    const Curve ninth =
        Spline{.Degree = 9, .Knots = knots, .Poles = nine, .Weights = std::vector<double>(10, 1.0)};
    for (const Curve &bezier : {flat, wild})
    {
      const auto [sure, hits] = Vaulting(bezier, rated, tolerance);
      REQUIRE(sure);
      CHECK(!*sure);
      CHECK(hits.Hits.empty());
    }
    const auto [ninthSure, ninthHits] = Vaulting(arch, ninth, tolerance);
    REQUIRE(ninthSure);
    CHECK(!*ninthSure);
    /**
     *  A rational span whose middle control point w (P - B) is 2^202, past the range, though every
     *  pole and weight is within it.
     */
    const double side = 0x1p41;
    const Curve  heavy =
        Spline{.Degree  = 2,
               .Knots   = {0, 0, 0, 1, 1, 1},
               .Poles   = {{.X = 0, .Y = 0}, {.X = side, .Y = side}, {.X = 0, .Y = 2 * side}},
               .Weights = {0x1p160, 0x1p161, 0x1p160}};
    const Curve across = Bezier::Quadratic{.Start   = {.X = -side, .Y = side},
                                           .Control = {.X = side, .Y = 2 * side},
                                           .End     = {.X = 2 * side, .Y = side}};
    const auto [heavySure, heavyHits] =
        Vaulting(heavy, across, {.Length = 0x1p-10, .Radians = 1e-10});
    REQUIRE(heavySure);
    CHECK(!*heavySure);
  }

  TEST_CASE("vault: crossings off the curve, at its end, shallow, tangent or coarse")
  {
    const Tolerance tolerance = {.Length = 2e-10, .Radians = 1e-10};
    /**
     *  y = 1 - x^2 from x = -1 to 1, its place s at x = 2 s - 1.
     */
    const Curve arch = Bezier::Quadratic{
        .Start = {.X = -1, .Y = 0}, .Control = {.X = 0, .Y = 2}, .End = {.X = 1, .Y = 0}};
    const auto line = [](Point from, Point onto) -> Curve
    {
      const Point middle = {.X = (from.X + onto.X) / 2.0, .Y = (from.Y + onto.Y) / 2.0};
      return Spline{.Degree  = 2,
                    .Knots   = {0, 0, 0, 1, 1, 1},
                    .Poles   = {from, middle, onto},
                    .Weights = {1, 2, 1}};
    };
    /**
     *  y = 1.75 - 2 x crosses the parabola at x = 0.5, on the curve, and at x = 1.5, off it.
     */
    const Curve off              = line({.X = 0, .Y = 1.75}, {.X = 2, .Y = -2.25});
    const auto [offSure, offMet] = Vaulting(arch, off, tolerance);
    REQUIRE(offSure);
    CHECK(*offSure);
    REQUIRE(offMet.Hits.size()       == 1);
    CHECK(offMet.Hits[0].At          == Point{.X = 0.5, .Y = 0.75});
    CHECK(offMet.Hits[0].First.Along == 0.75);
    CHECK(VaultBoth({off, arch}, tolerance) == std::array{true, true, true});
    /**
     *  y = x - 1 through the curve's end; y = 1 along its tangent at the top.
     */
    const Curve ending  = line({.X = 0.5, .Y = -0.5}, {.X = 1.5, .Y = 0.5});
    const Curve tangent = line({.X = -0.5, .Y = 1}, {.X = 0.5, .Y = 1});
    for (const Curve &spline : {ending, tangent})
    {
      const auto [sure, met] = Vaulting(spline, arch, tolerance);
      REQUIRE(sure);
      CHECK(!*sure);
      CHECK(met.Hits.empty());
    }
    /**
     *  y = 1 - x^2 + 0.2 (x - 0.5), raised to degree four, crosses at x = 0.5 at a sine near 0.11,
     *  shallow under an angle of 0.1.
     */
    const Curve shallow =
        Spline{.Degree  = 4,
               .Knots   = {0, 0, 0, 0, 0, 1, 1, 1, 1, 1},
               .Poles   = Raised({{.X = 0, .Y = 0.9}, {.X = 0.5, .Y = 1}, {.X = 1, .Y = 0.1}}),
               .Weights = {}};
    const auto [shallowSure, shallowMet] =
        Vaulting(arch, shallow, {.Length = 2e-10, .Radians = 0.1});
    REQUIRE(shallowSure);
    CHECK(!*shallowSure);
    CHECK(Vaulting(arch, shallow, tolerance).first == std::optional{true});
    /**
     *  The crossing on knots far from nothing, where a step of the spline's place is coarser than
     *  the tolerance allows.
     */
    const double far                   = 0x1p45;
    const Curve  coarse                = Spline{.Degree  = 2,
                                                .Knots   = {far, far, far, far + 1, far + 1, far + 1},
                                                .Poles   = std::get<Spline>(off).Poles,
                                                .Weights = {1, 2, 1}};
    const auto [coarseSure, coarseMet] = Vaulting(arch, coarse, tolerance);
    REQUIRE(coarseSure);
    CHECK(!*coarseSure);
  }

  TEST_CASE("vault: crossings the thorough way may take for one, or a place not settled")
  {
    const Curve arch = Bezier::Quadratic{
        .Start = {.X = -1, .Y = 0}, .Control = {.X = 0, .Y = 2}, .End = {.X = 1, .Y = 0}};
    const Tolerance tolerance = {.Length = 1e-9, .Radians = 1e-10};
    /**
     *  Three crossings, x evenly spread along the spline: over x from -0.8 to 0.9 the middle one is
     *  halfway between the outer two; from -0.7 to 0.8 the last is half the spline past halfway
     *  between the first two. A loop's double point on the parabola: two crossings a rounding
     *  apart. The graphs moved off the axes, where a coordinate all but zero would settle nothing.
     */
    const Curve loop = Spline{.Degree  = 4,
                              .Knots   = {0, 0, 0, 0, 0, 1, 1, 1, 1, 1},
                              .Poles   = Raised({{.X = -0.75, .Y = 0.125},
                                                 {.X =  2.25, .Y = 2.125},
                                                 {.X = -1.75, .Y = 2.125},
                                                 {.X = 1.25, .Y = 0.125}}),
                              .Weights = {}};
    double      low  = 0.1;
    double      high = 0.3;
    for (int step = 0; step < 60; ++step)
    {
      const double along                      = low + ((high - low) / 2.0);
      (At(loop, along).X < 0.25 ? low : high) = along;
    }
    const Point  twice   = At(loop, low);
    const Curve  crossed = Bezier::Quadratic{.Start   = {.X = twice.X - 0.5, .Y = twice.Y - 0.3},
                                             .Control = {.X = twice.X, .Y = twice.Y + 0.3},
                                             .End     = {.X = twice.X + 0.5, .Y = twice.Y - 0.3}};
    const Vector aside   = {.X = 0.375, .Y = 0.25};
    const Curve  moved   = VaultMoved(arch, 1.0, aside);
    const std::array<std::pair<Curve, Curve>, 3> pairs = {
        std::pair{moved, VaultMoved(Graphed(-0.8, 0.9), 1.0, aside)},
        std::pair{moved, VaultMoved(Graphed(-0.7, 0.8), 1.0, aside)}, std::pair{crossed, loop}};
    for (const auto &[bezier, spline] : pairs)
    {
      CHECK(Intersect(bezier, spline, tolerance).Hits.size() >= 2);
      const auto [sure, met] = Vaulting(bezier, spline, tolerance);
      REQUIRE(sure);
      CHECK(!*sure);
      CHECK(met.Hits.empty());
    }
    /**
     *  The comparison's case 3797: a rational cubic spline 1e6 off, whose crossing's place on the
     *  Bézier curve its point's bound does not settle.
     */
    const Curve far =
        Bezier::Quadratic{.Start   = {.X = -0x1.00001319bfe76p+20, .Y = -0x1.f8d87ee9d6e4bp+19},
                          .Control = {.X = -0x1.00000fae70ae3p+20, .Y = -0x1.f8d89fdc70168p+19},
                          .End     = {.X = -0x1.00000b45de8d6p+20, .Y = -0x1.f8d879e9579a2p+19}};
    const Curve afield =
        Spline{.Degree  = 3,
               .Knots   = {0, 0, 0, 0, 0x1.1133bc43f9cf6p-1, 0x1.2425b206f37c8p-1,
                           0x1.2425b206f37c8p-1, 1, 1, 1, 1},
               .Poles   = {{.X = -0x1.000007aa36859p+20, .Y = -0x1.f8d8a2045f23dp+19},
                           {.X = -0x1.ffffef4b2d38ap+19, .Y = -0x1.f8d89c06cf8a5p+19},
                           {.X = -0x1.0000136f7a15cp+20, .Y = -0x1.f8d88bca236d1p+19},
                           {.X = -0x1.fffff4491914cp+19, .Y = -0x1.f8d875e4d6c82p+19},
                           {.X = -0x1.000000fa658f2p+20, .Y = -0x1.f8d8a6796bb9fp+19},
                           {.X = -0x1.0000092b6a17ap+20, .Y = -0x1.f8d87dbc3a1d7p+19},
                           {.X = -0x1.00000cd87ac53p+20, .Y = -0x1.f8d88371d491fp+19}},
               .Weights = {0x1.c4f656385c06ap-1, 0x1.d2da721e13df2p+0, 0x1.a559f8f24923ap+0,
                           0x1.67440f927bda1p-1, 0x1.05b27cc5b3416p+1, 0x1.7abb892c9f1e7p+1,
                           0x1.49bd4d508a1ecp-2}};
    const auto [farSure, farMet] = Vaulting(far, afield, {.Length = 0x1p-27, .Radians = 0x1p-28});
    REQUIRE(farSure);
    CHECK(!*farSure);
  }
} // namespace Geoxy
