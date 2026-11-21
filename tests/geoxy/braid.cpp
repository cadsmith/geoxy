/***************************************************************************************************
 *  @file      tests/geoxy/braid.cpp
 *  @brief     Tests for the cubic's kernel: a cubic Bézier curve met with a rational spline or one
 *             past a cubic on the curve's implicit cubic, each crossing the double nearest the
 *             exact one, and every answer the thorough way's to the bit or none.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "swift/braid.hpp"

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
    auto BraidSame(const Intersections &one, const Intersections &other) -> bool
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
    auto BraidBoth(const Curves &curves, const Tolerance &tolerance) -> std::array<bool, 3>
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
      return {true, BraidSame(swiftly, thorough), !swiftly.Hits.empty()};
    }

    /**
     *  What the cubic's kernel answers for two curves, within the length Intersect would use.
     *
     *  @param  first     One curve.
     *  @param  second    The other.
     *  @param  tolerance The tolerance.
     *  @return           Whether it was sure, nothing where it does not take them, and the hits.
     */
    auto Braiding(const Curve &first, const Curve &second, const Tolerance &tolerance)
        -> std::pair<std::optional<bool>, Intersections>
    {
      Intersections             into;
      const Curves              curves = {first, second};
      const std::optional<bool> sure   = Swift::Braided(curves, Usable(curves, tolerance), into);
      return {sure, into};
    }

    /**
     *  A double from -1 to 1.
     *
     *  @param  draws Where it comes from.
     *  @return       The double.
     */
    auto BraidUnit(Draws &draws) -> double
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
    auto BraidPick(Draws &draws, std::uint64_t count) -> std::uint64_t
    {
      return draws.Next() % count;
    }

    /**
     *  A spline the kernel takes: rational of degree two to five, its weights from an eighth to
     *  eight times each other, or polynomial of degree four or five.
     *
     *  @param  draws Where the numbers come from.
     *  @return       The spline.
     */
    auto BraidSplined(Draws &draws) -> Spline
    {
      const bool rated  = BraidPick(draws, 3) != 0;
      Spline     spline = rated ? Tests::Splined(draws, 2, 3) : Tests::Splined(draws, 4, 1);
      spline.Weights.clear();
      if (rated)
      {
        for (std::size_t index = 0; index < spline.Poles.size(); ++index)
        {
          spline.Weights.push_back(std::ldexp(1.0 + (0.5 * BraidUnit(draws)),
                                              static_cast<int>(BraidPick(draws, 7)) - 3));
        }
      }
      return spline;
    }

    /**
     *  A cubic Bézier curve near a spline: at random, through one of its points at an angle, slight
     *  or not, along its tangent a few tolerances off, or through one of its ends.
     *
     *  @param  draws  Where the numbers come from.
     *  @param  spline The spline.
     *  @param  near   The tolerance's length.
     *  @return        The curve.
     */
    auto BraidTwisted(Draws &draws, const Spline &spline, double near) -> Curve
    {
      const auto kind = BraidPick(draws, 5);
      if (kind == 4)
      {
        return Bezier::Cubic{.Start = {.X = BraidUnit(draws), .Y = BraidUnit(draws)},
                             .Out   = {.X = BraidUnit(draws), .Y = BraidUnit(draws)},
                             .In    = {.X = BraidUnit(draws), .Y = BraidUnit(draws)},
                             .End   = {.X = BraidUnit(draws), .Y = BraidUnit(draws)}};
      }
      const Range  domain = Domain(spline);
      const double end    = BraidPick(draws, 2) == 0 ? 0.0 : 1.0;
      const double share  = kind == 2 ? end : 0.5 + (0.45 * BraidUnit(draws));
      const double along  = domain.Start + ((domain.End - domain.Start) * share);
      const Point  spot   = At(spline, along);
      Vector       way    = Derivative(spline, along);
      if (!(Length(way) > 0.0))
      {
        way = {.X = 1.0, .Y = 0.0};
      }
      const double angle = kind == 0
                             ? std::ldexp(BraidUnit(draws), -static_cast<int>(BraidPick(draws, 30)))
                             : 3.0 * BraidUnit(draws);
      const Vector turned = {.X = (way.X * std::cos(angle)) - (way.Y * std::sin(angle)),
                             .Y = (way.X * std::sin(angle)) + (way.Y * std::cos(angle))};
      const Vector toward = kind == 1 ? way : turned;
      const double off    = kind == 1 ? near * 16.0 * BraidUnit(draws) : 0.0;
      const Point  moved  = spot + (Perpendicular(way) * (off / Length(way)));
      const Vector reach  = toward * ((0.3 + (0.5 * std::abs(BraidUnit(draws)))) / Length(toward));
      const Vector side   = Perpendicular(reach) * (0.3 * BraidUnit(draws) / Length(reach));
      const Vector bend   = Perpendicular(reach) * (0.6 * BraidUnit(draws) / Length(reach));
      const Vector ending = kind == 2 ? reach + side : Vector{};
      return Bezier::Cubic{.Start = moved - reach + side - ending,
                           .Out   = moved - (reach * (1.0 / 3.0)) - side + bend - ending,
                           .In    = moved + (reach * (1.0 / 3.0)) - side - bend - ending,
                           .End   = moved + reach + side - ending};
    }

    /**
     *  A curve scaled by a power of two and moved.
     *
     *  @param  curve A cubic Bézier curve or a spline.
     *  @param  scale The scale.
     *  @param  by    How far it moves.
     *  @return       The curve, moved.
     */
    auto BraidMoved(const Curve &curve, double scale, Vector by) -> Curve
    {
      const auto shifted = [scale, by](Point point) -> Point
      {
        return {.X = (point.X * scale) + by.X, .Y = (point.Y * scale) + by.Y};
      };
      if (const auto *cubic = std::get_if<Bezier::Cubic>(&curve))
      {
        return Bezier::Cubic{.Start = shifted(cubic->Start),
                             .Out   = shifted(cubic->Out),
                             .In    = shifted(cubic->In),
                             .End   = shifted(cubic->End)};
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
    struct BraidPosed
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
    auto BraidPosing(Draws &draws) -> BraidPosed
    {
      const Spline spline = BraidSplined(draws);
      const bool   coarse = BraidPick(draws, 4) == 0;
      const double near =
          coarse ? std::pow(10.0, -3.0 - (9.0 * std::abs(BraidUnit(draws)))) : 2e-10;
      const bool   wide    = BraidPick(draws, 6) == 0;
      const double radians = wide ? std::pow(10.0, -6.0 * std::abs(BraidUnit(draws))) : 1e-10;
      const Curve  twisted = BraidTwisted(draws, spline, near);
      const auto   moving  = BraidPick(draws, 6);
      const double scale   = std::ldexp(
          moving == 0 ? 1.37 : 1.0, moving == 0 ? static_cast<int>(BraidPick(draws, 80)) - 40 : 0);
      const Vector by      = {.X = moving == 1 ? std::ldexp(BraidUnit(draws), 19) : 0.0,
                              .Y = moving == 1 ? 1e6 : 0.0};
      const Curve  one     = BraidMoved(spline, scale, by);
      const Curve  other   = BraidMoved(twisted, scale, by);
      const bool   swapped = BraidPick(draws, 2) == 0;
      return {.First     = swapped ? other : one,
              .Second    = swapped ? one : other,
              .Tolerance = {.Length = near * scale, .Radians = radians}};
    }

    /**
     *  Pairs as BraidPosing draws them: how many the kernels were sure of, how many they got wrong,
     *  and how many of those they were sure of cross.
     *
     *  @param  rounds How many pairs.
     *  @return        How many it was sure of, how many disagreed, how many crossing.
     */
    auto BraidDrawn(int rounds) -> std::array<int, 3>
    {
      Draws              draws{0xB4A1DU};
      std::array<int, 3> count = {};
      for (int round = 0; round < rounds; ++round)
      {
        const BraidPosed posed = BraidPosing(draws);
        const auto [sure, agrees, crossed] =
            BraidBoth({posed.First, posed.Second}, posed.Tolerance);
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
     *  A graph y = g(x) of a cubic g over x from a start to an end as a polynomial spline of degree
     *  four, one span from 0 to 1, x even along it: the cubic's control ordinates from its values
     *  and slopes at the ends, raised once, C_i = (i / 4) P_i-1 + (1 - i / 4) P_i.
     *
     *  @tparam GraphT The cubic's kind.
     *  @param  graph  The cubic, its value and its slope at x.
     *  @param  start  Where x starts.
     *  @param  end    Where it ends.
     *  @return        The spline.
     */
    template <typename GraphT>
    auto Charted(const GraphT &graph, double start, double end) -> Curve
    {
      const double                third = (end - start) / 3.0;
      const std::array<double, 2> first = graph(start);
      const std::array<double, 2> last  = graph(end);
      const std::array<Point, 4>  cubic = {
          Point{.X = start        , .Y = first[0]                     },
          Point{.X = start + third, .Y = first[0] + (third * first[1])},
          Point{.X = end - third, .Y = last[0] - (third * last[1])}, Point{.X = end, .Y = last[0]}};
      std::vector<Point> raised = {cubic[0]};
      for (std::size_t index = 1; index <= 3; ++index)
      {
        const double share = static_cast<double>(index) / 4.0;
        raised.push_back(
            {.X = (share * cubic.at(index - 1).X) + ((1.0 - share) * cubic.at(index).X),
             .Y = (share * cubic.at(index - 1).Y) + ((1.0 - share) * cubic.at(index).Y)});
      }
      raised.push_back(cubic[3]);
      return Spline{
          .Degree = 4, .Knots = {0, 0, 0, 0, 0, 1, 1, 1, 1, 1}, .Poles = raised, .Weights = {}};
    }

    /**
     *  A straight rational spline of degree two from one point to another, its middle pole halfway
     *  and of weight two.
     *
     *  @param  from Where it starts.
     *  @param  onto Where it ends.
     *  @return      The spline.
     */
    auto Ruled(Point from, Point onto) -> Curve
    {
      const Point middle = {.X = (from.X + onto.X) / 2.0, .Y = (from.Y + onto.Y) / 2.0};
      return Spline{.Degree  = 2,
                    .Knots   = {0, 0, 0, 1, 1, 1},
                    .Poles   = {from, middle, onto},
                    .Weights = {1, 2, 1}};
    }

    /**
     *  y = x^3 - 2.25 x from x = -1.5 to 1.5 as a cubic Bézier curve, its place s at x = 3 s - 1.5.
     *
     *  @return The curve.
     */
    auto Waved() -> Curve
    {
      return Bezier::Cubic{.Start = {.X = -1.5, .Y = 0},
                           .Out   = {.X = -0.5, .Y = 4.5},
                           .In    = {.X = 0.5, .Y = -4.5},
                           .End   = {.X = 1.5, .Y = 0}};
    }
  } // namespace

  TEST_CASE("braid: cubic Bézier curves with rational splines and splines past a cubic")
  {
    const int rounds                    = Rounds(20000);
    const auto [sure, differ, crossing] = BraidDrawn(rounds);
    MESSAGE("sure of " << sure << " of " << rounds << ", crossing in " << crossing);
    CHECK(differ == 0);
    CHECK(sure > rounds / 4);
  }

  TEST_CASE("braid: pairs the kernel does not take")
  {
    const Curve              wave      = Waved();
    const Tolerance          tolerance = {.Length = 2e-10, .Radians = 1e-10};
    const std::vector<Point> poles     = {{.X = -1, .Y = 1}, {.X = 0, .Y = -1}, {.X = 1, .Y = 1}};
    const Curve              plain =
        Spline{.Degree = 2, .Knots = {0, 0, 0, 1, 1, 1}, .Poles = poles, .Weights = {}};
    const Curve rated =
        Spline{.Degree = 2, .Knots = {0, 0, 0, 1, 1, 1}, .Poles = poles, .Weights = {1, 2, 1}};
    const Curve arch = Bezier::Quadratic{
        .Start = {.X = -1, .Y = 0}, .Control = {.X = 0, .Y = 2}, .End = {.X = 1, .Y = 0}};
    CHECK(!Braiding(wave, plain, tolerance).first);
    CHECK(!Braiding(arch, rated, tolerance).first);
    CHECK(!Braiding(wave, wave, tolerance).first);
    CHECK(!Braiding(rated, arch, tolerance).first);
    const auto [taken, met] = Braiding(rated, wave, tolerance);
    REQUIRE(taken);
    CHECK(*taken);
    CHECK(met.Hits.size() == 1);
    CHECK(BraidBoth({rated, wave}, tolerance) == std::array{true, true, true});
    /**
     *  A curve past the kernels' range; a spline of degree six, whose equation would be past the
     *  isolation's sixteenth degree; a spline with a weight too small to cut it by.
     */
    const Curve         wild = Bezier::Cubic{.Start = {.X = -1, .Y = 0},
                                             .Out   = {.X = 0, .Y = 1e300},
                                             .In    = {.X = 0.5, .Y = -1},
                                             .End   = {.X = 1, .Y = 0}};
    std::vector<double> knots(7, 0.0);
    knots.resize(14, 1.0);
    std::vector<Point> six;
    six.reserve(7);
    for (int index = 0; index < 7; ++index)
    {
      six.push_back({.X = (index / 3.0) - 1.0, .Y = index % 2 == 0 ? 1.0 : -1.0});
    }
    const Curve sixth =
        Spline{.Degree = 6, .Knots = knots, .Poles = six, .Weights = std::vector<double>(7, 1.0)};
    const Curve faint = Spline{
        .Degree = 2, .Knots = {0, 0, 0, 1, 1, 1}, .Poles = poles, .Weights = {1, 0x1p-210, 1}};
    for (const auto &[bezier, spline] :
         {std::pair{wild, rated}, std::pair{wave, sixth}, std::pair{wave, faint}})
    {
      const auto [sure, hits] = Braiding(bezier, spline, tolerance);
      REQUIRE(sure);
      CHECK(!*sure);
      CHECK(hits.Hits.empty());
    }
    /**
     *  A rational span whose middle control point w (P - B) is 2^202, past the range, though every
     *  pole and weight is within it; one whose entries would be past 2^300, their products past the
     *  doubles.
     */
    const double side = 0x1p41;
    const Curve  heavy =
        Spline{.Degree  = 2,
               .Knots   = {0, 0, 0, 1, 1, 1},
               .Poles   = {{.X = 0, .Y = 0}, {.X = side, .Y = side}, {.X = 0, .Y = 2 * side}},
               .Weights = {0x1p160, 0x1p161, 0x1p160}};
    const Curve across = Bezier::Cubic{.Start = {.X = -side, .Y = side},
                                       .Out   = {.X = 0, .Y = 2 * side},
                                       .In    = {.X = side, .Y = 2 * side},
                                       .End   = {.X = 2 * side, .Y = side}};
    const auto [heavySure, heavyHits] =
        Braiding(heavy, across, {.Length = 0x1p-10, .Radians = 1e-10});
    REQUIRE(heavySure);
    CHECK(!*heavySure);
    /**
     *  A curve 2^160 across through the origin, whose entries on a small spline there are past
     *  2^300, their products past the doubles.
     */
    const double vast = 0x1p160;
    const Curve  huge = BraidMoved(wave, vast, {.X = -0.75 * vast, .Y = 1.265625 * vast});
    const auto [hugeSure, hugeHits] = Braiding(huge, rated, tolerance);
    REQUIRE(hugeSure);
    CHECK(!*hugeSure);
  }

  TEST_CASE("braid: crossings off the curve, at its end, shallow, tangent or coarse")
  {
    const Curve     wave      = Waved();
    const Tolerance tolerance = {.Length = 2e-10, .Radians = 1e-10};
    /**
     *  y = 5.0625 (x - 1) crosses the cubic at x = 0.75, on the curve, and at x = 2.25, off it.
     */
    const Curve off              = Ruled({.X = 0.5, .Y = -2.53125}, {.X = 2.5, .Y = 7.59375});
    const auto [offSure, offMet] = Braiding(wave, off, tolerance);
    REQUIRE(offSure);
    CHECK(*offSure);
    REQUIRE(offMet.Hits.size()       == 1);
    CHECK(offMet.Hits[0].At          == Point{.X = 0.75, .Y = -1.265625});
    CHECK(offMet.Hits[0].First.Along == 0.75);
    CHECK(BraidBoth({off, wave}, tolerance) == std::array{true, true, true});
    /**
     *  y = x - 1.5 through the curve's end; y = -1.5 sqrt(0.75), rounded, along its tangent at the
     *  bottom of its trough.
     */
    const double low     = -1.5 * std::sqrt(0.75);
    const Curve  ending  = Ruled({.X = 1, .Y = -0.5}, {.X = 2, .Y = 0.5});
    const Curve  tangent = Ruled({.X = 0.4, .Y = low}, {.X = 1.3, .Y = low});
    for (const Curve &spline : {ending, tangent})
    {
      const auto [sure, met] = Braiding(spline, wave, tolerance);
      REQUIRE(sure);
      CHECK(!*sure);
      CHECK(met.Hits.empty());
    }
    /**
     *  y = x^3 - 2.15 x - 0.075 crosses it at x = 0.75 at a sine near 0.08, shallow under an angle
     *  of 0.1.
     */
    const Curve shallow =
        Charted([](double x) -> std::array<double, 2>
                { return {(x * x * x) - (2.15 * x) - 0.075, (3.0 * x * x) - 2.15}; }, 0.0, 1.4);
    const auto [shallowSure, shallowMet] =
        Braiding(wave, shallow, {.Length = 2e-10, .Radians = 0.1});
    REQUIRE(shallowSure);
    CHECK(!*shallowSure);
    CHECK(Braiding(wave, shallow, tolerance).first == std::optional{true});
    /**
     *  The crossing on knots far from nothing, where a step of the spline's place is coarser than
     *  the tolerance allows.
     */
    const double far                   = 0x1p45;
    const Curve  coarse                = Spline{.Degree  = 2,
                                                .Knots   = {far, far, far, far + 1, far + 1, far + 1},
                                                .Poles   = std::get<Spline>(off).Poles,
                                                .Weights = {1, 2, 1}};
    const auto [coarseSure, coarseMet] = Braiding(wave, coarse, tolerance);
    REQUIRE(coarseSure);
    CHECK(!*coarseSure);
  }

  TEST_CASE("braid: crossings the thorough way may take for one, near the double point, or a place "
            "not settled")
  {
    const Curve     wave      = Waved();
    const Tolerance tolerance = {.Length = 1e-9, .Radians = 1e-10};
    /**
     *  y = 0.5 x^3 - 2.125 x crosses the cubic at x = -0.5, 0 and 0.5, x even along the spline:
     *  over x from -0.8 to 0.9 the middle one is halfway between the outer two; from -0.7 to 0.8
     *  the last is half the spline past halfway between the first two. Moved off the axes, where a
     *  coordinate all but zero would settle nothing.
     */
    const auto graph = [](double x) -> std::array<double, 2>
    {
      return {(0.5 * x * x * x) - (2.125 * x), (1.5 * x * x) - 2.125};
    };
    const Vector aside = {.X = 0.375, .Y = 0.25};
    const Curve  moved = BraidMoved(wave, 1.0, aside);
    for (const auto &[start, end] : {std::pair{-0.8, 0.9}, std::pair{-0.7, 0.8}})
    {
      const Curve spline = BraidMoved(Charted(graph, start, end), 1.0, aside);
      CHECK(Intersect(moved, spline, tolerance).Hits.size() >= 2);
      const auto [sure, met] = Braiding(moved, spline, tolerance);
      REQUIRE(sure);
      CHECK(!*sure);
      CHECK(met.Hits.empty());
    }
    /**
     *  A cubic Bézier curve that loops, x(s) odd about s = 1/2: a line through its double point,
     *  where its place is two.
     */
    const Curve loop  = Bezier::Cubic{.Start = {.X = -1, .Y = 0},
                                      .Out   = {.X = 2, .Y = 2},
                                      .In    = {.X = -2, .Y = 2},
                                      .End   = {.X = 1, .Y = 0}};
    double      lower = 0.0;
    double      upper = 0.45;
    for (int step = 0; step < 60; ++step)
    {
      const double along                        = lower + ((upper - lower) / 2.0);
      (At(loop, along).X < 0.0 ? lower : upper) = along;
    }
    const double twice   = At(loop, lower).Y;
    const Curve  through = Ruled({.X = -0.5, .Y = twice}, {.X = 0.5, .Y = twice});
    CHECK(Intersect(loop, through, tolerance).Hits.size() >= 1);
    const auto [loopSure, loopMet] = Braiding(loop, through, tolerance);
    REQUIRE(loopSure);
    CHECK(!*loopSure);
    /**
     *  A quadratic curve raised to a cubic has no implicit cubic, its Bézout form nothing.
     */
    const Curve raised = Bezier::Cubic{.Start = {.X = -1.5, .Y = 0},
                                       .Out   = {.X = -0.5, .Y = 2},
                                       .In    = {.X = 0.5, .Y = 2},
                                       .End   = {.X = 1.5, .Y = 0}};
    const auto [raisedSure, raisedMet] =
        Braiding(raised, Ruled({.X = -1.5, .Y = 0.5}, {.X = 1.5, .Y = 0.75}), tolerance);
    REQUIRE(raisedSure);
    CHECK(!*raisedSure);
    /**
     *  The comparison's case 4184: a rational cubic spline 1e6 off, whose crossing's place on the
     *  Bézier curve its point's bound does not settle.
     */
    const Curve far =
        Bezier::Cubic{.Start = {.X = 0x1.000004e79ae43p+20, .Y = 0x1.f0dfe8a3287aap+18},
                      .Out   = {.X = 0x1.0000081deb168p+20, .Y = 0x1.f0dff6a4c7835p+18},
                      .In    = {.X = 0x1.fffffa0bc2643p+19, .Y = 0x1.f0e01862fe19fp+18},
                      .End   = {.X = 0x1.00000d5b5faa1p+20, .Y = 0x1.f0dfd306f8c28p+18}};
    const Curve afield =
        Spline{.Degree  = 3,
               .Knots   = {0, 0, 0, 0, 0x1.4026195ffc6d0p-2, 0x1.0a3cf042cac12p-1, 1, 1, 1, 1},
               .Poles   = {{.X = 0x1.0000015e83380p+20, .Y = 0x1.f0dff4939190cp+18},
                           {.X = 0x1.00000cae2132bp+20, .Y = 0x1.f0dfdd318931cp+18},
                           {.X = 0x1.0000102884c21p+20, .Y = 0x1.f0dffb772beb7p+18},
                           {.X = 0x1.0000093400368p+20, .Y = 0x1.f0e02de9ee951p+18},
                           {.X = 0x1.0000157da64e7p+20, .Y = 0x1.f0dfd91d6e290p+18},
                           {.X = 0x1.fffff3daa94c2p+19, .Y = 0x1.f0dff960a3ebdp+18}},
               .Weights = {0x1.b01dbcb0831c0p+0, 0x1.7c03f89a9a27cp+0, 0x1.c2ef4d1076436p+0,
                           0x1.a419361005787p+0, 0x1.7f44970aa98b5p+1, 0x1.66aab8fa73ccap-1}};
    const auto [farSure, farMet] = Braiding(far, afield, {.Length = 0x1p-27, .Radians = 0x1p-28});
    REQUIRE(farSure);
    CHECK(!*farSure);
  }
} // namespace Geoxy
