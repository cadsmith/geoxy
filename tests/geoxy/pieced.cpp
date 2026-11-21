/***************************************************************************************************
 *  @file      tests/geoxy/pieced.cpp
 *  @brief     Tests for polynomial splines met piece by piece by the swift kernels: every answer
 *             the thorough way's to the bit, or left to it.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "swift/pieced.hpp"

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
#include "geoxy/curve.hpp"
#include "geoxy/derivative.hpp"
#include "geoxy/intersect.hpp"
#include "geoxy/line.hpp"
#include "geoxy/point.hpp"
#include "geoxy/ray.hpp"
#include "geoxy/spline.hpp"
#include "geoxy/tolerance.hpp"
#include "geoxy/vector.hpp"
#include "geoxy/xline.hpp"

#include <doctest/doctest.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <numbers>
#include <optional>
#include <type_traits>
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
    auto PiecedSame(const Intersections &one, const Intersections &other) -> bool
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
    struct PiecedMet
    {
      bool Sure   = false;         /**< Whether the kernels were sure. */
      bool Agrees = true;          /**< Whether they found what the thorough
                                        way finds, to the bit.               */
      Intersections Swiftly  = {}; /**< What the kernels found.      */
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
    auto PiecedBoth(const Curves &curves, const Tolerance &tolerance) -> PiecedMet
    {
      PiecedMet                   met;
      const std::optional<double> length = Reachable(curves, tolerance);
      if (!length)
      {
        return met;
      }
      met.Sure = Swift::Met(curves, {.Length = *length, .Radians = tolerance.Radians}, met.Swiftly);
      if (met.Sure)
      {
        Thorough(curves, tolerance, met.Thorough);
        met.Agrees = PiecedSame(met.Swiftly, met.Thorough);
      }
      return met;
    }

    /**
     *  A number from -1 to 1.
     *
     *  @param  draws The draws.
     *  @return       The number.
     */
    auto PiecedUnit(Draws &draws) -> double
    {
      return draws.Near(0.0, 1.0);
    }

    /**
     *  A whole number below a count.
     *
     *  @param  draws The draws.
     *  @param  count The count.
     *  @return       The number.
     */
    auto PiecedPick(Draws &draws, std::uint64_t count) -> std::uint64_t
    {
      return draws.Next() % count;
    }

    /**
     *  A quadratic or cubic Bézier curve of control points anywhere.
     *
     *  @param  draws The draws.
     *  @return       The curve.
     */
    auto PiecedBent(Draws &draws) -> Curve
    {
      const auto point = [&draws]() -> Point
      {
        return {.X = PiecedUnit(draws), .Y = PiecedUnit(draws)};
      };
      if (PiecedPick(draws, 2) == 0)
      {
        const Point start = point();
        const Point bend  = point();
        return Bezier::Quadratic{.Start = start, .Control = bend, .End = point()};
      }
      const Point start  = point();
      const Point out    = point();
      const Point inward = point();
      return Bezier::Cubic{.Start = start, .Out = out, .In = inward, .End = point()};
    }

    /**
     *  A straight curve anywhere: a line, a ray or an endless line.
     *
     *  @param  draws The draws.
     *  @return       The curve.
     */
    auto PiecedStraight(Draws &draws) -> Curve
    {
      const Point  start = {.X = PiecedUnit(draws), .Y = PiecedUnit(draws)};
      const Vector way   = {.X = PiecedUnit(draws), .Y = PiecedUnit(draws)};
      switch (PiecedPick(draws, 3))
      {
      case 0:
        return Line{.Start = start, .End = start + (way * 2.0)};
      case 1:
        return Ray{.Start = start, .Way = way};
      default:
        return Xline{.Through = start, .Way = way};
      }
    }

    /**
     *  A curve scaled by a power of two and moved: every point as the curve's, so the same curve
     *  where both are exact.
     *
     *  @param  curve The curve.
     *  @param  scale The scale, a power of two or near it.
     *  @param  by    How far it is moved.
     *  @return       The curve moved.
     */
    auto PiecedMoved(const Curve &curve, double scale, Vector by) -> Curve
    {
      const auto shifted = [scale, by](Point point) -> Point
      {
        return {.X = (point.X * scale) + by.X, .Y = (point.Y * scale) + by.Y};
      };
      return std::visit(
          [&shifted, scale](auto shape) -> Curve
          {
            using Shape = std::decay_t<decltype(shape)>;
            if constexpr (std::is_same_v<Shape, Spline>)
            {
              for (Point &pole : shape.Poles)
              {
                pole = shifted(pole);
              }
            }
            else if constexpr (std::is_same_v<Shape, Bezier::Quadratic>)
            {
              shape = {.Start   = shifted(shape.Start),
                       .Control = shifted(shape.Control),
                       .End     = shifted(shape.End)};
            }
            else if constexpr (std::is_same_v<Shape, Bezier::Cubic>)
            {
              shape = {.Start = shifted(shape.Start),
                       .Out   = shifted(shape.Out),
                       .In    = shifted(shape.In),
                       .End   = shifted(shape.End)};
            }
            else if constexpr (std::is_same_v<Shape, Line>)
            {
              shape = {.Start = shifted(shape.Start), .End = shifted(shape.End)};
            }
            else if constexpr (std::is_same_v<Shape, Ray>)
            {
              shape = {.Start = shifted(shape.Start), .Way = shape.Way * scale};
            }
            else if constexpr (std::is_same_v<Shape, Xline>)
            {
              shape = {.Through = shifted(shape.Through), .Way = shape.Way * scale};
            }
            return shape;
          },
          curve);
    }

    /**
     *  A spline's poles moved by a little, square to its chord: a spline lying along it, a few
     *  tolerances off or within one.
     *
     *  @param  spline The spline.
     *  @param  by     How far.
     *  @return        The spline moved.
     */
    auto Beside(Spline spline, double by) -> Spline
    {
      const Vector chord  = spline.Poles.back() - spline.Poles.front();
      const double length = Length(chord);
      const Vector across = length > 0.0 ? Vector{.X = -chord.Y / length, .Y = chord.X / length}
                                         : Vector{.X = 0.0, .Y = 1.0};
      for (Point &pole : spline.Poles)
      {
        pole = pole + (across * by);
      }
      return spline;
    }

    /**
     *  Two curves and a tolerance, posed.
     */
    struct PiecedPosed
    {
      Curve            First     = Line{}; /**< The first curve. */
      Curve            Second    = Line{}; /**< The second.      */
      Geoxy::Tolerance Tolerance = {};     /**< The tolerance.   */
    };

    /**
     *  A random pair, a spline at least: with another spline, a Bézier curve or a straight curve,
     *  lying along it now and then, either way round, scaled and moved now and then.
     *
     *  @param  draws The draws.
     *  @return       The pair.
     */
    auto PiecedPosing(Draws &draws) -> PiecedPosed
    {
      const Spline one   = Tests::Splined(draws, 2, 1);
      const auto   kind  = PiecedPick(draws, 8);
      const double near  = PiecedPick(draws, 4) == 0
                             ? std::pow(10.0, -3.0 - (9.0 * std::abs(PiecedUnit(draws))))
                             : 2e-10;
      Curve        other = Tests::Splined(draws, 2, 1);
      if (kind == 0)
      {
        other = Beside(one, near * std::pow(10.0, 3.0 * PiecedUnit(draws)));
      }
      else if (kind <= 2)
      {
        other = PiecedBent(draws);
      }
      else if (kind <= 4)
      {
        other = PiecedStraight(draws);
      }
      const auto   moving = PiecedPick(draws, 6);
      const double scale  = std::ldexp(
          moving == 0 ? 1.37 : 1.0, moving == 0 ? static_cast<int>(PiecedPick(draws, 80)) - 40 : 0);
      const Vector by      = {.X = moving == 1 ? std::ldexp(PiecedUnit(draws), 26) : 0.0,
                              .Y = moving == 1 ? 1e8 : 0.0};
      const bool   swapped = PiecedPick(draws, 2) == 0;
      const Curve  first   = PiecedMoved(one, scale, by);
      const Curve  second  = PiecedMoved(other, scale, by);
      return {.First     = swapped ? second : first,
              .Second    = swapped ? first : second,
              .Tolerance = {.Length = near * scale, .Radians = 1e-10}};
    }
  } // namespace

  TEST_CASE("pieced: splines as met thoroughly, to the bit")
  {
    const int          rounds = Rounds(40000);
    Draws              draws{0x5B11E5U};
    std::array<int, 3> count = {};
    for (int round = 0; round < rounds; ++round)
    {
      const PiecedPosed posed  = PiecedPosing(draws);
      const Curves      curves = {posed.First, posed.Second};
      for (const Curve *curve : {&posed.First, &posed.Second})
      {
        const auto *spline = std::get_if<Spline>(curve);
        REQUIRE((spline == nullptr || Valid(*spline)));
      }
      const PiecedMet met  = PiecedBoth(curves, posed.Tolerance);
      count[0]            += static_cast<int>(met.Sure);
      count[1]            += static_cast<int>(!met.Agrees);
      count[2]            += static_cast<int>(met.Sure && !met.Swiftly.Hits.empty());
      INFO("round " << round);
      CHECK(met.Agrees);
    }
    MESSAGE("sure of "
            << count[0]
            << " of "
            << rounds
            << ", crossing in "
            << count[2]
            << ", unlike the thorough way in "
            << count[1]);
    CHECK(count[0] > rounds / 3);
    CHECK(count[2] > rounds / 8);
  }
} // namespace Geoxy

namespace Geoxy
{
  namespace
  {
    /**
     *  A cubic spline waving from left to right across eight poles, its knots clamped and even.
     *
     *  @return The spline.
     */
    auto Waved() -> Spline
    {
      return {.Degree  = 3,
              .Knots   = {0, 0, 0, 0, 1, 2, 3, 4, 5, 5, 5, 5},
              .Poles   = {{.X = 0.0, .Y = 0.0},
                          {.X = 1.0, .Y =  2.0},
                          {.X = 2.0, .Y = -1.0},
                          {.X = 3.0, .Y =  1.5},
                          {.X = 4.0, .Y = -2.0},
                          {.X = 5.0, .Y =  1.0},
                          {.X = 6.0, .Y = -0.5},
                          {.X = 7.0, .Y = 0.5}},
              .Weights = {}};
    }

    /**
     *  A line through a point at an angle, a few units long either way.
     *
     *  @param  through The point.
     *  @param  angle   The angle, in radians.
     *  @return         The line.
     */
    auto Across(Point through, double angle) -> Line
    {
      const Vector way = {.X = 3.0 * std::cos(angle), .Y = 3.0 * std::sin(angle)};
      return {.Start = through - way, .End = through + way};
    }
  } // namespace

  TEST_CASE("pieced: crossings at a knot or a sliver from one, as met thoroughly, or not sure")
  {
    const Tolerance tolerance = {.Length = 1e-9, .Radians = 1e-10};
    int             sure      = 0;
    int             asked     = 0;
    for (int index = 1; index < 5; ++index)
    {
      const auto knot = static_cast<double>(index);
      for (const double off : {0.0, 1e-13, 1e-11, 1e-9, 1e-7, 1e-4, 0.25})
      {
        const Point through = At(Waved(), knot + off);
        for (const double angle : {1.3, 2.1, 0.4})
        {
          for (const bool swapped : {false, true})
          {
            const Curve     line   = Across(through, angle);
            const Curve     wavy   = Waved();
            const Curves    curves = swapped ? Curves{wavy, line} : Curves{line, wavy};
            const PiecedMet met    = PiecedBoth(curves, tolerance);
            CAPTURE(knot);
            CAPTURE(off);
            CAPTURE(angle);
            CHECK(met.Agrees);
            CHECK((!met.Sure || !met.Swiftly.Hits.empty()));
            sure  += static_cast<int>(met.Sure);
            asked += 1;
            if (off < 1e-10)
            {
              CHECK(!met.Sure);
            }
          }
        }
      }
    }
    MESSAGE("sure of " << sure << " of " << asked);
    CHECK(sure > asked / 3);
  }

  TEST_CASE("pieced: a spline met with itself turned, through its knots, and closed ones")
  {
    const Tolerance tolerance = {.Length = 1e-9, .Radians = 1e-10};
    Spline          turned    = Waved();
    for (Point &pole : turned.Poles)
    {
      pole = {.X = pole.X + 0.37, .Y = (-pole.Y * 0.5) + 0.1};
    }
    Spline ring = {.Degree = 3, .Knots = {}, .Poles = {}, .Weights = {}};
    for (int index = 0; index < 9; ++index)
    {
      const double angle = 2.0 * std::numbers::pi * index / 8.0;
      ring.Poles.push_back({.X = 2.0 * std::cos(angle), .Y = std::sin(angle)});
    }
    ring.Poles.back() = ring.Poles.front();
    ring.Knots        = {0, 0, 0, 0, 1, 2, 3, 4, 5, 6, 6, 6, 6};
    REQUIRE(Valid(ring));
    REQUIRE(At(ring, 0.0) == At(ring, 6.0));
    Spline shifted = ring;
    for (Point &pole : shifted.Poles)
    {
      pole = {.X = pole.X + 0.75, .Y = pole.Y + 0.25};
    }
    const std::array<std::array<Curve, 2>, 6> pairs = {
        std::array<Curve, 2>{Waved(), turned},
        std::array<Curve, 2>{turned, Waved()},
        std::array<Curve, 2>{ring, shifted},
        std::array<Curve, 2>{ring, Line{.Start = {.X = -3, .Y = 0.1}, .End = {.X = 3, .Y = -0.2}}},
        std::array<Curve, 2>{Ray{.Start = {.X = 0, .Y = 0}, .Way = {.X = 1, .Y = 0.3}}, ring},
        std::array<Curve, 2>{Xline{.Through = {.X = 0.5, .Y = 0}, .Way = {.X = 0.1, .Y = 1}},
                             Waved()}};
    for (std::size_t index = 0; index < pairs.size(); ++index)
    {
      const PiecedMet met = PiecedBoth({pairs.at(index)[0], pairs.at(index)[1]}, tolerance);
      CAPTURE(index);
      CHECK(met.Agrees);
      CHECK(met.Sure);
      CHECK(!met.Swiftly.Hits.empty());
    }
  }

  TEST_CASE("pieced: splines the kernels do not cut, and curves they do not meet with one")
  {
    const Tolerance tolerance = {.Length = 1e-9, .Radians = 1e-10};
    Spline          heavy     = Waved();
    heavy.Weights.assign(heavy.Poles.size(), 1.0);
    heavy.Weights[3] = 2.0;
    Spline straight  = {.Degree  = 1,
                        .Knots   = {0, 0, 1, 1},
                        .Poles   = {{.X = 0, .Y = -3}, {.X = 1, .Y = 3}},
                        .Weights = {}};
    Spline steep     = Waved();
    steep.Degree     = 4;
    steep.Knots      = {0, 0, 0, 0, 0, 1, 2, 3, 4, 4, 4, 4, 4};
    Spline broken    = Waved();
    broken.Knots     = {0, 0, 0, 0, 1, 2, 3, 4, 5, 5, 5};
    Spline far       = Waved();
    far.Poles[2].Y   = 0x1p300;
    const Curve arc  = Arc{.Center = {.X = 3, .Y = 0}, .Radius = 1, .Start = {}, .Sweep = Turn};
    const Curve wavy = Waved();
    const Curve line = Line{.Start = {.X = 0.5, .Y = -5}, .End = {.X = 0.6, .Y = 5}};
    const std::array<Curve, 5> curves = {heavy, straight, steep, broken, far};
    for (const Curve &curve : curves)
    {
      Intersections meeting;
      CHECK(!Swift::Met({curve, wavy}, tolerance, meeting));
      CHECK(meeting.Hits.empty());
    }
    /**
     *  A line meets a rational spline, or one past a cubic, by the kernel that holds its equation
     *  as a spline (tests/geoxy/strand.cpp); one of degree 1, broken or past the kernels' numbers
     *  it meets with none.
     */
    for (const Curve &curve : {curves[1], curves[3], curves[4]})
    {
      Intersections meeting;
      CHECK(!Swift::Met({line, curve}, tolerance, meeting));
    }
    for (const Curve &curve : {curves[0], curves[2]})
    {
      const PiecedMet met = PiecedBoth({line, curve}, tolerance);
      CHECK(met.Sure);
      CHECK(met.Agrees);
    }
    /**
     *  An arc and a spline are the ring kernel's, span by span, as the thorough way meets them.
     */
    Intersections meeting;
    CHECK(Swift::Met({arc, wavy}, tolerance, meeting));
    CHECK(meeting.Hits.size() == Intersect(arc, wavy, tolerance).Hits.size());
    CHECK(Swift::Pieced({line, wavy}, tolerance, meeting) == std::nullopt);
  }
} // namespace Geoxy

namespace Geoxy
{
  TEST_CASE("pieced: two answers better than the thorough way's, each checked exactly")
  {
    /**
     *  An endless line across a spline a million units out, the crossing's y some 3e-6: the
     *  thorough way polishes it in double-double to 2^-64 of the size about it, short of the digits
     *  so small a coordinate keeps, and rounds it two steps off. The kernel proves the double
     *  nearest, as exact rational arithmetic finds it: x 0x1.e8480000046e5p+19, y
     *  -0x1.bab752e39985ep-19, the line's place 0x1.1608e17a02083p-17.
     */
    const Curve spline =
        Spline{.Degree  = 3,
               .Knots   = {-0x1.6c5c28820b794p+9, -0x1.6c5c28820b794p+9, -0x1.6c5c28820b794p+9,
                           -0x1.6c5c28820b794p+9, -0x1.6bdc28820b794p+9, -0x1.6bdc28820b794p+9,
                           -0x1.6bdc28820b794p+9, -0x1.6bdc28820b794p+9},
               .Poles   = {{.X = 0x1.e847fffffdf16p+19, .Y = -0x1.ea842abc6ca1p-18},
                           {.X = 0x1.e84800000efb2p+19, .Y = 0x1.31f94c907a59cp-18},
                           {.X = 0x1.e84800001beabp+19, .Y = 0x1.aee152508fd0dp-19},
                           {.X = 0x1.e84800002b548p+19, .Y = 0x1.c300ff1cc2db1p-19}},
               .Weights = {3, 3, 3, 3}};
    const Curve across = Xline{.Through = {.X = 0x1.e847fffff372ep+19, .Y = 0x1.298fe131a73f3p-18},
                               .Way     = {.X = 0x1.f45ccf3fb4756p-1, .Y = -0x1.ddcb7ae46dabcp-1}};
    const Intersections far =
        Intersect(spline, across, {.Length = 0x1.6ca62500c3cacp-36, .Radians = 1e-10});
    REQUIRE(far.Hits.size() == 1);
    CHECK(far.Hits[0].At == Point{.X = 0x1.e8480000046e5p+19, .Y = -0x1.bab752e39985ep-19});
    CHECK(far.Hits[0].Second.Along == 0x1.1608e17a02083p-17);

    /**
     *  Two splines of poles given twice and thrice, crossing steeply at a sine of 0.9: the thorough
     *  way finds no crossing, at any tolerance up to 1e-6; the kernel proves one, both curves'
     *  points there within 1e-15 of each other.
     */
    const Curve first =
        Spline{.Degree  = 3,
               .Knots   = {0, 0, 0, 0, 0x1.90252ee58cfep-2, 0x1.52b9dbf22a0e3p+0,
                           0x1.52b9dbf22a0e3p+0, 0x1.52b9dbf22a0e3p+0, 0x1.52b9dbf22a0e3p+0},
               .Poles   = {{.X = 0x1.e33cd456b4424p-2, .Y = -0x1.77f0595f53f94p-1},
                           {.X = 0x1.18bcbdbe0eda2p-1, .Y = 0x1.5b0618fe6478p-2},
                           {.X = 0x1.18bcbdbe0eda2p-1, .Y = 0x1.5b0618fe6478p-2},
                           {.X = 0x1.18bcbdbe0eda2p-1, .Y = 0x1.5b0618fe6478p-2},
                           {.X = -0x1.d91522f366808p-2, .Y = 0x1.a5ef79c58fc2p-2}},
               .Weights = {3, 3, 3, 3, 3}};
    const Curve second =
        Spline{.Degree  = 2,
               .Knots   = {0x1.3e843b5a4f136p+6, 0x1.3e843b5a4f136p+6, 0x1.3e843b5a4f136p+6,
                           0x1.41bb7a21fd741p+6, 0x1.46b5603f4fap+6, 0x1.497812ddf2acfp+6,
                           0x1.4cc7a2a781341p+6, 0x1.4f70e9ff06c2ap+6, 0x1.4f70e9ff06c2ap+6,
                           0x1.4f70e9ff06c2ap+6},
               .Poles   = {{.X = 0x1.4479edff12fbp-1, .Y = 0x1.e559930535668p-1},
                           {.X = -0x1.f5ba808f942ap-5, .Y =  0x1.9c9af6de2e1e4p-2},
                           {.X = -0x1.f5ba808f942ap-5, .Y =  0x1.9c9af6de2e1e4p-2},
                           {.X = 0x1.e8711725a548cp-2, .Y = -0x1.d979cbc829ba8p-1},
                           {.X =  0x1.e241ae953a36p-4, .Y = -0x1.33b7c15c270fcp-2},
                           {.X =  0x1.e241ae953a36p-4, .Y = -0x1.33b7c15c270fcp-2},
                           {.X = -0x1.9923e8765fdecp-2, .Y = -0x1.53f1e7e3f64a6p-1}},
               .Weights = {}};
    const Intersections missed =
        Intersect(first, second, {.Length = 0x1.b7cdfd9d7bdbbp-33, .Radians = 1e-10});
    REQUIRE(missed.Hits.size() == 1);
    const Hit &hit = missed.Hits[0];
    CHECK(hit.Kind == Contact::Cross);
    CHECK(Distance(At(first, hit.First.Along), At(second, hit.Second.Along)) < 1e-15);
    CHECK(Distance(hit.At, hit.First.At) < 1e-15);
    const Vector one   = Derivative(first, hit.First.Along);
    const Vector other = Derivative(second, hit.Second.Along);
    CHECK(std::abs(Cross(one, other)) > 0.8 * Length(one) * Length(other));
    Intersections thorough;
    Thorough({first, second}, {.Length = 1e-6, .Radians = 1e-10}, thorough);
    CHECK(thorough.Hits.empty());
  }
} // namespace Geoxy

namespace Geoxy
{
  TEST_CASE(
      "pieced: crossings of different pieces the thorough way may join, or one point: not sure")
  {
    const Tolerance tolerance = {.Length = 1e-9, .Radians = 1e-10};
    const auto      moved     = [](Point point) -> Point
    {
      return {.X = point.X + 0.25, .Y = point.Y + 0.375};
    };

    /**
     *  An S and a wave of three spans, each turned half about one point, crossing there and once
     *  either side, each crossing on a span of its own: the middle one halfway between the others
     *  along both curves, so that neither look halfway tells the outer two apart.
     */
    const Curve ess  = Bezier::Cubic{.Start = moved({.X = 0, .Y = 0}),
                                     .Out   = moved({.X = 1.2, .Y = 0.8}),
                                     .In    = moved({.X = 2.8, .Y = -0.8}),
                                     .End   = moved({.X = 4, .Y = 0})};
    Spline      wave = {.Degree  = 3,
                        .Knots   = {0, 0, 0, 0, 1, 2, 3, 3, 3, 3},
                        .Poles   = {{.X = 0.5, .Y = -1.875},
                                    {.X = 1.0, .Y =  2.5},
                                    {.X = 1.6, .Y =  1.0},
                                    {.X = 2.4, .Y = -1.0},
                                    {.X = 3.0, .Y = -2.5},
                                    {.X = 3.5, .Y = 1.875}},
                        .Weights = {}};
    for (Point &pole : wave.Poles)
    {
      pole = moved(pole);
    }
    const Curve   waved = wave;
    Intersections meeting;
    CHECK(Swift::Pieced({ess, waved}, tolerance, meeting) == std::optional<bool>{false});
    CHECK(meeting.Hits.empty());
    CHECK(PiecedBoth({waved, ess}, tolerance).Agrees);

    /**
     *  A loop of two spans, and a curve across its double point: two crossings at one point, on a
     *  span each.
     */
    Spline loop = {.Degree  = 3,
                   .Knots   = {0, 0, 0, 0, 0.5, 1, 1, 1, 1},
                   .Poles   = {{.X = -1, .Y = 0},
                               {.X =  0.5, .Y = 1},
                               {.X =    0, .Y = 2},
                               {.X = -0.5, .Y = 1},
                               {.X = 1, .Y = 0}},
                   .Weights = {}};
    for (Point &pole : loop.Poles)
    {
      pole = moved(pole);
    }
    const double twice  = 6.0 / 7.0;
    const Curve  looped = loop;
    const Curve  across = Bezier::Quadratic{.Start   = moved({.X = -0.8, .Y = twice - 0.3}),
                                            .Control = moved({.X = 0.0, .Y = twice + 0.05}),
                                            .End     = moved({.X = 0.8, .Y = twice + 0.2})};
    CHECK(Swift::Pieced({looped, across}, tolerance, meeting) == std::optional<bool>{false});
    CHECK(PiecedBoth({across, looped}, tolerance).Agrees);
  }

  TEST_CASE(
      "pieced: knots past the range the kernels take, or a fraction of them below the normals")
  {
    const Tolerance tolerance = {.Length = 1e-9, .Radians = 1e-10};
    const Curve     line      = Line{.Start = {.X = 0.5, .Y = -5}, .End = {.X = 0.6, .Y = 5}};
    Spline          far       = Waved();
    for (double &knot : far.Knots)
    {
      knot *= 1e70;
    }
    Spline faint = {
        .Degree  = 2,
        .Knots   = {0, 0, 0, 1e-310, 1, 1, 1},
        .Poles   = {{.X = 0, .Y = 0}, {.X = 0.4, .Y = 1}, {.X = 0.8, .Y = -1}, {.X = 1, .Y = 0}},
        .Weights = {}};
    REQUIRE(Valid(far));
    REQUIRE(Valid(faint));
    for (const Spline &spline : {far, faint})
    {
      CHECK(!Swift::Sectioned(spline));
      const Curve   curve = spline;
      Intersections meeting;
      CHECK(!Swift::Met({line, curve}, tolerance, meeting));
    }
  }
} // namespace Geoxy
