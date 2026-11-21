/***************************************************************************************************
 *  @file      tests/geoxy/rational.cpp
 *  @brief     Tests of a spline's strips, their roots and their integrals, and of the roots of a
 *             polynomial in Bernstein form.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "rational.hpp"

#include "aberth.hpp"
#include "eased.hpp"

#include "geoxy/area.hpp"
#include "geoxy/length.hpp"
#include "geoxy/point.hpp"
#include "geoxy/spline.hpp"
#include "geoxy/vector.hpp"

#include <doctest/doctest.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <span>
#include <vector>

namespace Geoxy
{
  namespace
  {
    /**
     *  Whether some roots hold one within a distance of a point.
     *
     *  @param  roots The roots.
     *  @param  root  The point.
     *  @param  reach The distance.
     *  @return       Whether they do.
     */
    auto Holds(const Aberth::Roots &roots, Vector root, double reach) -> bool
    {
      return std::ranges::any_of(std::span{roots.Values}.first(roots.Count),
                                 [root, reach](Vector value) -> bool
                                 { return Length(value - root) <= reach; });
    }

    /**
     *  A spline's reverse: the same curve run the other way.
     *
     *  @param  spline The spline.
     *  @return        Its reverse.
     */
    auto RationalReversed(Spline spline) -> Spline
    {
      std::ranges::reverse(spline.Poles);
      std::ranges::reverse(spline.Weights);
      std::ranges::reverse(spline.Knots);
      for (double &knot : spline.Knots)
      {
        knot = -knot;
      }
      return spline;
    }
  } // namespace

  TEST_CASE("rational: the roots of a polynomial in Bernstein form")
  {
    /**
     *  1 - 2 t, its root at the middle, where its power about the middle has no constant: the
     *  approximation starts there.
     */
    const std::array<Vector, 2> line = {Vector{.X = 1, .Y = 0}, Vector{.X = -1, .Y = 0}};
    const Aberth::Roots         half = Aberth::Solved(line, Aberth::Rounds);
    REQUIRE(half.Count == 1);
    CHECK(Holds(half, {.X = 0.5, .Y = 0}, 1e-15));

    /**
     *  (t - 1/4) (t - 3) and t^2 + 1 in Bernstein form, and a cubic with a root far past 1, where
     *  the sums are taken from 1.
     */
    const std::array<Vector, 3> real    = {Vector{.X = 0.75, .Y = 0}, Vector{.X = -0.875, .Y = 0},
                                           Vector{.X = -1.5, .Y = 0}};
    const Aberth::Roots         crossed = Aberth::Solved(real, Aberth::Rounds);
    REQUIRE(crossed.Count == 2);
    CHECK(Holds(crossed, {.X = 0.25, .Y = 0}, 1e-15));
    CHECK(Holds(crossed, {.X = 3, .Y = 0}, 1e-14));
    const std::array<Vector, 3> lifted = {Vector{.X = 1, .Y = 0}, Vector{.X = 1, .Y = 0},
                                          Vector{.X = 2, .Y = 0}};
    const Aberth::Roots         plane  = Aberth::Solved(lifted, Aberth::Rounds);
    REQUIRE(plane.Count == 2);
    CHECK(Holds(plane, {.X = 0, .Y = 1}, 1e-15));
    CHECK(Holds(plane, {.X = 0, .Y = -1}, 1e-15));
    const std::array<Vector, 4> far = {Vector{.X = 1, .Y = 0}, Vector{.X = 2.0 / 3.0, .Y = 1},
                                       Vector{.X = 1.0 / 3.0, .Y = -1}, Vector{.X = 1e-12, .Y = 0}};
    const Aberth::Roots         out = Aberth::Solved(far, Aberth::Rounds);
    REQUIRE(out.Count == 3);
    for (const Vector root : std::span{out.Values}.first(out.Count))
    {
      CHECK(std::isfinite(root.X));
      CHECK(std::isfinite(root.Y));
    }

    /**
     *  In closed form: t^2, its root at 0 twice; 1 - t^2, one root at 1 where b_2 is 0; (1 - t)^2,
     *  at 1 twice; 1 - 2 t, of degree 1 in a form of degree 2, its other root far off; and a
     *  constant in a form of degree 1, its root far off.
     */
    const auto solved = [](std::array<double, 3> values) -> Aberth::Roots
    {
      const std::array<Vector, 3> raised = {Vector{.X = values[0], .Y = 0},
                                            Vector{.X = values[1], .Y = 0},
                                            Vector{.X = values[2], .Y = 0}};
      return Aberth::Solved(raised, Aberth::Rounds);
    };
    CHECK(Holds(solved({0, 0, 1}), {.X = 0, .Y = 0}, 0));
    const Aberth::Roots square = solved({1, 1, 0});
    CHECK(Holds(square, {.X = 1, .Y = 0}, 0));
    CHECK(Holds(square, {.X = -1, .Y = 0}, 1e-15));
    CHECK(solved({1, 0, 0}).Count == 2);
    CHECK(Holds(solved({1, 0, 0}), {.X = 1, .Y = 0}, 0));
    const Aberth::Roots straight = solved({1, 0, -1});
    CHECK(Holds(straight, {.X = 0.5, .Y = 0}, 1e-15));
    CHECK(Holds(straight, {.X = 0x1p+500, .Y = 0}, 0));
    const std::array<Vector, 2> level = {Vector{.X = 1, .Y = 0}, Vector{.X = 1, .Y = 0}};
    CHECK(Holds(Aberth::Solved(level, Aberth::Rounds), {.X = 0x1p+500, .Y = 0}, 0));

    /**
     *  Coefficients so small, and so large, that a step's length squared leaves the doubles.
     */
    for (const double scale : {1e-200, 1e200})
    {
      std::array<Vector, 4> scaled = far;
      for (Vector &value : scaled)
      {
        value = value * scale;
      }
      const Aberth::Roots roots = Aberth::Solved(scaled, Aberth::Rounds);
      for (std::size_t index = 0; index < roots.Count; ++index)
      {
        CHECK(Length(roots.Values.at(index) - out.Values.at(index)) <= 1e-9);
      }
    }

    /**
     *  One round only: the approximations are where it leaves them.
     */
    const Aberth::Roots once = Aberth::Solved(far, 1);
    CHECK(once.Count == 3);
  }

  TEST_CASE("rational: a spline's strips")
  {
    /**
     *  A polynomial span is one strip; a rational span whose weight's roots crowd an end is cut so
     *  far from it, and its parts again, every strip's place on the spline told from both ends.
     */
    const Spline cubic = {
        .Degree  = 3,
        .Knots   = {0, 0, 0, 0, 1, 1, 1, 1},
        .Poles   = {{.X = 0, .Y = 0}, {.X = 1, .Y = 2}, {.X = 3, .Y = 2}, {.X = 4, .Y = 0}},
        .Weights = {1, 1e16, 1e16, 1}};
    const std::vector<Strip> strips = Stripped(cubic);
    REQUIRE(strips.size() > 1);
    double width = 0.0;
    for (const Strip &strip : strips)
    {
      CHECK_FALSE(strip.Plain);
      width += strip.Where.Width;
      CHECK(Local(strip.Where, Homed(strip.Where, 0.0)) == doctest::Approx(0.0).epsilon(1e-6));
      CHECK(Local(strip.Where, Homed(strip.Where, 0.5))
            == doctest::Approx(0.5).epsilon(0x1p-50 / strip.Where.Width));
    }
    CHECK(width == doctest::Approx(1.0).epsilon(1e-15));
    Spline even  = cubic;
    even.Weights = {};
    REQUIRE(Stripped(even).size() == 1);
    CHECK(Stripped(even).front().Plain);

    /**
     *  Weights some 2^1120 apart: cut as often as doubles allow, and then taken as they are.
     */
    const double light = std::ldexp(1.0, -560);
    const Spline steep = {.Degree  = 2,
                          .Knots   = {0, 0, 0, 1, 1, 1},
                          .Poles   = {{.X = 0, .Y = 0}, {.X = 1, .Y = 1}, {.X = 2, .Y = 0}},
                          .Weights = {light, std::ldexp(1.0, 560), light}};
    REQUIRE(Valid(steep));
    CHECK(std::isfinite(Length(steep, {})));
    CHECK(std::isfinite(Area(steep, {})));

    /**
     *  Poles so far apart that one less the span's middle is past every double: taken less nothing.
     */
    const Spline wide = {
        .Degree  = 2,
        .Knots   = {0, 0, 0, 1, 1, 1},
        .Poles   = {{.X = 1.7e308, .Y = 0}, {.X = -1.7e308, .Y = 1}, {.X = 1.7e308, .Y = 2}},
        .Weights = {1, 0.5, 1}};
    REQUIRE(Valid(wide));
    REQUIRE(Stripped(wide).size()       == 1);
    CHECK(Stripped(wide).front().Offset == Point{});
    CHECK_FALSE(std::isnan(Length(wide, {})));
  }

  TEST_CASE("rational: a strip's length at its edges")
  {
    /**
     *  A strip that does not move has no length; a stretch of no width none; one without end an
     *  endless one.
     */
    const Spline still = {.Degree  = 2,
                          .Knots   = {0, 0, 0, 1, 1, 1},
                          .Poles   = {{.X = 1, .Y = 2}, {.X = 1, .Y = 2}, {.X = 1, .Y = 2}},
                          .Weights = {1, 3, 1}};
    const Pace   point = Paced(Stripped(still).front());
    CHECK(Measured(point, {.Start = 0, .End = 1}, {}) == 0);
    CHECK(Length(still, {}) == 0);
    const Spline quartic = {.Degree  = 4,
                            .Knots   = {0, 0, 0, 0, 0, 1, 1, 1, 1, 1},
                            .Poles   = {{.X = 0, .Y = 0},
                                        {.X = 1, .Y =  2},
                                        {.X = 2, .Y = -1},
                                        {.X = 3, .Y =  2},
                                        {.X = 4, .Y = 0}},
                            .Weights = {}};
    const Pace   pace    = Paced(Stripped(quartic).front());
    CHECK(Measured(pace, {.Start = 0.5, .End = 0.5}, {}) == 0);
    CHECK(Measured(pace, {.Start = 0, .End = std::numeric_limits<double>::infinity()}, {})
          == std::numeric_limits<double>::infinity());
    CHECK(Length(quartic, {.Start = 0, .End = std::numeric_limits<double>::infinity()}, {})
          == std::numeric_limits<double>::infinity());

    /**
     *  With no halvings to spend, each piece takes the most points: near enough here.
     */
    const Spline bent  = {.Degree  = 2,
                          .Knots   = {0, 0, 0, 1, 1, 1},
                          .Poles   = {{.X = 0, .Y = 0}, {.X = 1, .Y = 1}, {.X = 2, .Y = 0}},
                          .Weights = {1, 3, 1}};
    const double whole = Length(bent, {});
    const Pace   stiff = Paced(Stripped(bent).front());
    CHECK(Measured(stiff, {.Start = 0, .End = 1}, {.Within = {}, .Halvings = 0})
          == doctest::Approx(whole).epsilon(1e-3));

    /**
     *  A distance no sum of doubles reaches to within nothing takes every step.
     */
    const double along = Walked(stiff, whole / 3, 0.0);
    CHECK(Measured(stiff, {.Start = 0, .End = along}, {}) == doctest::Approx(whole / 3));

    /**
     *  A rational spline of degree 1 runs straight from pole to pole: its speed's numerator is a
     *  constant, and the area it encloses the polygon's.
     */
    const Spline straight = {.Degree  = 1,
                             .Knots   = {0, 1, 2, 3, 4},
                             .Poles   = {{.X = 0, .Y = 0}, {.X = 3, .Y = 0}, {.X = 5, .Y = 0}},
                             .Weights = {1, 5, 2}};
    CHECK(Length(straight, {}) == doctest::Approx(5).epsilon(1e-15));
    CHECK(Area(straight, {})   == 0);
    const Spline triangle = {.Degree  = 1,
                             .Knots   = {0, 0, 1, 2, 2},
                             .Poles   = {{.X = 0, .Y = 0}, {.X = 3, .Y = 0}, {.X = 3, .Y = 4}},
                             .Weights = {1, 5, 2}};
    CHECK(Length(triangle, {}) == doctest::Approx(7).epsilon(1e-15));
    CHECK(Area(triangle, {})   == 6);
  }

  TEST_CASE("rational: what a strip encloses, and which way a spline is taken")
  {
    /**
     *  With no halvings to spend, each piece takes the most points.
     */
    const Spline bent  = {.Degree  = 2,
                          .Knots   = {0, 0, 0, 1, 1, 1},
                          .Poles   = {{.X = 0, .Y = 0}, {.X = 1, .Y = 1}, {.X = 2, .Y = 0}},
                          .Weights = {1, 3, 1}};
    const Strip  strip = Stripped(bent).front();
    const double area  = Enclosed(strip, {.Origin = {}, .Error = 0, .Moments = true}).Area;
    CHECK(Enclosed(strip, {.Origin = {}, .Error = 0, .Moments = false, .Halvings = 0}).Area
          == doctest::Approx(area).epsilon(1e-3));

    /**
     *  A spline whose poles read the same either way is taken by its weights, then its knots; one
     *  that reads the same either way in all of them runs out and back along itself, and encloses
     *  nothing.
     */
    const Spline weighed = {.Degree  = 2,
                            .Knots   = {0, 0, 0, 1, 1, 1},
                            .Poles   = {{.X = 0, .Y = 0}, {.X = 1, .Y = 1}, {.X = 0, .Y = 0}},
                            .Weights = {1, 2, 3}};
    CHECK(Area(RationalReversed(weighed), {}) == -Area(weighed, {}));
    const Spline knotted = {
        .Degree  = 2,
        .Knots   = {0, 0, 0, 0.25, 1, 1, 1},
        .Poles   = {{.X = 0, .Y = 0}, {.X = 1, .Y = 1}, {.X = 1, .Y = 1}, {.X = 0, .Y = 0}},
        .Weights = {1, 2, 2, 1}};
    CHECK(Area(RationalReversed(knotted), {}) == -Area(knotted, {}));
    const Spline folded = {
        .Degree  = 2,
        .Knots   = {-1, -1, -1, 0, 1, 1, 1},
        .Poles   = {{.X = 0, .Y = 0}, {.X = 1, .Y = 1}, {.X = 1, .Y = 1}, {.X = 0, .Y = 0}},
        .Weights = {1, 2, 2, 1}};
    CHECK(Area(folded, {}) == 0);
    CHECK_FALSE(Centroid(folded, {}).has_value());
  }
} // namespace Geoxy
