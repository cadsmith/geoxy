/***************************************************************************************************
 *  @file      tests/geoxy/dash.cpp
 *  @brief     Tests of dashes laid along curves: their lengths against the pattern's, their sum
 *             against the curve's, and PLINEGEN's two ways along a polyline.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "geoxy/dash.hpp"

#include "curves.hpp"
#include "real.hpp"

#include "geoxy/angle.hpp"
#include "geoxy/arc.hpp"
#include "geoxy/at.hpp"
#include "geoxy/bezier/cubic.hpp"
#include "geoxy/bezier/quadratic.hpp"
#include "geoxy/curve.hpp"
#include "geoxy/domain.hpp"
#include "geoxy/ellipse.hpp"
#include "geoxy/fit.hpp"
#include "geoxy/length.hpp"
#include "geoxy/line.hpp"
#include "geoxy/matrix.hpp"
#include "geoxy/nearest.hpp"
#include "geoxy/point.hpp"
#include "geoxy/polyline.hpp"
#include "geoxy/ray.hpp"
#include "geoxy/spline.hpp"
#include "geoxy/vector.hpp"
#include "geoxy/xline.hpp"

#include <doctest/doctest.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <vector>

namespace Geoxy
{
  namespace
  {
    /**
     *  The lengths a pattern's dashes have along a curve of a length from its start: each dash's,
     *  the last cut short by the end, a dot's nothing.
     *
     *  @param  elements The pattern's elements.
     *  @param  length   The curve's length.
     *  @return          The lengths, in order.
     */
    auto Expected(const std::vector<double> &elements, double length) -> std::vector<double>
    {
      std::vector<double> lengths;
      double              walked = 0.0;
      for (std::size_t index = 0; walked <= length; index = (index + 1) % elements.size())
      {
        const double element = elements[index];
        const bool   dot     = Zero(element);
        if (dot || (element > 0.0 && walked < length))
        {
          lengths.push_back(std::min(element, length - walked));
        }
        walked += std::abs(element);
      }
      return lengths;
    }

    /**
     *  A dash's length against the one the pattern lays, its kind, and its ends on the curve.
     *
     *  @param  curve    The curve.
     *  @param  dash     The dash.
     *  @param  expected The length it should have, and how far it may be off.
     *  @return          Its length.
     */
    auto Matched(const Curve &curve, const Curve &dash, std::array<double, 2> expected) -> double
    {
      const auto [wanted, slack] = expected;
      const double along         = Length(dash, {});
      const Point  start         = At(dash, Domain(dash).Start);
      const Point  end           = At(dash, Domain(dash).End);
      CHECK(std::abs(along - wanted) <= slack);
      CHECK(dash.index() == (wanted > 0.0 ? curve.index() : Curve{Line{}}.index()));
      CHECK(Distance(start, Nearest(curve, start).At) <= slack);
      CHECK(Distance(end, Nearest(curve, end).At) <= slack);
      return along;
    }

    /**
     *  Every dash's length against what the pattern lays, and their ends on the curve.
     *
     *  @param curve    The curve.
     *  @param elements The pattern's elements.
     *  @param slack    How far a length may be off.
     */
    void Checked(const Curve &curve, const std::vector<double> &elements, double slack)
    {
      const std::vector<Curve>  dashes   = Dashes(curve, {.Elements = elements}, {});
      const std::vector<double> expected = Expected(elements, Length(curve, {}));
      REQUIRE(dashes.size() == expected.size());
      double sum  = 0.0;
      double laid = 0.0;
      for (std::size_t index = 0; index < dashes.size(); ++index)
      {
        sum  += Matched(curve, dashes[index], {expected[index], slack});
        laid += expected[index];
      }
      CHECK(std::abs(sum - laid) <= slack * static_cast<double>(dashes.size()));
    }

    /**
     *  A polyline of a square's sides, 4 long, closed.
     *
     *  @return The polyline.
     */
    auto Square() -> Polyline
    {
      return {.Vertices = {{.Position = {.X = 0, .Y = 0}},
                           {.Position = {.X = 4, .Y = 0}},
                           {.Position = {.X = 4, .Y = 4}},
                           {.Position = {.X = 0, .Y = 4}}},
              .Closed   = true};
    }
  } // namespace

  TEST_CASE("dash: a line's dashes, gaps and dots")
  {
    const Line               line   = {.Start = {.X = 0, .Y = 0}, .End = {.X = 10, .Y = 0}};
    const std::vector<Curve> dashes = Dashes(line, {.Elements = {2, -1, 0, -1}}, {});
    const std::vector<Curve> wanted = {
        Line{.Start = {.X = 0, .Y = 0}, .End = {.X = 2, .Y = 0}},
        Line{.Start = {.X = 3, .Y = 0}, .End = {.X = 3, .Y = 0}},
        Line{.Start = {.X = 4, .Y = 0}, .End = {.X = 6, .Y = 0}},
        Line{.Start = {.X = 7, .Y = 0}, .End = {.X = 7, .Y = 0}},
        Line{.Start = {.X = 8, .Y = 0}, .End = {.X = 10, .Y = 0}},
    };
    CHECK(dashes == wanted);

    const std::vector<Curve> skipped = Dashes(line, {.Elements = {2, -1, 0, -1}, .Phase = 3.5}, {});
    REQUIRE(skipped.size() == 5);
    CHECK(skipped[0]       == Curve{Line{.Start = {.X = 0.5, .Y = 0}, .End = {.X = 2.5, .Y = 0}}});

    const std::vector<Curve> later = Dashes(line, {.Elements = {2, -2}, .Phase = 1}, {});
    CHECK(later
          == std::vector<Curve>{Line{.Start = {.X = 0, .Y = 0}, .End = {.X = 1, .Y = 0}},
                                Line{.Start = {.X = 3, .Y = 0}, .End = {.X = 5, .Y = 0}},
                                Line{.Start = {.X = 7, .Y = 0}, .End = {.X = 9, .Y = 0}}});
    const std::vector<Curve> back = Dashes(line, {.Elements = {2, -2}, .Phase = -1}, {});
    CHECK(back
          == std::vector<Curve>{Line{.Start = {.X = 1, .Y = 0}, .End = {.X = 3, .Y = 0}},
                                Line{.Start = {.X = 5, .Y = 0}, .End = {.X = 7, .Y = 0}},
                                Line{.Start = {.X = 9, .Y = 0}, .End = {.X = 10, .Y = 0}}});
    const double nan = std::numeric_limits<double>::quiet_NaN();
    CHECK(Dashes(line, {.Elements = {2, -2}, .Phase = nan}, {})
          == Dashes(line, {.Elements = {2, -2}}, {}));
  }

  TEST_CASE("dash: every kind laid by length")
  {
    const std::vector<double> pattern = {0.75, -0.25, 0, -0.25};
    for (const Vector offset : {Vector{}, Vector{.X = 5e6, .Y = -5e6}})
    {
      const Point    base   = Point{} + offset;
      const double   slack  = offset.X > 0 ? 1e-8 : 1e-12;
      const Polyline bulged = {
          .Vertices = {{.Position = base, .Bulge = 0.5},
                       {.Position = base + Vector{.X = 4, .Y = 0}, .Bulge = 0},
                       {.Position = base + Vector{.X = 4, .Y = 3}, .Bulge = -1},
                       {.Position = base + Vector{.X = 1, .Y = 5}}},
          .Closed   = false};
      Checked(Line{.Start = base, .End = base + Vector{.X = 3, .Y = 4}}, pattern, slack);
      Checked(Arc{.Center = base, .Radius = 3, .Start = Degrees(30), .Sweep = Degrees(250)},
              pattern, slack);
      Checked(Ellipse{.Map   = Translation(offset) * Scaling(4, 1.5),
                      .Start = Radians(0.5),
                      .Sweep = Radians(-4)},
              pattern, slack);
      Checked(bulged, pattern, slack);
      Checked(Bezier::Quadratic{.Start   = base,
                                .Control = base + Vector{.X = 2, .Y = 5},
                                .End     = base + Vector{.X = 5, .Y = 0}},
              pattern, slack);
      Checked(Bezier::Cubic{.Start = base + (Tests::Wave.Start - Point{}),
                            .Out   = base + (Tests::Wave.Out - Point{}),
                            .In    = base + (Tests::Wave.In - Point{}),
                            .End   = base + (Tests::Wave.End - Point{})},
              pattern, slack);
      const std::array<Point, 5> points = {
          base, base + Vector{.X = 2, .Y = 1}, base + Vector{.X = 3, .Y = -1},
          base + Vector{.X = 6, .Y = 2}, base + Vector{.X = 7, .Y = 0}};
      const auto fit = Fit(points, {});
      REQUIRE(fit.has_value());
      Checked(*fit, pattern, slack);
    }
    Checked(Tests::Uneven(), {1.5, -0.5}, 1e-9);
  }

  TEST_CASE("dash: a polyline's pattern across its vertices, or again at each")
  {
    const std::vector<Curve> across = Dashes(Square(), {.Elements = {5, -1}}, {});
    REQUIRE(across.size() == 3);
    CHECK(std::get<Polyline>(across[0]).Vertices.size()          == 3);
    CHECK(std::get<Polyline>(across[0]).Vertices.back().Position == Point{.X = 4, .Y = 1});
    CHECK(std::get<Polyline>(across[2]).Vertices.back().Position == Point{});

    const std::vector<Curve> each = Dashes(Square(), {.Elements = {5, -1}, .Across = false}, {});
    REQUIRE(each.size() == 4);
    CHECK(std::get<Polyline>(each[3]).Vertices.size()           == 2);
    CHECK(std::get<Polyline>(each[3]).Vertices.front().Position == Point{.X = 0, .Y = 4});

    const std::vector<Curve> again = Dashes(Square(), {.Elements = {1, -1}, .Across = false}, {});
    CHECK(again.size() == 8);
    CHECK(std::get<Polyline>(again[2]).Vertices.front().Position == Point{.X = 4, .Y = 0});
  }

  TEST_CASE("dash: what is drawn whole")
  {
    const Line               line  = {.Start = {.X = 0, .Y = 0}, .End = {.X = 10, .Y = 0}};
    const Curve              curve = line;
    const std::vector<Curve> whole = {curve};
    const double             nan   = std::numeric_limits<double>::quiet_NaN();
    const double             huge  = std::numeric_limits<double>::infinity();
    CHECK(Dashes(line, {}, {}) == whole);
    CHECK(Dashes(line, {.Elements = {0, 0}}, {})        == whole);
    CHECK(Dashes(line, {.Elements = {1, nan}}, {})      == whole);
    CHECK(Dashes(line, {.Elements = {1, -huge}}, {})    == whole);
    CHECK(Dashes(line, {.Elements = {8, -4}}, {})       == whole);
    CHECK(Dashes(line, {.Elements = {1e-6, -1e-6}}, {}) == whole);
    const Ray   ray   = {.Start = {}, .Way = {.X = 1, .Y = 0}};
    const Xline xline = {.Through = {}, .Way = {.X = 0, .Y = 1}};
    CHECK(Dashes(ray, {.Elements = {1, -1}}, {})   == std::vector<Curve>{ray});
    CHECK(Dashes(xline, {.Elements = {1, -1}}, {}) == std::vector<Curve>{xline});
    const Polyline lone = {.Vertices = {{.Position = {.X = 1, .Y = 1}}}, .Closed = false};
    CHECK(Dashes(lone, {.Elements = {1, -1}}, {})     == std::vector<Curve>{lone});
    CHECK(Dashes(Spline{}, {.Elements = {1, -1}}, {}) == std::vector<Curve>{Spline{}});

    const std::vector<Curve> fine =
        Dashes(Square(), {.Elements = {1e-6, -1e-6}, .Across = false}, {});
    REQUIRE(fine.size() == 4);
    CHECK(std::get<Polyline>(fine[1]).Vertices.front().Position == Point{.X = 4, .Y = 0});
    const std::vector<Curve> brief = Dashes(Square(), {.Elements = {6, -1}, .Across = false}, {});
    CHECK(brief.size() == 4);
    const Polyline one = {
        .Vertices = {{.Position = {.X = 0, .Y = 0}}, {.Position = {.X = 1, .Y = 0}}},
        .Closed   = false};
    CHECK(Dashes(one, {.Elements = {6, -1}, .Across = false}, {}) == std::vector<Curve>{one});
  }
} // namespace Geoxy
