/***************************************************************************************************
 *  @file      tests/geoxy/text.cpp
 *  @brief     Tests of curves written as text and read back.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "support/text.hpp"

#include "draws.hpp"

#include "geoxy/angle.hpp"
#include "geoxy/arc.hpp"
#include "geoxy/bezier/cubic.hpp"
#include "geoxy/bezier/quadratic.hpp"
#include "geoxy/curve.hpp"
#include "geoxy/ellipse.hpp"
#include "geoxy/line.hpp"
#include "geoxy/matrix.hpp"
#include "geoxy/polyline.hpp"
#include "geoxy/ray.hpp"
#include "geoxy/spline.hpp"
#include "geoxy/xline.hpp"

#include <doctest/doctest.h>

#include <array>
#include <bit>
#include <cmath>
#include <cstdint>
#include <limits>
#include <optional>
#include <string>
#include <variant>

namespace Geoxy
{
  namespace
  {
    /**
     *  Whether a curve reads back from its own text as the same curve.
     *
     *  @param  curve The curve.
     *  @return       Whether it does.
     */
    auto RoundTrips(const Curve &curve) -> bool
    {
      const std::optional<Curve> back = FromText(ToText(curve));
      return back.has_value() && *back == curve;
    }

    /**
     *  A spline of every list: knots, poles and weights.
     *
     *  @return The spline.
     */
    auto EveryList() -> Spline
    {
      return {
          .Degree  = 2,
          .Knots   = {0, 0, 0, 1, 1, 1},
          .Poles   = {{.X = 1, .Y = 0}, {.X = 1, .Y = 1}, {.X = 0, .Y = 1}},
          .Weights = {1, 0.7071067811865476, 1},
      };
    }
  } // namespace

  TEST_CASE("text: each kind writes as its words")
  {
    CHECK(ToText(Line{.Start = {.X = 1, .Y = 2}, .End = {.X = 3, .Y = 4}}) == "line 1 2 3 4");
    CHECK(ToText(Ray{.Start = {.X = 0, .Y = 0}, .Way = {.X = 1, .Y = -1}}) == "ray 0 0 1 -1");
    CHECK(ToText(Xline{.Through = {.X = 0.5, .Y = 0}, .Way = {.X = 0, .Y = 1}})
          == "xline 0.5 0 0 1");
    CHECK(ToText(Circle({.X = 0, .Y = 0}, 1)) == "arc 0 0 1 0 360d");
    CHECK(ToText(Polyline{.Vertices = {{.Position = {.X = 0, .Y = 0}, .Bulge = 1}}, .Closed = true})
          == "polyline closed 1 0 0 1");
    CHECK(ToText(EveryList()) == "spline 2 6 0 0 0 1 1 1 3 1 0 1 1 0 1 3 1 0.7071067811865476 1");
  }

  TEST_CASE("text: every kind reads back exactly")
  {
    const std::array<Curve, 10> curves = {
        Line{.Start = {.X = 0.1, .Y = 5e6 + 0.1}, .End = {.X = 1e-300, .Y = -0.0}},
        Ray{.Start = {.X = 1, .Y = 2}, .Way = {.X = 0x1.0p-1074, .Y = 1e300}},
        Xline{.Through = {.X = -1e8, .Y = 1e8}, .Way = {.X = 1.0 / 3, .Y = 2.0 / 3}},
        Arc{.Center = {.X = 1, .Y = 2}, .Radius = 3, .Start = Degrees(30), .Sweep = Degrees(-270)},
        Ellipse{
            .Map = Rotation(Radians(0.3)) * Scaling(2, 1e-9), .Start = Radians(1), .Sweep = -Turn},
        Polyline{.Vertices = {{.Position = {.X = 0, .Y = 0}, .Bulge = 0.5},
                              {.Position = {.X = 10, .Y = 0}, .Bulge = -1e-6}},
                 .Closed   = false},
        Polyline{.Vertices = {}, .Closed = true},
        EveryList(),
        Bezier::Quadratic{
            .Start = {.X = 0, .Y = 0}, .Control = {.X = 0.5, .Y = 1}, .End = {.X = 1, .Y = 0}},
        Bezier::Cubic{.Start = {.X = 0, .Y = 0},
                      .Out   = {.X = 0, .Y = 1},
                      .In    = {.X = 1, .Y = 1},
                      .End   = {.X = 1, .Y = 0}},
    };
    for (const Curve &curve : curves)
    {
      CHECK(RoundTrips(curve));
    }
  }

  TEST_CASE("text: random doubles read back to the bit")
  {
    Draws draws{0x74657874U};
    for (int index = 0; index < 2000; ++index)
    {
      const auto bits  = draws.Next();
      const auto value = std::bit_cast<double>(bits);
      if (std::isnan(value))
      {
        continue;
      }
      const Line                 line = {.Start = {.X = value, .Y = draws.Near(5e6, 1e3)},
                                         .End   = {.X = -value, .Y = value}};
      const std::optional<Curve> back = FromText(ToText(line));
      REQUIRE(back.has_value());
      const Line read = std::get<Line>(*back);
      CHECK(std::bit_cast<std::uint64_t>(read.Start.X) == bits);
    }
  }

  TEST_CASE("text: infinities and NaN read back")
  {
    constexpr double Infinity = std::numeric_limits<double>::infinity();
    CHECK(RoundTrips(Line{.Start = {.X = Infinity, .Y = -Infinity}, .End = {}}));
    const std::string text =
        ToText(Line{.Start = {.X = std::numeric_limits<double>::quiet_NaN(), .Y = 0}, .End = {}});
    const std::optional<Curve> back = FromText(text);
    REQUIRE(back.has_value());
    CHECK(ToText(*back) == text);
  }

  TEST_CASE("text: blanks of any kind and length between words")
  {
    CHECK(FromText("  line\t1  2\n3 4 \r\n")
          == std::optional<Curve>{Line{.Start = {.X = 1, .Y = 2}, .End = {.X = 3, .Y = 4}}});
  }

  TEST_CASE("text: what is not a curve reads as none")
  {
    CHECK_FALSE(FromText(""));
    CHECK_FALSE(FromText("   "));
    CHECK_FALSE(FromText("circle 0 0 1"));
    CHECK_FALSE(FromText("line 1 2 3"));
    CHECK_FALSE(FromText("line 1 2 3 x"));
    CHECK_FALSE(FromText("line 1 2 3 4x"));
    CHECK_FALSE(FromText("line 1 2 3 4 5"));
    CHECK_FALSE(FromText("polyline shut 0"));
    CHECK_FALSE(FromText("polyline open x"));
    CHECK_FALSE(FromText("polyline open 2 0 0 0"));
    CHECK_FALSE(FromText("polyline open 99999999999999 0 0 0"));
    CHECK_FALSE(FromText("spline 2.5 0 0 0"));
    CHECK_FALSE(FromText("spline 2 1 0 0 0 1"));
  }
} // namespace Geoxy
