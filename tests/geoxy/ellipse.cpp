/***************************************************************************************************
 *  @file      tests/geoxy/ellipse.cpp
 *  @brief     Tests of ellipses.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "geoxy/ellipse.hpp"

#include "elliptic.hpp"
#include "scratch.hpp"

#include "geoxy/angle.hpp"
#include "geoxy/matrix.hpp"

#include <doctest/doctest.h>

#include <algorithm>
#include <bit>
#include <cstdint>
#include <optional>

namespace Geoxy
{
  TEST_CASE("ellipse: equal by map, start and sweep")
  {
    const Ellipse ellipse = {.Map = Scaling(2, 1), .Start = Radians(0.5), .Sweep = Turn};
    Ellipse       other   = ellipse;
    CHECK(ellipse == other);
    other.Map = Scaling(2, 2);
    CHECK_FALSE(ellipse == other);
    other       = ellipse;
    other.Start = Radians(0.25);
    CHECK_FALSE(ellipse == other);
    other       = ellipse;
    other.Sweep = -Turn;
    CHECK_FALSE(ellipse == other);
  }

  TEST_CASE("ellipse: the first pass of an angle is the least of every pass, with no list")
  {
    std::uint64_t state = 0x9E3779B97F4A7C15U;
    const auto    draw  = [&state](double low, double high) -> double
    {
      state = (state * 6364136223846793005U) + 1442695040888963407U;
      return low + ((high - low) * static_cast<double>(state >> 11U) * 0x1p-53);
    };
    const auto first = [](const Ellipse &ellipse, double angle) -> std::optional<double>
    {
      const List<double> passes = Passes(ellipse, angle);
      if (passes.empty())
      {
        return std::nullopt;
      }
      return std::ranges::min(passes);
    };
    int differ = 0;
    for (int round = 0; round < 20000; ++round)
    {
      const double  reach = round % 4 == 0 ? 40.0 : 7.0;
      const Ellipse inDeg = {.Map   = Scaling(2, 1),
                             .Start = Degrees(draw(-720, 720)),
                             .Sweep = Degrees(draw(-360, 360) * (round % 5 == 0 ? 3.0 : 1.0))};
      const Ellipse inRad = {.Map   = Scaling(2, 1),
                             .Start = Radians(draw(-reach, reach)),
                             .Sweep = Radians(draw(-reach, reach))};
      const double  angle = draw(-reach, reach);
      for (const Ellipse &ellipse : {inDeg, inRad})
      {
        const std::optional<double> earliest = Earliest(ellipse, angle);
        const std::optional<double> least    = first(ellipse, angle);
        const bool                  same =
            earliest.has_value() == least.has_value()
            && (!earliest
                || std::bit_cast<std::uint64_t>(*earliest) == std::bit_cast<std::uint64_t>(*least));
        differ += same ? 0 : 1;
      }
    }
    CHECK(differ == 0);
    const Ellipse none = {.Map = Scaling(2, 1), .Start = {}, .Sweep = Radians(0)};
    CHECK(!Earliest(none, 1.0));
    const Ellipse past = {.Map = Scaling(2, 1), .Start = {}, .Sweep = Degrees(90)};
    CHECK(!Earliest(past, Pi));
    CHECK(Earliest(past, 0.0) == 0.0);
  }
} // namespace Geoxy
