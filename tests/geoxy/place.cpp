/***************************************************************************************************
 *  @file      tests/geoxy/place.cpp
 *  @brief     Tests of places on a curve.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "geoxy/place.hpp"

#include <doctest/doctest.h>

namespace Geoxy
{
  TEST_CASE("place: equal when both the point and the parameter are, exactly")
  {
    const Place place = {.At = {.X = 1, .Y = 2}, .Along = 0.5};
    CHECK(place       == Place{.At = {.X = 1, .Y = 2}, .Along = 0.5});
    CHECK_FALSE(place == Place{.At = {.X = 1, .Y = 3}, .Along = 0.5});
    CHECK_FALSE(place == Place{.At = {.X = 1, .Y = 2}, .Along = 0.25});
  }
} // namespace Geoxy
