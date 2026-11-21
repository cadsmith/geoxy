/***************************************************************************************************
 *  @file      tests/geoxy/version.cpp
 *  @brief     Tests of the version.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "geoxy/version.hpp"

#include <doctest/doctest.h>

namespace Geoxy
{
  TEST_CASE("version: the release is the project's")
  {
    CHECK(Version() == "0.1.0");
  }
} // namespace Geoxy
