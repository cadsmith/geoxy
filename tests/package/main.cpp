/***************************************************************************************************
 *  @file      tests/package/main.cpp
 *  @brief     A program built against the installed package.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "geoxy/geoxy.hpp"

#include <print>

auto main() -> int
{
  const Geoxy::Matrix turn   = Geoxy::Rotation(Geoxy::Degrees(90));
  const Geoxy::Point  turned = turn * Geoxy::Point{.X = 1, .Y = 0};
  std::println("geoxy {}: {} {}", Geoxy::Version(), turned.X, turned.Y);
  return 0;
}
