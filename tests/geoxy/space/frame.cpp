/***************************************************************************************************
 *  @file      tests/geoxy/space/frame.cpp
 *  @brief     Tests of frames and of AutoCAD's arbitrary axis.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "geoxy/space/frame.hpp"

#include "../draws.hpp"

#include "geoxy/matrix.hpp"
#include "geoxy/point.hpp"
#include "geoxy/space/vector.hpp"
#include "geoxy/vector.hpp"

#include <doctest/doctest.h>

#include <cmath>

namespace Geoxy::Space
{
  TEST_CASE("frame: equality axis by axis")
  {
    const Frame frame;
    CHECK_FALSE(frame == Frame{.X = North});
    CHECK_FALSE(frame == Frame{.Y = East});
    CHECK_FALSE(frame == Frame{.Z = East});
  }

  TEST_CASE("frame: facing up is the world's own axes")
  {
    CHECK(Facing(Up) == Frame{});
    CHECK(Facing({.X = 0, .Y = 0, .Z = 1e300}) == Frame{});
  }

  /**
   *  An extrusion of 0,0,-1, as a mirrored entity has: its x runs west, its y north, so its plan is
   *  mirrored.
   */
  TEST_CASE("frame: facing down mirrors the plan")
  {
    const Frame down = *Facing({.X = 0, .Y = 0, .Z = -1});
    CHECK(down.X == Vector{.X = -1, .Y = 0, .Z = 0});
    CHECK(down.Y == North);
    const Matrix plan = Plan(down, 0);
    CHECK(plan * Point{.X = 2, .Y = 3} == Point{.X = -2, .Y = 3});
    CHECK(Determinant(plan)            == -1.0);
  }

  /**
   *  What AutoCAD was measured to do through its console: the test is on the unit normal, both
   *  components strictly under 1/64.
   */
  TEST_CASE("frame: the arbitrary axis takes world y near the pole and world z elsewhere")
  {
    const Frame nearPole = *Facing({.X = 0.1, .Y = 0, .Z = 10});
    CHECK(nearPole.X.Y == 0.0);
    CHECK(nearPole.X.X > 0.99);
    const Frame elsewhere = *Facing({.X = 0.02, .Y = 0, .Z = 0.5});
    CHECK(elsewhere.X == Vector{.X = 0, .Y = 1, .Z = 0});
    const Frame under = *Facing({.X = 0.0156249, .Y = 0, .Z = 0.99987793});
    CHECK(under.X.Y == 0.0);
    const Frame over = *Facing({.X = 0.0156251, .Y = 0, .Z = 0.99987793});
    CHECK(over.X == Vector{.X = 0, .Y = 1, .Z = 0});
  }

  TEST_CASE("frame: every frame is unit, square and right-handed")
  {
    Draws draws{0x6672616DU};
    for (int index = 0; index < 2000; ++index)
    {
      const Vector normal = {.X = draws.Near(0, 1), .Y = draws.Near(0, 1), .Z = draws.Near(0, 1)};
      const Frame  frame  = *Facing(normal);
      CHECK(Length(frame.X) == doctest::Approx(1.0).epsilon(1e-15));
      CHECK(Length(frame.Y) == doctest::Approx(1.0).epsilon(1e-15));
      CHECK(std::abs(Dot(frame.X, frame.Y)) <= 1e-15);
      CHECK(std::abs(Dot(frame.X, frame.Z)) <= 1e-15);
      CHECK(Dot(Cross(frame.X, frame.Y), frame.Z) == doctest::Approx(1.0).epsilon(1e-15));
    }
  }

  TEST_CASE("frame: no plane without a normal")
  {
    CHECK_FALSE(Facing({}));
    CHECK_FALSE(Facing({.X = std::nan(""), .Y = 0, .Z = 1}));
  }

  TEST_CASE("frame: a plan at an elevation, and a plane seen edge on")
  {
    CHECK(Plan(Frame{}, 5) == Matrix{});
    const Frame  tilted = *Facing({.X = 1, .Y = 0, .Z = 1});
    const Matrix raised = Plan(tilted, 2);
    CHECK(raised.Origin.X == doctest::Approx(std::sqrt(2.0)));
    CHECK(raised.Origin.Y == 0.0);
    const Matrix edge = Plan(*Facing(East), 0);
    CHECK(edge.X        == Geoxy::Vector{.X = 0, .Y = 1});
    CHECK(edge.Y        == Geoxy::Vector{.X = 0, .Y = 0});
    CHECK(Stretch(edge) == 1.0);
  }
} // namespace Geoxy::Space
