/***************************************************************************************************
 *  @file      tests/geoxy/divide.cpp
 *  @brief     Tests of places at equal lengths: each place's length from the start against Length.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "geoxy/divide.hpp"

#include "curves.hpp"

#include "geoxy/angle.hpp"
#include "geoxy/arc.hpp"
#include "geoxy/at.hpp"
#include "geoxy/curve.hpp"
#include "geoxy/domain.hpp"
#include "geoxy/length.hpp"
#include "geoxy/line.hpp"
#include "geoxy/place.hpp"
#include "geoxy/point.hpp"
#include "geoxy/polyline.hpp"
#include "geoxy/ray.hpp"
#include "geoxy/spline.hpp"

#include <doctest/doctest.h>

#include <cmath>
#include <cstddef>
#include <limits>
#include <vector>

namespace Geoxy
{
  namespace
  {
    /**
     *  Every place's length from the curve's start against the one asked, and its point At's.
     *
     *  @param curve  The curve.
     *  @param places The places.
     *  @param step   The length between them.
     */
    void Spaced(const Curve &curve, const std::vector<Place> &places, double step)
    {
      const Range domain = Domain(curve);
      for (std::size_t index = 0; index < places.size(); ++index)
      {
        const Place &place = places[index];
        CHECK(place.At == At(curve, place.Along));
        const double length = Length(curve, {.Start = domain.Start, .End = place.Along}, {});
        CHECK(std::abs(length - (step * static_cast<double>(index + 1))) <= 1e-12 * (1 + length));
      }
    }
  } // namespace

  TEST_CASE("divide: equal parts of every kind")
  {
    const Polyline           path   = {.Vertices = {{.Position = {.X = 0, .Y = 0}, .Bulge = 0.4},
                                                    {.Position = {.X = 5, .Y = 0}},
                                                    {.Position = {.X = 5, .Y = 7}, .Bulge = -1},
                                                    {.Position = {.X = 1, .Y = 9}}},
                                       .Closed   = true};
    const std::vector<Curve> curves = {
        Line{.Start = {.X = 0, .Y = 0}, .End = {.X = 3, .Y = 4}},
        Arc{.Center = {}, .Radius = 2, .Start = Degrees(10), .Sweep = Degrees(-200)}, path,
        Tests::Wave, Tests::Uneven()};
    for (const Curve &curve : curves)
    {
      const double             length = Length(curve, {});
      const std::vector<Place> seven  = Divide(curve, 7, {});
      REQUIRE(seven.size() == 6);
      Spaced(curve, seven, length / 7);
      const std::vector<Place> measured = Measure(curve, 0.75, {});
      CHECK(measured.size() == static_cast<std::size_t>(std::floor(length / 0.75)));
      Spaced(curve, measured, 0.75);
    }
    const std::vector<Place> line = Divide(Line{.Start = {}, .End = {.X = 4, .Y = 0}}, 4, {});
    CHECK(line[1].At == Point{.X = 2, .Y = 0});
    CHECK(Measure(Line{.Start = {}, .End = {.X = 4, .Y = 0}}, 1, {}).size() == 4);
  }

  TEST_CASE("divide: what has no places")
  {
    const Line   line = {.Start = {}, .End = {.X = 4, .Y = 0}};
    const Ray    ray  = {.Start = {}, .Way = {.X = 1, .Y = 0}};
    const double nan  = std::numeric_limits<double>::quiet_NaN();
    CHECK(Divide(line, 1, {}).empty());
    CHECK(Divide(line, 0, {}).empty());
    CHECK(Divide(line, MostPlaces + 1, {}).empty());
    CHECK(Divide(ray, 4, {}).empty());
    CHECK(Divide(Spline{}, 4, {}).empty());
    CHECK(Measure(line, 0, {}).empty());
    CHECK(Measure(line, -1, {}).empty());
    CHECK(Measure(line, nan, {}).empty());
    CHECK(Measure(line, std::numeric_limits<double>::infinity(), {}).empty());
    CHECK(Measure(line, 1e-7, {}).empty());
    CHECK(Measure(line, 5, {}).empty());
    CHECK(Measure(ray, 1, {}).empty());
  }
} // namespace Geoxy
