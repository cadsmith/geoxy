/***************************************************************************************************
 *  @file      tests/geoxy/transform.cpp
 *  @brief     Tests of each curve under an affine map: its points mapped, its kind kept where the
 *             map allows, and two maps in turn their product.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "geoxy/transform.hpp"

#include "curves.hpp"

#include "geoxy/angle.hpp"
#include "geoxy/arc.hpp"
#include "geoxy/at.hpp"
#include "geoxy/bezier/quadratic.hpp"
#include "geoxy/curve.hpp"
#include "geoxy/domain.hpp"
#include "geoxy/ellipse.hpp"
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
#include <optional>
#include <variant>

namespace Geoxy
{
  namespace
  {
    /**
     *  A turn, a scale and a move, as a block's insertion composes them.
     *
     *  @return The map.
     */
    auto Placed() -> Matrix
    {
      return Translation({.X = 5e6, .Y = -2e6}) * Rotation(Degrees(30)) * Scaling(2.5);
    }

    /**
     *  A mirror across a slanted line.
     *
     *  @return The map.
     */
    auto Mirrored() -> Matrix
    {
      return Mirror({.X = 1, .Y = 2}, {.X = 4, .Y = 6}).value_or(Matrix{});
    }

    /**
     *  A map that shears and scales unevenly.
     */
    constexpr Matrix Sheared = {
        .X = {.X = 2, .Y = 0.5}, .Y = {.X = -0.3, .Y = 1.2}, .Origin = {.X = -7, .Y = 3}};

    /**
     *  Whether a curve mapped keeps its parameter: each of its points the curve's own mapped,
     *  within a part of the scale, an endless curve's from -2 to 2.
     *
     *  @param  curve  The curve.
     *  @param  mapped The curve mapped.
     *  @param  map    The map.
     *  @return        Whether it does.
     */
    auto Kept(const Curve &curve, const Curve &mapped, const Matrix &map) -> bool
    {
      const Range  domain = Domain(curve);
      const double start  = std::isinf(domain.Start) ? -2.0 : domain.Start;
      const double end    = std::isinf(domain.End) ? 2.0 : domain.End;
      bool         kept   = true;
      for (const double along : {0.0, 0.2, 0.5, 0.7, 1.0})
      {
        const double place = start + ((end - start) * along);
        const Point  image = map * At(curve, place);
        kept               = kept
                          && Distance(At(mapped, place), image) <= 1e-14 * std::max(1.0, Length(image - Point{}));
      }
      return kept;
    }
  } // namespace

  TEST_CASE("transform: lines, rays, xlines, Bézier curves and splines map their points")
  {
    const std::array<Curve, 7> curves = {
        Curve{Line{.Start = {.X = 1, .Y = 2}, .End = {.X = -3, .Y = 4}}},
        Curve{Ray{.Start = {.X = 1, .Y = 2}, .Way = {.X = 0.5, .Y = -1}}},
        Curve{Xline{.Through = {.X = 1, .Y = 2}, .Way = {.X = 0.5, .Y = -1}}},
        Curve{Bezier::Quadratic{
            .Start = {.X = 0, .Y = 0}, .Control = {.X = 1, .Y = 3}, .End = {.X = 4, .Y = 0}}},
        Curve{Tests::Wave},
        Curve{Tests::Uneven()},
        Curve{Tests::Quarter()},
    };
    for (const Matrix &map : {Placed(), Mirrored(), Sheared})
    {
      for (const Curve &curve : curves)
      {
        const Curve mapped = Transformed(curve, map);
        CHECK(mapped.index() == curve.index());
        CHECK(Kept(curve, mapped, map));
      }
    }
    CHECK(Transformed(Tests::Uneven(), Sheared).Weights == Tests::Uneven().Weights);
    CHECK(Transformed(Ray{.Start = {}, .Way = {.X = 1, .Y = 0}}, Translation({.X = 3, .Y = 0})).Way
          == Vector{.X = 1, .Y = 0});
  }

  TEST_CASE(
      "transform: an arc stays an arc under a similarity, in degrees exactly by quarter turns")
  {
    const Arc arc = {
        .Center = {.X = 3, .Y = 4}, .Radius = 2, .Start = Degrees(10), .Sweep = Degrees(135)};
    const Curve turned = Transformed(arc, Rotation(Degrees(90), {.X = 1, .Y = 1}));
    REQUIRE(std::holds_alternative<Arc>(turned));
    CHECK(std::get<Arc>(turned)
          == Arc{.Center = {.X = -2, .Y = 3},
                 .Radius = 2,
                 .Start  = Degrees(100),
                 .Sweep  = Degrees(135)});
    const Curve flipped = Transformed(arc, *Mirror({.X = 0, .Y = 0}, {.X = 1, .Y = 0}));
    REQUIRE(std::holds_alternative<Arc>(flipped));
    CHECK(std::get<Arc>(flipped)
          == Arc{.Center = {.X = 3, .Y = -4},
                 .Radius = 2,
                 .Start  = Degrees(-10),
                 .Sweep  = Degrees(-135)});
    for (const Matrix &map : {Placed(), Mirrored()})
    {
      const Curve mapped = Transformed(arc, map);
      CHECK(std::holds_alternative<Arc>(mapped));
      CHECK(Kept(arc, mapped, map));
      CHECK(Kept(arc, Transformed(Transformed(Curve{arc}, map), Placed()), Placed() * map));
    }
    const Arc backward = {
        .Center = {.X = -1, .Y = 0}, .Radius = -0.5, .Start = Radians(1), .Sweep = Radians(-4)};
    CHECK(Kept(backward, Transformed(backward, Mirrored() * Placed()), Mirrored() * Placed()));
  }

  TEST_CASE("transform: an arc becomes the ellipse of its circle under any other map")
  {
    const Arc arc = {
        .Center = {.X = 3, .Y = 4}, .Radius = 2, .Start = Degrees(10), .Sweep = Degrees(-300)};
    const Curve mapped = Transformed(arc, Sheared);
    REQUIRE(std::holds_alternative<Ellipse>(mapped));
    CHECK(std::get<Ellipse>(mapped).Start == arc.Start);
    CHECK(std::get<Ellipse>(mapped).Sweep == arc.Sweep);
    CHECK(Kept(arc, mapped, Sheared));
    CHECK(Kept(arc, Transformed(Transformed(Curve{arc}, Sheared), Placed()), Placed() * Sheared));
    const Ellipse ellipse = {.Map = Scaling(3, 1), .Start = Degrees(45), .Sweep = Degrees(90)};
    CHECK(Transformed(ellipse, Sheared).Map == Sheared * ellipse.Map);
    CHECK(Kept(ellipse, Transformed(Curve{ellipse}, Sheared), Sheared));
  }

  TEST_CASE("transform: a polyline keeps its bulges under a similarity, negated by a mirror")
  {
    const Polyline polyline = {.Vertices = {{.Position = {.X = 0, .Y = 0}, .Bulge = 0.4},
                                            {.Position = {.X = 4, .Y = 0}, .Bulge = 0},
                                            {.Position = {.X = 4, .Y = 3}, .Bulge = -1.5}},
                               .Closed   = true};
    const Curve    turned   = Transformed(polyline, Placed());
    REQUIRE(std::holds_alternative<Polyline>(turned));
    CHECK(std::get<Polyline>(turned).Vertices[0].Bulge == doctest::Approx(0.4));
    CHECK(Kept(polyline, turned, Placed()));
    const Curve flipped = Transformed(polyline, Mirrored());
    REQUIRE(std::holds_alternative<Polyline>(flipped));
    CHECK(std::get<Polyline>(flipped).Vertices[2].Bulge == doctest::Approx(1.5));
    CHECK(Kept(polyline, flipped, Mirrored()));

    const Polyline straight = {.Vertices = {{.Position = {.X = 0, .Y = 0}, .Bulge = 0x1p-61},
                                            {.Position = {.X = 4, .Y = 0}, .Bulge = 0},
                                            {.Position = {.X = 4, .Y = 3}, .Bulge = 5}},
                               .Closed   = false};
    const Curve    sheared  = Transformed(straight, Sheared);
    REQUIRE(std::holds_alternative<Polyline>(sheared));
    CHECK(Kept(straight, sheared, Sheared));
  }

  TEST_CASE("transform: a bulged polyline under any other map is the spline of its image exactly")
  {
    for (const double bulge : {0.2, 0.5, 1.0, 1.5, 3.0, -0.7})
    {
      const Polyline polyline = {.Vertices = {{.Position = {.X = 0, .Y = 0}, .Bulge = bulge},
                                              {.Position = {.X = 4, .Y = 0}, .Bulge = 0},
                                              {.Position = {.X = 4, .Y = 3}, .Bulge = -bulge}},
                                 .Closed   = true};
      const Curve    mapped   = Transformed(polyline, Sheared);
      REQUIRE(std::holds_alternative<Spline>(mapped));
      const auto &spline = std::get<Spline>(mapped);
      CHECK(Valid(spline));
      CHECK(spline.Degree == 2);
      CHECK(At(spline, Domain(spline).Start) == Sheared * polyline.Vertices[0].Position);
      CHECK(At(spline, Domain(spline).End)   == Sheared * polyline.Vertices[0].Position);
      const std::optional<Matrix> back = Inverse(Sheared);
      REQUIRE(back);
      const Range domain = Domain(spline);
      for (int step = 0; step <= 64; ++step)
      {
        const double along = domain.Start + ((domain.End - domain.Start) * step / 64.0);
        const Point  point = *back * At(spline, along);
        CHECK(Distance(Nearest(polyline, point).At, point) <= 1e-14);
      }
    }
  }
} // namespace Geoxy
