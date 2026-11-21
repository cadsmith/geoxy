/***************************************************************************************************
 *  @file      tests/geoxy/part.cpp
 *  @brief     Tests of a curve's stretch between two places: its ends At's own points to the bit,
 *             its points the curve's reparametrized, and the trims, extensions and splits made of
 *             it.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "geoxy/part.hpp"

#include "curves.hpp"
#include "real.hpp"

#include "geoxy/accuracy.hpp"
#include "geoxy/angle.hpp"
#include "geoxy/arc.hpp"
#include "geoxy/at.hpp"
#include "geoxy/bezier/cubic.hpp"
#include "geoxy/bezier/quadratic.hpp"
#include "geoxy/curve.hpp"
#include "geoxy/domain.hpp"
#include "geoxy/ellipse.hpp"
#include "geoxy/join.hpp"
#include "geoxy/length.hpp"
#include "geoxy/line.hpp"
#include "geoxy/matrix.hpp"
#include "geoxy/point.hpp"
#include "geoxy/polyline.hpp"
#include "geoxy/ray.hpp"
#include "geoxy/spline.hpp"
#include "geoxy/vector.hpp"
#include "geoxy/xline.hpp"

#include <doctest/doctest.h>

#include <array>
#include <cmath>
#include <cstdint>
#include <optional>
#include <variant>
#include <vector>

namespace Geoxy
{
  namespace
  {
    /**
     *  Splitmix64: numbers that look random, the same on every run.
     */
    struct Splitmix
    {
      std::uint64_t State = 0x9E3779B97F4A7C15; /**< Where the sequence is. */

      /**
       *  A double from 0 to 1.
       *
       *  @return The double.
       */
      auto Next() -> double
      {
        State             += 0x9E3779B97F4A7C15;
        std::uint64_t mix  = State;
        mix                = (mix ^ (mix >> 30U)) * 0xBF58476D1CE4E5B9;
        mix                = (mix ^ (mix >> 27U)) * 0x94D049BB133111EB;
        return static_cast<double>((mix ^ (mix >> 31U)) >> 11U) * 0x1p-53;
      }
    };

    /**
     *  Whether a part of a curve ends where the curve is at the stretch's ends, to the bit, and
     *  runs over its points between, reparametrized, within a distance.
     *
     *  @tparam ShapeT  The curve's kind.
     *  @param  shape   The curve.
     *  @param  stretch The stretch.
     *  @param  slack   How far apart the points between may be.
     *  @return         Whether it does.
     */
    template <typename ShapeT>
    auto Spanned(const ShapeT &shape, Range stretch, double slack) -> bool
    {
      const ShapeT part = Part(shape, stretch);
      bool kept = At(part, 0) == At(shape, stretch.Start) && At(part, 1) == At(shape, stretch.End);
      for (const double along : {0.125, 0.5, 0.875})
      {
        const double place = stretch.Start + ((stretch.End - stretch.Start) * along);
        kept               = kept && Distance(At(part, along), At(shape, place)) <= slack;
      }
      return kept;
    }

    /**
     *  How many of many random stretches of a curve end where the curve is, to the bit, and run
     *  over its points between.
     *
     *  @tparam ShapeT The curve's kind.
     *  @param  shape  The curve.
     *  @param  reach  How far past its domain the stretches may run.
     *  @param  slack  How far apart the points between may be.
     *  @return        How many do not.
     */
    template <typename ShapeT>
    auto Missed(const ShapeT &shape, double reach, double slack) -> int
    {
      Splitmix mixer;
      int      missed = 0;
      for (int index = 0; index < 2000; ++index)
      {
        const double start  = (mixer.Next() * (1 + (2 * reach))) - reach;
        const double end    = (mixer.Next() * (1 + (2 * reach))) - reach;
        missed             += Spanned(shape, {.Start = start, .End = end}, slack) ? 0 : 1;
      }
      return missed;
    }

    /**
     *  How many of many random stretches of an arc or an ellipse do not start where it is, to the
     *  bit, nor end there, nor run over its points between.
     *
     *  @tparam ShapeT An arc or an ellipse.
     *  @param  shape  The curve.
     *  @param  grid   Whether the stretches' ends are whole 1024ths, so that an arc's angles held
     *                 as binary fractions of degrees add exactly, and its part must end to the bit.
     *  @param  slack  How far apart the points may be, where they need not be one.
     *  @return        How many do not.
     */
    template <typename ShapeT>
    auto Swerved(const ShapeT &shape, bool grid, double slack) -> int
    {
      Splitmix mixer;
      int      missed = 0;
      for (int index = 0; index < 2000; ++index)
      {
        const double random  = (mixer.Next() * 3) - 1;
        const double other   = (mixer.Next() * 3) - 1;
        const double start   = grid ? std::round(random * 1024) / 1024 : random;
        const double end     = grid ? std::round(other * 1024) / 1024 : other;
        const ShapeT part    = Part(shape, {.Start = start, .End = end});
        const double off     = Distance(At(part, 1), At(shape, end));
        const double middle  = Distance(At(part, 0.5), At(shape, (start * 0.5) + (end * 0.5)));
        const bool   landed  = At(part, 0) == At(shape, start) && off <= (grid ? 0.0 : slack);
        missed              += landed && middle <= slack ? 0 : 1;
      }
      return missed;
    }

    /**
     *  Whether a spline's part ends where the spline is, to the bit, and its points between are the
     *  spline's at the same parameter within a distance.
     *
     *  @param  spline  The spline.
     *  @param  stretch The stretch.
     *  @param  slack   How far apart the points between may be.
     *  @return         Whether they are.
     */
    auto Kept(const Spline &spline, Range stretch, double slack) -> bool
    {
      const Spline part = Part(spline, stretch);
      bool         kept = Valid(part)
                       && At(part, stretch.Start) == At(spline, stretch.Start)
                       && At(part, stretch.End) == At(spline, stretch.End);
      for (const double along : {0.1, 0.37, 0.5, 0.81})
      {
        const double place = stretch.Start + ((stretch.End - stretch.Start) * along);
        kept               = kept && Distance(At(part, place), At(spline, place)) <= slack;
      }
      return kept;
    }
  } // namespace

  TEST_CASE("part: a line between two places, run on and the other way round")
  {
    const Line line = {.Start = {.X = 5e6, .Y = 1}, .End = {.X = 5e6 + 10, .Y = 3}};
    CHECK(Part(line, {.Start = 0.25, .End = 0.75})
          == Line{.Start = At(line, 0.25), .End = At(line, 0.75)});
    CHECK(Spanned(line, {.Start = 0.3, .End = 0.1}, 1e-9));
    CHECK(Spanned(line, {.Start = -2, .End = 3.5}, 1e-9));
    CHECK(Part(line, {.Start = 0, .End = 1}) == line);
  }

  TEST_CASE("part: a ray's and an xline's line, as a line, a ray or an xline")
  {
    const Ray   ray   = {.Start = {.X = 1, .Y = 2}, .Way = {.X = 3, .Y = 4}};
    const Xline xline = {.Through = {.X = 1, .Y = 2}, .Way = {.X = 3, .Y = 4}};
    CHECK(Part(ray, {.Start = 1, .End = 2})
          == Curve{Line{.Start = {.X = 4, .Y = 6}, .End = {.X = 7, .Y = 10}}});
    CHECK(Part(ray, {.Start = 1, .End = Infinite})
          == Curve{Ray{.Start = {.X = 4, .Y = 6}, .Way = ray.Way}});
    CHECK(Part(ray, {.Start = -1, .End = -Infinite})
          == Curve{Ray{.Start = {.X = -2, .Y = -2}, .Way = -ray.Way}});
    CHECK(Part(xline, {.Start = -Infinite, .End = 1})
          == Curve{Ray{.Start = {.X = 4, .Y = 6}, .Way = -xline.Way}});
    CHECK(Part(xline, {.Start = Infinite, .End = 1})
          == Curve{Ray{.Start = {.X = 4, .Y = 6}, .Way = xline.Way}});
    CHECK(Part(xline, {.Start = -Infinite, .End = Infinite}) == Curve{xline});
    CHECK(Part(ray, {.Start = Infinite, .End = -Infinite})
          == Curve{Xline{.Through = ray.Start, .Way = -ray.Way}});
  }

  TEST_CASE("part: an arc and an ellipse end where At puts them, to the bit")
  {
    const Arc degrees = {.Center = {.X = 5e6, .Y = -3e6},
                         .Radius = 7.25,
                         .Start  = Degrees(33.25),
                         .Sweep  = Degrees(-251.5)};
    CHECK(Swerved(degrees, true, 1e-8)  == 0);
    CHECK(Swerved(degrees, false, 1e-8) == 0);
    const Arc radians = {
        .Center = {.X = 0.5, .Y = 2}, .Radius = 3, .Start = Radians(-2.1), .Sweep = Radians(5.9)};
    CHECK(Swerved(radians, false, 4e-14) == 0);
    const Arc part = Part(radians, {.Start = 0.25, .End = 0.75});
    CHECK(part.Center            == radians.Center);
    CHECK(part.Sweep.InRadians() == doctest::Approx(2.95));
    CHECK(Spanned(Circle({.X = 1, .Y = 1}, 2), {.Start = 0.75, .End = 1.25}, 1e-15));

    const Ellipse ellipse = {.Map   = Rotation(Radians(0.4), {.X = 5e6, .Y = 1}) * Scaling(4, 0.5),
                             .Start = Radians(1.1),
                             .Sweep = Radians(-4.4)};
    CHECK(Swerved(ellipse, false, 1e-8) == 0);
    const Ellipse whole = {.Map = Scaling(3, 1), .Start = Degrees(10.5), .Sweep = Degrees(360)};
    CHECK(Swerved(whole, true, 4e-14) == 0);
  }

  TEST_CASE("part: Bézier curves by blossoming, their ends where At puts them")
  {
    const Bezier::Quadratic quadratic = {.Start   = {.X = 5e6, .Y = 0},
                                         .Control = {.X = 5e6 + 1, .Y = 3},
                                         .End     = {.X = 5e6 + 4, .Y = 0.5}};
    CHECK(Missed(quadratic, 0.5, 1e-8)    == 0);
    CHECK(Missed(Tests::Wave, 0.5, 1e-13) == 0);
    CHECK(Part(Tests::Wave, {.Start = 0, .End = 1}) == Tests::Wave);
    const Bezier::Cubic half = Part(Tests::Wave, {.Start = 0, .End = 0.5});
    CHECK(half.Out == Point{.X = 0.5, .Y = 1.5});
    CHECK(half.In  == Point{.X = 1.25, .Y = 1});
  }

  TEST_CASE("part: a polyline's vertices between, its cut bulges on their circles")
  {
    const Polyline polyline = {.Vertices = {{.Position = {.X = 0, .Y = 0}, .Bulge = 0.5},
                                            {.Position = {.X = 4, .Y = 0}, .Bulge = 0},
                                            {.Position = {.X = 4, .Y = 3}, .Bulge = -1},
                                            {.Position = {.X = 0, .Y = 3}, .Bulge = 1e-7}},
                               .Closed   = true};
    const std::optional<Polyline> part = Part(polyline, {.Start = 0.5, .End = 2.5});
    REQUIRE(part);
    REQUIRE(part->Vertices.size()    == 4);
    CHECK(part->Vertices[0].Position == At(polyline, 0.5));
    CHECK(part->Vertices[0].Bulge    == doctest::Approx(std::tan(std::atan(0.5) / 2)));
    CHECK(part->Vertices[1]          == polyline.Vertices[1]);
    CHECK(part->Vertices[2].Position == polyline.Vertices[2].Position);
    CHECK(part->Vertices[3].Position == At(polyline, 2.5));
    CHECK_FALSE(part->Closed);
    for (const double along : {0.0, 0.3, 0.5, 1.0})
    {
      CHECK(Distance(At(*part, along), At(polyline, 0.5 + (along / 2))) <= 1e-14);
      CHECK(Distance(At(*part, 1 + along), At(polyline, 1 + along)) <= 1e-14);
      CHECK(Distance(At(*part, 2 + along), At(polyline, 2 + (along / 2))) <= 1e-14);
    }

    const std::optional<Polyline> whole = Part(polyline, Domain(polyline));
    REQUIRE(whole);
    CHECK(whole->Vertices.size()      == 5);
    CHECK(whole->Vertices[3]          == polyline.Vertices[3]);
    CHECK(whole->Vertices[4].Position == polyline.Vertices[0].Position);

    const std::optional<Polyline> back = Part(polyline, {.Start = 3.5, .End = 1});
    REQUIRE(back);
    CHECK(back->Vertices.front().Position == At(polyline, 3.5));
    CHECK(back->Vertices.back().Position  == At(polyline, 1));
    CHECK(back->Vertices[1].Bulge         == doctest::Approx(1));

    const std::optional<Polyline> straight = Part(polyline, {.Start = 3.25, .End = 3.75});
    REQUIRE(straight);
    CHECK(straight->Vertices[0].Bulge == doctest::Approx(0));

    const std::optional<Polyline> point = Part(polyline, {.Start = 2, .End = 2});
    REQUIRE(point);
    CHECK(point->Vertices.size()      == 2);
    CHECK(point->Vertices[1].Position == polyline.Vertices[2].Position);
  }

  TEST_CASE("part: a polyline run on, but not through a whole turn of a bulge")
  {
    const Polyline polyline           = {.Vertices = {{.Position = {.X = 0, .Y = 0}, .Bulge = 0.25},
                                                      {.Position = {.X = 4, .Y = 0}, .Bulge = 0},
                                                      {.Position = {.X = 4, .Y = 3}, .Bulge = 0}},
                                         .Closed   = false};
    const std::optional<Polyline> run = Part(polyline, {.Start = -1, .End = 5});
    REQUIRE(run);
    CHECK(run->Vertices.front().Position == At(polyline, -1));
    CHECK(run->Vertices.back().Position  == At(polyline, 5));
    CHECK(Distance(At(*run, 0.25), At(polyline, -0.5)) <= 1e-14);
    CHECK_FALSE(Part(polyline, {.Start = -7, .End = 1}));
    CHECK(Part(polyline, {.Start = 2.5, .End = 4})->Vertices.size() == 2);
    CHECK(Part(polyline, {.Start = -3, .End = -2})->Vertices.size() == 2);
    CHECK(Part(Polyline{}, {.Start = 0, .End = 1}) == Polyline{});
  }

  TEST_CASE("part: a spline's knots between, its ends where At puts them, rational too")
  {
    Splitmix mixer;
    for (const Spline &spline : {Tests::Clamped(), Tests::Quarter(), Tests::Uneven()})
    {
      const Range domain = Domain(spline);
      int         missed = 0;
      for (int index = 0; index < 500; ++index)
      {
        const double start = domain.Start + ((domain.End - domain.Start) * mixer.Next());
        const double end   = domain.Start + ((domain.End - domain.Start) * mixer.Next());
        if (start < end)
        {
          missed += Kept(spline, {.Start = start, .End = end}, 1e-13) ? 0 : 1;
        }
      }
      CHECK(missed == 0);
      CHECK(Kept(spline, domain, 1e-13));
    }
    const Spline uneven = Tests::Uneven();
    CHECK(Kept(uneven, {.Start = 0.3, .End = 1.1}, 1e-13));
    CHECK(Part(uneven, {.Start = 0.3, .End = 1.1}).Knots
          == std::vector<double>{0.3, 0.3, 0.3, 0.3, 1.1, 1.1, 1.1, 1.1});
    CHECK(Kept(Tests::Clamped(), {.Start = -0.5, .End = 1.5}, 1e-13));
    CHECK(Kept(Tests::Quarter(), {.Start = -0.2, .End = 1.3}, 1e-13));
    CHECK_FALSE(Valid(Part(uneven, {.Start = -0.2, .End = 2.3})));
    CHECK(Part(Tests::Clamped(), {.Start = 0, .End = 1}).Weights.empty());
  }

  TEST_CASE("part: a spline the other way, of no length, not finite or not valid")
  {
    const Spline uneven = Tests::Uneven();
    const Spline back   = Part(uneven, {.Start = 1.5, .End = 0.2});
    CHECK(Domain(back)   == Range{.Start = -1.5, .End = -0.2});
    CHECK(At(back, -1.5) == At(uneven, 1.5));
    CHECK(At(back, -0.2) == At(uneven, 0.2));
    CHECK_FALSE(Valid(Part(uneven, {.Start = 1, .End = 1})));
    CHECK_FALSE(Valid(Part(uneven, {.Start = -Infinite, .End = 1})));
    CHECK_FALSE(Valid(Part(uneven, {.Start = 0, .End = Infinite})));
    CHECK_FALSE(Valid(Part(Spline{}, {.Start = 0, .End = 1})));
  }

  TEST_CASE("part: any curve, trimmed, extended and split")
  {
    const Arc arc = {
        .Center = {.X = 0, .Y = 0}, .Radius = 2, .Start = Degrees(0), .Sweep = Degrees(90)};
    CHECK(Part(Curve{arc}, {.Start = 0, .End = 0.5}) == Curve{Part(arc, {.Start = 0, .End = 0.5})});
    CHECK(Part(Curve{Ray{.Start = {}, .Way = {.X = 1, .Y = 0}}}, {.Start = 1, .End = 2})
          == Curve{Line{.Start = {.X = 1, .Y = 0}, .End = {.X = 2, .Y = 0}}});
    const Polyline polyline = {.Vertices = {{.Position = {.X = 0, .Y = 0}, .Bulge = 1},
                                            {.Position = {.X = 2, .Y = 0}, .Bulge = 0}},
                               .Closed   = false};
    CHECK(Part(Curve{polyline}, {.Start = 0, .End = 1}) == Curve{polyline});
    CHECK_FALSE(Part(Curve{polyline}, {.Start = -2, .End = 1}));

    CHECK(Trim(arc, 0.5, Keep::Start) == Curve{Part(arc, {.Start = 0, .End = 0.5})});
    CHECK(Trim(arc, 0.5, Keep::End)   == Curve{Part(arc, {.Start = 0.5, .End = 1})});
    const Xline xline = {.Through = {}, .Way = {.X = 1, .Y = 0}};
    CHECK(Trim(xline, 2, Keep::Start)
          == Curve{Ray{.Start = {.X = 2, .Y = 0}, .Way = {.X = -1, .Y = 0}}});

    const Line line = {.Start = {.X = 0, .Y = 0}, .End = {.X = 2, .Y = 0}};
    CHECK(Extend(line, 1.5) == Curve{Line{.Start = {.X = 0, .Y = 0}, .End = {.X = 3, .Y = 0}}});
    CHECK(Extend(line, -1)  == Curve{Line{.Start = {.X = -2, .Y = 0}, .End = {.X = 2, .Y = 0}}});
    CHECK(Extend(line, 0.5) == Curve{line});
    CHECK(Extend(arc, 2)    == Curve{Part(arc, {.Start = 0, .End = 2})});
    CHECK(Extend(Ray{.Start = {}, .Way = {.X = 1, .Y = 0}}, -1)
          == Curve{Ray{.Start = {.X = -1, .Y = 0}, .Way = {.X = 1, .Y = 0}}});
  }

  TEST_CASE("part: split at places, in order, once each, the lengths summing to the whole")
  {
    const std::array<double, 6> places = {0.7, 0.2, 0.7, -1, 0.45, 1};
    const Polyline polyline            = {.Vertices = {{.Position = {.X = 0, .Y = 0}, .Bulge = 0.5},
                                                       {.Position = {.X = 4, .Y = 0}, .Bulge = 0},
                                                       {.Position = {.X = 4, .Y = 3}, .Bulge = -1}},
                                          .Closed   = true};
    const std::array<Curve, 4> curves  = {
        Curve{Arc{.Center = {.X = 5e6, .Y = 1},
                  .Radius = 3,
                  .Start  = Degrees(20),
                  .Sweep  = Degrees(300)}},
        Curve{Tests::Wave},
        Curve{Tests::Uneven()},
        Curve{polyline},
    };
    for (const Curve &curve : curves)
    {
      const std::vector<Curve> parts = Split(curve, places);
      CHECK(parts.size() == (Domain(curve).End > 1 ? 5 : 4));
      double length = 0;
      for (const Curve &part : parts)
      {
        length += Length(part, Accuracy{});
      }
      CHECK(length == doctest::Approx(Length(curve, Accuracy{})).epsilon(1e-13));
      CHECK(At(parts.front(), Domain(parts.front()).Start) == At(curve, Domain(curve).Start));
    }
    const std::vector<Curve>   parts  = Split(Curve{polyline}, std::array{0.5, 2.0, 2.5});
    const std::optional<Curve> linked = Join(parts, {.Length = 1e-9, .Radians = 1e-9});
    REQUIRE(linked);
    REQUIRE(std::holds_alternative<Polyline>(*linked));
    CHECK(std::get<Polyline>(*linked).Vertices.size() == 5);
    CHECK(std::get<Polyline>(*linked).Closed);
    for (const auto &[along, place] : std::array<std::array<double, 2>, 5>{
             {{0.5, 0.25}, {1.5, 0.75}, {2.5, 1.5}, {3.5, 2.25}, {4.5, 2.75}}})
    {
      CHECK(Distance(At(*linked, along), At(polyline, place)) <= 1e-14);
    }
    const std::vector<Curve> rays =
        Split(Xline{.Through = {}, .Way = {.X = 1, .Y = 0}}, std::array{1.0});
    REQUIRE(rays.size() == 2);
    CHECK(rays[0]       == Curve{Ray{.Start = {.X = 1, .Y = 0}, .Way = {.X = -1, .Y = 0}}});
    CHECK(rays[1]       == Curve{Ray{.Start = {.X = 1, .Y = 0}, .Way = {.X = 1, .Y = 0}}});
  }
} // namespace Geoxy
