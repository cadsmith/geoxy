/***************************************************************************************************
 *  @file      tests/geoxy/extents.cpp
 *  @brief     Tests of the boxes of curves, against their own formulas and dense samples.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "geoxy/extents.hpp"

#include "curves.hpp"
#include "draws.hpp"
#include "real.hpp"

#include "geoxy/angle.hpp"
#include "geoxy/arc.hpp"
#include "geoxy/at.hpp"
#include "geoxy/bezier/cubic.hpp"
#include "geoxy/bezier/quadratic.hpp"
#include "geoxy/box.hpp"
#include "geoxy/curve.hpp"
#include "geoxy/domain.hpp"
#include "geoxy/ellipse.hpp"
#include "geoxy/line.hpp"
#include "geoxy/matrix.hpp"
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
#include <vector>

namespace Geoxy
{
  namespace
  {
    /**
     *  Whether a box is the box of a curve's dense samples over its domain, and of its points at
     *  whole numbers, where a polyline's vertices are: it holds them all, and reaches no further
     *  than a slack beyond them.
     *
     *  @tparam CurveT The curve's kind.
     *  @param  curve  The curve.
     *  @param  box    The box to check.
     *  @param  end    Where its domain ends, from 0.
     *  @return        Whether it is.
     */
    template <typename CurveT>
    auto Sampled(const CurveT &curve, const Box &box, double end) -> bool
    {
      constexpr int    Samples = 20'000;
      constexpr double Slack   = 1e-6;
      Box              samples;
      for (int index = 0; index <= Samples; ++index)
      {
        samples = Grown(samples, At(curve, end * index / Samples));
      }
      for (int whole = 0; whole <= static_cast<int>(end); ++whole)
      {
        samples = Grown(samples, At(curve, static_cast<double>(whole)));
      }
      return Contains(Padded(box, 1e-12), samples.Min)
          && Contains(Padded(box, 1e-12), samples.Max)
          && Contains(Padded(samples, Slack), box.Min)
          && Contains(Padded(samples, Slack), box.Max);
    }

    /**
     *  Where a curve's coordinate is greatest, near the best of some parameters, narrowed by
     *  thirds.
     *
     *  @tparam CurveT The curve's kind.
     *  @param  curve  The curve.
     *  @param  alongs The parameters, rising.
     *  @param  way    The way the coordinate is read: (1, 0) for x, (0, -1) for -y.
     *  @return        The parameter.
     */
    template <typename CurveT>
    auto Farthest(const CurveT &curve, const std::vector<double> &alongs, Vector way) -> double
    {
      const auto value = [&curve, way](double along) -> double
      {
        return Dot(At(curve, along) - Point{}, way);
      };
      std::size_t best = 0;
      for (std::size_t index = 0; index < alongs.size(); ++index)
      {
        best = value(alongs[index]) > value(alongs[best]) ? index : best;
      }
      double low  = alongs[best == 0 ? 0 : best - 1];
      double high = alongs[std::min(best + 1, alongs.size() - 1)];
      for (int round = 0; round < 200; ++round)
      {
        const double first  = low + ((high - low) / 3);
        const double second = high - ((high - low) / 3);
        if (value(first) < value(second))
        {
          low = first;
        }
        else
        {
          high = second;
        }
      }
      return low + ((high - low) / 2);
    }

    /**
     *  Whether a box holds every point At gives on a curve: at 2,001 even parameters, and about
     *  each extreme of each coordinate at 601 parameters a step of a double apart and 601 a
     *  billionth of the domain apart, where rounding pushes the points out furthest.
     *
     *  @tparam CurveT The curve's kind.
     *  @param  curve  The curve.
     *  @param  box    The box.
     *  @return        Whether it does.
     */
    template <typename CurveT>
    auto Holds(const CurveT &curve, const Box &box) -> bool
    {
      constexpr int Even   = 2000;
      const Range   domain = Domain(curve);
      const double  width  = domain.End - domain.Start;
      bool          held   = true;
      const auto    look   = [&curve, &box, &held, domain](double along) -> void
      {
        held =
            held && (along < domain.Start || along > domain.End || Contains(box, At(curve, along)));
      };
      std::vector<double> alongs;
      for (int index = 0; index <= Even; ++index)
      {
        alongs.push_back(domain.Start + (width * index / Even));
        look(alongs.back());
      }
      for (const Vector way : {Vector{.X = 1, .Y = 0}, Vector{.X = -1, .Y = 0},
                               Vector{.X = 0, .Y = 1}, Vector{.X = 0, .Y = -1}})
      {
        const double middle = Farthest(curve, alongs, way);
        const double step   = std::nextafter(middle, Infinite) - middle;
        for (int index = -300; index <= 300; ++index)
        {
          look(middle + (index * step));
          look(middle + (index * 1e-9 * width));
        }
      }
      return held;
    }
  } // namespace

  TEST_CASE("extents: a line's ends")
  {
    CHECK(Extents(Line{.Start = {.X = 3, .Y = -1}, .End = {.X = 1, .Y = 4}})
          == Box{.Min = {.X = 1, .Y = -1}, .Max = {.X = 3, .Y = 4}});
  }

  TEST_CASE("extents: a ray and an xline run without end where their way leaves an axis")
  {
    CHECK(Extents(Ray{.Start = {.X = 1, .Y = 2}, .Way = {.X = 1, .Y = -1}})
          == Box{.Min = {.X = 1, .Y = -Infinite}, .Max = {.X = Infinite, .Y = 2}});
    CHECK(Extents(Ray{.Start = {.X = 1, .Y = 2}, .Way = {.X = -1, .Y = 0}})
          == Box{.Min = {.X = -Infinite, .Y = 2}, .Max = {.X = 1, .Y = 2}});
    CHECK(Extents(Ray{.Start = {.X = 1, .Y = 2}, .Way = {}})
          == Box{.Min = {.X = 1, .Y = 2}, .Max = {.X = 1, .Y = 2}});
    CHECK(Extents(Xline{.Through = {.X = 1, .Y = 2}, .Way = {.X = 0, .Y = -3}})
          == Box{.Min = {.X = 1, .Y = -Infinite}, .Max = {.X = 1, .Y = Infinite}});
    CHECK(Extents(Xline{.Through = {.X = 1, .Y = 2}, .Way = {.X = 1, .Y = 1}})
          == Box{.Min = {.X = -Infinite, .Y = -Infinite}, .Max = {.X = Infinite, .Y = Infinite}});
  }

  TEST_CASE("extents: an arc's ends and the quarter turns it passes, exact")
  {
    CHECK(Extents(Circle({.X = 5e6, .Y = 1}, 2))
          == Box{.Min = {.X = 5e6 - 2, .Y = -1}, .Max = {.X = 5e6 + 2, .Y = 3}});
    const Arc quarter = {.Center = {}, .Radius = 1, .Start = {}, .Sweep = Degrees(90)};
    CHECK(Extents(quarter) == Box{.Min = {.X = 0, .Y = 0}, .Max = {.X = 1, .Y = 1}});
    const Arc top = {.Center = {}, .Radius = 1, .Start = Degrees(45), .Sweep = Degrees(90)};
    CHECK(Extents(top).Max.Y == 1.0);
    CHECK(Extents(top).Min.Y == doctest::Approx(std::sqrt(0.5)).epsilon(1e-15));
    const Arc back = {.Center = {}, .Radius = 1, .Start = Degrees(45), .Sweep = Degrees(-90)};
    CHECK(Extents(back).Max.X == 1.0);
    CHECK(Sampled(back, Extents(back), 1));
  }

  TEST_CASE("extents: an ellipse's extremes, turned and flat")
  {
    const Ellipse turned = {
        .Map   = Translation({.X = 5, .Y = 5}) * Rotation(Radians(0.6)) * Scaling(3, 1),
        .Start = {},
        .Sweep = Turn};
    const double across = std::hypot(3 * std::cos(0.6), std::sin(0.6));
    const double rising = std::hypot(3 * std::sin(0.6), std::cos(0.6));
    CHECK(Extents(turned).Max.X == doctest::Approx(5 + across).epsilon(1e-15));
    CHECK(Extents(turned).Min.Y == doctest::Approx(5 - rising).epsilon(1e-15));
    const Ellipse part = {
        .Map = Rotation(Radians(0.6)) * Scaling(3, 1), .Start = Radians(1), .Sweep = Radians(-2.5)};
    CHECK(Sampled(part, Extents(part), 1));
    const Ellipse flat = {.Map   = {.X = {.X = 2, .Y = 1}, .Y = {.X = 0, .Y = 0}, .Origin = {}},
                          .Start = {},
                          .Sweep = Turn};
    CHECK(Extents(flat) == Box{.Min = {.X = -2, .Y = -1}, .Max = {.X = 2, .Y = 1}});
  }

  TEST_CASE("extents: a polyline's vertices and its bulges' extremes")
  {
    const Polyline polyline = {
        .Vertices = {{.Position = {.X = 0, .Y = 0}, .Bulge = 1},
                     {.Position = {.X = 10, .Y = 0}, .Bulge = 0},
                     {.Position = {.X = 10, .Y = 10}, .Bulge = -0.3}},
        .Closed   = true,
    };
    CHECK(Extents(polyline).Min.Y == doctest::Approx(-5).epsilon(1e-15));
    CHECK(Sampled(polyline, Extents(polyline), 3));
    CHECK(IsEmpty(Extents(Polyline{})));
    CHECK(
        Extents(Polyline{.Vertices = {{.Position = {.X = 3, .Y = 4}, .Bulge = 1}}, .Closed = false})
        == Box{.Min = {.X = 3, .Y = 4}, .Max = {.X = 3, .Y = 4}});
  }

  TEST_CASE("extents: Bézier curves, where their derivatives are zero")
  {
    const Bezier::Quadratic quadratic = {
        .Start = {.X = 0, .Y = 0}, .Control = {.X = 1, .Y = 2}, .End = {.X = 2, .Y = 0}};
    const Box arch = Extents(quadratic);
    CHECK(arch.Min   == Point{.X = 0, .Y = 0});
    CHECK(arch.Max.X == 2);
    CHECK(arch.Max.Y >= 1);
    CHECK(arch.Max.Y <= 1 + 1e-14);
    CHECK(Sampled(Tests::Wave, Extents(Tests::Wave), 1));
    const Bezier::Cubic monotone = {.Start = {.X = 0, .Y = 0},
                                    .Out   = {.X = 1, .Y = 0.4},
                                    .In    = {.X = 2, .Y = 0.8},
                                    .End   = {.X = 3, .Y = 1.8}};
    CHECK(Extents(monotone) == Box{.Min = {.X = 0, .Y = 0}, .Max = {.X = 3, .Y = 1.8}});
    const Bezier::Cubic even = {.Start = {.X = 0, .Y = 0},
                                .Out   = {.X = 1, .Y = 1},
                                .In    = {.X = 2, .Y = 2},
                                .End   = {.X = 3, .Y = 3}};
    CHECK(Extents(even) == Box{.Min = {.X = 0, .Y = 0}, .Max = {.X = 3, .Y = 3}});
    const Bezier::Cubic flat = {.Start = {.X = 0, .Y = 0},
                                .Out   = {.X = 0, .Y = 0},
                                .In    = {.X = 0, .Y = 0},
                                .End   = {.X = 1, .Y = 1}};
    CHECK(Extents(flat) == Box{.Min = {.X = 0, .Y = 0}, .Max = {.X = 1, .Y = 1}});
    const Bezier::Cubic bent = {.Start = {.X = 0, .Y = 0},
                                .Out   = {.X = 1, .Y = 0},
                                .In    = {.X = 3, .Y = 0},
                                .End   = {.X = 6, .Y = 0}};
    CHECK(Extents(bent) == Box{.Min = {.X = 0, .Y = 0}, .Max = {.X = 6, .Y = 0}});
  }

  /**
   *  A spline's box is padded by what de Boor's fractions' rounding may move At's point, some tens
   *  of steps of a double of its poles' spread where its weights are near one another.
   */
  TEST_CASE("extents: a spline's, within some tens of steps of a double")
  {
    const Box wave    = Extents(Tests::Wave);
    const Box clamped = Extents(Tests::Clamped());
    CHECK(Distance(wave.Min, clamped.Min) <= 1e-13);
    CHECK(Distance(wave.Max, clamped.Max) <= 1e-13);
    CHECK(Sampled(Tests::Uneven(), Extents(Tests::Uneven()), 2));
    const Box quarter = Extents(Tests::Quarter());
    CHECK(quarter.Max.X == doctest::Approx(1.0).epsilon(1e-14));
    CHECK(quarter.Min.X == doctest::Approx(0.0).epsilon(1e-14));
  }

  /**
   *  Boxes were rounded points of the curve, and points At gives near an extreme poked past them:
   *  by up to 10 steps of a double on rational splines, 5 on others, 3 on bulges and 2 on ellipses.
   *  Every one now lies inside, near the origin, far from it, and over many turns.
   */
  TEST_CASE("extents: every point At gives lies in the box")
  {
    constexpr std::array<double, 3> Offsets = {0.0, 1e4, 1e8};
    Draws                           draws(1114);
    const auto                      uniform = [&draws](double low, double high) -> double
    {
      return draws.Near((low + high) / 2, (high - low) / 2);
    };
    for (int index = 0; index < Rounds(150); ++index)
    {
      CAPTURE(index);
      const double offset = Offsets.at(static_cast<std::size_t>(index % 3));
      const Arc    arc    = {
          .Center = {.X = offset + uniform(-1, 1), .Y = uniform(-1, 1)},
          .Radius = uniform(-3, 3),
          .Start  = index % 2 == 0 ? Degrees(uniform(-400, 400)) : Radians(uniform(-7, 7)),
          .Sweep  = index % 2 == 0 ? Degrees(uniform(-500, 500)) : Radians(uniform(-8, 8))};
      CHECK(Holds(arc, Extents(arc)));
      const Ellipse ellipse = {
          .Map   = {.X      = {.X = uniform(-3, 3), .Y = uniform(-3, 3)},
                    .Y      = {.X = uniform(-1, 1), .Y = uniform(-1, 1)},
                    .Origin = {.X = offset + uniform(-1, 1), .Y = uniform(-1, 1)}},
          .Start = Radians(uniform(-7, 7)),
          .Sweep = Radians(uniform(-8, 8))};
      CHECK(Holds(ellipse, Extents(ellipse)));
      const double   bulge = std::pow(10.0, uniform(-6, 6)) * (index % 2 == 0 ? 1 : -1);
      const double   chord = std::pow(10.0, uniform(-3, 3));
      const Polyline bowed = {
          .Vertices = {{.Position = {.X = (10 * offset) + uniform(-1, 1), .Y = uniform(-1, 1)},
                        .Bulge    = bulge},
                       {.Position = {.X = (10 * offset) + (chord * uniform(-1, 1)),
                                     .Y = chord * uniform(-1, 1)},
                        .Bulge    = 0}},
          .Closed   = false};
      CHECK(Holds(bowed, Extents(bowed)));
      const auto near = [&uniform, offset]() -> Point
      {
        return {.X = offset + uniform(-2, 2), .Y = uniform(-2, 2)};
      };
      const Bezier::Quadratic quadratic = {.Start = near(), .Control = near(), .End = near()};
      CHECK(Holds(quadratic, Extents(quadratic)));
      const Bezier::Cubic cubic = {.Start = near(), .Out = near(), .In = near(), .End = near()};
      CHECK(Holds(cubic, Extents(cubic)));
      Spline    spline = {.Degree = 1 + (index % 5), .Knots = {}, .Poles = {}, .Weights = {}};
      const int count  = spline.Degree + 2 + (index % 9);
      double    knot   = 0.0;
      for (int pole = 0; pole < count; ++pole)
      {
        spline.Poles.push_back({.X = offset + (uniform(-1, 1) * (1 + pole)), .Y = uniform(-1, 1)});
      }
      for (int at = 0; at < count + spline.Degree + 1; ++at)
      {
        knot += uniform(0.2, 1.5);
        spline.Knots.push_back(knot);
      }
      if (index % 2 == 0)
      {
        const auto degree = static_cast<std::size_t>(spline.Degree);
        const auto poles  = static_cast<std::size_t>(count);
        for (std::size_t at = 0; at < degree; ++at)
        {
          spline.Knots[at] = spline.Knots[degree];
        }
        for (std::size_t at = poles + 1; at < spline.Knots.size(); ++at)
        {
          spline.Knots[at] = spline.Knots[poles];
        }
      }
      REQUIRE(Valid(spline));
      CHECK(Holds(spline, Extents(spline)));
      for (int pole = 0; pole < count; ++pole)
      {
        spline.Weights.push_back(std::pow(10.0, uniform(-1, 1)));
      }
      CHECK(Holds(spline, Extents(spline)));
    }
  }

  /**
   *  A bulge's extremes were worked out along its chord by trigonometry: a semicircle's leftmost
   *  point came a step of a double inside -1.5, and a bulge of 1e3 sat 478 steps outside. From the
   *  bulge itself, rationally, they are where they are to a few steps.
   */
  TEST_CASE("extents: bulges from the bulge itself, with no trigonometry")
  {
    const Polyline closed = {.Vertices = {{.Position = {.X = 0, .Y = 0}, .Bulge = 0.4},
                                          {.Position = {.X = 4, .Y = 0}, .Bulge = 0},
                                          {.Position = {.X = 4, .Y = 3}, .Bulge = -0.3},
                                          {.Position = {.X = 0, .Y = 3}, .Bulge = 1}},
                             .Closed   = true};
    const Box      box    = Extents(closed);
    CHECK(box.Min.X <= -1.5);
    CHECK(box.Min.X >= -1.5 - 1e-14);
    CHECK(box.Min.Y <= -0.8);
    CHECK(box.Min.Y >= -0.8 - 1e-14);
    CHECK(Holds(closed, box));
    const Polyline big   = {.Vertices = {{.Position = {.X = 0, .Y = 0}, .Bulge = 1e3},
                                         {.Position = {.X = 1, .Y = 0}, .Bulge = 0}},
                            .Closed   = false};
    const Box      round = Extents(big);
    CHECK(round.Min.X <= -249.50025);
    CHECK(round.Min.X >= -249.50025 - 1e-12);
    CHECK(round.Max.X >= 250.50025);
    CHECK(round.Max.X <= 250.50025 + 1e-12);
    CHECK(round.Min.Y <= -500);
    CHECK(round.Min.Y >= -500 - 1e-12);
    CHECK(round.Max.Y >= 0);
    CHECK(round.Max.Y <= 1e-12);
    CHECK(Holds(big, round));
  }

  /**
   *  The pad a spline's box took for At's rounding left out how de Boor's fractions round: each
   *  mix's fraction off by its own few steps moves the point off the curve, by as much more as its
   *  weights differ. Here At's point came 345 steps of a double past the box.
   */
  TEST_CASE("extents: a spline of widely differing weights holds every point At gives")
  {
    const Spline spline = {.Degree  = 3,
                           .Knots   = {-71.74412736540022, -3.888060606003146, -0.27093200435802334,
                                       -0.010212948010281908, -0.005768818692651284,
                                       0.0749449935399895, 0.6879910458493327, 68.47742681972983},
                           .Poles   = {{.X = -1.6157020797409751, .Y = -7.065520695961191},
                                       {.X =  -3.820573335285868, .Y =  5.799337391642768},
                                       {.X = -3.2672687695523974, .Y = 3.3682978453244683},
                                       {.X = 2.555838438331825, .Y = -5.894382430316315}},
                           .Weights = {434.09492371917764, 0.08462361423600537, 0.00125580775998278,
                                       860.6455422309197}};
    REQUIRE(Valid(spline));
    const Box box = Extents(spline);
    CHECK(Contains(box, At(spline, -0x1.38c1643462c7fp-7)));
    CHECK(Holds(spline, box));
  }

  /**
   *  A cubic's turns are where a quadratic of its controls' differences is zero, and its
   *  discriminant squares them: below about 2^-540 the squares came to nothing, past 2^512 to
   *  infinity, no turn was found, and the box held the ends alone, missing the curve's bulge by
   *  most of its width.
   */
  TEST_CASE("extents: a cubic Bézier curve far below or above the normal doubles")
  {
    for (const int power : {-545, 520})
    {
      const double        unit  = std::ldexp(1.0, power);
      const Bezier::Cubic cubic = {.Start = {.X = 37 * unit, .Y = 35 * unit},
                                   .Out   = {.X = -39 * unit, .Y = 29 * unit},
                                   .In    = {.X = -39 * unit, .Y = 21 * unit},
                                   .End   = {.X = 36 * unit, .Y = -31 * unit}};
      const Box           box   = Extents(cubic);
      CHECK(box.Min.X < -20 * unit);
      CHECK(Holds(cubic, box));
    }
  }

  /**
   *  A spline of a degree past seven is worked out in cells for every degree: its box holds every
   *  point At gives, as a cubic's does.
   */
  TEST_CASE("extents: a spline of a high degree holds every point At gives")
  {
    Spline spline = {.Degree = 9, .Knots = {}, .Poles = {}, .Weights = {}};
    for (int index = 0; index < 12; ++index)
    {
      spline.Poles.push_back({.X = index * 0.5, .Y = index % 3 == 1 ? 1.5 : -0.75});
    }
    spline.Knots.assign(10, 0.0);
    spline.Knots.insert(spline.Knots.end(), {0.25, 0.5});
    spline.Knots.insert(spline.Knots.end(), 10, 1.0);
    REQUIRE(Valid(spline));
    const Box box = Extents(spline);
    CHECK(Holds(spline, box));
    CHECK(box.Max.Y < 1.5);
    CHECK(box.Max.Y > 0.0);
  }

  TEST_CASE("extents: any curve, by its kind")
  {
    CHECK(Extents(Curve{Circle({}, 1)}) == Box{.Min = {.X = -1, .Y = -1}, .Max = {.X = 1, .Y = 1}});
    CHECK(Extents(Curve{Tests::Clamped()}) == Extents(Tests::Clamped()));
  }
} // namespace Geoxy
