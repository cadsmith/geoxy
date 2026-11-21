/***************************************************************************************************
 *  @file      tests/geoxy/flatten.cpp
 *  @brief     Tests of curves as chords, against their own points.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "geoxy/flatten.hpp"

#include "curves.hpp"
#include "draws.hpp"
#include "kernel/piece.hpp"
#include "waves.hpp"

#include "geoxy/angle.hpp"
#include "geoxy/arc.hpp"
#include "geoxy/at.hpp"
#include "geoxy/bezier/cubic.hpp"
#include "geoxy/bezier/quadratic.hpp"
#include "geoxy/box.hpp"
#include "geoxy/curve.hpp"
#include "geoxy/derivative.hpp"
#include "geoxy/domain.hpp"
#include "geoxy/ellipse.hpp"
#include "geoxy/extents.hpp"
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
#include <cstddef>
#include <functional>
#include <limits>
#include <numbers>
#include <vector>

namespace Geoxy
{
  namespace
  {
    /**
     *  The distance from a point to the nearest of a run of chords.
     *
     *  @param  point  The point.
     *  @param  points The chords' ends.
     *  @return        The distance.
     */
    auto Stray(Point point, const std::vector<Point> &points) -> double
    {
      double nearest = std::numeric_limits<double>::infinity();
      for (std::size_t index = 0; index + 1 < points.size(); ++index)
      {
        const Vector chord  = points[index + 1] - points[index];
        const double square = Dot(chord, chord);
        const double along =
            square > 0.0 ? std::clamp(Dot(point - points[index], chord) / square, 0.0, 1.0) : 0.0;
        nearest = std::min(nearest, Distance(point, Lerp(points[index], points[index + 1], along)));
      }
      return nearest;
    }

    /**
     *  Whether a curve's chords hold it within a tolerance: every dense sample of it is near one.
     *
     *  @param  curve     The curve.
     *  @param  tolerance The tolerance.
     *  @return           Whether they do.
     */
    auto Holds(const Curve &curve, double tolerance) -> bool
    {
      const std::vector<Point> points  = Flatten(curve, {.Error = tolerance});
      const Range              domain  = Domain(curve);
      constexpr int            Samples = 4000;
      for (int index = 0; index <= Samples; ++index)
      {
        const double along = domain.Start + ((domain.End - domain.Start) * index / Samples);
        if (Stray(At(curve, along), points) > tolerance * (1 + 1e-9))
        {
          return false;
        }
      }
      return Distance(points.front(), At(curve, domain.Start)) <= 1e-12
          && Distance(points.back(), At(curve, domain.End))    <= 1e-12;
    }

    /**
     *  How far a point is from a chord.
     *
     *  @param  point The point.
     *  @param  start The chord's start.
     *  @param  end   Its end.
     *  @return       The distance.
     */
    auto Off(Point point, Point start, Point end) -> double
    {
      const Vector chord  = end - start;
      const double square = Dot(chord, chord);
      const double along =
          square > 0.0 ? std::clamp(Dot(point - start, chord) / square, 0.0, 1.0) : 0.0;
      return Distance(point, Lerp(start, end, along));
    }

    /**
     *  How far an ellipse's chords stray from it over the tolerance, each from its own stretch of
     *  the curve: each point's parameter found again on the unit circle through the inverse map,
     *  the turn nearest the point before.
     *
     *  @pre              The map has an inverse.
     *  @param  ellipse   The ellipse.
     *  @param  tolerance The tolerance.
     *  @return           The worst distance over the tolerance.
     */
    auto Strayed(const Ellipse &ellipse, double tolerance) -> double
    {
      const std::vector<Point> points = Flatten(ellipse, {.Error = tolerance});
      const Matrix             unit   = *Inverse(ellipse.Map);
      const double             start  = ellipse.Start.InRadians();
      const double             sweep  = ellipse.Sweep.InRadians();
      const double             whole  = 2.0 * std::numbers::pi;
      std::vector<double>      alongs = {0.0};
      double                   angle  = start;
      for (std::size_t index = 1; index + 1 < points.size(); ++index)
      {
        const Point  seen   = unit * points[index];
        const double turned = std::atan2(seen.Y, seen.X);
        angle               = turned + (whole * std::round((angle - turned) / whole));
        alongs.push_back((angle - start) / sweep);
      }
      alongs.push_back(1.0);
      double worst = 0.0;
      for (std::size_t index = 0; index + 1 < points.size(); ++index)
      {
        for (int step = 1; step < 32; ++step)
        {
          const double along = alongs[index] + ((alongs[index + 1] - alongs[index]) * step / 32.0);
          worst = std::max(worst, Off(At(ellipse, along), points[index], points[index + 1]));
        }
      }
      return worst / tolerance;
    }

    /**
     *  The fewest chords an ellipse could take within a tolerance, as their count tends to it: the
     *  integral of sqrt(curvature / (8 tolerance)) along it.
     *
     *  @param  ellipse   The ellipse.
     *  @param  tolerance The tolerance.
     *  @return           The count, not below one.
     */
    auto Fewest(const Ellipse &ellipse, double tolerance) -> double
    {
      constexpr int Samples = 20'000;
      const double  start   = ellipse.Start.InRadians();
      const double  sweep   = ellipse.Sweep.InRadians();
      double        sum     = 0.0;
      for (int index = 0; index < Samples; ++index)
      {
        const double angle = start + (sweep * (index + 0.5) / Samples);
        const Vector way   = (ellipse.Map.Y * std::cos(angle)) - (ellipse.Map.X * std::sin(angle));
        const Vector turn  = (ellipse.Map.X * -std::cos(angle)) - (ellipse.Map.Y * std::sin(angle));
        sum += std::sqrt(std::abs(Cross(way, turn)) / (8.0 * tolerance * Length(way)));
      }
      return std::max(1.0, sum * std::abs(sweep) / Samples);
    }

    /**
     *  The length of the diagonal of a curve's box.
     *
     *  @param  curve The curve.
     *  @return       The length.
     */
    auto Diagonal(const Curve &curve) -> double
    {
      const Box box = Extents(curve);
      return Distance(box.Min, box.Max);
    }

    /**
     *  How far a Bézier curve's chords stray from it over the tolerance, each from its own stretch
     *  of the curve: each point's parameter found again as its nearest place, the curve not
     *  crossing itself.
     *
     *  @tparam CurveT    The curve's kind.
     *  @param  curve     The curve.
     *  @param  tolerance The tolerance.
     *  @return           The worst distance over the tolerance.
     */
    template <typename CurveT>
    auto Wanders(const CurveT &curve, double tolerance) -> double
    {
      const std::vector<Point> points = Flatten(curve, {.Error = tolerance});
      std::vector<double>      alongs;
      alongs.reserve(points.size());
      for (const Point point : points)
      {
        alongs.push_back(Nearest(curve, point).Along);
      }
      double worst = 0.0;
      for (std::size_t index = 0; index + 1 < points.size(); ++index)
      {
        for (int step = 1; step < 32; ++step)
        {
          const double along = alongs[index] + ((alongs[index + 1] - alongs[index]) * step / 32.0);
          worst = std::max(worst, Off(At(curve, along), points[index], points[index + 1]));
        }
      }
      return worst / tolerance;
    }

    /**
     *  The fewest chords a Bézier curve could take within a tolerance, as their count tends to it.
     *
     *  @tparam CurveT    The curve's kind.
     *  @param  curve     The curve.
     *  @param  tolerance The tolerance.
     *  @return           The count, not below one.
     */
    template <typename CurveT>
    auto Least(const CurveT &curve, double tolerance) -> double
    {
      constexpr int Samples = 20'000;
      double        sum     = 0.0;
      for (int index = 0; index < Samples; ++index)
      {
        const double along = (index + 0.5) / Samples;
        const Vector way   = Derivative(curve, along);
        const Vector turn =
            (Derivative(curve, along + 1e-7) - Derivative(curve, along - 1e-7)) / 2e-7;
        sum += std::sqrt(std::abs(Cross(way, turn)) / (8.0 * tolerance * Length(way)));
      }
      return std::max(1.0, sum / Samples);
    }

    /**
     *  A circle seen from a point near its arcs' points, which they are moved by exactly.
     */
    struct Ring
    {
      Point  Origin = {};  /**< The point, within a factor of two of every coordinate far off. */
      Vector Center = {};  /**< The circle's center, less the point.                           */
      double Radius = 0.0; /**< The radius.                                                    */
    };

    /**
     *  How far a run of chords between points of one arc of a circle strays from it, over the
     *  tolerance, exactly enough however far from the origin: each chord's ends moved by the ring's
     *  point, which leaves them exact, and for each chord the circle's farthest from its line, the
     *  radius less the center's distance from it, by a cross product of products made exact with
     *  the fused multiply-add; and how far each point lies off the circle.
     *
     *  @param  points    The chords' ends.
     *  @param  ring      The circle, seen from its point.
     *  @param  tolerance The tolerance.
     *  @return           The worst over the tolerance.
     */
    auto Sagged(const std::vector<Point> &points, const Ring &ring, double tolerance) -> double
    {
      const auto cross = [](Vector lhs, Vector rhs) -> double
      {
        const double product = rhs.X * lhs.Y;
        const double error   = std::fma(-rhs.X, lhs.Y, product);
        return std::fma(lhs.X, rhs.Y, -product) + error;
      };
      double worst = 0.0;
      for (std::size_t index = 0; index < points.size(); ++index)
      {
        const Vector start = (points[index] - ring.Origin) - ring.Center;
        worst              = std::max(worst, std::abs(Length(start) - ring.Radius));
        if (index + 1 < points.size())
        {
          const Vector end   = (points[index + 1] - ring.Origin) - ring.Center;
          const Vector chord = end - start;
          worst = std::max(worst, ring.Radius - (std::abs(cross(start, end)) / Length(chord)));
        }
      }
      return worst / tolerance;
    }
  } // namespace

  TEST_CASE("flatten: a line is its two ends")
  {
    CHECK(Flatten(Line{.Start = {.X = 1, .Y = 2}, .End = {.X = 3, .Y = 4}}, {.Error = 1})
          == std::vector<Point>{{.X = 1, .Y = 2}, {.X = 3, .Y = 4}});
  }

  TEST_CASE(
      "flatten: an arc takes even chords, as few as the tolerance allows, a quarter turn at most")
  {
    const std::vector<Point> circle = Flatten(Circle({.X = 0, .Y = 0}, 1), {.Error = 1e-3});
    const double             step   = 4 * std::asin(std::sqrt(1e-3 / 2));
    CHECK(circle.size()  == static_cast<std::size_t>(std::ceil(2 * std::numbers::pi / step)) + 1);
    CHECK(circle.front() == Point{.X = 1, .Y = 0});
    CHECK(circle.back()  == Point{.X = 1, .Y = 0});
    CHECK(Flatten(Circle({}, 1), {.Error = 100}).size() == 5);
    CHECK(Flatten(Arc{.Center = {}, .Radius = 0, .Start = {}, .Sweep = Degrees(30)}, {.Error = 1})
              .size()
          == 2);
    CHECK(Flatten(Circle({}, 1), {.Error = 0}).size()      == MostChords + 1);
    CHECK(Flatten(Circle({}, 1), {.Error = 1e-300}).size() == MostChords + 1);
    CHECK(Flatten(Circle({}, 1), {.Error = 1e-12}).size()  == MostChords + 1);
    CHECK(Holds(Arc{.Center = {.X = 5e6, .Y = 1},
                    .Radius = 30,
                    .Start  = Degrees(10),
                    .Sweep  = Degrees(-200)},
                1e-4));
  }

  TEST_CASE("flatten: an ellipse, turned and part of one")
  {
    CHECK(
        Holds(Ellipse{.Map = Translation({.X = 3, .Y = 3}) * Rotation(Radians(0.4)) * Scaling(5, 1),
                      .Start = Radians(1),
                      .Sweep = Radians(-4)},
              1e-3));
  }

  /**
   *  Even steps sized for its largest radius took twice the chords an ellipse needs, and ninety
   *  times on a flat one: in each quarter between its axes' directions it now takes runs of even
   *  steps, each sized by the exact sagitta of an affine circle where the run is nearer the tip.
   */
  TEST_CASE("flatten: an ellipse's chords each stray no more than the tolerance, near the fewest")
  {
    const std::array ellipses{
        Ellipse{.Map = Rotation(Radians(0.3)) * Scaling(3, 1), .Start = {}, .Sweep = Turn},
        Ellipse{.Map = Translation({.X = 1, .Y = -2}) * Rotation(Radians(0.7)) * Scaling(10, 9.99),
                .Start = Radians(0.5),
                .Sweep = Radians(1.234)},
        Ellipse{.Map = Scaling(1, 1e-3), .Start = {}, .Sweep = Turn},
        Ellipse{.Map   = Rotation(Radians(0.2)) * Scaling(1e6, 1),
                .Start = Radians(0.1),
                .Sweep = Radians(2.5)},
        Ellipse{.Map = Scaling(2, 1), .Start = Radians(0.25), .Sweep = Radians(19.85)},
        Ellipse{.Map   = Translation({.X = 1e6, .Y = 1e6}) * Rotation(Radians(1)) * Scaling(5, 2),
                .Start = Radians(1),
                .Sweep = Radians(-4)},
        Ellipse{.Map   = {.X      = {.X = 2, .Y = 1},
                          .Y      = {.X = 0.5, .Y = 1},
                          .Origin = {.X = 0.5, .Y = 0.25}},
                .Start = Radians(0.3),
                .Sweep = Radians(4)},
        Ellipse{.Map = Scaling(-1.5, 0.8), .Start = Degrees(-40), .Sweep = Degrees(-300)},
        Ellipse{.Map = Scaling(1, 1e-9), .Start = Radians(-1), .Sweep = Radians(5)},
    };
    for (std::size_t index = 0; index < ellipses.size(); ++index)
    {
      const Ellipse &ellipse = ellipses.at(index);
      for (const double share : {1e-2, 1e-4, 1e-6, 1e-9})
      {
        CAPTURE(index);
        CAPTURE(share);
        const double             tolerance = share * Diagonal(ellipse);
        const std::vector<Point> points    = Flatten(ellipse, {.Error = tolerance});
        if (points.size() <= MostChords)
        {
          CHECK(Strayed(ellipse, tolerance) <= 1.0);
          CHECK(Holds(ellipse, tolerance));
        }
        if (index < 8 && share > 1e-8 && share < 1e-3)
        {
          CHECK(static_cast<double>(points.size() - 1) <= (1.6 * Fewest(ellipse, tolerance)) + 8.0);
        }
      }
    }
    const Ellipse flat = ellipses[3];
    CHECK(Flatten(flat, {.Error = 1e-6 * Diagonal(flat)}).size() == 4);
  }

  TEST_CASE("flatten: random ellipses take a fifth more chords than the fewest, or less")
  {
    Draws  draws{2026};
    double chords = 0.0;
    double fewest = 0.0;
    for (int round = 0; round < Rounds(64); ++round)
    {
      const double  turn      = draws.Near(0, 4);
      const double  ratio     = std::pow(10.0, draws.Near(-2, 2));
      const double  along     = draws.Near(0, 4);
      const Vector  major     = {.X = std::cos(turn), .Y = std::sin(turn)};
      const Vector  minor     = Vector{.X = -major.Y, .Y = major.X} * ratio;
      const Ellipse ellipse   = {.Map   = {.X = (major * std::cos(along)) + (minor * std::sin(along)),
                                           .Y = (minor * std::cos(along)) - (major * std::sin(along)),
                                           .Origin = {.X = draws.Near(0, 1), .Y = draws.Near(0, 1)}},
                                 .Start = Radians(draws.Near(0, 3)),
                                 .Sweep = round % 2 == 0 ? Turn : Radians(draws.Near(0, 6))};
      const double  tolerance = 1e-6 * Diagonal(ellipse);
      CAPTURE(round);
      CHECK(Strayed(ellipse, tolerance) <= 1.0);
      chords += static_cast<double>(Flatten(ellipse, {.Error = tolerance}).size() - 1);
      fewest += Fewest(ellipse, tolerance);
    }
    CHECK(chords <= 1.2 * fewest);
  }

  /**
   *  An ellipse whose start lies 6e-18 of a turn from a tip of its major axis took a chord from its
   *  start to that tip, of no length: the piece of the quarter before the tip, which the curve
   *  along is far within the tolerance, is now left to the next quarter.
   */
  TEST_CASE("flatten: an ellipse starting a hair from a tip takes no chord of no length")
  {
    const Ellipse ellipse   = {.Map   = {.X      = {.X = -3.5918387122349595, .Y = 0.5356727257601229},
                                         .Y      = {.X = -0.1599635315569264, .Y = -1.0726011939037554},
                                         .Origin = {.X = -2.713684879223665, .Y = 5.364095706285006}},
                               .Start = Radians(0.0),
                               .Sweep = Radians(6.283185307179586)};
    const double  tolerance = 0.0075800560041268646;
    const std::vector<Point> points = Flatten(ellipse, {.Error = tolerance});
    CHECK(points.size() == 49);
    CHECK(std::ranges::adjacent_find(points) == points.end());
    CHECK(points.front() == At(ellipse, 0.0));
    CHECK(points.back()  == At(ellipse, 1.0));
    CHECK(Strayed(ellipse, tolerance) <= 1.0);

    /**
     *  A piece as short that strays past the share, to a fine tolerance, keeps its chord.
     */
    const Ellipse before = {.Map = Scaling(2, 1), .Start = Radians(-1e-7), .Sweep = Turn};
    CHECK(Flatten(before, {.Error = 1e-8}).at(1) == Point{.X = 2, .Y = 0});
    CHECK(Flatten(before, {.Error = 0.1}).at(1) != Point{.X = 2, .Y = 0});
  }

  TEST_CASE("flatten: an ellipse past what chords can hold takes the most even ones")
  {
    const Ellipse wide = {.Map = Scaling(2, 1), .Start = {}, .Sweep = Turn};
    for (const double error : {0.0, -1.0, std::numeric_limits<double>::quiet_NaN(), 1e-20, 1e-11})
    {
      CAPTURE(error);
      CHECK(Flatten(wide, {.Error = error}).size() == MostChords + 1);
    }
    CHECK(Flatten(Ellipse{.Map = Scaling(2, 1), .Start = Degrees(90), .Sweep = Degrees(-90)},
                  {.Error = 1e-11})
              .size()
          == MostChords + 1);
    CHECK(Flatten(Ellipse{.Map = Scaling(2, 1), .Start = {}, .Sweep = Radians(1e6)}, {.Error = 1})
              .size()
          == MostChords + 1);
    const Ellipse nowhere = {.Map   = Scaling(2, 1),
                             .Start = Radians(std::numeric_limits<double>::infinity()),
                             .Sweep = Turn};
    CHECK(Flatten(nowhere, {.Error = 1e-3}).size() == MostChords + 1);
    const Ellipse broken = {
        .Map = Scaling(1, std::numeric_limits<double>::quiet_NaN()), .Start = {}, .Sweep = Turn};
    CHECK(Flatten(broken, {.Error = 1e-3}).size() == MostChords + 1);
  }

  TEST_CASE("flatten: an ellipse of no sweep, no size or no width")
  {
    const std::vector<Point> still =
        Flatten(Ellipse{.Map = Scaling(2, 1), .Start = Degrees(90), .Sweep = {}}, {.Error = 1e-3});
    CHECK(still == std::vector<Point>{{.X = 0, .Y = 1}, {.X = 0, .Y = 1}});
    CHECK(Flatten(Ellipse{.Map = Scaling(0), .Start = {}, .Sweep = Turn}, {.Error = 1e-3}).size()
          == 5);
    const Ellipse flat = {.Map = Scaling(2, 0), .Start = Degrees(10), .Sweep = Turn};
    CHECK(Holds(flat, 1e-3));
    CHECK(Flatten(flat, {.Error = 1e-3}).size() <= 16);
  }

  TEST_CASE("flatten: a polyline keeps its vertices exactly, its bulges as even chords")
  {
    const Polyline polyline = {
        .Vertices = {{.Position = {.X = 0, .Y = 0}, .Bulge = 1},
                     {.Position = {.X = 10, .Y = 0}, .Bulge = 0},
                     {.Position = {.X = 10, .Y = 10}, .Bulge = -0.4}},
        .Closed   = true,
    };
    const std::vector<Point> points = Flatten(polyline, {.Error = 1e-3});
    CHECK(std::ranges::find(points, Point{.X = 10, .Y = 0}) != points.end());
    CHECK(std::ranges::find(points, Point{.X = 10, .Y = 10}) != points.end());
    CHECK(points.back() == Point{.X = 0, .Y = 0});
    CHECK(Holds(polyline, 1e-3));
    CHECK(Flatten(Polyline{}, {.Error = 1}).empty());
    CHECK(
        Flatten(Polyline{.Vertices = {{.Position = {.X = 3, .Y = 4}, .Bulge = 1}}, .Closed = false},
                {.Error = 1})
        == std::vector<Point>{{.X = 3, .Y = 4}});
  }

  /**
   *  A polynomial spline of degree two or three was halved piece by piece until each half's control
   *  points were within the tolerance of its chord. Each piece now takes the greedy chords a Bézier
   *  curve takes, of its exact deviation, the span tried first carried from piece to piece: never
   *  over the tolerance, and fewer chords than halving gives.
   */
  TEST_CASE(
      "flatten: polynomial splines' chords stray no more than the tolerance, fewer than halved")
  {
    const auto halved = [](const Spline &spline, double tolerance) -> std::size_t
    {
      std::size_t                                   chords = 0;
      const std::function<void(Kernel::Piece, int)> halve  = [&](Kernel::Piece piece,
                                                                 int           depth) -> void
      {
        if (depth == 16 || Kernel::Flat(piece, tolerance))
        {
          ++chords;
          return;
        }
        const std::array<Kernel::Piece, 2> halves = Kernel::Split(piece, 0.5);
        halve(halves[0], depth + 1);
        halve(halves[1], depth + 1);
      };
      Kernel::Sections(spline, [&halve](const Kernel::Piece &piece, Range /*span*/) -> void
                       { halve(piece, 0); });
      return chords;
    };
    Draws  draws{0xF1A7};
    double chords = 0.0;
    double before = 0.0;
    for (int round = 0; round < Rounds(100); ++round)
    {
      Spline       spline = Tests::Splined(draws, 2, 1);
      const double scale =
          draws.Next() % 4 == 0 ? std::ldexp(1.0, static_cast<int>(draws.Next() % 40) - 20) : 1.0;
      for (Point &pole : spline.Poles)
      {
        pole = {.X = pole.X * scale, .Y = (pole.Y * scale) + (draws.Next() % 5 == 0 ? 1e4 : 0.0)};
      }
      const double tolerance = std::pow(10.0, -2.0 - (6.0 * std::abs(draws.Near(0, 1)))) * scale;
      CAPTURE(round);
      CHECK(Holds(spline, tolerance));
      const auto many  = static_cast<double>(Flatten(spline, {.Error = tolerance}).size() - 1);
      chords          += many;
      before          += static_cast<double>(halved(spline, tolerance));
    }
    CHECK(chords < 0.8 * before);

    /**
     *  A tolerance past rounding leaves each piece halved, as before.
     */
    const Spline wave = {
        .Degree  = 2,
        .Knots   = {0, 0, 0, 1, 2, 2, 2},
        .Poles   = {{.X = 0, .Y = 0}, {.X = 1, .Y = 1}, {.X = 2, .Y = -1}, {.X = 3, .Y = 0}},
        .Weights = {}};
    CHECK(Flatten(wave, {.Error = 1e-300}).size() == (2 * 65536) + 1);
    CHECK(Flatten(wave, {.Error = 0}).size()      == (2 * 65536) + 1);
    CHECK(Holds(wave, 1e-3));
  }

  TEST_CASE("flatten: splines halved until flat, and Bézier curves")
  {
    CHECK(Holds(Tests::Uneven(), 1e-4));
    CHECK(Holds(Tests::Quarter(), 1e-6));
    CHECK(Holds(Tests::Wave, 1e-5));
    CHECK(Holds(Bezier::Quadratic{.Start   = {.X = 0, .Y = 0},
                                  .Control = {.X = 1, .Y = 2},
                                  .End     = {.X = 2, .Y = 0}},
                1e-5));
    CHECK(Flatten(Tests::Wave, {.Error = 0}).size() == MostChords + 1);
    const Bezier::Cubic loop = {.Start = {.X = 0, .Y = 0},
                                .Out   = {.X = 2, .Y = 2},
                                .In    = {.X = -2, .Y = 2},
                                .End   = {.X = 0, .Y = 0}};
    CHECK(Holds(loop, 1e-4));
  }

  /**
   *  Bézier curves were halved until their control points were within the tolerance of their
   *  chords: two thirds more chords than they need. Each chord now is the longest of a few tries
   *  whose stretch strays no more than the tolerance, by the exact distance of the stretch, a cubic
   *  of its ends and their tangents.
   */
  TEST_CASE("flatten: Bézier curves' chords each stray no more than the tolerance, near the fewest")
  {
    Draws  draws{1114};
    double chords = 0.0;
    double fewest = 0.0;
    for (int round = 0; round < Rounds(48); ++round)
    {
      const Bezier::Cubic     cubic     = {.Start = {.X = draws.Near(0, 1), .Y = draws.Near(0, 1)},
                                           .Out   = {.X = draws.Near(0, 1), .Y = draws.Near(0, 1)},
                                           .In    = {.X = draws.Near(0, 1), .Y = draws.Near(0, 1)},
                                           .End   = {.X = draws.Near(0, 1), .Y = draws.Near(0, 1)}};
      const Bezier::Quadratic quadratic = {
          .Start = cubic.Start, .Control = cubic.Out, .End = cubic.End};
      const double tolerance = 1e-6 * Diagonal(cubic);
      CAPTURE(round);
      CHECK(Holds(cubic, tolerance));
      CHECK(Holds(quadratic, tolerance));
      CHECK(Wanders(quadratic, tolerance) <= 1.0);
      const double many = static_cast<double>(Flatten(cubic, {.Error = tolerance}).size() - 1);
      CHECK(many <= (1.6 * Least(cubic, tolerance)) + 3.0);
      chords += many + static_cast<double>(Flatten(quadratic, {.Error = tolerance}).size() - 1);
      fewest += Least(cubic, tolerance) + Least(quadratic, tolerance);
    }
    CHECK(chords <= 1.2 * fewest);
    for (const double share : {1e-2, 1e-4, 1e-7})
    {
      CAPTURE(share);
      CHECK(Wanders(Tests::Wave, share * Diagonal(Tests::Wave)) <= 1.0);
    }
  }

  TEST_CASE("flatten: Bézier curves through a cusp, back on themselves, or of no length")
  {
    const Bezier::Cubic cusp = {.Start = {.X = 0, .Y = 0},
                                .Out   = {.X = 1, .Y = 1},
                                .In    = {.X = 0, .Y = 1},
                                .End   = {.X = 1, .Y = 0}};
    CHECK(Holds(cusp, 1e-4));
    const std::vector<Point> sharp = Flatten(cusp, {.Error = 1e-4});
    CHECK(std::ranges::any_of(sharp, [&cusp](Point point) -> bool
                              { return Distance(point, At(cusp, 0.5)) < 1e-12; }));
    const Bezier::Cubic back = {.Start = {.X = 0, .Y = 0},
                                .Out   = {.X = 2, .Y = 0},
                                .In    = {.X = -1, .Y = 0},
                                .End   = {.X = 1, .Y = 0}};
    CHECK(Holds(back, 1e-6));
    CHECK(Flatten(back, {.Error = 1e-6}).size() == 4);
    CHECK(Flatten(Bezier::Quadratic{.Start   = {.X = 0, .Y = 0},
                                    .Control = {.X = 1, .Y = 1},
                                    .End     = {.X = 2, .Y = 2}},
                  {.Error = 1e-6})
              .size()
          == 2);
    const Bezier::Quadratic turned = {
        .Start = {.X = 0, .Y = 0}, .Control = {.X = 2, .Y = 0}, .End = {.X = 1, .Y = 0}};
    CHECK(Flatten(turned, {.Error = 1e-6}).size() == 3);
    const Bezier::Cubic bunched = {.Start = {.X = 0, .Y = 0},
                                   .Out   = {.X = 0, .Y = 0},
                                   .In    = {.X = 0, .Y = 0},
                                   .End   = {.X = 1, .Y = 1}};
    CHECK(Flatten(bunched, {.Error = 1e-6}).size() == 2);
    const Bezier::Cubic still = {.Start = {.X = 3, .Y = 4},
                                 .Out   = {.X = 3, .Y = 4},
                                 .In    = {.X = 3, .Y = 4},
                                 .End   = {.X = 3, .Y = 4}};
    CHECK(Flatten(still, {.Error = 1e-6})
          == std::vector<Point>{{.X = 3, .Y = 4}, {.X = 3, .Y = 4}});
    const std::vector<Point> fine = Flatten(Tests::Wave, {.Error = 1e-11});
    CHECK(fine.size() > MostChords / 2);
    CHECK(fine.size() <= MostChords + 1);
  }

  /**
   *  A rational spline of degree two or three is flattened greedy, its stretches' deviations bound
   *  by their control points over their weights: within the tolerance, weights near one or far
   *  apart, knots repeated or not, a span that closes on itself; any other degree, a tolerance
   *  below its rounding or one asking more chords than doubles hold is halved until flat.
   */
  TEST_CASE("flatten: rational splines greedy, their chords within the tolerance")
  {
    Draws draws{0x5EED'F1A7};
    for (int round = 0; round < Rounds(60); ++round)
    {
      const int    degree = 2 + (round % 2);
      const int    poles  = degree + 1 + static_cast<int>(draws.Next() % 6);
      const double spread = round % 3 == 2 ? 30.0 : 1.5;
      Spline       spline = {.Degree = degree, .Knots = {}, .Poles = {}, .Weights = {}};
      for (int index = 0; index < poles; ++index)
      {
        spline.Poles.push_back({.X = index + draws.Near(0, 0.6), .Y = draws.Near(0, 3)});
        spline.Weights.push_back(std::exp(draws.Near(0, std::log(spread))));
      }
      spline.Knots.assign(static_cast<std::size_t>(degree) + 1, 0.0);
      double knot = 0.0;
      for (int index = 0; index < poles - degree - 1; ++index)
      {
        knot += round % 4 == 3 && index == 0 ? 0.0 : draws.Near(1, 0.8);
        spline.Knots.push_back(knot);
      }
      spline.Knots.insert(spline.Knots.end(), static_cast<std::size_t>(degree) + 1, knot + 1);
      if (!Valid(spline))
      {
        continue;
      }
      CAPTURE(round);
      for (const double tolerance : {1e-2, 1e-4, 1e-7})
      {
        CHECK(Holds(spline, tolerance));
      }
    }
    const Spline loop = {.Degree  = 2,
                         .Knots   = {0, 0, 0, 1, 1, 1},
                         .Poles   = {{.X = 0, .Y = 0}, {.X = 1, .Y = 2}, {.X = 0, .Y = 0}},
                         .Weights = {1, 2, 1}};
    CHECK(Holds(loop, 1e-4));
    Spline plain = loop;
    plain.Degree = 1;
    plain.Knots  = {0, 0, 1, 2, 2};
    CHECK(Holds(plain, 1e-4));
    const Spline quartic = {.Degree  = 4,
                            .Knots   = {0, 0, 0, 0, 0, 1, 1, 1, 1, 1},
                            .Poles   = {{.X = 0, .Y = 0},
                                        {.X = 1, .Y =  2},
                                        {.X = 2, .Y = -1},
                                        {.X = 3, .Y =  1},
                                        {.X = 4, .Y = 0}},
                            .Weights = {1, 2, 1, 0.5, 1}};
    CHECK(Holds(quartic, 1e-4));
    Spline wave = quartic;
    wave.Degree = 3;
    wave.Knots  = {0, 0, 0, 0, 1, 2, 2, 2, 2};
    CHECK(Holds(wave, 1e-4));
    CHECK(Flatten(wave, {.Error = 1e-12}).size() > 2);
    CHECK(Flatten(wave, {.Error = 1e-20}).size() > 2);
  }

  /**
   *  A cubic Bézier curve of 1e-160 strayed a thousand times past a tolerance of 1e-6 of its size,
   *  the squares of its deviations below the normal doubles: a curve past 2^400 or below 2^-400 is
   *  flattened scaled by a power of two first, its chords those of the curve near 1 scaled back, to
   *  the bit.
   */
  TEST_CASE("flatten: Bézier curves far below and above the normal doubles, as near 1")
  {
    const Bezier::Cubic     cubic     = {.Start = {.X = 0, .Y = 0},
                                         .Out   = {.X = 1, .Y = 2},
                                         .In    = {.X = 3, .Y = -1},
                                         .End   = {.X = 4, .Y = 1}};
    const Bezier::Quadratic quadratic = {
        .Start = {.X = 0, .Y = 0}, .Control = {.X = 1, .Y = 2}, .End = {.X = 4, .Y = 1}};
    const auto scaled = [](Point point, int power) -> Point
    {
      return {.X = std::ldexp(point.X, power), .Y = std::ldexp(point.Y, power)};
    };
    for (const int power : {-530, -1000, 520, 1000})
    {
      CAPTURE(power);
      const Bezier::Cubic     far     = {.Start = scaled(cubic.Start, power),
                                         .Out   = scaled(cubic.Out, power),
                                         .In    = scaled(cubic.In, power),
                                         .End   = scaled(cubic.End, power)};
      const Bezier::Quadratic distant = {.Start   = scaled(quadratic.Start, power),
                                         .Control = scaled(quadratic.Control, power),
                                         .End     = scaled(quadratic.End, power)};
      for (const double share : {1e-3, 1e-6})
      {
        std::vector<Point> near = Flatten(cubic, {.Error = share});
        std::vector<Point> also = Flatten(quadratic, {.Error = share});
        for (Point &point : near)
        {
          point = scaled(point, power);
        }
        for (Point &point : also)
        {
          point = scaled(point, power);
        }
        CHECK(Flatten(far, {.Error = std::ldexp(share, power)})     == near);
        CHECK(Flatten(distant, {.Error = std::ldexp(share, power)}) == also);
      }
    }
    CHECK(Holds(cubic, 1e-6));
    CHECK(Holds(quadratic, 1e-6));
    CHECK(Flatten(Bezier::Cubic{}, {.Error = 1e-6}).size() == 2);
    const double        endless = std::numeric_limits<double>::infinity();
    const Bezier::Cubic beyond  = {.Start = {.X = 0, .Y = 0},
                                   .Out   = {.X = endless, .Y = 1},
                                   .In    = {.X = 2, .Y = 0},
                                   .End   = {.X = 3, .Y = 0}};
    CHECK(Flatten(beyond, {.Error = 1e-6}).size() == MostChords + 1);
  }

  /**
   *  An accuracy below nothing over an arc of no radius was a division by nothing, the count of
   *  chords read from what came of it.
   */
  /**
   *  An arc or a bulge far from the origin had its points rounded there by more than the chords
   *  were sized for: an arc of radius 3 centered at (1e8, -1e8) strayed 1.00075 times its tolerance
   *  of 1e-6 of its box's diagonal, and a bulge from (1e9, 1e9) to (1e9 + 1, 1e9) 1.035 times. The
   *  chords are sized to the tolerance less what rounding may move the points.
   */
  TEST_CASE("flatten: an arc or a bulge far from the origin stays within the tolerance, rounded")
  {
    const Arc                far       = {.Center = {.X = 1e8, .Y = -1e8},
                                          .Radius = 3,
                                          .Start  = Radians(0.1),
                                          .Sweep  = Radians(-1.0)};
    const double             tolerance = 2.8824227703074932e-06;
    const std::vector<Point> points    = Flatten(far, {.Error = tolerance});
    CHECK(Sagged(points, {.Origin = far.Center, .Center = {}, .Radius = 3}, tolerance) <= 1.0);
    const double step = 4.0 * std::asin(std::sqrt(tolerance / 6.0));
    CHECK(points.size() <= static_cast<std::size_t>(std::ceil(1.0 / step)) + 2);

    const Point    start = {.X = 1e9, .Y = 1e9};
    const Polyline bulge = {.Vertices = {{.Position = start, .Bulge = 0.5},
                                         {.Position = {.X = 1e9 + 1, .Y = 1e9}, .Bulge = 0}},
                            .Closed   = false};
    const double   fine  = 1.6007827349209112e-06;
    CHECK(Sagged(Flatten(bulge, {.Error = fine}),
                 {.Origin = start, .Center = {.X = 0.5, .Y = 0.375}, .Radius = 0.625}, fine)
          <= 1.0);

    /**
     *  Arcs and bulges at random, centered as far as 2^40 from the origin, to tolerances from their
     *  radius down to 2^-24 of it, or 2^-46 of their coordinates where that is more: no doubles
     *  hold a curve nearer than their rounding.
     */
    Draws draws{0xA2C5};
    for (int round = 0; round < Rounds(1600); ++round)
    {
      const double offset = std::ldexp(draws.Near(0, 1), static_cast<int>(draws.Next() % 41));
      const Point  origin = {.X = offset, .Y = offset * draws.Near(0, 1)};
      const double radius =
          std::ldexp(draws.Near(1.5, 0.5), static_cast<int>(draws.Next() % 9) - 4);
      const double within = std::max(radius * std::ldexp(1.0, -static_cast<int>(draws.Next() % 25)),
                                     0x1p-46 * (std::abs(origin.X) + std::abs(origin.Y) + radius));
      CAPTURE(round);
      if (round % 2 == 0)
      {
        const Arc arc = {.Center = origin,
                         .Radius = draws.Next() % 4 == 0 ? -radius : radius,
                         .Start  = Radians(draws.Near(0, 4)),
                         .Sweep  = Radians(draws.Near(0, 7))};
        CHECK(Sagged(Flatten(arc, {.Error = within}),
                     {.Origin = origin, .Center = {}, .Radius = radius}, within)
              <= 1.0);
        continue;
      }
      const double   sign   = draws.Next() % 2 == 0 ? -1.0 : 1.0;
      const double   height = sign * draws.Near(1.55, 1.45);
      const Vector   chord  = {.X = radius * draws.Near(0, 1), .Y = radius * draws.Near(0, 1)};
      const Point    end    = origin + chord;
      const Polyline bowed  = {
          .Vertices = {{.Position = origin, .Bulge = height}, {.Position = end, .Bulge = 0}},
          .Closed   = false};
      const Vector moved  = (end - origin) * 0.5;
      const Vector center = moved
                          + (Vector{.X = -2.0 * moved.Y, .Y = 2.0 * moved.X}
                             * ((1.0 - (height * height)) / (4.0 * height)));
      const double across =
          Length(end - origin) * (1.0 + (height * height)) / (4.0 * std::abs(height));
      CHECK(Sagged(Flatten(bowed, {.Error = within}),
                   {.Origin = origin, .Center = center, .Radius = across}, within)
            <= 1.0);
    }
  }

  TEST_CASE("flatten: an arc of no radius, to an accuracy that is not one")
  {
    const Arc dot = {.Center = {.X = 1, .Y = 2}, .Radius = 0, .Start = {}, .Sweep = Degrees(90)};
    for (const double error : {-1.0, 0.0, std::numeric_limits<double>::quiet_NaN()})
    {
      CAPTURE(error);
      const std::vector<Point> points = Flatten(dot, {.Error = error});
      REQUIRE(points.size() == MostChords + 1);
      CHECK(points.front()  == Point{.X = 1, .Y = 2});
      CHECK(points.back()   == Point{.X = 1, .Y = 2});
    }
    CHECK(Flatten(dot, {.Error = 1e-3}).size() == 2);
  }

  TEST_CASE("flatten: a ray and an xline never end, so give none")
  {
    CHECK(Flatten(Curve{Ray{.Start = {}, .Way = {.X = 1, .Y = 0}}}, {.Error = 1}).empty());
    CHECK(Flatten(Curve{Xline{.Through = {}, .Way = {.X = 1, .Y = 0}}}, {.Error = 1}).empty());
    CHECK(Flatten(Curve{Line{}}, {.Error = 1}).size() == 2);
  }
} // namespace Geoxy
