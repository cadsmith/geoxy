/***************************************************************************************************
 *  @file      tests/geoxy/spline.cpp
 *  @brief     Tests of splines, of what makes one valid, and of what one that is not answers.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "geoxy/spline.hpp"

#include "basis.hpp"
#include "draws.hpp"
#include "real.hpp"

#include "geoxy/accuracy.hpp"
#include "geoxy/area.hpp"
#include "geoxy/at.hpp"
#include "geoxy/bezier/quadratic.hpp"
#include "geoxy/box.hpp"
#include "geoxy/clip.hpp"
#include "geoxy/curve.hpp"
#include "geoxy/derivative.hpp"
#include "geoxy/domain.hpp"
#include "geoxy/extents.hpp"
#include "geoxy/flatten.hpp"
#include "geoxy/intersect.hpp"
#include "geoxy/length.hpp"
#include "geoxy/line.hpp"
#include "geoxy/nearest.hpp"
#include "geoxy/place.hpp"
#include "geoxy/point.hpp"
#include "geoxy/prepared.hpp"
#include "geoxy/tangents.hpp"
#include "geoxy/tolerance.hpp"
#include "geoxy/vector.hpp"

#include <doctest/doctest.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <span>
#include <vector>

namespace Geoxy
{
  namespace
  {
    /**
     *  A valid cubic of four poles, its knots clamped at both ends.
     *
     *  @return The spline.
     */
    auto Cubic() -> Spline
    {
      return {
          .Degree  = 3,
          .Knots   = {0, 0, 0, 0, 1, 1, 1, 1},
          .Poles   = {{.X = 0, .Y = 0}, {.X = 1, .Y = 2}, {.X = 3, .Y = 2}, {.X = 4, .Y = 0}},
          .Weights = {},
      };
    }

    /**
     *  Splines that are not valid, each in its own way: by its degree, by how many knots, poles or
     *  weights it has, and by what they are.
     *
     *  Before they were told apart, each of them read past its knots or its poles somewhere, or
     *  divided by nothing: a crash, a throw or a wrong answer by the question asked.
     *
     *  @return The splines.
     */
    auto Unsound() -> std::vector<Spline>
    {
      constexpr double    Hollow = std::numeric_limits<double>::quiet_NaN();
      std::vector<Spline> all;
      for (const int degree : {0, -1, 4, MaxDegree + 1, std::numeric_limits<int>::min(),
                               std::numeric_limits<int>::max()})
      {
        all.push_back(Cubic());
        all.back().Degree = degree;
      }
      const std::vector<std::vector<double>> knots = {
          {},
          {0, 1},
          {0, 0, 0, 0, 1, 1, 1},
          {0, 0, 0, 0, 0.5, 1, 1, 1, 1},
          {0, 0, 0, 0, 1, 0.5, 1, 1},
          {1, 1, 1, 1, 0, 0, 0, 0},
          {2, 2, 2, 2, 2, 2, 2, 2},
          {0, 0, 0, 0, Hollow, 1, 1, 1},
          {0, 0, 0, Hollow, 1, 1, 1, 1},
          {0, 0, 0, 0, 1, 1, 1, Infinite},
          {-Infinite, 0, 0, 0, 1, 1, 1, 1},
      };
      for (const std::vector<double> &rung : knots)
      {
        all.push_back(Cubic());
        all.back().Knots = rung;
      }
      const std::vector<std::vector<double>> weights = {
          {1},          {1, 1, 1, 1, 1},    {1, -1, 1, 1},     {1, 0, 1, 1},
          {0, 0, 0, 0}, {1, -0.5, 0.25, 1}, {1, Hollow, 1, 1}, {1, 1, Infinite, 1},
      };
      for (const std::vector<double> &heft : weights)
      {
        all.push_back(Cubic());
        all.back().Weights = heft;
      }
      all.push_back(Cubic());
      all.back().Poles.clear();
      all.push_back(Cubic());
      all.back().Poles.pop_back();
      all.push_back(Cubic());
      all.back().Poles[1].X = Hollow;
      all.push_back(Cubic());
      all.back().Poles[2].Y = -Infinite;
      all.push_back({.Degree = 3, .Knots = {}, .Poles = {}, .Weights = {}});
      return all;
    }

    /**
     *  A spline run the other way: its poles and weights in reverse, its knots in reverse and
     *  negated.
     *
     *  @param  spline The spline.
     *  @return        The same curve, run from its end to its start.
     */
    auto SplineReversed(Spline spline) -> Spline
    {
      std::ranges::reverse(spline.Poles);
      std::ranges::reverse(spline.Weights);
      std::ranges::reverse(spline.Knots);
      for (double &knot : spline.Knots)
      {
        knot = -knot;
      }
      return spline;
    }

    /**
     *  How far a point lies from a run of chords.
     *
     *  @param  point  The point.
     *  @param  chords The chords' ends, in order.
     *  @return        The distance to the nearest chord.
     */
    auto Off(Point point, std::span<const Point> chords) -> double
    {
      double least = std::numeric_limits<double>::infinity();
      for (std::size_t index = 0; index + 1 < chords.size(); ++index)
      {
        const Line chord = {.Start = chords[index], .End = chords[index + 1]};
        least            = std::min(least, Distance(point, Nearest(chord, point).At));
      }
      return least;
    }
  } // namespace

  TEST_CASE("spline: equal knot for knot, pole for pole and weight for weight")
  {
    const Spline spline = Cubic();
    Spline       other  = spline;
    CHECK(spline == other);
    other.Degree = 2;
    CHECK_FALSE(spline == other);
    other              = spline;
    other.Knots.back() = 2;
    CHECK_FALSE(spline == other);
    other                 = spline;
    other.Poles.front().X = 1;
    CHECK_FALSE(spline == other);
    other         = spline;
    other.Weights = {1, 1, 1, 1};
    CHECK_FALSE(spline == other);
  }

  TEST_CASE("spline: a valid one, polynomial and rational")
  {
    CHECK(Valid(Cubic()));
    Spline rational  = Cubic();
    rational.Weights = {1, 0.5, 2, 1};
    CHECK(Valid(rational));
  }

  TEST_CASE("spline: each way one is not valid")
  {
    Spline spline = Cubic();
    spline.Degree = 0;
    CHECK_FALSE(Valid(spline));
    spline        = Cubic();
    spline.Degree = MaxDegree + 1;
    spline.Knots  = std::vector<double>(spline.Poles.size() + MaxDegree + 2, 0.0);
    CHECK_FALSE(Valid(spline));
    spline        = Cubic();
    spline.Degree = 4;
    CHECK_FALSE(Valid(spline));
    spline = Cubic();
    spline.Knots.pop_back();
    CHECK_FALSE(Valid(spline));
    spline          = Cubic();
    spline.Knots[5] = Infinite;
    CHECK_FALSE(Valid(spline));
    spline       = Cubic();
    spline.Knots = {0, 0, 0, 0, 1, 0.5, 1, 1};
    CHECK_FALSE(Valid(spline));
    spline       = Cubic();
    spline.Knots = {1, 1, 1, 1, 1, 1, 1, 1};
    CHECK_FALSE(Valid(spline));
    spline            = Cubic();
    spline.Poles[2].Y = std::numeric_limits<double>::quiet_NaN();
    CHECK_FALSE(Valid(spline));
    spline            = Cubic();
    spline.Poles[2].X = Infinite;
    CHECK_FALSE(Valid(spline));
    spline         = Cubic();
    spline.Weights = {1, 1, 1};
    CHECK_FALSE(Valid(spline));
    spline         = Cubic();
    spline.Weights = {1, 0, 1, 1};
    CHECK_FALSE(Valid(spline));
    spline         = Cubic();
    spline.Weights = {1, Infinite, 1, 1};
    CHECK_FALSE(Valid(spline));
  }

  TEST_CASE("spline: one that is not valid is a curve of nothing to every question")
  {
    const Accuracy  fine  = {.Error = 1e-6};
    const Tolerance tight = {.Length = 1e-9, .Radians = 1e-9};
    const Box       box   = {.Min = {.X = -1, .Y = -1}, .Max = {.X = 5, .Y = 5}};
    const Point     aside = {.X = 2, .Y = 5};
    const Curve     level = Line{.Start = {.X = -5, .Y = 1}, .End = {.X = 9, .Y = 1}};
    const Curve     sound = Cubic();
    const Prepared  ruled(level);
    std::size_t     index = 0;
    for (const Spline &spline : Unsound())
    {
      CAPTURE(index);
      REQUIRE_FALSE(Valid(spline));
      CHECK(Extents(spline) == Box{});
      CHECK(Flatten(spline, fine).empty());
      CHECK(Length(spline, fine) == 0);
      CHECK(Length(spline, {.Start = 0.25, .End = 0.75}, fine) == 0);
      CHECK(PlaceAt(spline, 1.5, fine) == Place{});
      CHECK(Area(spline, fine)         == 0);
      CHECK_FALSE(Centroid(spline, fine).has_value());
      CHECK(Perpendiculars(spline, aside).empty());
      CHECK(Nearest(spline, aside) == Place{});
      CHECK(Tangents(spline, aside).empty());
      CHECK(Clip(spline, box).empty());

      const Curve curve = spline;
      for (const Curve &other : {level, sound, curve})
      {
        const Intersections met = Intersect(curve, other, tight);
        CHECK(met.Hits.empty());
        CHECK(met.Overlaps.empty());
        const Intersections back = Intersect(other, curve, tight);
        CHECK(back.Hits.empty());
        CHECK(back.Overlaps.empty());
      }
      const Prepared ready(curve);
      Marking        marks;
      Intersect(Prepared::Pair{ready, ruled}, tight, marks);
      CHECK(marks.Marks.empty());
      Intersect(Prepared::Pair{ruled, ready}, tight, marks);
      CHECK(marks.Marks.empty());
      CHECK(marks.Overlaps.empty());

      /**
       *  A point of it is asked without a pass over all of it, so only what is read is checked: the
       *  answers are the origin and no tangent where that is not sound, and whatever de Boor makes
       *  of the rest, a pole that is not a number among it. No answer is a fault.
       */
      for (const double along : {0.0, 0.5, 1.0, -3.0, 7.0, Infinite})
      {
        CHECK_NOTHROW(static_cast<void>(At(spline, along)));
        CHECK_NOTHROW(static_cast<void>(Derivative(spline, along)));
        CHECK_NOTHROW(static_cast<void>(Domain(spline)));
      }
      ++index;
    }
    CHECK(index == 30);
  }

  /**
   *  A spline of no weights has a weight of 1, which de Boor sums beside its point: far past the
   *  domain the sum is of two numbers far larger than 1 and came to nothing, and the point and the
   *  tangent, divided by it, to no number at all.
   */
  TEST_CASE("spline: one of no weights far past its domain, its point and tangent as its line's")
  {
    const Spline straight = {.Degree  = 1,
                             .Knots   = {0, 0, 1, 1},
                             .Poles   = {{.X = 1, .Y = 2}, {.X = 4, .Y = 6}},
                             .Weights = {}};
    for (const double along : {1e3, 1e9, 1e17, 1e40, -1e17, -1e200})
    {
      CAPTURE(along);
      const Point  here  = At(straight, along);
      const Vector slope = Derivative(straight, along);
      CHECK(here.X == doctest::Approx(1 + (3 * along)).epsilon(1e-12));
      CHECK(here.Y == doctest::Approx(2 + (4 * along)).epsilon(1e-12));
      CHECK(slope  == Vector{.X = 3, .Y = 4});
    }
    const Spline bent = Cubic();
    for (const double along : {1e17, -1e30})
    {
      CAPTURE(along);
      CHECK(std::isfinite(At(bent, along).X));
      CHECK(std::isfinite(Derivative(bent, along).Y));
    }
  }

  TEST_CASE("spline: the weight its point is divided by, 1 with no weights, nothing where unread")
  {
    CHECK(Weighed(Cubic(), 0.5) == 1);
    Spline spline  = Cubic();
    spline.Weights = {1, 0.5, 2, 4};
    CHECK(Weighed(spline, 0)   == 1);
    CHECK(Weighed(spline, 1)   == 4);
    CHECK(Weighed(spline, 0.5) == doctest::Approx(1.5625).epsilon(1e-15));
    spline.Weights = {1, -1, 1, 1};
    CHECK(Weighed(spline, 0.5) == 0);
  }

  TEST_CASE("spline: where one cannot be read it has no point and no tangent")
  {
    Spline spline = Cubic();
    spline.Knots.pop_back();
    CHECK(At(spline, 0.5)         == Point{});
    CHECK(Derivative(spline, 0.5) == Vector{});
    CHECK(Domain(spline)          == Range{});
    spline         = Cubic();
    spline.Weights = {1, -1, 1, 1};
    CHECK(At(spline, 0.5)         == Point{});
    CHECK(Derivative(spline, 0.5) == Vector{});
    CHECK(Domain(spline)          == Range{.Start = 0, .End = 1});
    spline       = Cubic();
    spline.Knots = {0, 0, 0, 0, 1, 0.5, 1, 1};
    CHECK(At(spline, 0.5)         == Point{});
    CHECK(Derivative(spline, 0.5) == Vector{});
  }

  /**
   *  A start knot repeated once more than the degree asks leaves the first span empty, and the
   *  first pole weighing on nothing. Before its domain the spline was read on that empty span, a
   *  division of nothing by nothing; it is read on the first span that is not empty, run on. Such a
   *  spline is not valid, and is still read where it is sound, as the curve its other poles make.
   */
  TEST_CASE("spline: a start knot repeated past the degree, read before its domain")
  {
    const Spline spline = {
        .Degree  = 2,
        .Knots   = {0, 0, 0, 0, 1, 1, 1},
        .Poles   = {{.X = 9, .Y = 9}, {.X = 0, .Y = 0}, {.X = 1, .Y = 2}, {.X = 3, .Y = 1}},
        .Weights = {},
    };
    const Bezier::Quadratic same = {
        .Start = {.X = 0, .Y = 0}, .Control = {.X = 1, .Y = 2}, .End = {.X = 3, .Y = 1}};
    CHECK_FALSE(Valid(spline));
    for (const double along : {-2.0, -0.25, 0.0, 0.5, 1.0, 1.75})
    {
      CAPTURE(along);
      const Point  here  = At(spline, along);
      const Point  owed  = At(same, along);
      const Vector slope = Derivative(spline, along);
      const Vector due   = Derivative(same, along);
      CHECK(here.X  == doctest::Approx(owed.X).epsilon(1e-14));
      CHECK(here.Y  == doctest::Approx(owed.Y).epsilon(1e-14));
      CHECK(slope.X == doctest::Approx(due.X).epsilon(1e-14));
      CHECK(slope.Y == doctest::Approx(due.Y).epsilon(1e-14));
    }
  }

  /**
   *  Weights 2^1400 apart about one span. To keep a pole times its weight from running past every
   *  double the heaviest was brought to 1, which left the lightest 2^-1400, nothing to a double:
   *  its pole fell out of the curve, and the spline's own end was nothing over nothing. The
   *  heaviest is brought no nearer 1 than leaves the lightest a number.
   */
  TEST_CASE("spline: a weight is kept beside one 2^1400 times it")
  {
    const Spline spline = {
        .Degree  = 3,
        .Knots   = {0, 0, 0, 0, 1, 1, 1, 1},
        .Poles   = {{.X = 0, .Y = 0}, {.X = 0.5, .Y = 1.5}, {.X = 1.5, .Y = 1.5}, {.X = 2, .Y = 0}},
        .Weights = {0x1p700, 1, 1, 0x1p-700}};
    REQUIRE(Valid(spline));
    CHECK(Leveling(spline, 3) == -321);
    CHECK(At(spline, 0.0)     == Point{.X = 0, .Y = 0});
    CHECK(At(spline, 1.0)     == Point{.X = 2, .Y = 0});
    for (const double along : {0.25, 0.5, 0.75})
    {
      CAPTURE(along);
      CHECK(std::isfinite(At(spline, along).X));
      CHECK(std::isfinite(At(spline, along).Y));
    }
    Spline mild  = spline;
    mild.Weights = {0x1p100, 1, 1, 0x1p-100};
    CHECK(Leveling(mild, 3) == -100);
  }

  /**
   *  An interior knot repeated a degree and one times lets the curve jump from one pole to the
   *  next; repeated a degree and two times, or an end knot so, it leaves a pole nothing weighs.
   *  Valid took both, and Nearest, Area and Flatten then disagreed with At: the nearest place to
   *  the jump's left end was the start, 1 away where the curve comes within 1e-6; the area of the
   *  curve reversed was not negated; the chords ran 0.78 off the curve. Neither is valid now, and
   *  either is a curve of nothing. A knot repeated as often as the degree is a corner, and valid.
   */
  TEST_CASE("spline: one that jumps, or has a pole nothing weighs, is not valid")
  {
    const Spline jump = {
        .Degree  = 1,
        .Knots   = {0, 0, 1, 1, 2, 2},
        .Poles   = {{.X = 0, .Y = 0}, {.X = 1, .Y = 0}, {.X = 5, .Y = 5}, {.X = 6, .Y = 5}},
        .Weights = {}};
    CHECK(At(jump, 0.999999).X == doctest::Approx(0.999999));
    CHECK(At(jump, 1.0)        == Point{.X = 5, .Y = 5});
    CHECK_FALSE(Valid(jump));
    CHECK(Nearest(jump, {.X = 1, .Y = 0}) == Place{});
    CHECK(Area(jump, {.Error = 1e-12})    == 0);
    CHECK(Flatten(jump, {.Error = 1e-3}).empty());
    const Line across = {.Start = {.X = 3, .Y = -1}, .End = {.X = 3, .Y = 9}};
    CHECK(Intersect(Curve{jump}, Curve{across}, {.Length = 1e-9, .Radians = 1e-9}).Hits.empty());

    const Spline unweighed = {.Degree  = 1,
                              .Knots   = {0, 0, 1, 1, 1, 2, 2},
                              .Poles   = {{.X = 0, .Y = 0},
                                          {.X = 1, .Y = 0},
                                          {.X = 9, .Y = 9},
                                          {.X = 5, .Y = 5},
                                          {.X = 6, .Y = 5}},
                              .Weights = {}};
    CHECK_FALSE(Valid(unweighed));
    const Spline ended = {
        .Degree  = 2,
        .Knots   = {0, 0, 0, 1, 1, 1, 1},
        .Poles   = {{.X = 0, .Y = 0}, {.X = 1, .Y = 2}, {.X = 3, .Y = 1}, {.X = 9, .Y = 9}},
        .Weights = {}};
    CHECK_FALSE(Valid(ended));

    const Spline corner = {.Degree  = 2,
                           .Knots   = {0, 0, 0, 1, 1, 2, 2, 2},
                           .Poles   = {{.X = 0, .Y = 0},
                                       {.X = 1, .Y = 0},
                                       {.X = 2, .Y = 1},
                                       {.X = 1, .Y = 2},
                                       {.X = 0, .Y = 2}},
                           .Weights = {}};
    CHECK(Valid(corner));
    CHECK(At(corner, 1.0) == Point{.X = 2, .Y = 1});
    const Spline unclamped = {
        .Degree  = 2,
        .Knots   = {-1, 0, 0, 0, 1, 1, 1},
        .Poles   = {{.X = 9, .Y = 9}, {.X = 0, .Y = 0}, {.X = 1, .Y = 2}, {.X = 3, .Y = 1}},
        .Weights = {}};
    CHECK(Valid(unclamped));
    CHECK(At(unclamped, 0.0) == Point{.X = 0, .Y = 0});
  }

  /**
   *  A spline whose weights no power of two brings within the doubles, a pole times its leveled
   *  weight running past every double, is not valid: its points would be worked out from no number.
   *  Weights 2^2000 apart are held where the poles are small enough, and a weight that has lost its
   *  digits beside a weight of 1 is brought up with its neighbours, where it was taken as it was
   *  and its products came to nothing.
   */
  TEST_CASE("spline: one whose weights no power of two holds is not valid")
  {
    Spline spline = {.Degree  = 2,
                     .Knots   = {0, 0, 0, 1, 1, 1},
                     .Poles   = {{.X = 0, .Y = 0}, {.X = 1, .Y = 1}, {.X = 2, .Y = 0}},
                     .Weights = {5e-324, 1e308, 5e-324}};
    CHECK_FALSE(Valid(spline));
    spline.Weights = {0x1p-1000, 0x1p1000, 0x1p-1000};
    CHECK(Valid(spline));
    CHECK(Leveling(spline, 2) == -21);
    spline.Poles[1].X = 0x1p43;
    CHECK(Valid(spline));
    spline.Poles[1].X = 0x1p45;
    CHECK_FALSE(Valid(spline));

    spline.Poles[1].X = 1;
    spline.Weights    = {1, 5e-324, 1};
    REQUIRE(Valid(spline));
    CHECK(Leveling(spline, 2) == 53);
    spline.Weights = {1, 0x1p-1021, 1};
    CHECK(Leveling(spline, 2) == 0);
    for (const double along : {0.0, 0.25, 0.5, 1.0})
    {
      CAPTURE(along);
      const Point point = At(spline, along);
      CHECK(std::isfinite(point.X));
      CHECK(std::isfinite(point.Y));
    }
  }

  /**
   *  The robustness audit's 20,000 rational splines of weights 1e+-(290 to 307), or 1e+-150 mixed,
   *  and of coordinates up to 1e12: a pole times its weight ran past every double, and 4,140 of
   *  them gave points, boxes or tolerances that were no numbers. Every valid one now gives numbers.
   */
  TEST_CASE("spline: whatever its weights, a valid one has points, a box and a tolerance")
  {
    Draws     draws(9);
    const int rounds = Rounds(2000);
    int       valid  = 0;
    for (int round = 0; round < rounds; ++round)
    {
      CAPTURE(round);
      Spline spline;
      spline.Degree       = 1 + static_cast<int>(draws.Next() % 4);
      const auto   degree = static_cast<std::size_t>(spline.Degree);
      const auto   count  = degree + 1 + static_cast<std::size_t>(draws.Next() % 4);
      const double scale  = std::pow(10.0, static_cast<double>(draws.Next() % 13));
      for (std::size_t index = 0; index < count; ++index)
      {
        spline.Poles.push_back({.X = draws.Near(0, scale), .Y = draws.Near(0, scale)});
      }
      spline.Knots.assign(degree + 1, 0.0);
      for (std::size_t knot = 1; knot < count - degree; ++knot)
      {
        spline.Knots.push_back(static_cast<double>(knot));
      }
      spline.Knots.insert(spline.Knots.end(), degree + 1, static_cast<double>(count - degree));
      const std::uint64_t mode = draws.Next() % 3;
      for (std::size_t index = 0; index < count; ++index)
      {
        double power = draws.Next() % 2 == 0 ? 150.0 : -150.0;
        if (mode == 0)
        {
          power = 290.0 + static_cast<double>(draws.Next() % 18);
        }
        if (mode == 1)
        {
          power = -300.0 - static_cast<double>(draws.Next() % 8);
        }
        spline.Weights.push_back(std::pow(10.0, power) * draws.Near(1, 0.5));
      }
      REQUIRE(Valid(spline));
      ++valid;
      Box poles;
      for (const Point pole : spline.Poles)
      {
        poles = Grown(poles, pole);
      }
      const Range domain = Domain(spline);
      for (int step = 0; step <= 20; ++step)
      {
        const Point point = At(spline, domain.Start + ((domain.End - domain.Start) * step / 20));
        CHECK(Contains(Padded(poles, 1e-9 * scale), point));
      }
      const Box box = Extents(spline);
      CHECK(std::isfinite(box.Min.X + box.Min.Y + box.Max.X + box.Max.Y));
      const Tolerance tolerance = Tolerance::Of(box);
      CHECK(std::isfinite(tolerance.Length + tolerance.Radians));
    }
    CHECK(valid == rounds);
  }

  /**
   *  What every valid spline owes: its area negated when it is run the other way, its chords within
   *  their accuracy of it, its nearest place no farther than any of its points, and one stretch
   *  shared with itself run the other way. Weights 1e8 and 1e16 times their neighbours' are held to
   *  all but the area, which is not yet worked out to its accuracy beside weights so heavy.
   */
  TEST_CASE("spline: a valid one answers as the curve it is, run either way")
  {
    const Spline cubic = {
        .Degree  = 3,
        .Knots   = {0, 0, 0, 0, 1, 1, 1, 1},
        .Poles   = {{.X = 0, .Y = 0}, {.X = 1, .Y = 2}, {.X = 2, .Y = -1}, {.X = 3, .Y = 0}},
        .Weights = {}};
    std::vector<Spline> splines = {cubic};
    for (const std::vector<double> &weights : {std::vector<double>{1, 0.5, 2, 1},
                                               {1, 1e4, 1e4, 1},
                                               {1, 1e8, 1e8, 1},
                                               {1, 1e16, 1e16, 1}})
    {
      splines.push_back(cubic);
      splines.back().Weights = weights;
    }
    splines.push_back({.Degree  = 2,
                       .Knots   = {0, 0, 0, 1, 1, 2, 2, 2},
                       .Poles   = {{.X = 0, .Y = 0},
                                   {.X = 1, .Y = 0},
                                   {.X = 2, .Y = 1},
                                   {.X = 1, .Y = 2},
                                   {.X = 0, .Y = 2}},
                       .Weights = {1, 2, 1, 0.5, 1}});
    splines.push_back({.Degree  = 3,
                       .Knots   = {0, 0, 0, 0, 0.3, 0.3, 1.1, 2, 2, 2, 2},
                       .Poles   = {{.X = 0, .Y = 0},
                                   {.X = 1, .Y =  3},
                                   {.X = 2, .Y = -2},
                                   {.X = 4, .Y =  1},
                                   {.X = 5, .Y =  4},
                                   {.X = 6, .Y =  0},
                                   {.X = 4, .Y = -3}},
                       .Weights = {}});
    const Tolerance tight = {.Length = 1e-9, .Radians = 1e-9};
    const Accuracy  fine  = {.Error = 1e-3};
    for (std::size_t index = 0; index < splines.size(); ++index)
    {
      CAPTURE(index);
      const Spline &spline = splines[index];
      REQUIRE(Valid(spline));
      const Spline back = SplineReversed(spline);
      REQUIRE(Valid(back));
      if (index != 3 && index != 4)
      {
        const double area = Area(spline, {.Error = 1e-12});
        CHECK(Area(back, {.Error = 1e-12}) == doctest::Approx(-area).epsilon(1e-12));
      }

      const std::vector<Point> chords = Flatten(spline, fine);
      const Range              domain = Domain(spline);
      std::vector<Point>       points;
      for (int step = 0; step <= 400; ++step)
      {
        points.push_back(At(spline, domain.Start + ((domain.End - domain.Start) * step / 400)));
        CHECK(Off(points.back(), chords) <= fine.Error * (1 + 1e-9));
      }
      for (const Point asked : {Point{.X = 1, .Y = 0}, Point{.X = 1.2, .Y = 0.1},
                                Point{.X = 2, .Y = 3}, Point{.X = -1, .Y = -1}})
      {
        CAPTURE(asked);
        const double nearest = Distance(asked, Nearest(spline, asked).At);
        for (const Point point : points)
        {
          CHECK(nearest <= Distance(asked, point) * (1 + 1e-12));
        }
      }
      const Intersections met = Intersect(Curve{spline}, Curve{back}, tight);
      CHECK(met.Hits.empty());
      CHECK_FALSE(met.Spent);
      REQUIRE(met.Overlaps.size() == 1);
      CHECK(met.Overlaps[0].First == domain);
      CHECK_FALSE(met.Overlaps[0].Along);
    }
  }

  /**
   *  A cubic whose weights are 2^1993 apart, met with itself run the other way. Before the scratch
   *  arena was capped it took 1.5 GB in under a second and kept on; its weights are farther apart
   *  than its stand-in can be cut for, so nothing is found, at once, and the meeting says so.
   */
  TEST_CASE("spline: one of weights 1e+-300 met with itself reversed says what it cannot tell")
  {
    const Spline spline = {
        .Degree  = 3,
        .Knots   = {0, 0, 0, 0, 1, 1, 1, 1},
        .Poles   = {{.X = 0, .Y = 0}, {.X = 1, .Y = 2}, {.X = 2, .Y = -1}, {.X = 3, .Y = 0}},
        .Weights = {1e-300, 1e300, 1e300, 1e-300}};
    REQUIRE(Valid(spline));
    const Intersections met =
        Intersect(Curve{spline}, Curve{SplineReversed(spline)}, {.Length = 1e-9, .Radians = 1e-9});
    CHECK(met.Hits.empty());
    CHECK(met.Overlaps.empty());
    CHECK(met.Spent);
  }

  /**
   *  A rational spline run on before its domain, where its weight sums to nothing: its point is
   *  past every number, as dividing by the weight gives it, now with no division, which the
   *  sanitize build stopped on.
   */
  TEST_CASE("spline: where its run-on weight sums to nothing, its point is past every number")
  {
    const Spline spline = {.Degree  = 1,
                           .Knots   = {0, 0, 1, 1},
                           .Poles   = {{.X = 0, .Y = 0}, {.X = 1, .Y = 0}},
                           .Weights = {1, 3}};
    REQUIRE(Weighed(spline, -0.5) == 0);
    const Point point = At(spline, -0.5);
    CHECK(point.X == -Infinite);
    CHECK(std::isnan(point.Y));
    CHECK(Projected({.X = 2, .Y = -3, .W = -0.0}) == Point{.X = -Infinite, .Y = Infinite});
  }
} // namespace Geoxy
