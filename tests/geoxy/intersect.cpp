/***************************************************************************************************
 *  @file      tests/geoxy/intersect.cpp
 *  @brief     Tests of where two curves meet: closed forms, tangents near and exact, shared
 *             carriers, ends, and the same answer either way round.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "geoxy/intersect.hpp"

#include "../suites/respelled.hpp"
#include "curves.hpp"
#include "draws.hpp"

#include "geoxy/angle.hpp"
#include "geoxy/arc.hpp"
#include "geoxy/at.hpp"
#include "geoxy/bezier/cubic.hpp"
#include "geoxy/bezier/quadratic.hpp"
#include "geoxy/curve.hpp"
#include "geoxy/domain.hpp"
#include "geoxy/ellipse.hpp"
#include "geoxy/line.hpp"
#include "geoxy/matrix.hpp"
#include "geoxy/point.hpp"
#include "geoxy/polyline.hpp"
#include "geoxy/prepared.hpp"
#include "geoxy/ray.hpp"
#include "geoxy/spline.hpp"
#include "geoxy/tolerance.hpp"
#include "geoxy/vector.hpp"
#include "geoxy/xline.hpp"

#include <doctest/doctest.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <ctime>
#include <functional>
#include <limits>
#include <numbers>
#include <optional>
#include <tuple>
#include <utility>
#include <vector>

namespace Geoxy
{
  namespace
  {
    /**
     *  A tight tolerance.
     */
    constexpr Tolerance TightMeet = {.Length = 1e-9, .Radians = 1e-9};

    /**
     *  The points where two curves meet.
     *
     *  @param  meeting What meeting them found.
     *  @return         The hits' points, in order along the first.
     */
    auto Points(const Intersections &meeting) -> std::vector<Point>
    {
      std::vector<Point> points;
      points.reserve(meeting.Hits.size());
      for (const Hit &hit : meeting.Hits)
      {
        points.push_back(hit.At);
      }
      return points;
    }

    /**
     *  Whether every hit lies on both curves and the places agree with them.
     *
     *  @param  first   One curve.
     *  @param  second  The other.
     *  @param  meeting What meeting them found.
     *  @return         Whether they do.
     */
    auto OnBoth(const Curve &first, const Curve &second, const Intersections &meeting) -> bool
    {
      return std::ranges::all_of(meeting.Hits,
                                 [&first, &second](const Hit &hit) -> bool
                                 {
                                   return hit.First.At == At(first, hit.First.Along)
                                       && hit.Second.At == At(second, hit.Second.Along)
                                       && Distance(hit.First.At, hit.Second.At) <= TightMeet.Length;
                                 });
    }

    /**
     *  Whether meeting two curves either way round finds the same points.
     *
     *  @param  one The one curve.
     *  @param  two The other.
     *  @return     Whether it does.
     */
    auto Symmetric(const Curve &one, const Curve &two) -> bool
    {
      const Intersections forward  = Intersect(one, two, TightMeet);
      const Intersections backward = Intersect(two, one, TightMeet);
      if (forward.Hits.size() != backward.Hits.size())
      {
        return false;
      }
      for (const Hit &hit : forward.Hits)
      {
        bool found = false;
        for (const Hit &back : backward.Hits)
        {
          found = found || Distance(hit.At, back.At) <= 1e-9;
        }
        if (!found)
        {
          return false;
        }
      }
      return true;
    }

    /**
     *  Whether a run of points matches another to within a tolerance, in order.
     *
     *  @param  points   The points found.
     *  @param  expected The points expected.
     *  @return          Whether they match.
     */
    auto Match(const std::vector<Point> &points, const std::vector<Point> &expected) -> bool
    {
      if (points.size() != expected.size())
      {
        return false;
      }
      for (std::size_t index = 0; index < points.size(); ++index)
      {
        if (Distance(points[index], expected[index]) > 1e-12)
        {
          return false;
        }
      }
      return true;
    }

    /**
     *  Whether two meetings are one, to the bit: their hits, their stretches and whether either ran
     *  out.
     *
     *  @param  lhs One meeting.
     *  @param  rhs The other.
     *  @return     Whether they are.
     */
    auto Alike(const Intersections &lhs, const Intersections &rhs) -> bool
    {
      const auto hit = [](const Hit &one, const Hit &other) -> bool
      {
        return one.At     == other.At
            && one.First  == other.First
            && one.Second == other.Second
            && one.Kind   == other.Kind;
      };
      const auto overlap = [](const Overlap &one, const Overlap &other) -> bool
      {
        return one.First == other.First && one.Second == other.Second && one.Along == other.Along;
      };
      return std::ranges::equal(lhs.Hits, rhs.Hits, hit)
          && std::ranges::equal(lhs.Overlaps, rhs.Overlaps, overlap)
          && lhs.Spent == rhs.Spent;
    }

    /**
     *  A point either way round.
     */
    using Points2 = std::array<std::optional<Point>, 2>;

    /**
     *  Where a meeting touches, where it finds that and nothing else.
     *
     *  @param  meeting The meeting.
     *  @return         The touch's point, where it has one hit, a touch, and no stretch; else
     *                  nothing.
     */
    auto Touched(const Intersections &meeting) -> std::optional<Point>
    {
      const bool touch = meeting.Hits.size() == 1
                      && meeting.Overlaps.empty()
                      && meeting.Hits[0].Kind == Contact::Touch;
      return touch ? std::optional<Point>{meeting.Hits[0].At} : std::nullopt;
    }
  } // namespace

  TEST_CASE("intersect: lines cross, share a carrier exactly, or meet at their ends")
  {
    const Line          rising  = {.Start = {.X = 0, .Y = 0}, .End = {.X = 2, .Y = 2}};
    const Line          falling = {.Start = {.X = 0, .Y = 2}, .End = {.X = 2, .Y = 0}};
    const Intersections crossed = Intersect(rising, falling, TightMeet);
    REQUIRE(crossed.Hits.size()       == 1);
    CHECK(crossed.Hits[0].At          == Point{.X = 1, .Y = 1});
    CHECK(crossed.Hits[0].First.Along == 0.5);
    CHECK(crossed.Hits[0].Kind        == Contact::Cross);

    const Line          base = {.Start = {.X = 0, .Y = 0}, .End = {.X = 4, .Y = 0}};
    const Intersections shared =
        Intersect(base, Line{.Start = {.X = 2, .Y = 0}, .End = {.X = 6, .Y = 0}}, TightMeet);
    REQUIRE(shared.Overlaps.size() == 1);
    CHECK(shared.Hits.empty());
    CHECK(shared.Overlaps[0].First  == Range{.Start = 0.5, .End = 1});
    CHECK(shared.Overlaps[0].Second == Range{.Start = 0, .End = 0.5});
    CHECK(shared.Overlaps[0].Along);
    const Intersections against =
        Intersect(base, Line{.Start = {.X = 6, .Y = 0}, .End = {.X = 2, .Y = 0}}, TightMeet);
    REQUIRE(against.Overlaps.size() == 1);
    CHECK(!against.Overlaps[0].Along);
    CHECK(against.Overlaps[0].Second == Range{.Start = 0.5, .End = 1});

    const Intersections ends =
        Intersect(base, Line{.Start = {.X = 4, .Y = 0}, .End = {.X = 9, .Y = 0}}, TightMeet);
    REQUIRE(ends.Hits.size() == 1);
    CHECK(ends.Hits[0].Kind  == Contact::End);
    CHECK(ends.Overlaps.empty());

    const Intersections tee =
        Intersect(base, Line{.Start = {.X = 1, .Y = 0}, .End = {.X = 1, .Y = 3}}, TightMeet);
    REQUIRE(tee.Hits.size()       == 1);
    CHECK(tee.Hits[0].Kind        == Contact::End);
    CHECK(tee.Hits[0].First.Along == 0.25);

    CHECK(Intersect(base, Line{.Start = {.X = 0, .Y = 1}, .End = {.X = 4, .Y = 1}}, TightMeet)
              .Hits.empty());
    CHECK(Intersect(base, Line{.Start = {.X = 5, .Y = -1}, .End = {.X = 5, .Y = 1}}, TightMeet)
              .Hits.empty());

    const Intersections endless =
        Intersect(Ray{.Start = {.X = 0, .Y = 0}, .Way = {.X = 1, .Y = 1}},
                  Xline{.Through = {.X = 0, .Y = 4}, .Way = {.X = 1, .Y = -1}}, TightMeet);
    REQUIRE(endless.Hits.size() == 1);
    CHECK(endless.Hits[0].At    == Point{.X = 2, .Y = 2});
    CHECK(Intersect(Xline{.Through = {}, .Way = {.X = 1, .Y = 0}},
                    Xline{.Through = {.X = 5, .Y = 0}, .Way = {.X = -2, .Y = 0}}, TightMeet)
              .Overlaps.size()
          == 1);
  }

  TEST_CASE("intersect: a line and a circle, across, tangent, near tangent and stopping short")
  {
    const Arc           circle = Circle({}, 1);
    const Intersections across =
        Intersect(Line{.Start = {.X = -2, .Y = 0}, .End = {.X = 2, .Y = 0}}, circle, TightMeet);
    CHECK(Match(Points(across), {{.X = -1, .Y = 0}, {.X = 1, .Y = 0}}));
    CHECK(across.Hits[0].Kind == Contact::Cross);

    const Intersections tangent =
        Intersect(Line{.Start = {.X = -2, .Y = 1}, .End = {.X = 2, .Y = 1}}, circle, TightMeet);
    REQUIRE(tangent.Hits.size() == 1);
    CHECK(tangent.Hits[0].Kind  == Contact::Touch);
    CHECK(Distance(tangent.Hits[0].At, {.X = 0, .Y = 1}) < 1e-9);

    const Line near = {.Start = {.X = -2, .Y = 1 + 1e-12}, .End = {.X = 2, .Y = 1 + 1e-12}};
    const Intersections close = Intersect(near, circle, {.Length = 1e-8, .Radians = 1e-9});
    REQUIRE(close.Hits.size() == 1);
    CHECK(close.Hits[0].Kind  == Contact::Touch);
    CHECK(Intersect(near, circle, {}).Hits.empty());

    const Line          stopping = {.Start = {.X = 0, .Y = 0}, .End = {.X = 0, .Y = 1 - 1e-12}};
    const Intersections stops    = Intersect(stopping, circle, {.Length = 1e-8, .Radians = 1e-9});
    REQUIRE(stops.Hits.size() == 1);
    CHECK(stops.Hits[0].Kind  == Contact::End);

    const Vector        far   = {.X = 5e6, .Y = -3e6};
    const Intersections moved = Intersect(
        Line{.Start = Point{.X = -2, .Y = 0.5} + far, .End = Point{.X = 2, .Y = 0.5} + far},
        Circle(Point{} + far, 1), {.Length = 1e-8, .Radians = 1e-9});
    CHECK(moved.Hits.size() == 2);
  }

  /**
   *  An arc swept round more than once passes each crossing on every turn. Each hit was placed on
   *  the arc's first pass, whatever turn found it, and two passes four turns apart were taken for
   *  one hit, every place looked at between them being the same point again: no more than four hits
   *  a crossing, and all of them at one place of the arc.
   */
  TEST_CASE("intersect: an arc swept round many times, a hit on each turn it passes")
  {
    const Line line = {.Start = {.X = -3, .Y = -1}, .End = {.X = 4, .Y = 2}};
    const std::array<std::pair<double, std::size_t>, 5> cases = {
        {{20.0, 7}, {26.0, 9}, {30.0, 10}, {-30.0, 9}, {100.0, 32}}};
    for (const auto &[sweep, count] : cases)
    {
      CAPTURE(sweep);
      const Arc coil = {
          .Center = {.X = 1, .Y = 1}, .Radius = 2, .Start = Degrees(10), .Sweep = Radians(sweep)};
      const Arc           once  = {.Center = {.X = 1, .Y = 1},
                                   .Radius = 2,
                                   .Start  = Degrees(10),
                                   .Sweep  = sweep > 0 ? Turn : -Turn};
      const Intersections whole = Intersect(once, line, TightMeet);
      REQUIRE(whole.Hits.size() == 2);

      /**
       *  What is owed, from the one turn alone: each of its two crossings again a turn on, for as
       *  long as the sweep lasts.
       */
      const double        turns = std::abs(sweep) / (2 * std::numbers::pi);
      std::vector<double> owed;
      for (const Hit &hit : whole.Hits)
      {
        const auto passes = static_cast<std::size_t>(std::floor(turns - hit.First.Along)) + 1;
        for (std::size_t pass = 0; pass < passes; ++pass)
        {
          owed.push_back((hit.First.Along + static_cast<double>(pass)) / turns);
        }
      }
      std::ranges::sort(owed);
      REQUIRE(owed.size() == count);

      const Intersections met = Intersect(coil, line, TightMeet);
      REQUIRE(met.Hits.size() == count);
      for (std::size_t index = 0; index < count; ++index)
      {
        CAPTURE(index);
        const Hit &hit = met.Hits[index];
        CHECK(hit.First.Along == doctest::Approx(owed[index]).epsilon(1e-12));
        CHECK(Distance(hit.At, At(coil, hit.First.Along)) < 1e-9);
        CHECK(Distance(hit.At, At(line, hit.Second.Along)) < 1e-9);
        CHECK(hit.Kind == Contact::Cross);
      }

      const Intersections back = Intersect(line, coil, TightMeet);
      std::vector<double> passes;
      passes.reserve(back.Hits.size());
      for (const Hit &hit : back.Hits)
      {
        passes.push_back(hit.Second.Along);
      }
      std::ranges::sort(passes);
      REQUIRE(passes.size() == count);
      for (std::size_t index = 0; index < count; ++index)
      {
        CHECK(passes[index] == doctest::Approx(owed[index]).epsilon(1e-12));
      }
    }
  }

  /**
   *  Every turn of a curve swept round many times that passes a crossing is a hit, and two such
   *  curves meet once for every pair of their turns: an arc of 1e308 radians never came back from a
   *  meeting, two arcs of 1000 radians took four seconds over fifty thousand hits, and one cut off
   *  at its 65,536th turn did not say so. A meeting looks at 1024 turns of such a curve and 2048
   *  pairs of turns at a point, and is spent where there were more.
   */
  TEST_CASE(
      "intersect: curves swept round more turns than a meeting looks at, cut short and said so")
  {
    const Line line = {.Start = {.X = -3, .Y = -1}, .End = {.X = 4, .Y = 2}};
    const auto coil = [](double sweep) -> Arc
    {
      return {
          .Center = {.X = 1, .Y = 1}, .Radius = 2, .Start = Degrees(10), .Sweep = Radians(sweep)};
    };
    const Intersections under = Intersect(coil(6000), line, TightMeet);
    CHECK(under.Hits.size() == 1910);
    CHECK_FALSE(under.Spent);
    for (const double sweep : {7000.0, -7000.0, 1e9, 1e308})
    {
      CAPTURE(sweep);
      /**
       *  The arc's own end is weighed as any curve's is, and may add the crossing beside it, on its
       *  last turn, to the first 1024 turns' two each.
       */
      const Intersections over = Intersect(coil(sweep), line, TightMeet);
      CHECK(over.Hits.size() >= 2048);
      CHECK(over.Hits.size() <= 2050);
      CHECK(over.Spent);
      CHECK(Intersect(line, coil(sweep), TightMeet).Spent);
    }
    const Line          aside = {.Start = {.X = 9, .Y = 9}, .End = {.X = 12, .Y = 9}};
    const Intersections none  = Intersect(coil(1e308), aside, TightMeet);
    CHECK(none.Hits.empty());
    CHECK_FALSE(none.Spent);

    const auto other = [](double sweep) -> Arc
    {
      return {.Center = {.X = 2.5, .Y = 1.5},
              .Radius = 2,
              .Start  = Degrees(40),
              .Sweep  = Radians(sweep)};
    };
    const Intersections few = Intersect(coil(30), other(30), TightMeet);
    CHECK(few.Hits.size() == 45);
    CHECK_FALSE(few.Spent);
    for (const double sweep : {300.0, 1e308})
    {
      CAPTURE(sweep);
      const Intersections many = Intersect(coil(sweep), other(sweep), TightMeet);
      CHECK(many.Hits.size() == 4096);
      CHECK(many.Spent);
    }

    const auto oval = [](double sweep) -> Ellipse
    {
      return {.Map   = {.X = {.X = 3, .Y = 0}, .Y = {.X = 0, .Y = 1}, .Origin = {}},
              .Start = Degrees(10),
              .Sweep = Radians(sweep)};
    };
    const Line across = {.Start = {.X = -4, .Y = -1}, .End = {.X = 4, .Y = 1.5}};
    CHECK(Intersect(oval(30), across, TightMeet).Hits.size() == 10);
    CHECK_FALSE(Intersect(oval(30), across, TightMeet).Spent);
    const Intersections wound = Intersect(oval(1e308), across, TightMeet);
    CHECK(wound.Hits.size() == 2048);
    CHECK(wound.Spent);

    const Intersections same = Intersect(coil(7000), coil(-7000), TightMeet);
    CHECK(same.Spent);
    CHECK_FALSE(same.Overlaps.empty());
  }

  namespace
  {
    /**
     *  Curves of every kind, their lengths times a power of two: the same curves at another size,
     *  to the bit of every number.
     *
     *  @param  scale The power of two.
     *  @return       A line, a ray, an endless line, an arc, an ellipse, a polyline, a spline, and
     *                a quadratic and a cubic Bézier curve.
     */
    auto Kinds(double scale) -> std::vector<Curve>
    {
      const auto placed = [scale](double x, double y) -> Point
      {
        return {.X = x * scale, .Y = y * scale};
      };
      const auto headed = [scale](double x, double y) -> Vector
      {
        return {.X = x * scale, .Y = y * scale};
      };
      return {
          Line{.Start = placed(-3, -1), .End = placed(4, 2)},
          Ray{.Start = placed(-2, 3), .Way = headed(1, -1)},
          Xline{.Through = placed(0, 0.5), .Way = headed(1, 0.125)},
          Arc{.Center = placed(0.5, 0.25),
              .Radius = 2 * scale,
              .Start  = Degrees(10),
              .Sweep  = Degrees(250)},
          Ellipse{.Map   = {.X = headed(3, 0.5), .Y = headed(-0.25, 1), .Origin = placed(0.25, 0)},
                  .Start = {},
                  .Sweep = Turn},
          Polyline{.Vertices = {{.Position = placed(-2, -2), .Bulge = 0.5},
                                {.Position = placed(2, -1), .Bulge = 0},
                                {.Position = placed(1, 3), .Bulge = -0.75}},
                   .Closed   = true},
          Spline{.Degree  = 2,
                 .Knots   = {0, 0, 0, 0.5, 1, 1, 1},
                 .Poles   = {placed(-3, 0), placed(-1, 3), placed(1, -3), placed(3, 0.5)},
                 .Weights = {1, 2, 0.5, 1}},
          Bezier::Quadratic{
              .Start = placed(-3, 1), .Control = placed(0.5, -4), .End = placed(3, 1.5)},
          Bezier::Cubic{.Start = placed(-3, -2),
                        .Out   = placed(0, 5),
                        .In    = placed(1, -5),
                        .End   = placed(3, 2)},
      };
    }
  } // namespace

  /**
   *  A power of two changes no bit of a number but its exponent, so curves scaled by one meet where
   *  they did, scaled. Past 2^100 either way they did not: lengths multiplied ran past every double
   *  or came to nothing, and of the eighty-one pairs of kinds seven met right at 2^600, eleven at
   *  2^-600. Such curves are met at a size near 1 and their hits scaled back.
   */
  TEST_CASE("intersect: curves of any size, as the same curves at a size near 1, to the bit")
  {
    constexpr Tolerance      Near  = {.Length = 0x1p-30, .Radians = 1e-9};
    const std::vector<Curve> small = Kinds(1);
    for (const int power : {-1060, -300, -150, -100, 100, 150, 300, 1020})
    {
      CAPTURE(power);
      const double             scale = std::ldexp(1.0, power);
      const std::vector<Curve> sized = Kinds(scale);
      for (std::size_t one = 0; one < small.size(); ++one)
      {
        for (std::size_t other = 0; other < small.size(); ++other)
        {
          CAPTURE(one);
          CAPTURE(other);
          const Intersections owed = Intersect(small[one], small[other], Near);
          const Intersections met  = Intersect(
              sized[one], sized[other], {.Length = Near.Length * scale, .Radians = Near.Radians});
          REQUIRE(met.Hits.size()     == owed.Hits.size());
          REQUIRE(met.Overlaps.size() == owed.Overlaps.size());
          for (std::size_t index = 0; index < owed.Overlaps.size(); ++index)
          {
            CHECK(met.Overlaps[index].First  == owed.Overlaps[index].First);
            CHECK(met.Overlaps[index].Second == owed.Overlaps[index].Second);
            CHECK(met.Overlaps[index].Along  == owed.Overlaps[index].Along);
          }
          CHECK(met.Spent == owed.Spent);
          for (std::size_t index = 0; index < owed.Hits.size(); ++index)
          {
            const Hit &hit  = met.Hits[index];
            const Hit &same = owed.Hits[index];
            CHECK(hit.At           == Point{.X = same.At.X * scale, .Y = same.At.Y * scale});
            CHECK(hit.First.Along  == same.First.Along);
            CHECK(hit.Second.Along == same.Second.Along);
            CHECK(hit.Kind         == same.Kind);
          }
        }
      }
    }
  }

  /**
   *  Curves too large or too small to meet as they are, met at a size near 1 where scaling loses a
   *  bit: an end's coordinate 2^-1880 of the others' size, lost scaling down, two lines sharing it
   *  exactly; and a ray and a line, and an arc and a line, crossing below the normal doubles, whose
   *  point loses a bit scaled back. Each crossing is worked out again, the double nearest the true
   *  one, met plain or made ready, either first: the lines' from the lines as given, where a
   *  coordinate was lost, and each line's point at its place taken on it as given; the arc's from
   *  the curves as met, rounded at the size given.
   */
  TEST_CASE("intersect: curves met scaled cross where they do as given, where scaling loses a bit")
  {
    const Point shared = {.X = 0x1.0c791f2b3e9a3p-923, .Y = -0x1.3977e06e61f85p+957};
    const Curve into =
        Line{.Start = {.X = 0x1.d93c25355d81dp+963, .Y = 0x1.a829ee9f60c6fp+962}, .End = shared};
    const Curve out =
        Line{.Start = shared, .End = {.X = -0x1.6a0011d78c18ap+957, .Y = 0x1.e26f1ca107fd9p+957}};
    const Curve ray  = Ray{.Start = {.X = -0x1.83a251cf4c4bfp-1017, .Y = -0x1.65046447a6531p-1015},
                           .Way   = {.X = 0x1.6bc7e17001547p-1014, .Y = 0x1.19c3f85799b86p-1012}};
    const Point lead = {.X = 0x0.000000004b8a5p-1022, .Y = -0x1.e5f50dd468bd2p-1017};
    const Curve line =
        Line{.Start = lead, .End = {.X = -0x1.8bdcc2620fb15p-1016, .Y = 0x1.7c65bb8c2fa1ap-1011}};
    const Point below = {.X = -0x0.f05c6640e0167p-1022, .Y = -0x1.f4c07f4ff36b8p-1018};
    const Curve arc   = Arc{.Center = {.X = -0x1.bf404e2fed381p-1020, .Y = 0x1.250ff27155be2p-1018},
                            .Radius = 0x1.e729b0c78c0a0p-1020,
                            .Start  = Radians(-0x1.772417da0c299p+1),
                            .Sweep  = Radians(-0x1.78782d04d7389p+2)};
    const Curve chord = Line{.Start = {.X = -0x1.cdd09bc772b6ep-1019, .Y = 0x1.0b93fdcae0bb2p-1018},
                             .End   = {.X = 0x1.461bcc9fc2dc2p-1018, .Y = 0x1.d2fadf3ad72a2p-1019}};
    const Point across = {.X = 0x0.3491d33be25a7p-1022, .Y = 0x1.fa7df5678b36dp-1019};

    struct Pair
    {
      const Curve *First  = nullptr; /**< The first curve.             */
      const Curve *Second = nullptr; /**< The second.                  */
      Point        At     = {};      /**< Where they cross.            */
      bool         Ends   = false;   /**< Whether at both lines' ends. */
    };

    for (const Pair &pair : {Pair{.First = &into, .Second = &out, .At = shared, .Ends = true},
                             Pair{.First = &out , .Second = &into , .At = shared, .Ends = true },
                             Pair{.First = &ray , .Second = &line , .At = below , .Ends = false},
                             Pair{.First = &line, .Second = &ray  , .At = below , .Ends = false},
                             Pair{.First = &arc , .Second = &chord, .At = across, .Ends = false},
                             Pair{.First = &chord, .Second = &arc, .At = across, .Ends = false}})
    {
      const Tolerance tolerance = {.Length = 0, .Radians = 1e-12};
      const Prepared  one{*pair.First};
      const Prepared  other{*pair.Second};
      Intersections   ready;
      Intersect(Prepared::Pair{one, other}, tolerance, ready);
      const Intersections plain = Intersect(*pair.First, *pair.Second, tolerance);
      for (const Intersections *meeting : std::array<const Intersections *, 2>{&plain, &ready})
      {
        REQUIRE(!meeting->Hits.empty());
        const Hit &hit = meeting->Hits.back();
        CHECK(hit.At == pair.At);
        if (pair.Ends)
        {
          CHECK(hit.First.At  == shared);
          CHECK(hit.Second.At == shared);
        }
      }
    }
  }

  /**
   *  An arc too large to meet as it is, its centre's coordinate 2^-1880 of its radius, lost scaled
   *  down: no crossing of a circle is worked out from curves as given at that size, and the
   *  meeting's is the scaled one's, here the double nearest the true one all the same.
   */
  TEST_CASE("intersect: an arc that loses a number scaled down meets as scaled")
  {
    const Curve arc = Arc{.Center = {.X = 0x1p-923, .Y = 0}, .Radius = 0x1p957, .Sweep = Turn};
    const Curve chord =
        Line{.Start = {.X = -0x1p958, .Y = 0x1p950}, .End = {.X = 0x1p958, .Y = 0x1p950}};
    const Intersections met  = Intersect(arc, chord, {.Length = 0, .Radians = 1e-12});
    const Intersections back = Intersect(chord, arc, {.Length = 0, .Radians = 1e-12});
    REQUIRE(met.Hits.size()  == 2);
    REQUIRE(back.Hits.size() == 2);
    CHECK(met.Hits[0].At     == Point{.X = 0x1.fffbfffbfff8p+956, .Y = 0x1p950});
    CHECK(met.Hits[1].At     == Point{.X = -0x1.fffbfffbfff8p+956, .Y = 0x1p950});
    CHECK(back.Hits[0].At    == met.Hits[1].At);
    CHECK(back.Hits[1].At    == met.Hits[0].At);
  }

  /**
   *  A ray's or an endless line's way is of any length its maker chose; one so long or so short
   *  that its length squared ran past every double or came to nothing met wrongly, a division by
   *  nothing among the ways. It meets as the same curve with a way near 1, each place along it as
   *  many times less as the way is longer.
   */
  TEST_CASE("intersect: a ray or an endless line of any way, as with a way near 1, to the bit")
  {
    constexpr Tolerance      Near   = {.Length = 0x1p-30, .Radians = 1e-9};
    const std::vector<Curve> others = Kinds(1);
    for (const int power : {-1060, -600, -200, 200, 600, 1020})
    {
      CAPTURE(power);
      const double                                 scale = std::ldexp(1.0, power);
      const std::array<std::pair<Curve, Curve>, 3> runs  = {{
          {Ray{.Start = {.X = -2, .Y = 3}, .Way = {.X = 1, .Y = -1}},
           Ray{.Start = {.X = -2, .Y = 3}, .Way = {.X = scale, .Y = -scale}}},
          {Xline{.Through = {.X = 0, .Y = 0.5}, .Way = {.X = 1, .Y = 0.125}},
           Xline{.Through = {.X = 0, .Y = 0.5}, .Way = {.X = scale, .Y = 0.125 * scale}}},
          {Ray{.Start = {.X = -4, .Y = 5}, .Way = {.X = 2, .Y = -2}},
           Ray{.Start = {.X = -4, .Y = 5}, .Way = {.X = 2 * scale, .Y = -2 * scale}}},
      }};
      for (const auto &[plain, odd] : runs)
      {
        for (std::size_t other = 0; other < others.size(); ++other)
        {
          CAPTURE(other);
          const Intersections owed = Intersect(plain, others[other], Near);
          const Intersections met  = Intersect(odd, others[other], Near);
          const Intersections back = Intersect(others[other], odd, Near);
          REQUIRE(met.Hits.size()     == owed.Hits.size());
          REQUIRE(met.Overlaps.size() == owed.Overlaps.size());
          CHECK(back.Hits.size()      == owed.Hits.size());
          CHECK(back.Overlaps.size()  == owed.Overlaps.size());
          for (std::size_t index = 0; index < owed.Hits.size(); ++index)
          {
            CHECK(met.Hits[index].At           == owed.Hits[index].At);
            CHECK(met.Hits[index].First.Along  == owed.Hits[index].First.Along / scale);
            CHECK(met.Hits[index].Second.Along == owed.Hits[index].Second.Along);
          }
          for (std::size_t index = 0; index < owed.Overlaps.size(); ++index)
          {
            CHECK(met.Overlaps[index].First.Start == owed.Overlaps[index].First.Start / scale);
            CHECK(met.Overlaps[index].First.End   == owed.Overlaps[index].First.End / scale);
            CHECK(met.Overlaps[index].Second      == owed.Overlaps[index].Second);
          }
        }
      }
    }
  }

  /**
   *  Beside a weight far heavier than its neighbours a spline rushes along its poles' legs within a
   *  sliver of its parameter, and one step of a double of the parameter moves it by far more than a
   *  step of its coordinates: 2^k of them for a weight 2^k times the rest. A crossing there was
   *  found, then dropped, no place of the spline lying within the tolerance of it: from a weight of
   *  64 with no tolerance asked, of 16,384 with 1e-12. A place is now held to what a step of its
   *  parameter moves its curve, and the crossing kept, its point to the bit. The points owed are
   *  worked out at 120 digits.
   */
  TEST_CASE(
      "intersect: a spline with one weight far heavier than the rest, its crossings to the bit")
  {
    struct Case
    {
      int                  Power  = 0;
      Line                 Across = {};
      std::array<Point, 2> Owed   = {};
    };

    const std::array<Case, 3> conics = {{
        {.Power  = 8,
         .Across = {.Start = {.X = -0.5025140721928083, .Y = 0.5536027710735691},
                    .End   = {.X = 1.7468397103216273, .Y = 0.742519375605899}},
         .Owed   = {{{.X = 0.6504457033340373, .Y = 0.6504364695553646},
                     {.X = 1.2953842196749232, .Y = 0.7046029575424346}}}},
        {.Power  = 20,
         .Across = {.Start = {.X = -0.5959892544217351, .Y = 1.024913898990881},
                    .End   = {.X = 1.9372061048571447, .Y = 0.47509347213260256}},
         .Owed   = {{{.X = 0.7358444535114849, .Y = 0.7358444535105527},
                     {.X = 1.410611078373992, .Y = 0.5893889216256232}}}},
        {.Power  = 30,
         .Across = {.Start = {.X = -0.6445459428702307, .Y = 0.3626647297821429},
                    .End   = {.X = 2.7150354994188453, .Y = 0.16686812065544876}},
         .Owed   = {{{.X = 0.3071970884990951, .Y = 0.3071970884990951},
                     {.X = 1.7785536613263222, .Y = 0.22144633867367788}}}},
    }};
    const auto                conic  = [](int power) -> Spline
    {
      return {.Degree  = 2,
              .Knots   = {0, 0, 0, 1, 1, 1},
              .Poles   = {{.X = 0, .Y = 0}, {.X = 1, .Y = 1}, {.X = 2, .Y = 0}},
              .Weights = {1, std::ldexp(1.0, power), 1}};
    };
    const std::array<Case, 2> cubics = {{
        {.Power  = 12,
         .Across = {.Start = {.X = -0.997260768129417, .Y = 0.05797146626039218},
                    .End   = {.X = 1.9852943707760615, .Y = 1.185332185492061}},
         .Owed   = {{{.X = 0.16588811777510776, .Y = 0.4976241477756277},
                     {.X = 1.6473335250453778, .Y = 1.057588098837923}}}},
        {.Power  = 28,
         .Across = {.Start = {.X = 0.0007223205902138208, .Y = 1.1114487934219606},
                    .End   = {.X = 2.5722987627976783, .Y = -0.11995835182453701}},
         .Owed   = {{{.X = 0.31958656392698664, .Y = 0.9587596875630408},
                     {.X = 1.9388815327690905, .Y = 0.18335540162931452}}}},
    }};
    const auto                cubic  = [](int power) -> Spline
    {
      const double heavy = std::ldexp(1.0, power);
      return {
          .Degree = 3,
          .Knots  = {0, 0, 0, 0, 1, 1, 1, 1},
          .Poles = {{.X = 0, .Y = 0}, {.X = 0.5, .Y = 1.5}, {.X = 1.5, .Y = 1.5}, {.X = 2, .Y = 0}},
          .Weights = {1, heavy, heavy, 1}};
    };
    const auto holds = [](const Curve &curve, const Case &held) -> bool
    {
      for (const Tolerance &tolerance :
           {Tolerance{}, Tolerance{.Length = 1e-12, .Radians = 1e-9}, TightMeet})
      {
        for (const bool swapped : {false, true})
        {
          const Intersections met     = swapped ? Intersect(Curve{held.Across}, curve, tolerance)
                                                : Intersect(curve, Curve{held.Across}, tolerance);
          const bool          counted = met.Hits.size() == 2 && !met.Spent;
          for (const Point owed : held.Owed)
          {
            const bool found = std::ranges::any_of(met.Hits, [owed](const Hit &hit) -> bool
                                                   { return hit.At == owed; });
            if (!counted || !found)
            {
              return false;
            }
          }
        }
      }
      return true;
    };
    for (const Case &held : conics)
    {
      CAPTURE(held.Power);
      CHECK(holds(Curve{conic(held.Power)}, held));
    }
    for (const Case &held : cubics)
    {
      CAPTURE(held.Power);
      CHECK(holds(Curve{cubic(held.Power)}, held));
    }

    /**
     *  Past a weight of 2^32 or so one step of the parameter is too far for a crossing to be told
     *  by it, and past 2^52 no double lies on the leg at all: the spline is met through one that
     *  stands for it, cut until its weights are even, and each crossing polished from the spline's
     *  own numbers, the place told from the end of the span it lies by. To the bit, at any weight
     *  doubles hold.
     */
    const std::array<Case, 3> far = {{
        {.Power  = 60,
         .Across = {.Start = {.X = -0.9584593461557238, .Y = 0.4705701413978771},
                    .End   = {.X = 2.97593582170314, .Y = 0.6647300940468062}},
         .Owed   = {{{.X = 0.5447527212818957, .Y = 0.5447527212818957},
                     {.X = 1.4124280417214141, .Y = 0.5875719582785859}}}},
        {.Power  = 300,
         .Across = {.Start = {.X = 0.2859792714070186, .Y = 1.0978908998369188},
                    .End   = {.X = 2.8787254168904037, .Y = 0.24458759577094824}},
         .Owed   = {{{.X = 0.8968471001783226, .Y = 0.8968471001783226},
                     {.X = 1.2043585366451308, .Y = 0.7956414633548691}}}},
        {.Power  = 1000,
         .Across = {.Start = {.X = -0.3833934145647919, .Y = 0.8141496191582065},
                    .End   = {.X = 2.7011902645420767, .Y = 0.8449241230388949}},
         .Owed   = {{{.X = 0.8262177578741333, .Y = 0.8262177578741333},
                     {.X = 1.1703488890598197, .Y = 0.8296511109401802}}}},
    }};
    for (const Case &held : far)
    {
      CAPTURE(held.Power);
      CHECK(holds(Curve{conic(held.Power)}, held));
    }
    const std::array<Case, 2> farther = {{
        {.Power  = 500,
         .Across = {.Start = {.X = -0.13373841802257824, .Y = 0.12793698585892777},
                    .End   = {.X = 2.262105850510716, .Y = 0.5606456202482955}},
         .Owed   = {{{.X = 0.05394468505701764, .Y = 0.16183405517105293},
                     {.X = 1.8386134952553483, .Y = 0.48415951423395515}}}},
        {.Power  = 1000,
         .Across = {.Start = {.X = -0.3352117726848791, .Y = 1.122968687336656},
                    .End   = {.X = 2.8438270842020863, .Y = -0.40465088609114824}},
         .Owed   = {{{.X = 0.2763631221832064, .Y = 0.8290893665496193},
                     {.X = 1.9996696571508759, .Y = 0.0009910285473724713}}}},
    }};
    for (const Case &held : farther)
    {
      CAPTURE(held.Power);
      CHECK(holds(Curve{cubic(held.Power)}, held));
    }

    /**
     *  Two such splines, legs across legs: each all but its poles' two legs, they cross where the
     *  legs do, a hair from it that no double tells. And one met with itself shares all of itself.
     */
    Spline turned = conic(50);
    turned.Poles  = {{.X = 0, .Y = 1}, {.X = 1, .Y = 0}, {.X = 2, .Y = 1}};
    for (const Tolerance &tolerance : {Tolerance{}, TightMeet})
    {
      const Intersections legs = Intersect(Curve{conic(60)}, Curve{turned}, tolerance);
      CHECK_FALSE(legs.Spent);
      REQUIRE(legs.Hits.size() == 2);
      CHECK(legs.Hits[0].At    == Point{.X = 0.5, .Y = 0.5});
      CHECK(legs.Hits[1].At    == Point{.X = 1.5, .Y = 0.5});
    }
    const Intersections self = Intersect(Curve{conic(60)}, Curve{conic(60)}, TightMeet);
    REQUIRE(self.Overlaps.size()  == 1);
    CHECK(self.Overlaps[0].First  == Range{.Start = 0, .End = 1});
    CHECK(self.Overlaps[0].Second == Range{.Start = 0, .End = 1});

    /**
     *  A line that ends on such a spline's own start meets it there, at an end; one far off meets
     *  it in nothing, told from the boxes before anything is cut.
     */
    const Intersections ended =
        Intersect(Curve{conic(60)},
                  Curve{Line{.Start = {.X = 0, .Y = 0}, .End = {.X = -1, .Y = 5}}}, TightMeet);
    REQUIRE(ended.Hits.size() == 1);
    CHECK(ended.Hits[0].Kind  == Contact::End);
    CHECK(ended.Hits[0].At    == Point{.X = 0, .Y = 0});
    const Intersections apart =
        Intersect(Curve{conic(60)},
                  Curve{Line{.Start = {.X = 90, .Y = 90}, .End = {.X = 95, .Y = 91}}}, TightMeet);
    CHECK(apart.Hits.empty());
    CHECK_FALSE(apart.Spent);

    /**
     *  Weights farther apart than doubles hold beside each other, 2^2000 apart, cannot be stood
     *  for: nothing is found, and the meeting says so.
     */
    Spline beyond            = conic(0);
    beyond.Weights           = {0x1p-1000, 0x1p1000, 0x1p-1000};
    const Intersections none = Intersect(Curve{beyond}, Curve{far[0].Across}, TightMeet);
    CHECK(none.Hits.empty());
    CHECK(none.Spent);
    CHECK(Intersect(Curve{far[0].Across}, Curve{beyond}, TightMeet).Spent);

    /**
     *  Weights all as heavy, or all as light, weigh as weights of 1: the same curve, the same
     *  crossings.
     */
    const Intersections plain = Intersect(Curve{conic(3)}, Curve{conics[0].Across}, TightMeet);
    for (const int power : {-1000, -300, 300, 1000})
    {
      CAPTURE(power);
      Spline scaled = conic(3);
      for (double &weight : scaled.Weights)
      {
        weight = std::ldexp(weight, power);
      }
      const Intersections met = Intersect(Curve{scaled}, Curve{conics[0].Across}, TightMeet);
      REQUIRE(met.Hits.size() == plain.Hits.size());
      for (std::size_t index = 0; index < plain.Hits.size(); ++index)
      {
        CHECK(met.Hits[index].At == plain.Hits[index].At);
      }
    }
  }

  TEST_CASE(
      "intersect: a ray with no way, or a way that is no number, as before; curves of no size")
  {
    constexpr double Endless = std::numeric_limits<double>::infinity();
    const Line       line    = {.Start = {.X = -3, .Y = -1}, .End = {.X = 4, .Y = 2}};
    const Ray        still   = {.Start = {.X = 0.5, .Y = 0.5}, .Way = {}};
    const Ray        wild    = {.Start = {.X = 0.5, .Y = 0.5}, .Way = {.X = Endless, .Y = 1}};
    CHECK(Intersect(still, line, TightMeet).Hits.size() == 1);
    CHECK(Intersect(line, still, TightMeet).Hits.size() == 1);
    CHECK_NOTHROW(static_cast<void>(Intersect(wild, line, TightMeet)));
    const Ray           level = {.Start = {}, .Way = {.X = 1, .Y = 0}};
    const Intersections dot   = Intersect(level, Line{}, {});
    CHECK(dot.Hits.size() == 1);
    CHECK(dot.Hits[0].At  == Point{});
  }

  TEST_CASE("intersect: circles and ellipses, closed form and implicit")
  {
    const Arc           unit = Circle({}, 1);
    const Intersections two  = Intersect(unit, Circle({.X = 1, .Y = 0}, 1), TightMeet);
    REQUIRE(two.Hits.size() == 2);
    CHECK(Distance(two.Hits[0].At, {.X = 0.5, .Y = std::sqrt(0.75)}) < 1e-12);
    CHECK(Distance(two.Hits[1].At, {.X = 0.5, .Y = -std::sqrt(0.75)}) < 1e-12);

    const Intersections kiss = Intersect(unit, Circle({.X = 2, .Y = 0}, 1), TightMeet);
    REQUIRE(kiss.Hits.size() == 1);
    CHECK(kiss.Hits[0].Kind  == Contact::Touch);

    const Intersections same =
        Intersect(unit, Arc{.Center = {}, .Radius = 1, .Start = Degrees(90), .Sweep = Degrees(90)},
                  TightMeet);
    REQUIRE(same.Overlaps.size()       == 1);
    CHECK(same.Overlaps[0].First.Start == doctest::Approx(0.25));
    CHECK(same.Overlaps[0].First.End   == doctest::Approx(0.5));

    const Ellipse       wide = {.Map = Scaling(3, 1), .Start = {}, .Sweep = Turn};
    const Ellipse       tall = {.Map = Scaling(1, 3), .Start = {}, .Sweep = Turn};
    const Intersections four = Intersect(wide, tall, TightMeet);
    const double        side = std::sqrt(0.9);
    REQUIRE(four.Hits.size() == 4);
    for (const Hit &hit : four.Hits)
    {
      CHECK(std::abs(std::abs(hit.At.X) - side) < 1e-12);
      CHECK(std::abs(std::abs(hit.At.Y) - side) < 1e-12);
    }
    CHECK(Intersect(wide, wide, TightMeet).Overlaps.size() == 1);
    CHECK(Match(Points(Intersect(Line{.Start = {.X = -5, .Y = 0}, .End = {.X = 5, .Y = 0}}, wide,
                                 TightMeet)),
                {{.X = -3, .Y = 0}, {.X = 3, .Y = 0}}));
    CHECK(Symmetric(wide, Circle({.X = 1, .Y = 0.5}, 1.2)));
  }

  TEST_CASE("intersect: polylines segment by segment, a bulge on a circle sharing it")
  {
    const Polyline      square = {.Vertices = {{.Position = {.X = 0, .Y = 0}, .Bulge = 0},
                                               {.Position = {.X = 2, .Y = 0}, .Bulge = 0},
                                               {.Position = {.X = 2, .Y = 2}, .Bulge = 0},
                                               {.Position = {.X = 0, .Y = 2}, .Bulge = 0}},
                                  .Closed   = true};
    const Intersections diagonal =
        Intersect(square, Line{.Start = {.X = -1, .Y = 0.5}, .End = {.X = 3, .Y = 1.5}}, TightMeet);
    CHECK(Match(Points(diagonal), {{.X = 2, .Y = 1.25}, {.X = 0, .Y = 0.75}}));

    const Polyline      bulged = {.Vertices = {{.Position = {.X = -1, .Y = 0}, .Bulge = 1},
                                               {.Position = {.X = 1, .Y = 0}, .Bulge = 0}},
                                  .Closed   = false};
    const Intersections shared = Intersect(bulged, Circle({}, 1), TightMeet);
    REQUIRE(shared.Overlaps.size() == 1);
    CHECK(shared.Overlaps[0].First == Range{.Start = 0, .End = 1});
    CHECK(Symmetric(square, Circle({.X = 1, .Y = 1}, 1.2)));

    /**
     *  The edges either side of the bulge are tangent to the circle where it starts and ends: the
     *  touches there go into the stretch, which stays the bulge's own, as its ends are no curve's.
     *  It is cut in two where the circle's parameter starts over.
     */
    const Polyline      capsule = {.Vertices = {{.Position = {.X = -2, .Y = -1}, .Bulge = 0},
                                                {.Position = {.X = 0, .Y = -1}, .Bulge = 1},
                                                {.Position = {.X = 0, .Y = 1}, .Bulge = 0},
                                                {.Position = {.X = -2, .Y = 1}, .Bulge = 0}},
                                   .Closed   = true};
    const Intersections rounded = Intersect(capsule, Circle({}, 1), TightMeet);
    REQUIRE(rounded.Overlaps.size() == 2);
    CHECK(rounded.Overlaps[0].First == Range{.Start = 1, .End = 1.5});
    CHECK(rounded.Overlaps[1].First == Range{.Start = 1.5, .End = 2});
    CHECK(rounded.Hits.empty());
  }

  TEST_CASE("intersect: splines and Bézier curves by their pieces")
  {
    const Bezier::Quadratic arch = {
        .Start = {.X = 0, .Y = 0}, .Control = {.X = 1, .Y = 2}, .End = {.X = 2, .Y = 0}};
    const Line          level = {.Start = {.X = -1, .Y = 0.5}, .End = {.X = 3, .Y = 0.5}};
    const Intersections cut   = Intersect(arch, level, TightMeet);
    REQUIRE(cut.Hits.size()       == 2);
    CHECK(cut.Hits[0].First.Along == doctest::Approx((1 - std::sqrt(0.5)) / 2).epsilon(1e-12));

    const Intersections top =
        Intersect(arch, Line{.Start = {.X = -1, .Y = 1}, .End = {.X = 3, .Y = 1}}, TightMeet);
    REQUIRE(top.Hits.size() == 1);
    CHECK(top.Hits[0].Kind  == Contact::Touch);

    const Bezier::Cubic      hump   = {.Start = {.X = 0, .Y = 0},
                                       .Out   = {.X = 1, .Y = 3},
                                       .In    = {.X = 3, .Y = 2},
                                       .End   = {.X = 4, .Y = 0}};
    const std::vector<Curve> curves = {Tests::Wave,
                                       Tests::Uneven(),
                                       Tests::Quarter(),
                                       arch,
                                       hump,
                                       Circle({.X = 2, .Y = 0.5}, 1.5),
                                       Line{.Start = {.X = -1, .Y = 0.3}, .End = {.X = 9, .Y = 1}}};
    for (std::size_t index = 0; index < curves.size(); ++index)
    {
      for (std::size_t other = index + 1; other < curves.size(); ++other)
      {
        const Intersections meeting = Intersect(curves[index], curves[other], TightMeet);
        CHECK(OnBoth(curves[index], curves[other], meeting));
        CHECK(!meeting.Spent);
        CHECK(Symmetric(curves[index], curves[other]));
      }
    }
    CHECK(Intersect(Tests::Uneven(), Tests::Uneven(), TightMeet).Overlaps.size() == 1);
  }

  TEST_CASE("intersect: degenerate parts: points, lone vertices, flat ellipses, pieces on a line")
  {
    const Line base = {.Start = {.X = 0, .Y = 0}, .End = {.X = 4, .Y = 0}};
    const Line dot  = {.Start = {.X = 1, .Y = 0}, .End = {.X = 1, .Y = 0}};
    CHECK(Intersect(dot, base, TightMeet).Hits.size() == 1);
    const Bezier::Quadratic arch = {
        .Start = {.X = 0, .Y = 0}, .Control = {.X = 1, .Y = 2}, .End = {.X = 2, .Y = 0}};
    CHECK(Intersect(Line{.Start = {.X = 1, .Y = 1}, .End = {.X = 1, .Y = 1}}, arch, TightMeet)
              .Hits.size()
          == 1);
    CHECK(Intersect(
              Polyline{.Vertices = {{.Position = {.X = 2, .Y = 0}, .Bulge = 0}}, .Closed = false},
              base, TightMeet)
              .Hits.size()
          == 1);
    CHECK(Intersect(Line{.Start = {.X = 0, .Y = 0}, .End = {.X = 1, .Y = 0}},
                    Line{.Start = {.X = 2, .Y = 0}, .End = {.X = 3, .Y = 0}}, TightMeet)
              .Hits.empty());

    const Ellipse flat = {.Map = Scaling(2, 0), .Start = {}, .Sweep = Turn};
    CHECK(!Intersect(flat, Line{.Start = {.X = 1, .Y = -1}, .End = {.X = 1, .Y = 1}}, TightMeet)
               .Hits.empty());
    CHECK(Intersect(Line{.Start = {.X = 9, .Y = 9}, .End = {.X = 9, .Y = 9}}, flat, TightMeet)
              .Hits.empty());
    CHECK(!Intersect(flat, arch, TightMeet).Hits.empty());
    CHECK(!Intersect(flat, Circle({.X = 2, .Y = 0}, 1), TightMeet).Hits.empty());
    CHECK(Intersect(flat, Ellipse{.Map = Scaling(0, 2), .Start = {}, .Sweep = Turn}, TightMeet)
              .Spent);

    const Bezier::Cubic     loop = {.Start = {.X = 0, .Y = 0},
                                    .Out   = {.X = 3, .Y = 3},
                                    .In    = {.X = -3, .Y = 3},
                                    .End   = {.X = 0, .Y = 0}};
    const Bezier::Quadratic tall = {
        .Start = {.X = -2, .Y = 0}, .Control = {.X = 0, .Y = 4}, .End = {.X = 2, .Y = 0}};
    const Intersections looped = Intersect(loop, tall, TightMeet);
    REQUIRE(looped.Hits.size() == 2);
    CHECK(OnBoth(loop, tall, looped));
    CHECK(looped.Hits[0].At.X == doctest::Approx(-looped.Hits[1].At.X));
    CHECK(looped.Hits[0].At.Y == doctest::Approx(looped.Hits[1].At.Y));

    const Bezier::Quadratic lying = {
        .Start = {.X = 1, .Y = 0}, .Control = {.X = 2, .Y = 0}, .End = {.X = 3, .Y = 0}};
    const Intersections onLine = Intersect(lying, base, TightMeet);
    REQUIRE(onLine.Overlaps.size() == 1);
    CHECK(onLine.Hits.empty());
    CHECK(onLine.Overlaps[0].First        == Range{.Start = 0, .End = 1});
    CHECK(onLine.Overlaps[0].Second.Start == doctest::Approx(0.25));
    CHECK(onLine.Overlaps[0].Second.End   == doctest::Approx(0.75));
  }

  TEST_CASE("intersect: a cubic and its own half share a stretch, found by sampling")
  {
    const Bezier::Cubic wave = Tests::Wave;
    const Bezier::Cubic half = {
        .Start = wave.Start,
        .Out   = Lerp(wave.Start, wave.Out, 0.5),
        .In    = Lerp(Lerp(wave.Start, wave.Out, 0.5), Lerp(wave.Out, wave.In, 0.5), 0.5),
        .End   = At(wave, 0.5)};
    const Intersections shared = Intersect(wave, half, TightMeet);
    REQUIRE(shared.Overlaps.size() == 1);
    CHECK(shared.Hits.empty());
    CHECK(!shared.Spent);
    CHECK(shared.Overlaps[0].First.Start == 0);
    CHECK(shared.Overlaps[0].First.End   == doctest::Approx(0.5).epsilon(1e-12));
    CHECK(shared.Overlaps[0].Second.End  == doctest::Approx(1).epsilon(1e-12));
    CHECK(shared.Overlaps[0].Along);

    /**
     *  Two cubics leaving one start along one tangent, at 1e8: their first 3e-5 of parameter stays
     *  within 1e-8 of each other, finer than a double there, yet their directions part by 8e-5
     *  radians. Read from control points split off at 1e8 the directions are one, and the start a
     *  stretch; read from the cubics whole it is an end.
     */
    const Point         far     = {.X = 1e8, .Y = 1e8};
    const Bezier::Cubic rising  = {.Start = far,
                                   .Out   = far + Vector{.X = 1, .Y = 3},
                                   .In    = far + Vector{.X = 3, .Y = 2},
                                   .End   = far + Vector{.X = 4, .Y = 0}};
    const Bezier::Cubic parting = {.Start = far,
                                   .Out   = far + Vector{.X = 1, .Y = 3},
                                   .In    = far + Vector{.X = 3, .Y = -2},
                                   .End   = far + Vector{.X = 4, .Y = 1}};
    const Intersections apart   = Intersect(rising, parting, {.Length = 0x1p-20, .Radians = 1e-9});
    CHECK(apart.Overlaps.empty());
    REQUIRE(apart.Hits.size() == 2);
    CHECK(apart.Hits[0].Kind  == Contact::End);
    CHECK(apart.Hits[1].Kind  == Contact::Cross);
  }

  TEST_CASE("intersect: every turn of the rules, each side")
  {
    const Line base = {.Start = {.X = 0, .Y = 0}, .End = {.X = 4, .Y = 0}};
    CHECK(Intersect(Polyline{}, base, TightMeet).Hits.empty());
    CHECK(Intersect(base, Polyline{}, TightMeet).Hits.empty());
    CHECK(Intersect(base, Line{.Start = {.X = 1, .Y = 0}, .End = {.X = 1, .Y = 0}}, TightMeet)
              .Hits.size()
          == 1);
    CHECK(Intersect(base, Line{.Start = {.X = 1, .Y = 1}, .End = {.X = 1, .Y = 2}}, TightMeet)
              .Hits.empty());
    CHECK(Intersect(base, Line{.Start = {.X = 1, .Y = -2}, .End = {.X = 1, .Y = -1}}, TightMeet)
              .Hits.empty());

    const Arc           circle = Circle({}, 1);
    const Intersections grazing =
        Intersect(Line{.Start = {.X = 0, .Y = 1}, .End = {.X = 2, .Y = 1}}, circle, TightMeet);
    REQUIRE(grazing.Hits.size() == 1);
    CHECK(grazing.Hits[0].Kind  == Contact::End);

    CHECK(Intersect(Arc{.Center = {}, .Radius = 1, .Start = {}, .Sweep = Degrees(90)},
                    Arc{.Center = {}, .Radius = 1, .Start = Degrees(180), .Sweep = Degrees(90)},
                    TightMeet)
              .Overlaps.empty());
    const Ellipse turned  = {.Map = Rotation(Degrees(30)), .Start = {}, .Sweep = Turn};
    double        covered = 0.0;
    for (const Overlap &overlap : Intersect(circle, turned, TightMeet).Overlaps)
    {
      covered += overlap.First.End - overlap.First.Start;
    }
    CHECK(covered == doctest::Approx(1));

    /**
     *  Two circles rounded from one, centers and radii a hair apart: one circle within the
     *  tolerance, all of it shared; twice the tolerance apart, not.
     */
    const Tolerance     loose  = {.Length = 1e-9, .Radians = 1e-9};
    const Intersections hair   = Intersect(circle, Circle({.X = 3e-10, .Y = 0}, 1 + 4e-10), loose);
    const Intersections parted = Intersect(circle, Circle({.X = 1e-9, .Y = 0}, 1 + 1e-9), loose);
    REQUIRE(hair.Overlaps.size() == 1);
    CHECK(hair.Overlaps[0].First.End - hair.Overlaps[0].First.Start == doctest::Approx(1));
    CHECK(parted.Overlaps.empty());

    /**
     *  Two arcs of one circle, a rounding apart, meeting end to end: a point, not a stretch.
     */
    const Intersections ended =
        Intersect(Arc{.Center = {}, .Radius = 2, .Start = Radians(Pi), .Sweep = Radians(-Pi)},
                  Arc{.Center = {.X = 5.360630966617862e-16, .Y = -4.440892098500626e-16},
                      .Radius = 1.9999999999999996,
                      .Start  = Radians(-1.5707963267948968),
                      .Sweep  = Radians(1.570796326794897)},
                  {.Length = 1e-8, .Radians = 1e-9});
    CHECK(ended.Overlaps.empty());
    REQUIRE(ended.Hits.size() == 1);
    CHECK(Distance(ended.Hits[0].At, Point{.X = 2, .Y = 0}) <= 1e-8);
    CHECK(Intersect(circle, Circle({}, 2), TightMeet).Hits.empty());
    CHECK(Intersect(circle, Circle({.X = 5, .Y = 0}, 1), TightMeet).Hits.empty());
    const Intersections near =
        Intersect(circle, Circle({.X = 2 + 1e-12, .Y = 0}, 1), {.Length = 1e-8, .Radians = 1e-9});
    REQUIRE(near.Hits.size() == 1);
    CHECK(near.Hits[0].Kind  == Contact::Touch);
    const Ellipse mirrored = {
        .Map = Scaling(-1, 1) * Translation({.X = -1, .Y = 0}), .Start = {}, .Sweep = Turn};
    CHECK(Intersect(circle, mirrored, TightMeet).Hits.size() == 2);
    const Ellipse skewed = {
        .Map = Rotation(Degrees(20)) * Scaling(1.5, 0.5), .Start = {}, .Sweep = Turn};
    CHECK(Intersect(circle, skewed, TightMeet).Hits.size() == 4);

    const Ellipse sheared = {.Map = {.X = {.X = 1, .Y = 0}, .Y = {.X = 0.5, .Y = 1}, .Origin = {}},
                             .Start = {},
                             .Sweep = Turn};
    CHECK(Intersect(circle, sheared, TightMeet).Hits.size() == 4);
    const Intersections secant =
        Intersect(Line{.Start = {.X = -2, .Y = 1 - 1e-10}, .End = {.X = 2, .Y = 1 - 1e-10}}, circle,
                  {.Length = 1e-8, .Radians = 1e-9});
    REQUIRE(secant.Hits.size() == 1);
    CHECK(secant.Hits[0].Kind  == Contact::Touch);
    const Spline bent = {.Degree  = 1,
                         .Knots   = {0, 0, 1, 2, 2},
                         .Poles   = {{.X = -1, .Y = -1}, {.X = 0, .Y = 1}, {.X = 1, .Y = -1}},
                         .Weights = {}};
    CHECK(Intersect(bent, base, TightMeet).Hits.size() == 1);
    const Ellipse flat = {.Map = Scaling(2, 0), .Start = {}, .Sweep = Turn};
    CHECK(!Intersect(Circle({.X = 2, .Y = 0}, 1), flat, TightMeet).Hits.empty());
    const Bezier::Quadratic lying = {
        .Start = {.X = 1, .Y = 0}, .Control = {.X = 2, .Y = 0}, .End = {.X = 3, .Y = 0}};
    CHECK(!Intersect(lying, flat, TightMeet).Hits.empty());

    const Polyline      zigzag = {.Vertices = {{.Position = {.X = 2, .Y = 1}, .Bulge = 0},
                                               {.Position = {.X = 2, .Y = -1}, .Bulge = 0},
                                               {.Position = {.X = 2.5, .Y = 1}, .Bulge = 0},
                                               {.Position = {.X = 3, .Y = 0}, .Bulge = 0},
                                               {.Position = {.X = 5, .Y = 0}, .Bulge = 0},
                                               {.Position = {.X = 5.5, .Y = 1}, .Bulge = 0},
                                               {.Position = {.X = 5.5, .Y = -1}, .Bulge = 0},
                                               {.Position = {.X = 5.8, .Y = 1}, .Bulge = 0}},
                                  .Closed   = false};
    const Intersections crossed =
        Intersect(Line{.Start = {.X = 0, .Y = 0}, .End = {.X = 6, .Y = 0}}, zigzag, TightMeet);
    REQUIRE(crossed.Overlaps.size() == 1);
    CHECK(Match(Points(crossed),
                {{.X = 2, .Y = 0}, {.X = 2.25, .Y = 0}, {.X = 5.5, .Y = 0}, {.X = 5.65, .Y = 0}}));
  }

  TEST_CASE("intersect: a hit is a place on both curves, a crossing passes to the other side")
  {
    /**
     *  A loop through (1, 0) twice: a line along y = 0 crosses it there on both of its passes, two
     *  hits at one point, whichever curve comes first; and at (0, 0), its nose, once more.
     */
    const Bezier::Cubic loop = {.Start = {.X = 2.25, .Y = -1.875},
                                .Out   = {.X = -0.75, .Y = 3.875},
                                .In    = {.X = -0.75, .Y = -3.875},
                                .End   = {.X = 2.25, .Y = 1.875}};
    const Line          axis = {.Start = {.X = -1, .Y = 0}, .End = {.X = 3, .Y = 0}};
    for (const Intersections &meeting :
         {Intersect(loop, axis, TightMeet), Intersect(axis, loop, TightMeet)})
    {
      REQUIRE(meeting.Hits.size() == 3);
      CHECK(std::ranges::all_of(meeting.Hits,
                                [](const Hit &hit) -> bool { return hit.Kind == Contact::Cross; }));
      CHECK(std::ranges::count_if(meeting.Hits, [](const Hit &hit) -> bool
                                  { return Distance(hit.At, {.X = 1, .Y = 0}) < 1e-9; })
            == 2);
    }

    /**
     *  A line that meets a square only at a corner, from outside, touches it: the square comes up
     *  to the line and turns back, and passes to neither side.
     */
    const Polyline square = {.Vertices = {{.Position = {.X = 0, .Y = 0}, .Bulge = 0},
                                          {.Position = {.X = 4, .Y = 0}, .Bulge = 0},
                                          {.Position = {.X = 4, .Y = 4}, .Bulge = 0},
                                          {.Position = {.X = 0, .Y = 4}, .Bulge = 0}},
                             .Closed   = true};
    const Line     corner = {.Start = {.X = 2, .Y = 6}, .End = {.X = 6, .Y = 2}};
    for (const Intersections &meeting :
         {Intersect(square, corner, TightMeet), Intersect(corner, square, TightMeet)})
    {
      REQUIRE(meeting.Hits.size() == 1);
      CHECK(meeting.Hits[0].Kind  == Contact::Touch);
    }

    /**
     *  A circle that osculates an ellipse at both their seams, their contact of order four: the
     *  touch is found on either side of the circle's seam, and is one.
     */
    const Ellipse   oval   = {.Map = Scaling(4, 2), .Start = {}, .Sweep = Turn};
    const Arc       inner  = Circle({.X = 3, .Y = 0}, 1);
    const Tolerance coarse = {.Length = 0x1p-24, .Radians = 1e-9};
    for (const Intersections &meeting :
         {Intersect(oval, inner, coarse), Intersect(inner, oval, coarse)})
    {
      REQUIRE(meeting.Hits.size() == 1);
      CHECK(meeting.Hits[0].Kind  == Contact::Touch);
    }

    /**
     *  A line four tolerances inside a circle's tangent crosses it twice at a shallow angle: where
     *  along each crossing its candidates were polished does not make it a touch.
     */
    const Intersections shallow =
        Intersect(Line{.Start = {.X = -8, .Y = 5 - 4e-9}, .End = {.X = 8, .Y = 5 - 4e-9}},
                  Circle({}, 5), TightMeet);
    REQUIRE(shallow.Hits.size() == 2);
    CHECK(shallow.Hits[0].Kind  == Contact::Cross);
    CHECK(shallow.Hits[1].Kind  == Contact::Cross);
  }

  TEST_CASE("intersect: two crossings within the tolerance are one touch; a closed curve first")
  {
    const Polyline      dip  = {.Vertices = {{.Position = {.X = -1, .Y = 1}, .Bulge = 0},
                                             {.Position = {.X = 0, .Y = -1e-10}, .Bulge = 0},
                                             {.Position = {.X = 1, .Y = 1}, .Bulge = 0}},
                                .Closed   = false};
    const Intersections kiss = Intersect(Line{.Start = {.X = -2, .Y = 0}, .End = {.X = 2, .Y = 0}},
                                         dip, {.Length = 1e-8, .Radians = 1e-9});
    REQUIRE(kiss.Hits.size() == 1);
    CHECK(kiss.Hits[0].Kind  == Contact::Touch);

    const Intersections round = Intersect(
        Circle({}, 1), Line{.Start = {.X = -2, .Y = 0}, .End = {.X = 2, .Y = 0}}, TightMeet);
    CHECK(Match(Points(round), {{.X = 1, .Y = 0}, {.X = -1, .Y = 0}}));
    const Intersections upright = Intersect(
        Circle({}, 1), Line{.Start = {.X = 0, .Y = -2}, .End = {.X = 0, .Y = 2}}, TightMeet);
    CHECK(upright.Hits.size() == 2);
  }

  TEST_CASE("intersect: pieces on one carrier with a line, an arc or a spline of other spans")
  {
    const Spline        quarter = Tests::Quarter();
    const Intersections round   = Intersect(quarter, Circle({}, 1), TightMeet);
    REQUIRE(round.Overlaps.size() == 1);
    CHECK(round.Hits.empty());
    CHECK(round.Overlaps[0].Second.End == doctest::Approx(0.25).epsilon(1e-12));
    const Intersections part = Intersect(
        quarter, Arc{.Center = {}, .Radius = 1, .Start = Degrees(45), .Sweep = Degrees(-90)},
        TightMeet);
    REQUIRE(part.Overlaps.size()     == 1);
    CHECK(part.Overlaps[0].First.End == doctest::Approx(0.5).epsilon(1e-12));
    CHECK(!part.Overlaps[0].Along);
    CHECK(Intersect(quarter,
                    Arc{.Center = {}, .Radius = 1, .Start = Degrees(180), .Sweep = Degrees(45)},
                    TightMeet)
              .Overlaps.empty());

    const Bezier::Cubic wave  = Tests::Wave;
    const Point         left  = Lerp(wave.Start, wave.Out, 0.5);
    const Point         mid   = Lerp(wave.Out, wave.In, 0.5);
    const Point         right = Lerp(wave.In, wave.End, 0.5);
    const Point         near  = Lerp(left, mid, 0.5);
    const Point         far   = Lerp(mid, right, 0.5);
    const Spline spans = {.Degree  = 3,
                          .Knots   = {0, 0, 0, 0, 0.5, 0.5, 0.5, 1, 1, 1, 1},
                          .Poles   = {wave.Start, left, near, At(wave, 0.5), far, right, wave.End},
                          .Weights = {}};
    const Intersections same = Intersect(spans, wave, TightMeet);
    REQUIRE(same.Overlaps.size() == 1);
    CHECK(same.Hits.empty());
    CHECK(same.Overlaps[0].First == Range{.Start = 0, .End = 1});

    const Polyline      straight = {.Vertices = {{.Position = {.X = 1, .Y = 0}, .Bulge = 0},
                                                 {.Position = {.X = 2, .Y = 0}, .Bulge = 0},
                                                 {.Position = {.X = 3, .Y = 0}, .Bulge = 0}},
                                    .Closed   = false};
    const Intersections along =
        Intersect(Line{.Start = {.X = 0, .Y = 0}, .End = {.X = 4, .Y = 0}}, straight, TightMeet);
    REQUIRE(along.Overlaps.size()  == 1);
    CHECK(along.Overlaps[0].First  == Range{.Start = 0.25, .End = 0.75});
    CHECK(along.Overlaps[0].Second == Range{.Start = 0, .End = 2});
  }

  TEST_CASE("intersect: shared stretches partly on both, the other way round, and out of reach")
  {
    const Bezier::Cubic wave  = Tests::Wave;
    const auto          after = [&wave](double along) -> Bezier::Cubic
    {
      const Point a = Lerp(wave.Start, wave.Out, along);
      const Point b = Lerp(wave.Out, wave.In, along);
      const Point c = Lerp(wave.In, wave.End, along);
      const Point d = Lerp(a, b, along);
      const Point q = Lerp(b, c, along);
      return {.Start = Lerp(d, q, along), .Out = q, .In = c, .End = wave.End};
    };
    const auto before = [&wave](double along) -> Bezier::Cubic
    {
      const Point a = Lerp(wave.Start, wave.Out, along);
      const Point b = Lerp(wave.Out, wave.In, along);
      const Point c = Lerp(wave.In, wave.End, along);
      const Point d = Lerp(a, b, along);
      const Point q = Lerp(b, c, along);
      return {.Start = wave.Start, .Out = a, .In = d, .End = Lerp(d, q, along)};
    };
    const Intersections middle = Intersect(after(0.25), before(0.75), TightMeet);
    REQUIRE(middle.Overlaps.size()        == 1);
    CHECK(middle.Overlaps[0].First.Start  == doctest::Approx(0).epsilon(1e-12));
    CHECK(middle.Overlaps[0].First.End    == doctest::Approx(2.0 / 3.0).epsilon(1e-12));
    CHECK(middle.Overlaps[0].Second.Start == doctest::Approx(1.0 / 3.0).epsilon(1e-12));
    CHECK(middle.Overlaps[0].Second.End   == doctest::Approx(1).epsilon(1e-12));
    CHECK(middle.Hits.empty());
    const Intersections turned = Intersect(before(0.75), after(0.25), TightMeet);
    REQUIRE(turned.Overlaps.size() == 1);
    CHECK(turned.Hits.empty());
    CHECK(turned.Overlaps[0].First.Start == doctest::Approx(1.0 / 3.0).epsilon(1e-12));
    CHECK(turned.Overlaps[0].First.End   == doctest::Approx(1).epsilon(1e-12));

    const Spline        backward = {.Degree = 2,
                                    .Knots  = {0, 0, 0, 1, 1, 1},
                                    .Poles = {{.X = 0, .Y = 1}, {.X = 1, .Y = 1}, {.X = 1, .Y = 0}},
                                    .Weights = {1, std::sqrt(0.5), 1}};
    const Intersections reversed = Intersect(backward, Circle({}, 1), TightMeet);
    REQUIRE(reversed.Overlaps.size() == 1);
    CHECK(!reversed.Overlaps[0].Along);

    const Bezier::Quadratic lying = {
        .Start = {.X = 1, .Y = 0}, .Control = {.X = 2, .Y = 0}, .End = {.X = 3, .Y = 0}};
    CHECK(Intersect(lying, Line{.Start = {.X = 5, .Y = 0}, .End = {.X = 8, .Y = 0}}, TightMeet)
              .Overlaps.empty());

    const Bezier::Cubic hump  = {.Start = {.X = 0, .Y = 0},
                                 .Out   = {.X = 1, .Y = 3},
                                 .In    = {.X = 3, .Y = 2},
                                 .End   = {.X = 4, .Y = 0}};
    const Intersections loose = Intersect(wave, hump, {.Length = 10, .Radians = 1e-9});
    CHECK(!loose.Hits.empty());

    Bezier::Cubic lifted = wave;
    for (Point *pole : {&lifted.Start, &lifted.Out, &lifted.In, &lifted.End})
    {
      *pole = *pole + Vector{.X = 0, .Y = 1e-10};
    }
    /**
     *  Lifted 1e-10 all along, within the tolerance, it shares the cubic end to end: closeness
     *  decides, not the angle, here none at all.
     */
    const Intersections raised = Intersect(wave, lifted, {.Length = 1e-9, .Radians = 0});
    CHECK(!raised.Spent);
    REQUIRE(raised.Overlaps.size() == 1);
    CHECK(raised.Overlaps[0].First == Range{.Start = 0, .End = 1});
  }

  TEST_CASE("intersect: shared stretches joined where they run on, apart where they do not")
  {
    const Line   axis   = {.Start = {.X = 0, .Y = 0}, .End = {.X = 4, .Y = 0}};
    const Spline uneven = {
        .Degree  = 1,
        .Knots   = {0, 0, 0.1, 0.7, 1, 1},
        .Poles   = {{.X = 0.5, .Y = 0}, {.X = 1, .Y = 0}, {.X = 2, .Y = 0}, {.X = 3, .Y = 0}},
        .Weights = {}};
    const Intersections forward  = Intersect(axis, uneven, TightMeet);
    const Intersections backward = Intersect(uneven, axis, TightMeet);
    REQUIRE(forward.Overlaps.size()       == 1);
    REQUIRE(backward.Overlaps.size()      == 1);
    CHECK(forward.Overlaps[0].First.Start == doctest::Approx(0.125));
    CHECK(forward.Overlaps[0].First.End   == doctest::Approx(0.75));
    CHECK(backward.Overlaps[0].First      == Range{.Start = 0, .End = 1});

    const Polyline twice = {.Vertices = {{.Position = {.X = 0, .Y = 0}, .Bulge = 0},
                                         {.Position = {.X = 1, .Y = 0}, .Bulge = 0},
                                         {.Position = {.X = 1, .Y = 1}, .Bulge = 0},
                                         {.Position = {.X = 2, .Y = 1}, .Bulge = 0},
                                         {.Position = {.X = 2, .Y = 0}, .Bulge = 0},
                                         {.Position = {.X = 3, .Y = 0}, .Bulge = 0}},
                            .Closed   = false};
    CHECK(Intersect(axis, twice, TightMeet).Overlaps.size() == 2);

    const Polyline back = {.Vertices = {{.Position = {.X = 0, .Y = 0}, .Bulge = 0},
                                        {.Position = {.X = 2, .Y = 0}, .Bulge = 0},
                                        {.Position = {.X = 1, .Y = 0}, .Bulge = 0}},
                           .Closed   = false};
    CHECK(Intersect(axis, back, TightMeet).Overlaps.size() == 2);

    const Polyline loop = {.Vertices = {{.Position = {.X = 1, .Y = 0}, .Bulge = 0},
                                        {.Position = {.X = 2, .Y = 0}, .Bulge = 0},
                                        {.Position = {.X = 2, .Y = 1}, .Bulge = 0},
                                        {.Position = {.X = 0, .Y = 1}, .Bulge = 0},
                                        {.Position = {.X = 0, .Y = 0}, .Bulge = 0},
                                        {.Position = {.X = 1, .Y = 0}, .Bulge = 0}},
                           .Closed   = false};
    CHECK(Intersect(axis, loop, TightMeet).Overlaps.size() == 2);

    const Polyline      reverse = {.Vertices = {{.Position = {.X = 3, .Y = 0}, .Bulge = 0},
                                                {.Position = {.X = 2, .Y = 0}, .Bulge = 0},
                                                {.Position = {.X = 1, .Y = 0}, .Bulge = 0}},
                                   .Closed   = false};
    const Intersections against = Intersect(axis, reverse, TightMeet);
    REQUIRE(against.Overlaps.size() == 1);
    CHECK(!against.Overlaps[0].Along);
    CHECK(against.Overlaps[0].Second == Range{.Start = 0, .End = 2});
  }

  /**
   *  The processor time is measured, not the wall clock's, after a first call: a machine loaded
   *  with other work, or a process that has not touched its memory yet, does not slow the check
   *  down.
   */
  TEST_CASE("intersect: curves leaving one tangent are halved no deeper than the tolerance needs")
  {
    const Bezier::Cubic hump   = {.Start = {.X = 0, .Y = 0},
                                  .Out   = {.X = 1, .Y = 3},
                                  .In    = {.X = 3, .Y = 2},
                                  .End   = {.X = 4, .Y = 0}};
    const Bezier::Cubic parted = {.Start = {.X = 0, .Y = 0},
                                  .Out   = {.X = 1, .Y = 3},
                                  .In    = {.X = 3, .Y = -2},
                                  .End   = {.X = 4, .Y = 1}};
    for (const Tolerance tolerance : {TightMeet, Tolerance{.Length = 1e-3, .Radians = 1e-9}})
    {
      static_cast<void>(Intersect(hump, parted, tolerance));
      const std::clock_t  start   = std::clock();
      const Intersections meeting = Intersect(hump, parted, tolerance);
      const std::clock_t  took    = std::clock() - start;
      CHECK(!meeting.Spent);
      CHECK(meeting.Overlaps.empty());
      REQUIRE(meeting.Hits.size() == 2);
      CHECK(meeting.Hits[0].Kind  == Contact::End);
      CHECK(Distance(meeting.Hits[0].At, {}) <= tolerance.Length);
      CHECK(meeting.Hits[1].Kind         == Contact::Cross);
      CHECK(meeting.Hits[1].First.Along  == doctest::Approx(12.0 / 13).epsilon(1e-12));
      CHECK(meeting.Hits[1].Second.Along == doctest::Approx(12.0 / 13).epsilon(1e-12));
      CHECK(took < CLOCKS_PER_SEC / 20);
    }
  }
} // namespace Geoxy

namespace Geoxy
{
  TEST_CASE("intersect: one geometry, one answer, whatever its spelling")
  {
    const Tolerance tight = {.Length = 1e-9, .Radians = 1e-9};

    /**
     *  Two lines 1e-10 apart all along share it, as lines and as straight Bézier curves.
     */
    const Line under = {.Start = {.X = 0, .Y = 0}, .End = {.X = 10, .Y = 0}};
    const Line above = {.Start = {.X = 0, .Y = 1e-10}, .End = {.X = 10, .Y = 1e-10}};
    const auto bent  = [](const Line &line) -> Bezier::Quadratic
    {
      return {.Start = line.Start, .Control = Midpoint(line.Start, line.End), .End = line.End};
    };
    for (const Intersections &meeting :
         {Intersect(under, above, tight), Intersect(bent(under), bent(above), tight),
          Intersect(under, bent(above), tight)})
    {
      CHECK(meeting.Hits.empty());
      REQUIRE(meeting.Overlaps.size()       == 1);
      CHECK(meeting.Overlaps[0].First.Start == doctest::Approx(0));
      CHECK(meeting.Overlaps[0].First.End   == doctest::Approx(1));
    }

    /**
     *  Straight and parallel, 1e-3 apart: nothing. Collinear but overlapping by less than the
     *  tolerance: an end, not a stretch.
     */
    const Line          away  = {.Start = {.X = 0, .Y = 1e-3}, .End = {.X = 10, .Y = 1e-3}};
    const Intersections apart = Intersect(bent(under), bent(away), tight);
    CHECK(apart.Hits.empty());
    CHECK(apart.Overlaps.empty());
    const Intersections kissing =
        Intersect(under, Line{.Start = {.X = 10 - 5e-10, .Y = 0}, .End = {.X = 20, .Y = 0}}, tight);
    CHECK(kissing.Overlaps.empty());
    REQUIRE(kissing.Hits.size() == 1);
    CHECK(kissing.Hits[0].Kind  == Contact::End);

    /**
     *  Endless: an xline and a ray 1e-10 off it, exactly parallel, share; tilted by 1e-12 they do
     *  not, however near they run where they start.
     */
    const Xline         axis = {.Through = {.X = 0, .Y = 0}, .Way = {.X = 1, .Y = 0}};
    const Intersections along =
        Intersect(axis, Ray{.Start = {.X = 5, .Y = 1e-10}, .Way = {.X = 1, .Y = 0}}, tight);
    REQUIRE(along.Overlaps.size() == 1);
    CHECK(!std::isfinite(along.Overlaps[0].First.End));
    CHECK(Intersect(axis, Ray{.Start = {.X = 5, .Y = 1e-10}, .Way = {.X = 1, .Y = 1e-12}}, tight)
              .Overlaps.empty());
    CHECK(Intersect(axis, Ray{.Start = {.X = 5, .Y = 1e-3}, .Way = {.X = 1, .Y = 0}}, tight)
              .Overlaps.empty());

    /**
     *  A cubic on the line that runs past its ends, or out and back: a stretch for each pass along
     *  the line, between where it turns back, one way at a time.
     */
    const Line          line    = {.Start = {.X = -5, .Y = 0}, .End = {.X = 5, .Y = 0}};
    const Bezier::Cubic past    = {.Start = {.X = 0, .Y = 0},
                                   .Out   = {.X = 4, .Y = 0},
                                   .In    = {.X = -1, .Y = 0},
                                   .End   = {.X = 3, .Y = 0}};
    const Bezier::Cubic back    = {.Start = {.X = 0, .Y = 0},
                                   .Out   = {.X = 3, .Y = 0},
                                   .In    = {.X = 3, .Y = 0},
                                   .End   = {.X = 0, .Y = 0}};
    const Intersections beyond  = Intersect(line, past, tight);
    const Intersections turning = Intersect(line, back, tight);
    REQUIRE(beyond.Overlaps.size()      == 3);
    CHECK(beyond.Overlaps[0].Second     == Range{.Start = 0, .End = 1.0 / 3});
    CHECK(beyond.Overlaps[1].Second.End == doctest::Approx(2.0 / 3));
    CHECK(!beyond.Overlaps[1].Along);
    CHECK(beyond.Overlaps[2].First.End  == doctest::Approx(0.8));
    REQUIRE(turning.Overlaps.size()     == 2);
    CHECK(turning.Overlaps[0].First.End == doctest::Approx(0.725));
    CHECK(turning.Overlaps[1].First.End == doctest::Approx(0.725));
    CHECK(turning.Hits.empty());

    /**
     *  Lines whose carriers cross past an end, farther than the tolerance: nothing.
     */
    CHECK(Intersect(under, Line{.Start = {.X = 11, .Y = -1}, .End = {.X = 12, .Y = 1}}, tight)
              .Hits.empty());

    /**
     *  A segment stopping just short of another's line still meets it there, at its end.
     */
    const Intersections stopping = Intersect(
        under, Line{.Start = {.X = 10 + 1e-10, .Y = -1}, .End = {.X = 10 + 1e-10, .Y = 1}}, tight);
    REQUIRE(stopping.Hits.size() == 1);
    CHECK(stopping.Hits[0].Kind  == Contact::End);

    /**
     *  A cubic along its inflection's tangent, a line spelled as a Bézier curve: the band where
     *  they stay within the tolerance is a touch, not a stretch; so are three crossings closer than
     *  the tolerance parts them; and a cusp along its tangent, the line first.
     */
    const Line          flat       = {.Start = {.X = -2, .Y = 0}, .End = {.X = 2, .Y = 0}};
    const Bezier::Cubic inflection = {.Start = {.X = -1, .Y = -1},
                                      .Out   = {.X = -1.0 / 3, .Y = 1},
                                      .In    = {.X = 1.0 / 3, .Y = -1},
                                      .End   = {.X = 1, .Y = 1}};
    const Bezier::Cubic cluster    = {.Start = {.X = -1, .Y = -1 + 0x1p-40},
                                      .Out   = {.X = -1.0 / 3, .Y = 1 - (0x1p-40 / 3)},
                                      .In    = {.X = 1.0 / 3, .Y = -1 + (0x1p-40 / 3)},
                                      .End   = {.X = 1, .Y = 1 - 0x1p-40}};
    const Bezier::Cubic cusp       = {.Start = {.X = 1, .Y = -1},
                                      .Out   = {.X = -1.0 / 3, .Y = 1},
                                      .In    = {.X = -1.0 / 3, .Y = -1},
                                      .End   = {.X = 1, .Y = 1}};
    const Tolerance     loose      = {.Length = 0x1p-20, .Radians = 1e-9};
    for (const Intersections &meeting :
         {Intersect(inflection, bent(flat), loose), Intersect(cluster, bent(flat), loose),
          Intersect(bent(flat), cusp, loose), Intersect(inflection, bent(flat), tight)})
    {
      CHECK(meeting.Overlaps.empty());
      REQUIRE(meeting.Hits.size() == 1);
      CHECK(meeting.Hits[0].Kind  == Contact::Touch);
    }

    /**
     *  Two straight pieces crossing at 1e-8 radians stay within the tolerance over a long band, and
     *  cross: the band is no stretch, and its hit a crossing.
     */
    const Line          tilted = {.Start = {.X = -10, .Y = -1e-7}, .End = {.X = 10, .Y = 1e-7}};
    const Intersections shallow =
        Intersect(bent(Line{.Start = {.X = -10, .Y = 0}, .End = {.X = 10, .Y = 0}}),
                  Bezier::Cubic{.Start = tilted.Start,
                                .Out   = Lerp(tilted.Start, tilted.End, 0.25),
                                .In    = Lerp(tilted.Start, tilted.End, 0.5),
                                .End   = tilted.End},
                  tight);
    CHECK(shallow.Overlaps.empty());
    REQUIRE(shallow.Hits.size() == 1);
    CHECK(shallow.Hits[0].Kind  == Contact::Cross);

    /**
     *  Arcs of one circle as rational splines, their knots apart: one stretch, its ends where the
     *  arcs end, across the knots where it was found piece by piece; and circles touching, as
     *  splines, one touch.
     */
    const auto quarters = [](Point center, double start) -> Spline
    {
      Spline spline = {.Degree = 2, .Knots = {0, 0, 0}, .Poles = {}, .Weights = {}};
      for (int quarter = 0; quarter < 2; ++quarter)
      {
        const double from = start + (quarter * std::numbers::pi / 2);
        const Point  head = center + Direction(Radians(from));
        if (quarter == 0)
        {
          spline.Poles.push_back(head);
          spline.Weights.push_back(1);
        }
        else
        {
          spline.Knots.push_back(spline.Knots.back());
        }
        spline.Poles.push_back(
            center + (Direction(Radians(from + (std::numbers::pi / 4))) * std::numbers::sqrt2));
        spline.Weights.push_back(std::sqrt(0.5));
        spline.Poles.push_back(center + Direction(Radians(from + (std::numbers::pi / 2))));
        spline.Weights.push_back(1);
        spline.Knots.push_back(spline.Knots.back() + 1);
      }
      spline.Knots.push_back(spline.Knots.back());
      spline.Knots.push_back(spline.Knots.back());
      return spline;
    };
    const Intersections arcs =
        Intersect(quarters({}, 0), quarters({}, std::numbers::pi / 4), tight);
    CHECK(arcs.Hits.empty());
    REQUIRE(arcs.Overlaps.size() == 1);
    CHECK(Distance(At(quarters({}, 0), arcs.Overlaps[0].First.Start),
                   Point{.X = std::sqrt(0.5), .Y = std::sqrt(0.5)})
          < 1e-9);
    const Intersections kiss = Intersect(quarters({}, -std::numbers::pi / 2),
                                         quarters({.X = 2, .Y = 0}, std::numbers::pi / 2), tight);
    CHECK(kiss.Overlaps.empty());
    REQUIRE(kiss.Hits.size() == 1);
    CHECK(Distance(kiss.Hits[0].At, {.X = 1, .Y = 0}) < 1e-4);

    /**
     *  Far from the origin, as splines, a band's kind does not turn on where its hit is polished
     *  to: three crossings closer than the tolerance, the curves parallel only between samples, are
     *  one touch; a cusp along its tangent, the line first, all but stalled at it, a touch; a
     *  secant four tolerances inside a circle, a band and no tangent, two crossings.
     */
    const Vector        shift  = {.X = 5'000'000, .Y = -3'000'000};
    const double        sliver = 0x1p-28;
    const Line          level  = {.Start = Point{.X = -2, .Y = 0} + shift,
                                  .End   = Point{.X = 2, .Y = 0} + shift};
    const Bezier::Cubic three  = {.Start = Point{.X = -1.5, .Y = -3.375 + (1.5 * sliver)} + shift,
                                  .Out   = Point{.X = -0.5, .Y = 3.375 + (0.5 * sliver)} + shift,
                                  .In    = Point{.X = 0.5, .Y = -3.375 - (0.5 * sliver)} + shift,
                                  .End   = Point{.X = 1.5, .Y = 3.375 - (1.5 * sliver)} + shift};
    const Intersections clustered = Intersect(*Suites::Respelled(three), *Suites::Respelled(level),
                                              {.Length = 0x1p-24, .Radians = 1e-9});
    REQUIRE(clustered.Hits.size() == 1);
    CHECK(clustered.Hits[0].Kind  == Contact::Touch);

    const Vector        far    = {.X = 1e8, .Y = 1e8};
    const Bezier::Cubic cusped = {.Start = Point{.X = 2.25, .Y = -3.375} + far,
                                  .Out   = Point{.X = -0.75, .Y = 3.375} + far,
                                  .In    = Point{.X = -0.75, .Y = -3.375} + far,
                                  .End   = Point{.X = 2.25, .Y = 3.375} + far};
    const Line ground = {.Start = Point{.X = -1, .Y = 0} + far, .End = Point{.X = 3, .Y = 0} + far};
    const Intersections cusping = Intersect(*Suites::Respelled(ground), *Suites::Respelled(cusped),
                                            {.Length = 0x1p-20, .Radians = 1e-9});
    REQUIRE(cusping.Hits.size() == 1);
    CHECK(cusping.Hits[0].Kind  == Contact::Touch);

    const double        inside = 4 * 0x1p-30;
    const Intersections secant =
        Intersect(*Suites::Respelled(
                      Line{.Start = {.X = -8, .Y = 5 - inside}, .End = {.X = 8, .Y = 5 - inside}}),
                  *Suites::Respelled(Circle({}, 5)), {.Length = 0x1p-30, .Radians = 1e-9});
    REQUIRE(secant.Hits.size() == 2);
    CHECK(secant.Hits[0].Kind  == Contact::Cross);
    CHECK(secant.Hits[1].Kind  == Contact::Cross);

    /**
     *  Straight pieces crossing: no stretch, a crossing. A segment whose carrier crosses another's
     *  past its end, beyond the tolerance: nothing. A crossing past a shared stretch stays.
     */
    const Intersections crossed =
        Intersect(bent(Line{.Start = {.X = 0, .Y = 0}, .End = {.X = 2, .Y = 2}}),
                  bent(Line{.Start = {.X = 0, .Y = 2}, .End = {.X = 2, .Y = 0}}), tight);
    CHECK(crossed.Overlaps.empty());
    REQUIRE(crossed.Hits.size() == 1);
    CHECK(crossed.Hits[0].Kind  == Contact::Cross);
    CHECK(Intersect(under, Line{.Start = {.X = 5, .Y = -2}, .End = {.X = 6, .Y = -1}}, tight)
              .Hits.empty());
    const Polyline      runs = {.Vertices = {{.Position = {.X = 0, .Y = 0}, .Bulge = 0},
                                             {.Position = {.X = 4, .Y = 0}, .Bulge = 0},
                                             {.Position = {.X = 6, .Y = 2}, .Bulge = 0},
                                             {.Position = {.X = 8, .Y = -2}, .Bulge = 0}},
                                .Closed   = false};
    const Intersections runOn =
        Intersect(runs, Line{.Start = {.X = -1, .Y = 0}, .End = {.X = 10, .Y = 0}}, tight);
    REQUIRE(runOn.Overlaps.size() == 1);
    REQUIRE(runOn.Hits.size()     == 1);
    CHECK(Distance(runOn.Hits[0].At, {.X = 7, .Y = 0}) < 1e-9);

    /**
     *  The cubic out and back on the line, met first: a stretch each way.
     */
    const Intersections backFirst = Intersect(back, line, tight);
    REQUIRE(backFirst.Overlaps.size() == 2);
    CHECK(backFirst.Hits.empty());

    /**
     *  Out and back, met by the line spelled as a Bézier curve: the same two stretches.
     */
    const Intersections backBent = Intersect(back, bent(line), tight);
    REQUIRE(backBent.Overlaps.size() == 2);
    CHECK(backBent.Hits.empty());

    /**
     *  A line shorter than the cubic's passes: each pass shares only the line, its ends found
     *  within the pass.
     */
    const Intersections shorter =
        Intersect(Line{.Start = {.X = 1, .Y = 0}, .End = {.X = 2, .Y = 0}}, back, tight);
    REQUIRE(shorter.Overlaps.size() == 2);
    CHECK(At(back, shorter.Overlaps[0].Second.Start).X == doctest::Approx(1));
    CHECK(At(back, shorter.Overlaps[0].Second.End).X   == doctest::Approx(2));

    /**
     *  A square bulged on one side, its edges tangent to the bulge's circle either side: the
     *  touches where the edges leave the circle go into the stretch.
     */
    const Polyline      bulged = {.Vertices = {{.Position = {.X = 0, .Y = 0}, .Bulge = 0},
                                               {.Position = {.X = 4, .Y = 0}, .Bulge = 1},
                                               {.Position = {.X = 4, .Y = 4}, .Bulge = 0},
                                               {.Position = {.X = 0, .Y = 4}, .Bulge = 0}},
                                  .Closed   = true};
    const Intersections onIt   = Intersect(bulged, Circle({.X = 4, .Y = 2}, 2), tight);
    CHECK(onIt.Hits.empty());
    CHECK(!onIt.Overlaps.empty());

    /**
     *  The bulge first, then the edge leaving its circle: the touch past the stretch's end goes
     *  too.
     */
    const Polyline      after    = {.Vertices = {{.Position = {.X = 4, .Y = 0}, .Bulge = 1},
                                                 {.Position = {.X = 4, .Y = 4}, .Bulge = 0},
                                                 {.Position = {.X = 0, .Y = 4}, .Bulge = 0}},
                                    .Closed   = false};
    const Intersections beyondIt = Intersect(after, Circle({.X = 4, .Y = 2}, 2), tight);
    CHECK(beyondIt.Hits.size() <= 1);
    CHECK(!beyondIt.Overlaps.empty());

    /**
     *  Two circles of one radius 2^-10 apart cross at 2e-4 radians: the curves stay within the
     *  tolerance about each crossing over a band, and it is no tangent; two crossings.
     */
    const Intersections apartBy =
        Intersect(*Suites::Respelled(Circle({}, 5)),
                  *Suites::Respelled(Circle({.X = 0x1p-10, .Y = 0}, 5)), tight);
    CHECK(apartBy.Overlaps.empty());
    REQUIRE(apartBy.Hits.size() == 2);
    CHECK(apartBy.Hits[0].Kind  == Contact::Cross);
    CHECK(apartBy.Hits[1].Kind  == Contact::Cross);

    /**
     *  Arcs of radius 1000 and 1100 through one point at 1e-5 radians: they part, and cross again
     *  0.22 along, each crossing found within a band the curves stay within the tolerance over, no
     *  tangent in it: two crossings.
     */
    const double        tilt  = 1e-5;
    const Point         other = {.X = 1100 * std::sin(tilt), .Y = -1100 * std::cos(tilt)};
    const double        top   = std::atan2(-other.Y, -other.X);
    const Intersections twice = Intersect(
        *Suites::Respelled(Arc{.Center = {.X = 0, .Y = -1000},
                               .Radius = 1000,
                               .Start  = Radians((std::numbers::pi / 2) - 0.01),
                               .Sweep  = Radians(0.02)}),
        *Suites::Respelled(Arc{
            .Center = other, .Radius = 1100, .Start = Radians(top - 0.01), .Sweep = Radians(0.02)}),
        tight);
    REQUIRE(twice.Hits.size() == 2);
    CHECK(twice.Hits[0].Kind  == Contact::Cross);
    CHECK(twice.Hits[1].Kind  == Contact::Cross);

    /**
     *  At 1e-7 radians they cross again only 2.2e-3 along, the curves within 6e-11 between: each
     *  crossing in a band with no tangent, the two one touch.
     */
    const double        slight = 1e-7;
    const Point         near   = {.X = 1100 * std::sin(slight), .Y = -1100 * std::cos(slight)};
    const double        peak   = std::atan2(-near.Y, -near.X);
    const Intersections once   = Intersect(
        *Suites::Respelled(Arc{.Center = {.X = 0, .Y = -1000},
                               .Radius = 1000,
                               .Start  = Radians((std::numbers::pi / 2) - 0.01),
                               .Sweep  = Radians(0.02)}),
        *Suites::Respelled(Arc{
            .Center = near, .Radius = 1100, .Start = Radians(peak - 0.01), .Sweep = Radians(0.02)}),
        tight);
    REQUIRE(once.Hits.size() == 1);
    CHECK(once.Hits[0].Kind  == Contact::Touch);

    /**
     *  The cusp along its tangent with the line the other way: a touch still.
     */
    const Intersections backward =
        Intersect(*Suites::Respelled(Line{.Start = ground.End, .End = ground.Start}),
                  *Suites::Respelled(cusped), {.Length = 0x1p-20, .Radians = 1e-9});
    REQUIRE(backward.Hits.size() == 1);
    CHECK(backward.Hits[0].Kind  == Contact::Touch);

    /**
     *  The cusp turned a little, its tip rounded by its numbers, as a spline cut at the tip: no
     *  sample along the band parallel nor turning either way, but the curve all but stalled at one,
     *  a touch.
     */
    const Matrix        turned = Rotation(Degrees(13));
    const Intersections stalled =
        Intersect(*Suites::Respelled(Line{.Start = turned * Point{.X = -1, .Y = 0},
                                          .End   = turned * Point{.X = 3, .Y = 0}}),
                  *Suites::Respelled(Bezier::Cubic{.Start = turned * Point{.X = 2.25, .Y = -3.375},
                                                   .Out   = turned * Point{.X = -0.75, .Y = 3.375},
                                                   .In    = turned * Point{.X = -0.75, .Y = -3.375},
                                                   .End   = turned * Point{.X = 2.25, .Y = 3.375}}),
                  {.Length = 0x1p-12, .Radians = 1e-9});
    REQUIRE(stalled.Hits.size() == 1);
    CHECK(stalled.Hits[0].Kind  == Contact::Touch);

    /**
     *  A line that ends within the band where it touches a circle: one hit, an end.
     */
    const Intersections endsIn = Intersect(
        Line{.Start = {.X = -8, .Y = 5}, .End = {.X = 5e-5, .Y = 5}}, Circle({}, 5), tight);
    REQUIRE(endsIn.Hits.size() == 1);
    CHECK(endsIn.Hits[0].Kind  == Contact::End);

    /**
     *  With no tolerance at all, two cubics 1e-12 apart all along never lie along each other: the
     *  budget runs out, and says so.
     */
    const Bezier::Cubic wave = Tests::Wave;
    Bezier::Cubic       lift = wave;
    for (Point *pole : {&lift.Start, &lift.Out, &lift.In, &lift.End})
    {
      *pole = *pole + Vector{.X = 0, .Y = 1e-12};
    }
    CHECK(Intersect(wave, lift, {.Length = 0, .Radians = 0}).Spent);
  }

  TEST_CASE(
      "intersect: a cusp along a line is one touch, however turned, moved, spelled or ordered")
  {
    /**
     *  Turned, the cusp's control points round and its tip is a hairpin far narrower than the
     *  tolerance: a band about the tip, the curve all but stalled in it, a touch. As one Bézier
     *  piece against a run, the band was once not seen, and the hit a crossing.
     */
    for (const double degrees : {13.0, 30.0, 45.0, 137.0})
    {
      for (const double away : {0.0, 1e5})
      {
        for (const double length : {0x1p-10, 0x1p-16, 0x1p-22})
        {
          const Matrix map = Translation({.X = away, .Y = away}) * Rotation(Degrees(degrees));
          const Curve  ground =
              Line{.Start = map * Point{.X = -1, .Y = 0}, .End = map * Point{.X = 3, .Y = 0}};
          const Curve     cusped    = Bezier::Cubic{.Start = map * Point{.X = 2.25, .Y = -3.375},
                                                    .Out   = map * Point{.X = -0.75, .Y = 3.375},
                                                    .In    = map * Point{.X = -0.75, .Y = -3.375},
                                                    .End   = map * Point{.X = 2.25, .Y = 3.375}};
          const Tolerance tolerance = {.Length = length, .Radians = 1e-9};
          for (const auto &[one, other] :
               {std::pair{ground, cusped}, std::pair{cusped, ground},
                std::pair{*Suites::Respelled(ground), *Suites::Respelled(cusped)},
                std::pair{*Suites::Respelled(cusped), *Suites::Respelled(ground)}})
          {
            const Intersections meeting = Intersect(one, other, tolerance);
            REQUIRE(meeting.Hits.size() == 1);
            CHECK(meeting.Hits[0].Kind  == Contact::Touch);
          }
        }
      }
    }
  }

  TEST_CASE("intersect: hits a far end of a part tells apart are told so still")
  {
    /**
     *  Cases of the curves audit, met within 2e-10, where a candidate from an end of a part far
     *  from the other curve, polished from a long way off, lands on the far side of a cusp or of
     *  where the curves come nearest, and so tells a touch from a crossing. A line that passes a
     *  cubic's cusp within the tolerance, its exact reference crossing the cubic only elsewhere:
     *  the crossing found at the cusp without the far ends is not found exactly, and they are met
     *  with.
     */
    const Tolerance tolerance = {.Length = 0x1.b7cdfd9d7bdbbp-33, .Radians = 0x1.b7cdfd9d7bdbbp-34};
    const Curve     line = Line{.Start = {.X = 0x1.bfcb195f57ed8p-6, .Y = -0x1.4f0e0cfbe72e0p-6},
                                .End   = {.X = -0x1.1e511a7c43c7bp-3, .Y = -0x1.1b54ccda9150cp+0}};
    const Curve     cusped =
        Bezier::Cubic{.Start = {.X = -0x1.27b6809610988p-1, .Y = -0x1.ce6f8dee1d482p-3},
                      .Out   = {.X = 0x1.5bf8a387fd416p-1, .Y = -0x1.5d14d15474838p-1},
                      .In    = {.X = -0x1.01b351110b279p+0, .Y = -0x1.547e8e15935ccp-1},
                      .End   = {.X = 0x1.1bd4628a017c0p+0, .Y = -0x1.f0c89ae9a1e32p-3}};
    const Intersections past = Intersect(line, cusped, tolerance);
    REQUIRE(past.Hits.size() == 2);
    CHECK(past.Hits[0].Kind  == Contact::Cross);
    CHECK(Distance(past.Hits[0].At, {.X = -0.05096837056473041834, .Y = -0.52936874186740395472})
          < 1e-16);
    CHECK(past.Hits[1].Kind == Contact::Touch);

    /**
     *  An arc and a polyline's arc that cross twice 9e-9 apart at an angle of 1e-8, within the
     *  tolerance of each other all the way between: one touch. Only a candidate from the arc's far
     *  end lands beyond where they come nearest; the crossing found without it is not clear, and
     *  the far ends are met with.
     */
    const Curve arc = Arc{.Center = {.X = -0x1.7e4f739de7086p-2, .Y = 0x1.413440ef31106p-2},
                          .Radius = 0x1.4ffaaebf2d058p-2,
                          .Start  = Degrees(0x1.98p+5),
                          .Sweep  = Degrees(0x1.6p+4)};
    const Curve arced =
        Polyline{.Vertices = {{.Position = {.X = 0x1.21b71acd7911cp-2, .Y = 0x1.36acad9f39913p-3},
                               .Bulge    = 0},
                              {.Position = {.X = 0x1.3bb9e33e86b50p-4, .Y = 0x1.7ad0812236282p-2},
                               .Bulge    = 0x1.5f9f0b0ba1f1ap-3},
                              {.Position = {.X = -0x1.241bfa3a038dep-1, .Y = 0x1.71f1a67b5f5c4p-1},
                               .Bulge    = 0}},
                 .Closed   = false};
    const Intersections twice = Intersect(arc, arced, tolerance);
    REQUIRE(twice.Hits.size() == 1);
    CHECK(twice.Hits[0].Kind  == Contact::Touch);
  }

  TEST_CASE("intersect: a tolerance finer than the curves' numbers hold is taken at their grain")
  {
    /**
     *  At 1e8 a step of a double is 2^-26: two lines 2^-22 apart there are one line to 16 steps, as
     *  they are to a tolerance of 16 steps asked for.
     */
    const Vector far = {.X = 1e8, .Y = 1e8};
    const Curve  one =
        Line{.Start = Point{.X = 0, .Y = 0} + far, .End = Point{.X = 1, .Y = 0} + far};
    const Curve other =
        Line{.Start = Point{.X = 0, .Y = 0x1p-22} + far, .End = Point{.X = 1, .Y = 0x1p-22} + far};
    const Intersections fine  = Intersect(one, other, {.Length = 0x1p-30, .Radians = 1e-9});
    const Intersections grain = Intersect(one, other, {.Length = 0x1p-22, .Radians = 1e-9});
    CHECK(fine.Overlaps.size() == grain.Overlaps.size());
    CHECK(fine.Overlaps.size() == 1);
  }

  TEST_CASE("intersect: a band along a line that ends inside it is the line's end")
  {
    /**
     *  A hump within 1e-3 of the line y = 1 a long way either side of its top: a line that ends
     *  within that stretch meets it at the end, one whole across it touches it at the top.
     */
    const Curve hump = Bezier::Quadratic{
        .Start = {.X = 0, .Y = 0}, .Control = {.X = 1, .Y = 2}, .End = {.X = 2, .Y = 0}};
    const Tolerance loose = {.Length = 1e-3, .Radians = 1e-9};
    for (const Line line : {Line{.Start = {.X = 0.99, .Y = 1}, .End = {.X = 1.01, .Y = 1}},
                            Line{.Start = {.X = -1, .Y = 1}, .End = {.X = 1.01, .Y = 1}},
                            Line{.Start = {.X = 0.99, .Y = 1}, .End = {.X = 3, .Y = 1}}})
    {
      const Intersections meeting = Intersect(line, hump, loose);
      REQUIRE(meeting.Hits.size() == 1);
      CHECK(meeting.Hits[0].Kind  == Contact::End);
    }
    const Intersections whole =
        Intersect(Line{.Start = {.X = -1, .Y = 1}, .End = {.X = 3, .Y = 1}}, hump, loose);
    REQUIRE(whole.Hits.size() == 1);
    CHECK(whole.Hits[0].Kind  == Contact::Touch);
    CHECK(whole.Hits[0].At    == Point{.X = 1, .Y = 1});

    /**
     *  A curve with a number past every double has no size to take the grain from, and meets
     *  nothing.
     */
    const Curve endless = Arc{.Center = {},
                              .Radius = std::numeric_limits<double>::infinity(),
                              .Start  = {},
                              .Sweep  = Turn};
    CHECK(Intersect(Line{.Start = {.X = -1, .Y = 1}, .End = {.X = 3, .Y = 1}}, endless, loose)
              .Hits.empty());
  }

  /**
   *  A tolerance whose length was not a number failed every test of a distance: a cubic and a line
   *  that cross twice met in nothing, and said nothing of it. One without end had every hit near
   *  every other, joined each to each through a search for the nearest point, and a quadratic and
   *  an ellipse with an entry past every double took gigabytes. Each is now met as the tolerance
   *  fit for it, to the bit: not a number or below 0 as 0, past every double as the largest double.
   */
  TEST_CASE("intersect: a tolerance not a number, below 0 or past every double, met as one fit")
  {
    constexpr double Most    = std::numeric_limits<double>::max();
    constexpr double None    = std::numeric_limits<double>::quiet_NaN();
    constexpr double Endless = std::numeric_limits<double>::infinity();
    const Curve      cubic   = Bezier::Cubic{.Start = {.X = 0, .Y = 0},
                                             .Out   = {.X = 1, .Y = 2},
                                             .In    = {.X = 2, .Y = -2},
                                             .End   = {.X = 3, .Y = 0}};
    const Curve      line    = Line{.Start = {.X = -1, .Y = 0.1}, .End = {.X = 4, .Y = 0.1}};
    CHECK(Intersect(cubic, line, {.Length = None, .Radians = 1e-9}).Hits.size() == 2);
    const std::array<std::array<Tolerance, 2>, 7> fits = {{
        {{{.Length = None, .Radians = None}, {.Length = 0, .Radians = 0}}},
        {{{.Length = None, .Radians = 1e-9}, {.Length = 0, .Radians = 1e-9}}},
        {{{.Length = 1e-9, .Radians = None}, {.Length = 1e-9, .Radians = 0}}},
        {{{.Length = -1, .Radians = -1}, {.Length = 0, .Radians = 0}}},
        {{{.Length = -Endless, .Radians = 1e-9}, {.Length = 0, .Radians = 1e-9}}},
        {{{.Length = Endless, .Radians = 1e-9}, {.Length = Most, .Radians = 1e-9}}},
        {{{.Length = 1e-9, .Radians = Endless}, {.Length = 1e-9, .Radians = Most}}},
    }};
    for (const auto &[asked, fit] : fits)
    {
      CAPTURE(asked.Length);
      CAPTURE(asked.Radians);
      CHECK(Alike(Intersect(cubic, line, asked), Intersect(cubic, line, fit)));
      CHECK(Alike(Intersect(line, cubic, asked), Intersect(line, cubic, fit)));
    }

    /**
     *  The quadratic and the ellipse with an entry past every double, found by fuzzing.
     */
    const Curve quadratic =
        Bezier::Quadratic{.Start   = {.X = 0x1.3ca81ed26cb0cp+3, .Y = 0x1.e37792a794394p+2},
                          .Control = {.X = 0x1.249f6bd3cd044p-37, .Y = 0x1.b42159e9eeba8p+19},
                          .End     = {.X = 0x1.3ca81ed26dd56p+2, .Y = 0x1.b4224ba5b80e4p+18}};
    const Curve ellipse = Ellipse{.Map   = {.X      = {.X = 0x1p+4, .Y = -Endless},
                                            .Y      = {.X = 0x1.2f2871bbbefep-13, .Y = -Endless},
                                            .Origin = {.X = 0x1.27536f3e0b3c3p+32, .Y = -0x1.6p+3}},
                                  .Start = Degrees(0x1.f8p+5),
                                  .Sweep = Radians(0x1.4ce527a661c13p+27)};
    CHECK(Alike(Intersect(quadratic, ellipse, {.Length = Endless, .Radians = None}),
                Intersect(quadratic, ellipse, {.Length = Most, .Radians = 0})));
  }

  /**
   *  A line 1e-170 long has a way whose square is nothing, and three of Intersect's steps divided
   *  by it; two such lines, by the product of their lengths. Two lines parallel exactly, met within
   *  an angle whose sine is below 0, were crossed where their carriers cross, at no place; and an
   *  arc that starts at an angle past every double made a count of turns that is no number, cast to
   *  a whole one. The sanitizers stopped at each; the answers are those they gave before.
   */
  TEST_CASE(
      "intersect: a line too short to square and an angle past every double divide by nothing")
  {
    const Curve tiny    = Line{.Start = {.X = 0, .Y = 0}, .End = {.X = 1e-170, .Y = 0}};
    const Curve upright = Line{.Start = {.X = 0, .Y = -1}, .End = {.X = 0, .Y = 1}};
    const Curve dot     = Line{.Start = {.X = 0, .Y = 0}, .End = {.X = 0, .Y = 1e-170}};
    for (const auto &[one, other] : {std::pair{&tiny, &upright}, std::pair{&upright, &tiny}})
    {
      const Intersections met = Intersect(*one, *other, TightMeet);
      REQUIRE(met.Hits.size() == 1);
      CHECK(met.Hits[0].At    == Point{.X = 0, .Y = 0});
      CHECK(met.Hits[0].Kind  == Contact::End);
      CHECK(met.Hits[0].First.Along + met.Hits[0].Second.Along == 0.5);
    }
    const Intersections across = Intersect(tiny, dot, {});
    REQUIRE(across.Hits.size() == 1);
    CHECK(across.Hits[0].At    == Point{.X = 0, .Y = 0});
    CHECK(across.Hits[0].Kind  == Contact::End);
    const Intersections itself = Intersect(tiny, tiny, TightMeet);
    CHECK(itself.Hits.empty());
    CHECK(itself.Overlaps.size() == 1);

    const Tolerance turned = {.Length = 1e-9, .Radians = 4};
    CHECK(Intersect(Line{.Start = {.X = 0, .Y = 0}, .End = {.X = 1, .Y = 0}},
                    Line{.Start = {.X = 0, .Y = 1}, .End = {.X = 1, .Y = 1}}, turned)
              .Hits.empty());
    CHECK(Intersect(Xline{.Through = {.X = 0, .Y = 0}, .Way = {.X = 1, .Y = 0}},
                    Xline{.Through = {.X = 0, .Y = 1}, .Way = {.X = 1, .Y = 0}}, turned)
              .Hits.empty());

    const double endless = std::numeric_limits<double>::infinity();
    const Curve  lost =
        Arc{.Center = {}, .Radius = 1, .Start = Radians(endless), .Sweep = Radians(1)};
    const Curve found = Arc{.Center = {}, .Radius = 1, .Start = Radians(0.5), .Sweep = Radians(2)};
    for (const Tolerance &tolerance : {TightMeet, Tolerance{}})
    {
      const Intersections forth = Intersect(lost, found, tolerance);
      const Intersections back  = Intersect(found, lost, tolerance);
      CHECK(forth.Hits.empty());
      CHECK(forth.Overlaps.empty());
      CHECK(back.Hits.empty());
      CHECK(back.Overlaps.empty());
    }
  }

  /**
   *  Two endless lines crossing at 4e-13 rad, well inside the drawing, stay within the tolerance of
   *  each other for some 500 units either way and then part: the band about a shallow crossing, one
   *  touch, at the crossing of their carriers. The thorough way worked that crossing out in
   *  doubles, its two places too far apart along lines so near parallel to lie within the
   *  tolerance, and met them in nothing. Rays and lines that end past the band meet as they do;
   *  ones that end within it meet at that end, or share it, as before. Each point owed is the exact
   *  crossing of the carriers, rounded, worked out in rationals.
   */
  TEST_CASE("intersect: two lines at a tiny angle touch once, where their carriers cross")
  {
    const Tolerance tolerance = {.Length = 0x1.b7cdfd9d7bdbbp-33, .Radians = 0x1.b7cdfd9d7bdbbp-34};
    const Point     p         = {.X = -0x1.c25037bc7227cp-2, .Y = -0x1.98fc0bcfd416cp-2};
    const Point     q         = {.X = -0x1.4e8239e329a0ep-2, .Y = -0x1.e0217d23c6b92p-3};
    const Vector    u         = {.X = 0x1.2181c7d6d9ad6p-1, .Y = 0x1.a64a819d555bap-1};
    const Vector    v         = {.X = 0x1.2181c7d6da673p-1, .Y = 0x1.a64a819d54dc4p-1};
    const auto      touches   = [&tolerance](const Curve &one, const Curve &other) -> Points2
    {
      return {Touched(Intersect(one, other, tolerance)), Touched(Intersect(other, one, tolerance))};
    };
    const Point crossing = {.X = -0x1.df3dcd586f0ecp-2, .Y = -0x1.c32e3eee8dcf3p-2};
    CHECK(touches(Xline{.Through = p, .Way = u}, Xline{.Through = q, .Way = v})
          == Points2{crossing, crossing});
    const Point back  = {.X = -0x1.218fda589790fp+11, .Y = -0x1.a657497db3da4p+11};
    const Point other = {.X = -0x1.218c3be8a9808p+11, .Y = -0x1.a6520223496b6p+11};
    const Point owed  = {.X = -0x1.d7937ede7f7f3p-2, .Y = -0x1.b7ffea6d04358p-2};
    CHECK(touches(Ray{.Start = back, .Way = u}, Ray{.Start = other, .Way = v})
          == Points2{owed, owed});
    CHECK(touches(
              Line{.Start = back, .End = {.X = 0x1.2173b5551bc9dp+11, .Y = 0x1.a63db9bcf6dd0p+11}},
              Line{.Start = other, .End = {.X = 0x1.217753c50b4dep+11, .Y = 0x1.a6430117604d2p+11}})
          == Points2{owed, owed});

    /**
     *  Rays from the same starts running away from the crossing never come within the tolerance;
     *  endless lines whose crossing lies so far out that the doubles there are coarser than the
     *  tolerance have no places within it of each other there, and meet in nothing, as before.
     */
    CHECK(Intersect(Ray{.Start = back, .Way = -u}, Ray{.Start = other, .Way = -v}, tolerance)
              .Hits.empty());
    const Point from    = {.X = 0x1.682aa62147d0ap+2, .Y = -0x1.07eb82081babp+3};
    const Point onto    = {.X = -0x1.203a74db9ab6p-2, .Y = -0x1.bfb27589cbb44p+2};
    const Point start   = {.X = 0x1.57003cb71dc99p+0, .Y = -0x1.caa913216de52p+2};
    const Point end     = {.X = -0x1.bb43b40529934p-2, .Y = -0x1.b29e1b92b4016p+2};
    const Curve distant = Xline{.Through = from, .Way = onto - from};
    const Curve farther = Xline{.Through = start, .Way = end - start};
    CHECK(Intersect(distant, farther, {.Length = 0.0625, .Radians = 1e-9}).Hits.empty());

    const Curve         near   = Ray{.Start = p - (u * 4.0), .Way = u};
    const Curve         beside = Ray{.Start = q - (v * 4.0), .Way = v};
    const Intersections ended  = Intersect(near, beside, tolerance);
    REQUIRE(ended.Hits.size() == 1);
    CHECK(ended.Hits[0].Kind  == Contact::End);
    const Curve shorter = Line{.Start = q - (v * 4.0), .End = q + (v * 4.0)};
    for (const Curve &one : {Curve{Line{.Start = p - (u * 4.0), .End = p + (u * 4.0)}},
                             Curve{Xline{.Through = p, .Way = u}}})
    {
      const Intersections shared = Intersect(one, shorter, tolerance);
      CHECK(shared.Hits.empty());
      CHECK(shared.Overlaps.size() == 1);
    }
  }

  /**
   *  Two hits found apart and each put where the curves meet exactly, to the same point and places,
   *  are one hit there: a bulge starting where a polyline passes twice, tangent there to the
   *  polyline's second pass, met one end, a crossing and the other end, an end if either is.
   */
  TEST_CASE("intersect: hits put where another was put are one")
  {
    const Tolerance tolerance = {.Length = 0x1.b7cdfd9d7bdbbp-34, .Radians = 0x1.12e0be826d695p-30};
    const Curve     bulge = Polyline{.Vertices = {{.Position = {.X = 1, .Y = -1}, .Bulge = -0.125},
                                                  {.Position = {.X = 3, .Y = -3}}},
                                     .Closed   = false};
    const Curve     twice = Polyline{.Vertices = {{.Position = {.X = -4, .Y = 1}},
                                                  {.Position = {.X = 1, .Y = -1}, .Bulge = -1},
                                                  {.Position = {.X = 2, .Y = -2}},
                                                  {.Position = {.X = -4, .Y = -1}},
                                                  {.Position = {.X = 1, .Y = -1}, .Bulge = -0.125},
                                                  {.Position = {.X = 2, .Y = -2}}},
                                     .Closed   = false};
    const Point     start = {.X = 1, .Y = -1};
    const Point     crossing = {.X = 0x1.0df6fc8240df7p+1, .Y = -0x1.dbf209037dbf2p+0};
    const Prepared  one      = Prepared{bulge};
    const Prepared  other    = Prepared{twice};
    for (const bool flip : {false, true})
    {
      const Curves         curves = flip ? Curves{twice, bulge} : Curves{bulge, twice};
      const Prepared::Pair pair   = flip ? Prepared::Pair{std::cref(other), std::cref(one)}
                                         : Prepared::Pair{std::cref(one), std::cref(other)};
      Intersections        ready;
      Intersect(pair, tolerance, ready);
      for (const Intersections &met : {Intersect(curves[0], curves[1], tolerance), ready})
      {
        REQUIRE(met.Hits.size() == 3);
        const auto passes = [&met, &start, flip](double along) -> bool
        {
          return std::ranges::count_if(met.Hits,
                                       [&start, flip, along](const Hit &hit) -> bool
                                       {
                                         const double place =
                                             flip ? hit.First.Along : hit.Second.Along;
                                         return hit.At == start
                                             && hit.Kind == Contact::End
                                             && std::equal_to<double>{}(place, along);
                                       })
              == 1;
        };
        CHECK(passes(1.0));
        CHECK(passes(4.0));
        CHECK(std::ranges::count_if(met.Hits, [&crossing](const Hit &hit) -> bool
                                    { return hit.At == crossing && hit.Kind == Contact::Cross; })
              == 1);
        CHECK(std::ranges::is_sorted(met.Hits, {},
                                     [](const Hit &hit) -> double { return hit.First.Along; }));
      }
    }
  }

  /**
   *  Hits the thorough way put past the next one along the first curve, two random pairs shrunk: a
   *  bulge ending where another polyline passes twice, 1e8 off, its end met at the first pass and
   *  put there exactly, past the second pass's hit left where it was found a step of a double short
   *  of the end; and a polyline touching another's bulge, the run's touch where the curves come
   *  nearest lying past the next run's at a vertex. Each meeting's hits in order along the first
   *  curve, as many as before, plain and ready, either first.
   */
  TEST_CASE("intersect: hits put past the next one along the first curve are put in order")
  {
    const Curve ending = Polyline{
        .Vertices = {{.Position = {.X = 0x1.7d7840cp+26, .Y = -0x1.8p+1}, .Bulge = -0x1p-3},
                     {.Position = {.X = 0x1.7d783fcp+26, .Y = 0x1p+2}, .Bulge = -0x1p-1}},
        .Closed   = false};
    const Curve twice =
        Polyline{.Vertices = {{.Position = {.X = 0x1.7d783fcp+26, .Y = 0x1p+2}},
                              {.Position = {.X = 0x1.7d784p+26, .Y = -0x1p+0}, .Bulge = -0x1p+1},
                              {.Position = {.X = 0x1.7d783fcp+26, .Y = 0x1p+2}, .Bulge = 0x1p+1}},
                 .Closed   = false};
    const Curve touching =
        Polyline{.Vertices = {{.Position = {.X = -0x1.f8552790ecd6p-2, .Y = -0x1.87855fcf46714p-1}},
                              {.Position = {.X = 0x1.df2a7a1af3a94p-2, .Y = -0x1.480909780ff2p-2},
                               .Bulge    = 0x1.b1f77b3e10e08p-2},
                              {.Position = {.X = 0x1.992342a96bc5cp-2, .Y = -0x1.850f8f111e278p-2},
                               .Bulge    = 0x1.03d42d385c8d8p-1}},
                 .Closed   = false};
    const Curve touched = Polyline{
        .Vertices = {{.Position = {.X = 0x1.7dde28e45e874p-3, .Y = -0x1.44da26ee41ecp-6},
                      .Bulge    = -0x1.8287575117cccp-2},
                     {.Position = {.X = 0x1.f156ff0a7df08p-2, .Y = -0x1.78e37678ae4eep-1}},
                     {.Position = {.X = 0x1.bc26de622fb78p-2, .Y = -0x1.668c4c44970ccp-2},
                      .Bulge    = 0x1.18fcf5111c794p+1},
                     {.Position = {.X = -0x1.e6d83c9ad128p-6, .Y = -0x1.67a7e5a9b22e2p-1}}},
        .Closed   = false};

    struct Pairing
    {
      const Curve *One   = nullptr;
      const Curve *Other = nullptr;
      Tolerance    Within;
      std::size_t  Hits = 0;
    };

    for (const Pairing &pairing :
         {Pairing{.One    = &ending,
                  .Other  = &twice,
                  .Within = {.Length = 0x1p-9, .Radians = 0x1.0624dd2f1a9fcp-10},
                  .Hits   = 2},
          Pairing{.One    = &touching,
                  .Other  = &touched,
                  .Within = {.Length = 0x1p-6, .Radians = 0x1.12e0be826d695p-30},
                  .Hits   = 3}})
    {
      const Prepared one   = Prepared{*pairing.One};
      const Prepared other = Prepared{*pairing.Other};
      for (const bool flip : {false, true})
      {
        const Curves curves =
            flip ? Curves{*pairing.Other, *pairing.One} : Curves{*pairing.One, *pairing.Other};
        const Prepared::Pair pair = flip ? Prepared::Pair{std::cref(other), std::cref(one)}
                                         : Prepared::Pair{std::cref(one), std::cref(other)};
        Intersections        ready;
        Intersect(pair, pairing.Within, ready);
        for (const Intersections &met : {Intersect(curves[0], curves[1], pairing.Within), ready})
        {
          CHECK(met.Hits.size() == pairing.Hits);
          CHECK(std::ranges::is_sorted(met.Hits, {},
                                       [](const Hit &hit) -> double { return hit.First.Along; }));
        }
      }
    }
  }

  /**
   *  A curve's end within the tolerance of two passes of a polyline, a random pair shrunk, 1e8 off:
   *  a bulge starting on a polyline's straight segment, where the polyline's bulge before it comes
   *  back within the tolerance, the bulge's circle crossing it a little before the bulge starts.
   *  The thorough way took one candidate for the end, the nearest place, and nothing else found the
   *  farther pass but a far end's candidate polished there, a crossing whose place on the bulge lay
   *  past the tolerance. An end at each pass, both points of every hit within the tolerance of each
   *  other, plain and ready, either first.
   */
  /**
   *  A bulge starting a quarter along a straight segment, a random pair shrunk: the crossing is its
   *  start but for the rounding, and its place on the bulge was told by the angle from the chord to
   *  the way from the start to the crossing, a few roundings long and its angle noise, which put it
   *  0.466 along the bulge, its point 1.94 from the hit. Now at the start, both points of the hit
   *  within the tolerance of each other, plain and ready, either first.
   */
  TEST_CASE("intersect: a crossing a rounding from a bulge's start is placed at the start")
  {
    const Curve straight =
        Polyline{.Vertices = {{.Position = {.X = -0x1.50a38e14e3314p-2, .Y = 0x1.d0f4d86395058p-1}},
                              {.Position = {.X = 0x1.d0dadef4085cp-5, .Y = 0x1.a8b797dc43572p-1}}},
                 .Closed   = false};
    const Curve bow =
        Polyline{.Vertices = {{.Position = {.X = -0x1.dbe7a73014442p-3, .Y = 0x1.c6e58841c099ep-1},
                               .Bulge    = -0x1.6e00aaa341a48p+1},
                              {.Position = {.X = 0x1.8eb8fd162cb2p-4, .Y = -0x1.9204d9d61912p-2}}},
                 .Closed   = false};
    const Tolerance  within = {.Length = 0x1.19799812dea11p-40, .Radians = 0x1.0624dd2f1a9fcp-10};
    const std::array ready  = {Prepared{straight}, Prepared{bow}};
    for (const bool flip : {false, true})
    {
      const Curves         curves = flip ? Curves{bow, straight} : Curves{straight, bow};
      const Prepared::Pair pair   = flip ? Prepared::Pair{std::cref(ready[1]), std::cref(ready[0])}
                                         : Prepared::Pair{std::cref(ready[0]), std::cref(ready[1])};
      Intersections        readied;
      Intersect(pair, within, readied);
      for (const Intersections &met : {Intersect(curves[0], curves[1], within), readied})
      {
        REQUIRE(met.Hits.size() == 1);
        const Hit &hit = met.Hits.front();
        CHECK(hit.Kind == Contact::End);
        CHECK((flip ? hit.First.Along : hit.Second.Along) < 1e-15);
        CHECK(Distance(hit.First.At, hit.Second.At) <= within.Length);
      }
    }
  }

  TEST_CASE(
      "intersect: a curve's end within the tolerance of two passes of a polyline is met at each")
  {
    const Curve spike = Polyline{
        .Vertices = {{.Position = {.X = 0x1.7d783ffef71a5p+26, .Y = 0x1.51b8b2bc1de04p-1},
                      .Bulge    = 0x1.691a4a084a178p-2},
                     {.Position = {.X = 0x1.7d78400538fc2p+26, .Y = 0x1.76a3c3eb43934p-2}},
                     {.Position = {.X = 0x1.7d783ff2bb4afp+26, .Y = 0x1.4cf16eb20d0c6p-1}}},
        .Closed   = false};
    const Curve bow = Polyline{
        .Vertices = {{.Position = {.X = 0x1.7d784000998fdp+26, .Y = 0x1.bf738a497934ap-2},
                      .Bulge    = 0x1.ec4d47d8f9694p-2},
                     {.Position = {.X = 0x1.7d783fd48137ap+26, .Y = 0x1.b4785913511acp-2}}},
        .Closed   = false};
    const Tolerance  within = {.Length = 0x1p-9, .Radians = 0x1.0624dd2f1a9fcp-10};
    const std::array ready  = {Prepared{spike}, Prepared{bow}};
    for (const bool flip : {false, true})
    {
      const Curves         curves = flip ? Curves{bow, spike} : Curves{spike, bow};
      const Prepared::Pair pair   = flip ? Prepared::Pair{std::cref(ready[1]), std::cref(ready[0])}
                                         : Prepared::Pair{std::cref(ready[0]), std::cref(ready[1])};
      Intersections        readied;
      Intersect(pair, within, readied);
      for (const Intersections &met : {Intersect(curves[0], curves[1], within), readied})
      {
        const auto places = [flip](const Hit &hit) -> std::array<double, 2>
        {
          return flip ? std::array{hit.Second.Along, hit.First.Along}
                      : std::array{hit.First.Along, hit.Second.Along};
        };
        const auto passed = [&places](const Hit &hit, Range spiked) -> bool
        {
          const auto [along, bowed] = places(hit);
          return hit.Kind     == Contact::End
              && spiked.Start < along
              && along        < spiked.End
              && bowed        < 1e-6;
        };
        CHECK(std::ranges::count_if(met.Hits, [&passed](const Hit &hit) -> bool
                                    { return passed(hit, {.Start = 0.69, .End = 0.7}); })
              == 1);
        CHECK(std::ranges::count_if(met.Hits, [&passed](const Hit &hit) -> bool
                                    { return passed(hit, {.Start = 1.24, .End = 1.26}); })
              == 1);
        for (const Hit &hit : met.Hits)
        {
          CHECK(Distance(hit.First.At, hit.Second.At) <= within.Length);
        }
      }
    }
  }

  /**
   *  A random pair (6466) shrunk by 2^40: a line from a closed polyline's first vertex crosses its
   *  segment from the third vertex nearly square, at 88.8 degrees, the segment fifty tolerances
   *  long and its bulge too flat to bend it. The line's stretch about the crossing that the
   *  segment's shadow covers is a little longer than the tolerance and lies within it of the
   *  segment's line, and the two were taken to share it, though the segment runs far from the line:
   *  the line first, the crossing was lost. Either first, a crossing, plain and ready.
   */
  TEST_CASE("intersect: a line crossing a long segment nearly square crosses it, either first")
  {
    const Curve line = Line{.Start = {.X = -0x1.e46246d8642eap-1, .Y = -0x1.6e9fdf326ca9p-3},
                            .End   = {.X = -0x1.3281ad407123p-3, .Y = 0x1.ad243009ce18p-7}};
    const Curve ring =
        Polyline{.Vertices = {{.Position = {.X = -0x1.e46246d8642eap-1, .Y = -0x1.6e9fdf326ca9p-3}},
                              {.Position = {.X = -0x1.a3175aeb31b04p-2, .Y = -0x1.e410ec2c83504p-2},
                               .Bulge    = 0x1.1bcba3b5f4845p+0},
                              {.Position = {.X = -0x1.94e77ec446a28p-1, .Y = 0x1.330a7acd57aa8p-3},
                               .Bulge    = 0x1.57505665e768cp-22},
                              {.Position = {.X = -0x1.2ea0539f24ecp-1, .Y = -0x1.375060c97106ap-1}},
                              {.Position = {.X = 0x1.a92e9191fb2e8p-3, .Y = -0x1.0e51d1ceac3a8p-1},
                               .Bulge    = -0x1.04165f78fb36ap-1},
                              {.Position = {.X = 0x1.6562c8c86878cp-1, .Y = -0x1.0568f025158d8p-1}},
                              {.Position = {.X = 0x1.b1ae152bd686ap-1, .Y = 0x1.a3acf2b7f5a6cp-1},
                               .Bulge    = -0x1.2462803911a1p+0},
                              {.Position = {.X = 0x1.60e578be96c04p-2, .Y = 0x1.cefe19445c7cp-6}}},
                 .Closed   = true};
    const Tolerance  within = {.Length = 0x1p-6, .Radians = 0x1.0624dd2f1a9fcp-10};
    const std::array ready  = {Prepared{line}, Prepared{ring}};
    for (const bool flip : {false, true})
    {
      const Curves         curves = flip ? Curves{ring, line} : Curves{line, ring};
      const Prepared::Pair pair   = flip ? Prepared::Pair{std::cref(ready[1]), std::cref(ready[0])}
                                         : Prepared::Pair{std::cref(ready[0]), std::cref(ready[1])};
      Intersections        readied;
      Intersect(pair, within, readied);
      for (const Intersections &met : {Intersect(curves[0], curves[1], within), readied})
      {
        CHECK(met.Hits.size() == 3);
        CHECK(met.Overlaps.empty());
        CHECK(std::ranges::count_if(met.Hits,
                                    [flip](const Hit &hit) -> bool
                                    {
                                      const double lined =
                                          flip ? hit.Second.Along : hit.First.Along;
                                      const double ringed =
                                          flip ? hit.First.Along : hit.Second.Along;
                                      return hit.Kind                  == Contact::Cross
                                          && std::abs(lined - 0.2856)  < 1e-4
                                          && std::abs(ringed - 2.3616) < 1e-4;
                                    })
              == 1);
        /**
         *  The segment's bulge is slight: the crossing's point on it is told from how far along its
         *  chord the crossing lies, and was 8.7e-11 off, told from an angle of 1e-6.
         */
        for (const Hit &hit : met.Hits)
        {
          CHECK(Distance(hit.First.At, hit.Second.At) <= 0x1p-50);
        }
      }
    }
  }

  /**
   *  A line along a polyline's first segment, the polyline coming back across it: the crossing of
   *  the later pass lies in the stretch they share along the line, and outside it along the
   *  polyline. Judged by the line's place alone, it was the stretch's and dropped with the line
   *  first; a candidate is the stretch's only where both its places are. Either first, the stretch
   *  and the crossing, plain and ready.
   */
  TEST_CASE("intersect: another pass's crossing in a shared stretch is a hit, either first")
  {
    const Curve      line   = Line{.Start = {.X = 0, .Y = 0}, .End = {.X = 4, .Y = 0}};
    const Curve      back   = Polyline{.Vertices = {{.Position = {.X = 0, .Y = 0}},
                                                    {.Position = {.X = 4, .Y = 0}},
                                                    {.Position = {.X = 2, .Y = 1}},
                                                    {.Position = {.X = 1, .Y = -1}}},
                                       .Closed   = false};
    const Tolerance  within = {.Length = 1e-9, .Radians = 1e-9};
    const std::array ready  = {Prepared{line}, Prepared{back}};
    for (const bool flip : {false, true})
    {
      const Curves         curves = flip ? Curves{back, line} : Curves{line, back};
      const Prepared::Pair pair   = flip ? Prepared::Pair{std::cref(ready[1]), std::cref(ready[0])}
                                         : Prepared::Pair{std::cref(ready[0]), std::cref(ready[1])};
      Intersections        readied;
      Intersect(pair, within, readied);
      for (const Intersections &met : {Intersect(curves[0], curves[1], within), readied})
      {
        REQUIRE(met.Hits.size()   == 1);
        CHECK(met.Overlaps.size() == 1);
        const Hit &hit = met.Hits[0];
        CHECK(hit.Kind == Contact::Cross);
        CHECK(hit.At   == Point{.X = 1.5, .Y = 0});
        CHECK((flip ? hit.First.Along : hit.Second.Along) == 2.5);
        CHECK((flip ? hit.Second.Along : hit.First.Along) == 0.375);
      }
    }
  }

  /**
   *  A random pair (82939, 1e8 off): the polylines share a stretch, the first's first segment along
   *  the second's fourth, and the first's second segment crosses that stretch. With the second
   *  first, that crossing, inside the stretch along the curve met first, and a crossing past the
   *  stretch's end, the curves together between them only along the stretch, were one hit, the
   *  farther lost; hits are never one across an end of a shared stretch. Either first, ten hits.
   */
  TEST_CASE("intersect: hits either side of a shared stretch's end are two, either first")
  {
    const Curve one = Polyline{
        .Vertices = {{.Position = {.X = 0x1.7d78400b53764p+26, .Y = -0x1.0708a49ce7328p-1}},
                     {.Position = {.X = 0x1.7d783fc912f9cp+26, .Y = 0x1.0897a1d6e07b2p-1},
                      .Bulge    = -0x1.a5ad79aab884cp-3},
                     {.Position = {.X = 0x1.7d783feb8da11p+26, .Y = -0x1.850ab0c563f22p-1},
                      .Bulge    = 0x1.9bfa403ae15fp+1},
                     {.Position = {.X = 0x1.7d783ffffdadap+26, .Y = -0x1.7ef7b7e225d4ap-1},
                      .Bulge    = -0x1.2848bd77c838p-5},
                     {.Position = {.X = 0x1.7d78401e40b0cp+26, .Y = 0x1.f5af501a4a3f8p-1}}},
        .Closed   = false};
    const Curve other = Polyline{
        .Vertices = {{.Position = {.X = 0x1.7d783fddabebap+26, .Y = 0x1.3e1b1be93938ep-1},
                      .Bulge    = 0x1.fe7f5f782b14ep-1},
                     {.Position = {.X = 0x1.7d78402690f08p+26, .Y = -0x1.b2d90b7729936p-1},
                      .Bulge    = -0x1.97532aff0c42cp-2},
                     {.Position = {.X = 0x1.7d783fff1193fp+26, .Y = 0x1.3c45c83e08a82p-1},
                      .Bulge    = 0x1.2ed808e8c87cfp+0},
                     {.Position = {.X = 0x1.7d783fc912f9cp+26, .Y = 0x1.0897a1d6e07b2p-1},
                      .Bulge    = 0},
                     {.Position = {.X = 0x1.7d783ffac3572p+26, .Y = -0x1.064125ffea8e3p-2},
                      .Bulge    = 0x1.7316829c38a6cp-1},
                     {.Position = {.X = 0x1.7d783fcc8f283p+26, .Y = -0x1.a8f2ae27b40b8p-2},
                      .Bulge    = 0x1.496c6e9a465p-8},
                     {.Position = {.X = 0x1.7d783feb8da11p+26, .Y = -0x1.850ab0c563f22p-1},
                      .Bulge    = 0},
                     {.Position = {.X = 0x1.7d783fc08349cp+26, .Y = -0x1.e2c97b0d3aa26p-1}}},
        .Closed   = false};
    const Tolerance  within = {.Length = 0x1.19799812dea11p-40, .Radians = 0x1.12e0be826d695p-30};
    const std::array ready  = {Prepared{one}, Prepared{other}};
    for (const bool flip : {false, true})
    {
      const Curves         curves = flip ? Curves{other, one} : Curves{one, other};
      const Prepared::Pair pair   = flip ? Prepared::Pair{std::cref(ready[1]), std::cref(ready[0])}
                                         : Prepared::Pair{std::cref(ready[0]), std::cref(ready[1])};
      Intersections        readied;
      Intersect(pair, within, readied);
      for (const Intersections &met : {Intersect(curves[0], curves[1], within), readied})
      {
        CHECK(met.Hits.size()     == 10);
        CHECK(met.Overlaps.size() == 1);
      }
    }
  }

  /**
   *  Endless lines through points 3 and 5 units either side of their crossing, at angles from 1e-9
   *  down to 1e-15 rad, the tolerance 2e-10 and its angle 1e-10: one hit each, at the exact
   *  crossing of their carriers, a crossing where they part by more than the angle and a touch
   *  where they part by less. Lines parallel exactly share all of themselves within the tolerance,
   *  and nothing past it.
   */
  TEST_CASE("intersect: endless lines at angles down to 1e-15 rad cross or touch once")
  {
    struct Row
    {
      Vector  Way     = {};             /**< The second line's way.      */
      Point   Through = {};             /**< A point of the second line. */
      Point   At      = {};             /**< The crossing owed.          */
      Contact Kind    = Contact::Cross; /**< How they meet there.        */
    };

    const std::array rows = {
        Row{.Way     = {.X = 0x1.3333332c53fb4p-1, .Y = 0x1.9999999ec103ap-1},
            .Through = {.X = -0x1.4ffffff768fa1p+1, .Y = -0x1.5000000338a24p+2},
            .At      = {.X = 0x1.8000039387005p-2, .Y = -0x1.3ffffeced2ffdp+0},
            .Kind    = Contact::Cross},
        Row{.Way     = {.X = 0x1.333333322b513p-1, .Y = 0x1.9999999a5f832p-1},
            .Through = {.X = -0x1.4ffffffeb6258p+1, .Y = -0x1.500000007bb1fp+2},
            .At      = {.X = 0x1.7fffffffffffep-2, .Y = -0x1.3ffffffffffffp+0},
            .Kind    = Contact::Cross},
        Row{.Way     = {.X = 0x1.33333333219b9p-1, .Y = 0x1.99999999a6cb5p-1},
            .Through = {.X = -0x1.4fffffffea027p+1, .Y = -0x1.50000000083f1p+2},
            .At      = {.X = 0x1.8003d779708d7p-2, .Y = -0x1.3ffeb82cda7b7p+0},
            .Kind    = Contact::Touch},
        Row{.Way     = {.X = 0x1.333333333170dp-1, .Y = 0x1.999999999aeb6p-1},
            .Through = {.X = -0x1.4ffffffffdcd0p+1, .Y = -0x1.5000000000d32p+2},
            .At      = {.X = 0x1.8045d9645709ap-2, .Y = -0x1.3fe8b78938520p+0},
            .Kind    = Contact::Touch},
        Row{.Way     = {.X = 0x1.3333333333062p-1, .Y = 0x1.9999999999bb6p-1},
            .Through = {.X = -0x1.4fffffffffc7ap+1, .Y = -0x1.5000000000152p+2},
            .At      = {.X = 0x1.83461da59ab0cp-2, .Y = -0x1.3ee8a0c8cc6fap+0},
            .Kind    = Contact::Touch},
        Row{.Way     = {.X = 0x1.33333333332ebp-1, .Y = 0x1.99999999999d0p-1},
            .Through = {.X = -0x1.4ffffffffffa6p+1, .Y = -0x1.5000000000022p+2},
            .At      = {.X = 0x1.883126e978d4cp-2, .Y = -0x1.3d44f3078263ap+0},
            .Kind    = Contact::Touch},
        Row{.Way     = {.X = 0x1.333333333332cp-1, .Y = 0x1.999999999999fp-1},
            .Through = {.X = -0x1.4fffffffffff7p+1, .Y = -0x1.5000000000003p+2},
            .At      = {.X = 0x1.8e49d2c0e49d2p-2, .Y = -0x1.3b3cb9bfb3cb8p+0},
            .Kind    = Contact::Touch},
    };
    const Tolerance tolerance = {.Length = 2e-10, .Radians = 1e-10};
    const Curve     level = Xline{.Through = {.X = 0x1.1666666666666p+1, .Y = 0x1.2666666666668p+0},
                                  .Way     = {.X = 0.6, .Y = 0.8}};
    for (const Row &row : rows)
    {
      CAPTURE(row.Way.Y);
      const Curve tilted = Xline{.Through = row.Through, .Way = row.Way};
      for (const Intersections &met :
           {Intersect(level, tilted, tolerance), Intersect(tilted, level, tolerance)})
      {
        REQUIRE(met.Hits.size() == 1);
        CHECK(met.Overlaps.empty());
        CHECK(met.Hits[0].At   == row.At);
        CHECK(met.Hits[0].Kind == row.Kind);
      }
    }

    const Vector way = {.X = 3, .Y = 4};
    const Curve  one = Xline{.Through = {.X = 0, .Y = 0}, .Way = way};
    for (const double gap : {1e-10, 4e-10})
    {
      const Curve parallel    = Xline{.Through = {.X = -0.8 * gap, .Y = 0.6 * gap}, .Way = way};
      const Intersections met = Intersect(one, parallel, tolerance);
      CHECK(met.Hits.empty());
      CHECK(met.Overlaps.size() == (gap < tolerance.Length ? 1 : 0));
    }
  }

  /**
   *  A quarter of an ellipse and the same quarter of one 1e-13 to 1e-11 larger and moved as far, a
   *  tolerance of 2e-10: as ellipses both they met in one hit at an end, as an ellipse and a spline
   *  likewise, and only as splines both did they share the stretch. Two ellipses shared stretches
   *  only where their maps were equal to the bit, and a piece lay on an ellipse only where the
   *  ellipse's implicit equation was nothing on it to 64 steps of a double. Every spelling now
   *  shares the one stretch; four tolerances apart, none does.
   */
  TEST_CASE("intersect: ellipses a hair apart share their stretch, however spelled")
  {
    const Tolerance tolerance = {.Length = 2e-10, .Radians = 1e-10};
    const auto      quarters  = [](double gap) -> std::array<Curve, 2>
    {
      const double size = 1.0 + gap;
      const Matrix map  = {
          .X = {.X = 2 * size, .Y = 0}, .Y = {.X = 0, .Y = size}, .Origin = {.X = gap, .Y = -gap}};
      return {Ellipse{.Map = map, .Start = Degrees(0), .Sweep = Degrees(90)},
              Spline{.Degree  = 2,
                     .Knots   = {0, 0, 0, 1, 1, 1},
                     .Poles   = {{.X = (2 * size) + gap, .Y = -gap},
                                 {.X = (2 * size) + gap, .Y = size - gap},
                                 {.X = gap, .Y = size - gap}},
                     .Weights = {1, std::sqrt(0.5), 1}}};
    };
    const std::array<Curve, 2> exact = quarters(0.0);
    const Curve whole = Ellipse{.Map = std::get<Ellipse>(exact[0]).Map, .Start = {}, .Sweep = Turn};
    for (const double gap : {1e-13, 1e-12, 1e-11, 8e-10})
    {
      CAPTURE(gap);
      const std::size_t          owed = gap < tolerance.Length ? 1 : 0;
      const std::array<Curve, 2> near = quarters(gap);
      for (const Curve &one : exact)
      {
        for (const Curve &other : near)
        {
          CHECK(Intersect(one, other, tolerance).Overlaps.size() == owed);
          CHECK(Intersect(other, one, tolerance).Overlaps.size() == owed);
        }
      }
      const double        size  = 1.0 + gap;
      const Curve         arc   = Ellipse{.Map   = {.X      = {.X = 2 * size, .Y = 0},
                                                    .Y      = {.X = 0, .Y = size},
                                                    .Origin = {.X = gap, .Y = -gap}},
                                          .Start = Degrees(10),
                                          .Sweep = Degrees(60)};
      const Intersections along = Intersect(whole, arc, tolerance);
      REQUIRE(along.Overlaps.size() == owed);
      CHECK((owed == 0 || along.Hits.empty()));
      CHECK((owed == 0 || along.Overlaps[0].Second == Range{.Start = 0, .End = 1}));
      CHECK(Intersect(whole, near[1], tolerance).Overlaps.size() == owed);
    }

    /**
     *  The quarter 1.4e-10 larger: within the tolerance at its middle, where the nearest point is
     *  nearer than the one toward the center, and 2.8e-10 off at its end on the long axis, so along
     *  it nowhere, as two splines tell. A piece through the center is along it nowhere too.
     */
    const std::array<Curve, 2> plain   = quarters(0.0);
    const double               grown   = 1.0 + 1.4e-10;
    const Curve                swollen = Spline{
        .Degree  = 2,
        .Knots   = {0, 0, 0, 1, 1, 1},
        .Poles   = {{.X = 2 * grown, .Y = 0}, {.X = 2 * grown, .Y = grown}, {.X = 0, .Y = grown}},
        .Weights = {1, std::sqrt(0.5), 1}};
    CHECK(Intersect(whole, swollen, tolerance).Overlaps.empty());
    CHECK(Intersect(plain[1], swollen, tolerance).Overlaps.empty());
    /**
     *  An arc 1e-170 across has a map whose entries square to nothing: how much it stretches is
     *  worked out at a size near 1, and a cubic through it meets it at its end, as before.
     */
    const Curve speck = Arc{.Center = {}, .Radius = 1e-170, .Start = {}, .Sweep = Degrees(90)};
    const Curve wave  = Bezier::Cubic{.Start = {.X = -1, .Y = 0},
                                      .Out   = {.X = 0, .Y = 1},
                                      .In    = {.X = 0, .Y = -1},
                                      .End   = {.X = 1, .Y = 0}};
    for (const Intersections &met :
         {Intersect(speck, wave, tolerance), Intersect(wave, speck, tolerance)})
    {
      REQUIRE(met.Hits.size() == 1);
      CHECK(met.Hits[0].Kind  == Contact::End);
      CHECK(met.Overlaps.empty());
    }
    const Curve through = Bezier::Quadratic{
        .Start = {.X = -1, .Y = 0}, .Control = {.X = 0, .Y = 0}, .End = {.X = 1, .Y = 0}};
    const Intersections none = Intersect(through, whole, tolerance);
    CHECK(none.Hits.empty());
    CHECK(none.Overlaps.empty());
  }

  /**
   *  One ellipse spelled by maps its parameter runs turned a quarter in, or the other way round:
   *  the same carrier, which only an equal map was found to share. Now any spelling of it shares,
   *  and so does one a hair off it.
   */
  TEST_CASE("intersect: an ellipse spelled turned or mirrored shares its stretch")
  {
    const Tolerance tolerance = {.Length = 2e-10, .Radians = 1e-10};
    const Curve     whole     = Ellipse{.Map   = {.X = {.X = 2, .Y = 0}, .Y = {.X = 0, .Y = 1}},
                                        .Start = Degrees(0),
                                        .Sweep = Degrees(360)};
    for (const double gap : {0.0, 1e-12})
    {
      CAPTURE(gap);
      const Curve turned   = Ellipse{.Map   = {.X      = {.X = 0, .Y = 1 + gap},
                                               .Y      = {.X = -2 - gap, .Y = 0},
                                               .Origin = {.X = gap, .Y = 0}},
                                     .Start = Degrees(0),
                                     .Sweep = Degrees(90)};
      const Curve mirrored = Ellipse{.Map   = {.X      = {.X = 2 + gap, .Y = 0},
                                               .Y      = {.X = 0, .Y = -1 - gap},
                                               .Origin = {.X = 0, .Y = gap}},
                                     .Start = Degrees(0),
                                     .Sweep = Degrees(-90)};
      const std::array<std::pair<const Curve *, Range>, 2> owed = {
          {{&turned, {.Start = 0.25, .End = 0.5}}, {&mirrored, {.Start = 0, .End = 0.25}}}};
      for (const auto &[other, stretch] : owed)
      {
        const Intersections met = Intersect(whole, *other, tolerance);
        CHECK(met.Hits.empty());
        REQUIRE(met.Overlaps.size()        == 1);
        CHECK(met.Overlaps[0].First.Start  == doctest::Approx(stretch.Start).epsilon(1e-12));
        CHECK(met.Overlaps[0].First.End    == doctest::Approx(stretch.End).epsilon(1e-12));
        CHECK(met.Overlaps[0].Second.Start == doctest::Approx(0).epsilon(1e-12));
        CHECK(met.Overlaps[0].Second.End   == doctest::Approx(1).epsilon(1e-12));
        CHECK(met.Overlaps[0].Along);
        CHECK(Intersect(*other, whole, tolerance).Overlaps.size() == 1);
      }
    }
  }

  /**
   *  Two arcs crossing within twice the tolerance of internally tangent, deeper than the direct way
   *  decides alone, so that the thorough way does: one touch, where they cross deepest, halfway
   *  between their points on the line of centers, to the bit, either first, plain or ready. A
   *  circle of radius 1/4 crossing the unit circle's inside 1.8 tolerances deep touches it halfway
   *  between (1, 0) and its own point on the x axis; the arcs of the offset case 21982, of radii
   *  0.011 and 0.0034, cross 1.8 tolerances deep at two crossings 2.7e-4 apart. Each touch is the
   *  double nearest the exact one. Each was put at one of its two crossings, which of them hanging
   *  on the order of the arcs.
   */
  TEST_CASE("intersect: arcs crossing within twice the tolerance of tangent touch where deepest")
  {
    /**
     *  Two arcs crossing near tangent, the tolerance and where they touch.
     */
    struct Lens
    {
      Curve     One    = {}; /**< One arc.                   */
      Curve     Other  = {}; /**< The other.                 */
      Tolerance Within = {}; /**< The tolerance.             */
      Point     Touch  = {}; /**< Where they touch, exactly. */
    };

    const std::array lenses = {
        Lens{.One    = Arc{.Center = {.X = 0, .Y = 0},
                           .Radius = 1,
                           .Start  = Degrees(-60),
                           .Sweep  = Degrees(120)},
             .Other  = Arc{.Center = {.X = 0x1.80003c65e1dc3p-1, .Y = 0},
                           .Radius = 0.25,
                           .Start  = Degrees(-90),
                           .Sweep  = Degrees(180)},
             .Within = {.Length = 1e-6, .Radians = 1e-9},
             .Touch  = {.X = 0x1.00000f1978771p+0, .Y = 0}},
        Lens{.One    = Arc{.Center = {.X = -0x1.26e978d4fdf3ap-9, .Y = 0x1.1a9fbe76c8b43p-7},
                           .Radius = 0x1.6904361b6543fp-7,
                           .Start  = Radians(-0x1.45a5ea77ad263p-4),
                           .Sweep  = Radians(-0x1.d2b018e001425p+0)},
             .Other  = Arc{.Center = {.X = -0x1.47ae147ae147cp-9, .Y = 0x1.0624dd2f1aap-10},
                           .Radius = 0x1.bc575e3cdcec2p-9,
                           .Start  = Radians(-0x1.58eea0b8fbe77p+1),
                           .Sweep  = Radians(0x1.78354c7488563p+1)},
             .Within = {.Length = 1e-6, .Radians = 1e-10},
             .Touch  = {.X = -0x1.563ca354ef5aap-9, .Y = -0x1.38e996631fecep-9}},
    };
    for (const Lens &lens : lenses)
    {
      const std::array ready = {Prepared{lens.One}, Prepared{lens.Other}};
      for (const bool flip : {false, true})
      {
        const Curves curves = flip ? Curves{lens.Other, lens.One} : Curves{lens.One, lens.Other};
        const Prepared::Pair pair = flip ? Prepared::Pair{std::cref(ready[1]), std::cref(ready[0])}
                                         : Prepared::Pair{std::cref(ready[0]), std::cref(ready[1])};
        Intersections        readied;
        Intersect(pair, lens.Within, readied);
        for (const Intersections &met : {Intersect(curves[0], curves[1], lens.Within), readied})
        {
          REQUIRE(met.Hits.size() == 1);
          CHECK(met.Hits[0].Kind  == Contact::Touch);
          CHECK(met.Hits[0].At    == lens.Touch);
        }
      }
    }
  }

  /**
   *  A touch of a line, a ray or an endless line near tangent to an arc, or of two arcs, the direct
   *  way's closed form, is the double nearest the true point halfway between where they come
   *  nearest, as where the same line as a polyline's segment touches the arc: either first, plain
   *  or ready, and mirrored. An endless line touching an arc of radius 1 (case 1530 of the
   *  comparison) was put 173 steps of a double off it in x, the line as a polyline's segment not; a
   *  circle of radius 78,000 and an arc of 3.4 million that come within 0.19 tolerances, 108 steps.
   *  Each was worked out in doubles from the centers or the foot, not rounded once.
   */
  TEST_CASE("intersect: a touch in closed form is the double nearest the true one, however spelled")
  {
    /**
     *  Two curves near tangent, the tolerance, and the double nearest where they touch.
     */
    struct Graze
    {
      Curve     One    = {}; /**< One curve.                 */
      Curve     Other  = {}; /**< The other, an arc.         */
      Tolerance Within = {}; /**< The tolerance.             */
      Point     Touch  = {}; /**< Where they touch, rounded. */
    };

    const Tolerance  within  = {.Length = 0x1.b7cdfd9d7bdbbp-33, .Radians = 0x1.b7cdfd9d7bdbbp-34};
    const Point      through = {.X = -0x1.4a70e75c75b76p-2, .Y = -0x1.1f8db092c2a50p-3};
    const Vector     way     = {.X = 0x1.eeada2f538635p-1, .Y = 0x1.081a0143a84c3p-2};
    const Point      touch   = {.X = -0x1.8e1a9c7401353p-9, .Y = -0x1.c32f607c42459p-5};
    const Arc        arc     = {.Center = {.X = 0x1.05e2554f0dfe3p-2, .Y = -0x1.0646547302276p+0},
                                .Radius = 0x1.00dd868f1285fp+0,
                                .Start  = Degrees(0x1.03c8e4f300c4ap+6),
                                .Sweep  = Degrees(80)};
    const Arc        mirror  = {.Center = {.X = -arc.Center.X, .Y = arc.Center.Y},
                                .Radius = arc.Radius,
                                .Start  = Degrees(180.0 - 0x1.03c8e4f300c4ap+6),
                                .Sweep  = Degrees(-80)};
    const std::array grazes  = {
        Graze{.One    = Xline{.Through = through, .Way = way},
              .Other  = arc,
              .Within = within,
              .Touch  = touch},
        Graze{.One    = Xline{.Through = {.X = -through.X, .Y = through.Y},
                              .Way     = {.X = -way.X, .Y = way.Y}},
              .Other  = mirror,
              .Within = within,
              .Touch  = {.X = -touch.X, .Y = touch.Y}},
        Graze{.One    = Ray{.Start = through, .Way = way},
              .Other  = arc,
              .Within = within,
              .Touch  = touch},
        Graze{.One    = Line{.Start = through, .End = through + way},
              .Other  = arc,
              .Within = within,
              .Touch  = touch},
        Graze{.One    = Polyline{.Vertices = {{.Position = through}, {.Position = through + way}}},
              .Other  = arc,
              .Within = within,
              .Touch  = touch},
        Graze{.One    = Circle({.X = 0x1.ea23fe94fa1bcp+15, .Y = 0x1.01983a3d4a078p+15},
                               0x1.31f57bc1924b9p+16),
              .Other  = Arc{.Center = {.X = -0x1.4c7c353fa784cp+21, .Y = -0x1.eaecabc085afdp+20},
                            .Radius = 0x1.9c451bc6b0422p+21,
                            .Start  = Radians(0x1.baa01cf6d2f63p+2),
                            .Sweep  = Radians(0x1.fd1008f6289a5p-17)},
              .Within = {.Length = 0x1.c54bad0ea4432p-19, .Radians = 1e-10},
              .Touch  = {.X = -0x1.a3c4ebb4e220fp+8, .Y = -0x1.a144b1f80eab6p+13}},
    };
    for (const Graze &graze : grazes)
    {
      const std::array ready = {Prepared{graze.One}, Prepared{graze.Other}};
      for (const bool flip : {false, true})
      {
        const Curves curves =
            flip ? Curves{graze.Other, graze.One} : Curves{graze.One, graze.Other};
        const Prepared::Pair pair = flip ? Prepared::Pair{std::cref(ready[1]), std::cref(ready[0])}
                                         : Prepared::Pair{std::cref(ready[0]), std::cref(ready[1])};
        Intersections        readied;
        Intersect(pair, graze.Within, readied);
        for (const Intersections &met : {Intersect(curves[0], curves[1], graze.Within), readied})
        {
          REQUIRE(met.Hits.size() == 1);
          CHECK(met.Hits[0].Kind  == Contact::Touch);
          CHECK(met.Hits[0].At    == graze.Touch);
        }
      }
    }
  }

  /**
   *  An arc lying within the tolerance of an arc of another circle all along, both its ends inside
   *  the band where the two stay within the tolerance, shares a stretch with it: the stretch's ends
   *  are where a curve ends, not where the curves part. Two short arcs of radii 466 and 65, the
   *  second 3e-6 long and within 1e-14 of the first all along, against a tolerance of 4.7e-10, were
   *  one end hit, at the second's end one way and at its start the other, as two circles that do
   *  not lie within the tolerance all the way round shared nothing. The same as bulges too.
   */
  TEST_CASE("intersect: an arc within the tolerance of another circle's arc all along shares it")
  {
    const Arc        wide   = {.Center = {.X = 0x1.f37a5036058p+4, .Y = 0x1.8430561148c2ep+8},
                               .Radius = 0x1.d2945d42078p+8,
                               .Start  = Radians(0x1.cbcd767aa978ep+1),
                               .Sweep  = Radians(0x1.acd982f946feap-25)};
    const Arc        narrow = {.Center = {.X = -0x1.bf4b3950ad668p+8, .Y = 0x1.395b5bf1e2934p+7},
                               .Radius = 0x1.03f97655c69d7p+6,
                               .Start  = Radians(0x1.aef695e22025dp+2),
                               .Sweep  = Radians(0x1.93d76e28019ddp-25)};
    const Tolerance  within = {.Length = 0x1.00813bd405008p-31, .Radians = 1e-10};
    const Curve      one    = wide;
    const Curve      other  = narrow;
    const std::array ready  = {Prepared{one}, Prepared{other}};
    for (const bool flip : {false, true})
    {
      const Curves         curves = flip ? Curves{other, one} : Curves{one, other};
      const Prepared::Pair pair   = flip ? Prepared::Pair{std::cref(ready[1]), std::cref(ready[0])}
                                         : Prepared::Pair{std::cref(ready[0]), std::cref(ready[1])};
      Intersections        readied;
      Intersect(pair, within, readied);
      for (const Intersections &met : {Intersect(curves[0], curves[1], within), readied})
      {
        CHECK(met.Hits.empty());
        REQUIRE(met.Overlaps.size() == 1);
        const Overlap &shared  = met.Overlaps[0];
        const Range   &longer  = flip ? shared.Second : shared.First;
        const Range   &shorter = flip ? shared.First : shared.Second;
        CHECK(longer.Start  == 0.29125829934685415);
        CHECK(longer.End    == 0.42243325279015431);
        CHECK(shorter.Start == 0.0);
        CHECK(shorter.End   == 1.0);
        CHECK_FALSE(shared.Along);
      }
    }
    const auto bulged = [](const Arc &arc) -> Curve
    {
      return Polyline{.Vertices = {{.Position = At(Curve{arc}, 0.0),
                                    .Bulge    = std::tan(arc.Sweep.InRadians() / 4.0)},
                                   {.Position = At(Curve{arc}, 1.0)}}};
    };
    for (const bool flip : {false, true})
    {
      const Curve         first  = bulged(flip ? narrow : wide);
      const Curve         second = bulged(flip ? wide : narrow);
      const Intersections met    = Intersect(first, second, within);
      CHECK(met.Hits.empty());
      CHECK(met.Overlaps.size() == 1);
    }
  }

  /**
   *  A line, a ray or an endless line lying within the tolerance of an arc of a circle shares a
   *  stretch with it, between ends of either within the tolerance of the other, as two arcs do. A
   *  line 2^-14 long 2^-31 above the bottom of the unit circle, the arc running 0.2 either side,
   *  was one end hit at the line's end; the arc's foot inside the stretch or outside it, the arc
   *  ending past the line's ends or within them, the same either first, plain, ready and as a
   *  polyline's segment. An arc swept round more than once is left to its passes.
   */
  TEST_CASE("intersect: a line within the tolerance of an arc all along shares it")
  {
    const double    low    = -1.0 + 0x1p-31;
    const double    half   = 0x1p-15;
    const Tolerance within = {.Length = 1e-9, .Radians = 1e-9};
    const auto      arc    = [](double from, double sweep) -> Arc
    {
      return {.Center = {},
              .Radius = 1,
              .Start  = Radians((-std::numbers::pi / 2) + from),
              .Sweep  = Radians(sweep)};
    };
    const auto line = [low](double from, double until) -> Curve
    {
      return Line{.Start = {.X = from, .Y = low}, .End = {.X = until, .Y = low}};
    };
    const Curve endless  = Xline{.Through = {.X = 0, .Y = low}, .Way = {.X = 1, .Y = 0}};
    const Curve polyline = Polyline{
        .Vertices = {{.Position = {.X = -half, .Y = low}}, {.Position = {.X = half, .Y = low}}}};
    /**
     *  The curve along, the arc, and each one's stretch.
     */
    const std::vector<std::tuple<Curve, Curve, Range, Range>> cases = {
        {line(-half, half),
         arc(-0.2, 0.4),
         {.Start = 0, .End = 1},
         {.Start = 0.49992370605467518, .End = 0.50007629394532405}},
        {polyline,
         arc(-0.2, 0.4),
         {.Start = 0, .End = 1},
         {.Start = 0.49992370605467518, .End = 0.50007629394532405}},
        {line(half / 4, half),
         arc(-0.2, 0.4),
         {.Start = 0, .End = 1},
         {.Start = 0.50001907348633634, .End = 0.50007629394532405}},
        {line(-half / 2, half / 2),
         arc(-half, 2 * half),
         {.Start = 0, .End = 1},
         {.Start = 0.24999999990177457, .End = 0.75000000009458745}},
        {Ray{.Start = {.X = -half / 2, .Y = low}, .Way = {.X = 1, .Y = 0}},
         arc(-0.2, 0.2 + half),
         {.Start = 0, .End = 4.5776367182824283e-05},
         {.Start = 0.99977115308330056, .End = 1}},
        {endless,
         arc(-half, 2 * half),
         {.Start = -3.0517578120201816e-05, .End = 3.051757812032428e-05},
         {.Start = 0, .End = 1}},
        {endless,
         arc(half / 4, half),
         {.Start = 7.6293945312372174e-06, .End = 3.8146972647059374e-05},
         {.Start = 0, .End = 1}},
        {endless,
         arc(-half / 4, -half),
         {.Start = -3.8146972646936907e-05, .End = -7.6293945311147526e-06},
         {.Start = 0, .End = 1}},
    };
    for (const auto &[along, round, onAlong, onRound] : cases)
    {
      const std::array ready = {Prepared{along}, Prepared{round}};
      for (const bool flip : {false, true})
      {
        const Prepared::Pair pair = flip ? Prepared::Pair{std::cref(ready[1]), std::cref(ready[0])}
                                         : Prepared::Pair{std::cref(ready[0]), std::cref(ready[1])};
        Intersections        readied;
        Intersect(pair, within, readied);
        for (const Intersections &met :
             {flip ? Intersect(round, along, within) : Intersect(along, round, within), readied})
        {
          CHECK(met.Hits.empty());
          REQUIRE(met.Overlaps.size() == 1);
          const Overlap &shared = met.Overlaps[0];
          const Range   &first  = flip ? shared.Second : shared.First;
          const Range   &second = flip ? shared.First : shared.Second;
          CHECK(first.Start  == onAlong.Start);
          CHECK(first.End    == onAlong.End);
          CHECK(second.Start == onRound.Start);
          CHECK(second.End   == onRound.End);
          CHECK(shared.Along == (std::get<Arc>(round).Sweep.InRadians() > 0.0));
        }
      }
    }
    const Curve coiled = arc(-0.2, (4 * std::numbers::pi) + 0.4);
    CHECK(Intersect(line(-half, half), coiled, within).Overlaps.empty());
    CHECK(Intersect(coiled, line(-half, half), within).Overlaps.empty());
  }

  /**
   *  An end hit stands at its end, the same at any scale and whichever curve is first. Two arcs of
   *  radii 6.3e-11 and 1.1e-11 crossing 2.3 tolerances deep, the second's start within the
   *  tolerance of the first: their circles' other crossing lies past that start, off the second arc
   *  by 3,000 tolerances, and was taken for on it, held to a few steps of a double at 1, not at the
   *  arcs' own size; the end hit went there, and scaled by 2^11 stayed at the start. And an arc's
   *  start 0.8 tolerances from another, in the band about where they come nearest 0.72 apart: the
   *  run's hit was a candidate polished from a far end to 0.0125 along the arc, as near the other
   *  as the start is, and so neither at the end nor where they come nearest.
   */
  TEST_CASE("intersect: an end hit stands at the end, at any scale, either first")
  {
    for (const int power : {0, 11})
    {
      const auto scaled = [power](double value) -> double
      {
        return std::ldexp(value, power);
      };
      const Curve one =
          Arc{.Center = {.X = scaled(0x1.0c820340c07f8p-34), .Y = scaled(-0x1.0de33831babfp-38)},
              .Radius = scaled(0x1.14918fa5cca57p-34),
              .Start  = Radians(0x1.75773d5d8b96p+2),
              .Sweep  = Radians(0x1.7a11cf1eb274fp-15)};
      const Curve other =
          Arc{.Center = {.X = scaled(0x1.181f05e35e482p-33), .Y = scaled(-0x1.39f8af755c78dp-35)},
              .Radius = scaled(0x1.788f8d07b2faep-37),
              .Start  = Radians(0x1.1f43c9b74fc92p+3),
              .Sweep  = Radians(-0x1.8b3c2e0f194abp-17)};
      const Tolerance within = {.Length = scaled(0x1.301723c4b30a2p-74), .Radians = 1e-10};
      for (const bool flip : {false, true})
      {
        const Intersections met =
            flip ? Intersect(other, one, within) : Intersect(one, other, within);
        REQUIRE(met.Hits.size() == 2);
        const Hit &end = met.Hits[0];
        CHECK(end.Kind == Contact::End);
        CHECK(end.At
              == Point{.X = scaled(0x1.02e7fad033e43p-33), .Y = scaled(-0x1.1136d9ee4d7b6p-35)});
        CHECK((flip ? end.First : end.Second).Along == 0.0);
        CHECK(met.Hits[1].Kind == Contact::Cross);
      }
    }
    const Curve     one   = Arc{.Center = {.X = 0x1.0c20a0f066ed8p+18, .Y = -0x1.0bcb9becc47d2p+19},
                                .Radius = 0x1.3f592564cc60bp+20,
                                .Start  = Radians(0x1.8bbcacb388b4ep+1),
                                .Sweep  = Radians(0x1.cb2e046f8bb7ep-6)};
    const Curve     other = Arc{.Center = {.X = -0x1.93c019902d8bp+18, .Y = -0x1.fc22ecaf85374p+18},
                                .Radius = 0x1.2e3313b158481p+19,
                                .Start  = Radians(0x1.a3115939fa79dp+1),
                                .Sweep  = Radians(-0x1.39c7362ae1579p-2)};
    const Tolerance within = {.Length = 0x1.8b052e04714cfp+9, .Radians = 1e-10};
    for (const bool flip : {false, true})
    {
      const Intersections met =
          flip ? Intersect(other, one, within) : Intersect(one, other, within);
      REQUIRE(met.Hits.size() == 1);
      CHECK(met.Hits[0].Kind  == Contact::End);
      CHECK(met.Hits[0].At    == Point{.X = -0x1.f7af4703ca93p+19, .Y = -0x1.d7e5a99c73f38p+18});
      CHECK((flip ? met.Hits[0].Second : met.Hits[0].First).Along == 0.0);
    }
  }

  /**
   *  An end within the tolerance of a curve that crosses it at a sine of 2e-12 to 2e-9 far along is
   *  one end hit with that crossing, the curves staying within twice the tolerance between: the hit
   *  stands at the end. Comparison cases 729, 771 and 4737, a tolerance of 2e-10: two rays, a ray
   *  and an endless line, two polyline segments; the end hit stood at the crossing, 1e8 to 1e9
   *  tolerances from every end.
   */
  TEST_CASE(
      "intersect: an end in the band of a far shallow crossing is met at the end, either first")
  {
    const Tolerance within  = {.Length = 0x1.b7cdfd9d7bdbbp-33, .Radians = 0x1.b7cdfd9d7bdbbp-34};
    const auto      segment = [](Point start, Point end) -> Curve
    {
      return Polyline{.Vertices = {{.Position = start}, {.Position = end}}};
    };
    const std::vector<std::tuple<Curve, Curve, Point>> cases = {
        {Ray{.Start = {.X = -0x1.fe75ce3a04c08p-1, .Y = 0x1.ece5f052baae0p-1},
             .Way   = {.X = 0x1.4465499f9393bp-1, .Y = -0x1.62516462c6802p+0}},
         Ray{.Start = {.X = -0x1.cef8ac9edd8a0p-1, .Y = 0x1.8528ed30c3262p-1},
             .Way   = {.X = 0x1.aa398130f5aeap-2, .Y = -0x1.d18a24d010709p-1}},
         {.X = -0x1.cef8ac9edd8a0p-1, .Y = 0x1.8528ed30c3262p-1}},
        {Ray{.Start = {.X = -0x1.b29f763734b40p-2, .Y = -0x1.86ae412a82fc3p-1},
             .Way   = {.X = 0x1.9e2c63cb24dbcp-2, .Y = 0x1.65e968edeaecfp-1}},
         Xline{.Through = {.X = -0x1.53f8f0b38be58p-2, .Y = -0x1.34e339e3fa754p-1},
               .Way     = {.X = 0x1.006a0f5dc4a80p-1, .Y = 0x1.bb2a61c04b090p-1}},
         {.X = -0x1.b29f763734b40p-2, .Y = -0x1.86ae412a82fc3p-1}},
        {segment({.X = -0x1.039f63ca3a0e0p-1, .Y = 0x1.68fe507f26b68p-1},
                 {.X = -0x1.1734f32f8fdc4p-2, .Y = 0x1.488a156bcbc88p-3}),
         segment({.X = -0x1.039f63c8118a7p-1, .Y = 0x1.68fe507935159p-1},
                 {.X = -0x1.1734f320287dap-2, .Y = 0x1.488a158d68895p-3}),
         {.X = -0x1.039f63c8118a7p-1, .Y = 0x1.68fe507935159p-1}},
    };
    for (const auto &[one, other, end] : cases)
    {
      for (const bool flip : {false, true})
      {
        const Curve        &first  = flip ? other : one;
        const Curve        &second = flip ? one : other;
        const Intersections met    = Intersect(first, second, within);
        REQUIRE(met.Hits.size() == 1);
        CHECK(met.Hits[0].Kind  == Contact::End);
        CHECK(Distance(met.Hits[0].At, end) < within.Length);
        const Prepared ready(first);
        const Prepared readied(second);
        Marking        marking;
        Intersect(Prepared::Pair{ready, readied}, within, marking);
        REQUIRE(marking.Marks.size() == 1);
        CHECK(marking.Marks[0]       == Mark{.At = met.Hits[0].At, .Kind = Contact::End});
      }
    }
  }

  namespace
  {
    /**
     *  Two curves met plain and ready, the first first and then the second.
     *
     *  @param  one    The one curve.
     *  @param  other  The other.
     *  @param  within The tolerance.
     *  @return        Plain and ready with the one first, then plain and ready with the other.
     */
    auto Fourfold(const Curve &one, const Curve &other, const Tolerance &within)
        -> std::array<Intersections, 4>
    {
      const std::array             ready = {Prepared{one}, Prepared{other}};
      std::array<Intersections, 4> met   = {
          Intersect(one, other, within), {}, Intersect(other, one, within), {}};
      Intersect(Prepared::Pair{std::cref(ready[0]), std::cref(ready[1])}, within, met[1]);
      Intersect(Prepared::Pair{std::cref(ready[1]), std::cref(ready[0])}, within, met[3]);
      return met;
    }

    /**
     *  How two curves meet at each hit.
     *
     *  @param  meeting What meeting them found.
     *  @return         The hits' kinds, in order along the first.
     */
    auto Kinds(const Intersections &meeting) -> std::vector<Contact>
    {
      std::vector<Contact> kinds;
      kinds.reserve(meeting.Hits.size());
      for (const Hit &hit : meeting.Hits)
      {
        kinds.push_back(hit.Kind);
      }
      return kinds;
    }
  } // namespace

  /**
   *  Two hits the curves part more than twice the tolerance between stay two. The hunt's RA-85389:
   *  a ray starting 0.35 tolerances off an arc runs along it nearly as a chord and crosses it 0.7
   *  along, 480 tolerances, at a sine of 0.016; between, the two part by 2.004 tolerances (MPFR) at
   *  0.48 of the way, where the samples at a quarter, a half and three quarters fell short of twice
   *  the tolerance, and the crossing was merged into the end hit. And BB-5328: two segments of
   *  bulge -14.1, nearly whole circles, the first's ends each within the tolerance of the second,
   *  which part by 386,590 tolerances between them; the samples found the second near, and one end
   *  hit was missing met the first way round.
   */
  TEST_CASE(
      "intersect: hits the curves part more than twice the tolerance between are two, either first")
  {
    const Curve arc = Arc{.Center = {.X = 0x1.59062c9e9f55dp+5, .Y = 0x1.b8e369ddbd741p+4},
                          .Radius = 0x1.79883f4ac523fp+4,
                          .Start  = Radians(0x1.16114bff0aep-2),
                          .Sweep  = Radians(0x1.c10f818ab3a2cp+1)};
    const Curve ray = Ray{.Start = {.X = 0x1.3b0b96dcaa90ap+4, .Y = 0x1.8d550d5074734p+4},
                          .Way   = {.X = 0x1.0967a3111297ap-3, .Y = -0x1.fbaeb8164f38bp-1}};
    const std::vector<Point> hits = {{.X = 0x1.3b0a8aeaa95ccp+4, .Y = 0x1.8d54ee315a0bap+4},
                                     {.X = 0x1.3c849af086716p+4, .Y = 0x1.8210560865151p+4}};
    for (const Intersections &met :
         Fourfold(arc, ray, {.Length = 0x1.8016144bfad99p-10, .Radians = 1e-9}))
    {
      CHECK(Points(met) == hits);
      CHECK(Kinds(met)  == std::vector{Contact::End, Contact::Cross});
    }
    const Curve one = Polyline{
        .Vertices = {{.Position = {.X = 0x1.b1549b33a4d7cp+4, .Y = -0x1.c32ccf048af04p+4},
                      .Bulge    = -0x1.c33ed50b88778p+3},
                     {.Position = {.X = 0x1.9b9c055db31bcp+4, .Y = -0x1.96f889d9076a5p+4}}},
        .Closed   = false};
    const Curve other = Polyline{
        .Vertices = {{.Position = {.X = 0x1.539fd0dfaf591p+4, .Y = -0x1.4c1079ca4cd59p+5},
                      .Bulge    = -0x1.c33ed50b88778p+3},
                     {.Position = {.X = 0x1.7d235aa4e421ap+4, .Y = -0x1.3ed01410e89b6p+5}}},
        .Closed   = false};
    const std::vector<Point> ends = {{.X = 0x1.b154997d11108p+4, .Y = -0x1.c32ccf9384ad4p+4},
                                     {.X = 0x1.9b9c0422ea31bp+4, .Y = -0x1.96f88aaf9c023p+4}};
    const std::array<Intersections, 4> met =
        Fourfold(one, other, {.Length = 0x1.0ab0422e98e14p-18, .Radians = 1e-9});
    for (std::size_t index = 0; index < met.size(); ++index)
    {
      CHECK(Points(met.at(index)) == (index < 2 ? ends : std::vector{ends[1], ends[0]}));
      CHECK(Kinds(met.at(index))  == std::vector{Contact::End, Contact::End});
    }
  }

  /**
   *  An end hit where the curves cross within the tolerance of a bulge's end stands at the double
   *  nearest that crossing: no exact way places a crossing on a bulge, and the hit stayed where its
   *  candidates found it, 1 to 10,642 steps of a double off, in 10,503 of the hunt's million pairs.
   *  The hunt's BB-1, two segments crossing at a sine of 0.3, 0.97 tolerances from the second's
   *  start, 2 steps off; and AB-6, a segment and an arc crossing at a sine of 0.36, 0.06 tolerances
   *  from the arc's end, 5 steps off.
   */
  TEST_CASE("intersect: an end hit at a crossing near a bulge's end is the double nearest it, "
            "either first")
  {
    const Curve bow = Polyline{
        .Vertices = {{.Position = {.X = 0x1.78e787c5c4fb2p+16, .Y = -0x1.30337f79800acp+16},
                      .Bulge    = -0x1.ab41f88729db3p-7},
                     {.Position = {.X = 0x1.62620e04d339p+16, .Y = -0x1.07ce4577c398ep+16}}},
        .Closed   = false};
    const Curve bent = Polyline{
        .Vertices = {{.Position = {.X = 0x1.6c55287c7b419p+16, .Y = -0x1.1a44041922098p+16},
                      .Bulge    = -0x1.56e9a90fbe989p-1},
                     {.Position = {.X = 0x1.7c9783ae942b1p+16, .Y = -0x1.0f378907b3f42p+16}}},
        .Closed   = false};
    const Curve segment = Polyline{
        .Vertices = {{.Position = {.X = -0x1.501565884d9p+17, .Y = -0x1.c6617425cb0ap+18},
                      .Bulge    = -0x1.5cb13735c2d89p-6},
                     {.Position = {.X = 0x1.f3c9838322864p+19, .Y = -0x1.3b6d9824526p+16}}},
        .Closed   = false};
    const Curve arc = Arc{.Center = {.X = 0x1.178b9e435141fp+19, .Y = 0x1.af13f782541e8p+20},
                          .Radius = 0x1.eddc70d911736p+20,
                          .Start  = Radians(0x1.5286218054518p+3),
                          .Sweep  = Radians(0x1.69acec441835cp-2)};
    for (const Intersections &met :
         Fourfold(bow, bent, {.Length = 0x1.3142473cfae2dp-7, .Radians = 1e-9}))
    {
      CHECK(Points(met)
            == std::vector<Point>{{.X = 0x1.6c5528039fd6bp+16, .Y = -0x1.1a4401d7bb4cdp+16}});
      CHECK(Kinds(met) == std::vector{Contact::End});
    }
    for (const Intersections &met :
         Fourfold(segment, arc, {.Length = 0x1.01aec163793d1p-4, .Radians = 1e-9}))
    {
      CHECK(Points(met)
            == std::vector<Point>{{.X = 0x1.b1bae35ff87fdp+18, .Y = -0x1.ee4d2c2f5b3ebp+17}});
      CHECK(Kinds(met) == std::vector{Contact::End});
    }
  }

  TEST_CASE("intersect: two curves a ten million millionth from parallel cross where they cross")
  {
    /**
     *  A level quadratic, and a line and a cubic that rise across it by a slant a side: each
     *  crosses it once, and at 1e-13 was found not to cross at all, taken for parallel under 1e-12
     *  radians. The heights are doubles, and a step of a double in height moves the crossing some
     *  ten thousandths along at such a slant: the hit is where the curves as their doubles tell
     *  them cross, worked out exactly, the line's to the bit and the cubic's within a step of a
     *  double there; it was left where it was first found, as far as 1e-4 from the crossing.
     */
    struct Case
    {
      double Slant = 0.0; /**< How far the curves rise a side.          */
      double Line  = 0.0; /**< Where the line crosses, along the level. */
      double Cubic = 0.0; /**< Where the cubic does.                    */
    };

    const Bezier::Quadratic level = {.Start   = {.X = -1, .Y = 0.25},
                                     .Control = {.X = 0.1, .Y = 0.25},
                                     .End     = {.X = 1, .Y = 0.25}};
    const Tolerance         fine  = {.Length = 1e-14, .Radians = 1e-12};
    for (const Case &row :
         {Case{.Slant = 1e-12, .Line = 1.3877902216300983e-05, .Cubic = 2.3282317879477894e-06},
          Case{.Slant = 1e-13, .Line = 0.00013879250520471896, .Cubic = -3.8311961381332961e-05}})
    {
      const double slant = row.Slant;
      const Line line = {.Start = {.X = -1, .Y = 0.25 - slant}, .End = {.X = 1, .Y = 0.25 + slant}};
      const Bezier::Cubic rising = {.Start = {.X = -1, .Y = 0.25 - slant},
                                    .Out   = {.X = -0.3, .Y = 0.25 - (0.3 * slant)},
                                    .In    = {.X = 0.4, .Y = 0.25 + (0.4 * slant)},
                                    .End   = {.X = 1, .Y = 0.25 + slant}};
      for (const Intersections &met :
           {Intersect(Curve{level}, Curve{line}, fine), Intersect(Curve{line}, Curve{level}, fine)})
      {
        REQUIRE(met.Hits.size() == 1);
        CHECK(met.Overlaps.empty());
        CHECK(met.Hits[0].At.Y == 0.25);
        CHECK(met.Hits[0].At.X == row.Line);
      }
      for (const Intersections &met : {Intersect(Curve{level}, Curve{rising}, fine),
                                       Intersect(Curve{rising}, Curve{level}, fine)})
      {
        REQUIRE(met.Hits.size() == 1);
        CHECK(met.Overlaps.empty());
        CHECK(met.Hits[0].At.Y == 0.25);
        CHECK(std::abs(met.Hits[0].At.X - row.Cubic) <= 0x1p-52 * std::abs(row.Cubic));
      }
    }
  }

  TEST_CASE("intersect: a polyline's slight arc is met as the arc it is, however slight")
  {
    /**
     *  A bulge under 1e-4 puts its circle's center thousands of chords off, and one of 1e-17 a
     *  hundred thousand million million: every way of meeting worked the arc out through that
     *  center, and the polyline was sized by it, which made the tolerance 128 where 1e-10 was
     *  asked. Far from an arc and from a line, the arc met both at an end.
     */
    const auto bowed = [](double chord, double bulge) -> Curve
    {
      return Polyline{.Vertices = {{.Position = {.X = 0, .Y = 0}, .Bulge = bulge},
                                   {.Position = {.X = chord, .Y = 0}}}};
    };
    const Tolerance usual = {.Length = 1e-10, .Radians = 1e-10};
    const Curve     round =
        Arc{.Center = {.X = 5, .Y = 5}, .Radius = 1, .Start = Degrees(0), .Sweep = Degrees(90)};
    const Curve far = Line{.Start = {.X = 3, .Y = 3}, .End = {.X = 4, .Y = 5}};
    for (const double bulge : {0x1p-60, 1e-17, 1e-12, 1e-6})
    {
      for (const Curve &other : {round, far})
      {
        for (const Intersections &met :
             {Intersect(bowed(1, bulge), other, usual), Intersect(other, bowed(1, bulge), usual)})
        {
          CHECK(met.Hits.empty());
          CHECK(met.Overlaps.empty());
        }
      }
    }

    /**
     *  A line across the arc's middle crosses it on the arc, half the bulge below the chord, which
     *  the chord stood for while the crossing was found: (1 / 2, -5e-18).
     */
    const Curve         upright = Line{.Start = {.X = 0.5, .Y = -1}, .End = {.X = 0.5, .Y = 1}};
    const Intersections across  = Intersect(bowed(1, 1e-17), upright, usual);
    REQUIRE(across.Hits.size()       == 1);
    CHECK(across.Hits[0].Kind        == Contact::Cross);
    CHECK(across.Hits[0].At          == Point{.X = 0.5, .Y = -5.0000000000000004e-18});
    CHECK(across.Hits[0].First.At    == across.Hits[0].At);
    CHECK(across.Hits[0].First.Along == 0.5);

    /**
     *  A level line halfway between the chord and the arc's middle crosses the arc twice, at 1 / 2
     *  -+ 1 / sqrt(8) of the chord, the tolerance a thousandth of how far the arc stands off: it
     *  was told to share the arc's whole length. The places are the doubles nearest the true ones,
     *  worked out at 50 digits; a polyline made ready meets the same.
     */
    struct Case
    {
      double Bulge = 0.0; /**< The arc's bulge.                 */
      double Left  = 0.0; /**< Where the line first crosses it. */
      double Right = 0.0; /**< Where it crosses it again.       */
    };

    for (const Case &row :
         {Case{.Bulge = 1e-8, .Left = 146.44660940672625, .Right = 853.5533905932738},
          Case{.Bulge = 1e-11, .Left = 146.44660940672625, .Right = 853.5533905932738},
          Case{.Bulge = 1e-13, .Left = 146.44660940672622, .Right = 853.5533905932738}})
    {
      const double sagitta = row.Bulge * 1000 / 2;
      const Curve  bow     = bowed(1000, row.Bulge);
      const Curve  level =
          Line{.Start = {.X = -1000, .Y = -sagitta / 2}, .End = {.X = 2000, .Y = -sagitta / 2}};
      const Tolerance  fine  = {.Length = sagitta / 1e3, .Radians = 1e-14};
      const std::array ready = {Prepared{bow}, Prepared{level}};
      Intersections    readied;
      Intersect(Prepared::Pair{std::cref(ready[0]), std::cref(ready[1])}, fine, readied);
      Marking marked;
      Intersect(Prepared::Pair{std::cref(ready[0]), std::cref(ready[1])}, fine, marked);
      CHECK(marked.Marks.size() == 2);
      for (const Intersections &met :
           {Intersect(bow, level, fine), Intersect(level, bow, fine), readied})
      {
        REQUIRE(met.Hits.size() == 2);
        CHECK(met.Overlaps.empty());
        const bool flipped = met.Hits[0].At.X > met.Hits[1].At.X;
        const Hit &left    = met.Hits[flipped ? 1 : 0];
        const Hit &right   = met.Hits[flipped ? 0 : 1];
        CHECK(left.Kind  == Contact::Cross);
        CHECK(right.Kind == Contact::Cross);
        CHECK(left.At    == Point{.X = row.Left, .Y = -sagitta / 2});
        CHECK(right.At   == Point{.X = row.Right, .Y = -sagitta / 2});
        CHECK(Distance(left.First.At, left.Second.At) <= 0x1p-42);
        CHECK(Distance(right.First.At, right.Second.At) <= 0x1p-42);
      }
    }
  }

  TEST_CASE("intersect: a line that touches a slight arc touches it once, at one place")
  {
    /**
     *  A line nearly along a chord of a million, the arc standing 1.3e-8 off it and the tolerance
     *  1.1e-8: the touch was polished through the arc's center, 9e18 away, and came back as three
     *  touches at one point with three places on the arc, each hundreds from it.
     */
    const Curve bow = Polyline{
        .Vertices = {{.Position = {.X = 0x1.13b56ecc6ae3ap+18, .Y = -0x1.8414ca608908cp+19},
                      .Bulge    = 0x1.dfbdba2065d58p-46},
                     {.Position = {.X = 0x1.52209b53cf40cp+19, .Y = 0x1.0451d4d8b1abp+16}}}};
    const Curve     line = Line{.Start = {.X = 0x1.3716ff4883b67p+17, .Y = -0x1.01223f8936027p+20},
                                .End   = {.X = 0x1.8e3592e7e3cf6p+19, .Y = 0x1.3d73de99f2594p+18}};
    const Tolerance within = {.Length = 0x1.6aeef7e04421ep-27, .Radians = 0x1.b7cdfd9d7bdbbp-34};
    for (const Intersections &met : {Intersect(bow, line, within), Intersect(line, bow, within)})
    {
      REQUIRE(met.Hits.size() == 2);
      CHECK(met.Overlaps.empty());
      CHECK(std::ranges::count(met.Hits, Contact::Touch, &Hit::Kind) == 1);
      CHECK(std::ranges::count(met.Hits, Contact::End, &Hit::Kind)   == 1);
      for (const Hit &hit : met.Hits)
      {
        CHECK(Distance(hit.First.At, hit.At) <= within.Length);
        CHECK(Distance(hit.Second.At, hit.At) <= within.Length);
      }
    }
  }

  TEST_CASE("intersect: a crossing just past an arc's start is placed on the arc there")
  {
    /**
     *  A line crosses a polyline 5e-6 of the way into its second segment, bulged 0.6: the place was
     *  told from the tangent at the start, as though the arc ran straight on, and the arc's point
     *  there lay 1e-8 from the crossing. It is told from the angle at the arc's far end.
     */
    const Curve strung = Polyline{
        .Vertices = {{.Position = {.X = -0x1.5d622031a7c73p+12, .Y = 0x1.3317a0dbc2283p+10},
                      .Bulge    = -0x1.cf33f03876c7ep-18},
                     {.Position = {.X = -0x1.73deb002f87c8p+12, .Y = -0x1.7f093c766205cp+5},
                      .Bulge    = -0x1.3333333333333p-1},
                     {.Position = {.X = -0x1.75401e58a04e5p+12, .Y = -0x1.6b96834af865p+9}}}};
    const Curve     line = Line{.Start = {.X = -0x1.56a32987e580cp+12, .Y = 0x1.92cf9bedb9ab2p+10},
                                .End   = {.X = -0x1.7a9d95addad42p+12, .Y = -0x1.aec1607f6dd07p+8}};
    const Tolerance within = {.Length = 0x1.6e969ab7c47fbp-16, .Radians = 1e-10};
    for (const Intersections &met :
         {Intersect(strung, line, within), Intersect(line, strung, within)})
    {
      REQUIRE(met.Hits.size() == 3);
      for (const Hit &hit : met.Hits)
      {
        CHECK(hit.Kind == Contact::Cross);
        CHECK(Distance(hit.First.At, hit.Second.At) <= 0x1p-38);
      }
    }
  }

  TEST_CASE(
      "intersect: a crossing past a joint is not lost in a long shallow approach, either first")
  {
    /**
     *  A line runs within about twice the tolerance of a polyline's first segment, a slight arc a
     *  million long, for all its length, crossing it at a sine of 1.5e-8, and then crosses the next
     *  segment twice: at the joint, and 160,000 farther on at a sine of 0.19. The approach, its end
     *  and the crossing at the joint are one hit, an end. The far crossing was taken into it too,
     *  three places between the two hits being looked at and all three on the first segment: the
     *  line first, it was lost; the polyline first, the crossing at the joint was.
     */
    const Curve line   = Line{.Start = {.X = -0x1.63a6511f945f1p+19, .Y = 0x1.152cc27d5d8d2p+20},
                              .End   = {.X = 0x1.83cfa60ff8830p+14, .Y = -0x1.b3f0b443e9c09p+19}};
    const Curve strung = Polyline{
        .Vertices = {{.Position = {.X = -0x1.1eb16a82d9ce4p+19, .Y = 0x1.70ab9a3831b58p+19},
                      .Bulge    = -0x1.15f56bc787b78p-27},
                     {.Position = {.X = -0x1.c6b34be45650cp+16, .Y = -0x1.f485930ede6e8p+18},
                      .Bulge    = -0x1.3333333333333p-1},
                     {.Position = {.X = 0x1.fd0a172de05acp+17, .Y = 0x1.f4ab9d08c2a8ap+16},
                      .Bulge    = -0x1.87c49b0243b4cp-47}}};
    const Tolerance within = {.Length = 0x1.45c74d464b798p-9, .Radians = 1e-10};
    for (const Intersections &met :
         {Intersect(line, strung, within), Intersect(strung, line, within)})
    {
      REQUIRE(met.Hits.size() == 2);
      CHECK(met.Overlaps.empty());
      CHECK(std::ranges::count(met.Hits, Contact::End, &Hit::Kind) == 1);
      const auto crossing = std::ranges::find(met.Hits, Contact::Cross, &Hit::Kind);
      REQUIRE(crossing != met.Hits.end());
      CHECK(crossing->At == Point{.X = -172766.20692255773, .Y = -360766.04877798544});
      CHECK(Distance(crossing->First.At, crossing->Second.At) <= 1e-9);
    }
  }

  TEST_CASE("intersect: a line along a segment into a joint shares it and no more, either first")
  {
    /**
     *  A line runs along a polyline's first segment and a step past its end, where the polyline
     *  turns back along a slight arc whose circle the line crosses 0.8 of a tolerance before the
     *  arc starts. The line first, that crossing was told beside the shared stretch, its place
     *  within the stretch along the line and a step past it along the polyline; the polyline first,
     *  it was the stretch's.
     */
    const Curve line   = Line{.Start = {.X = -0x1.ab01ffdb91df3p+8, .Y = -0x1.0905670d8f93cp+13},
                              .End   = {.X = 0x1.a8b7c553ad3a2p+9, .Y = -0x1.076ebb1891b31p+13}};
    const Curve strung = Polyline{
        .Vertices = {{.Position = {.X = -0x1.77596bc61d378p+7, .Y = -0x1.08b926cf98e02p+13},
                      .Bulge    = 0x1.deed7c123e0cep-33},
                     {.Position = {.X = 0x1.310d205762570p+9, .Y = -0x1.07bafb5679e04p+13},
                      .Bulge    = -0x1.17e053ca8aafep-14},
                     {.Position = {.X = -0x1.95c0d5c42229cp+9, .Y = -0x1.068b7d676b85bp+13},
                      .Bulge    = 0x1.3333333333333p-2}}};
    const Tolerance within = {.Length = 0x1.c3008ce94f813p-21, .Radians = 1e-10};
    for (const Intersections &met :
         {Intersect(line, strung, within), Intersect(strung, line, within)})
    {
      CHECK(met.Hits.empty());
      CHECK(met.Overlaps.size() == 1);
    }
  }

  TEST_CASE("intersect: polylines with remote arcs meet at points of both, made ready or not")
  {
    /**
     *  Polylines of a few segments, most bulged from 1e-17 to 1e-3, at sizes from a thousandth to a
     *  thousand and as far from the origin, met with a line anywhere, a line nearly along the first
     *  segment, between its chord and its arc or just off them, a circle, and another such
     *  polyline, within tolerances from the numbers' grain up. Every hit's point on each curve is
     *  the curve's own at its place, the two within twice the length used of each other and the
     *  hit's point between them; and the curves made ready meet in the same hits, to the bit.
     */
    Draws       draws{0x5EED'5117};
    std::size_t hits   = 0;
    std::size_t uneven = 0;
    for (int round = 0; round < Rounds(4000); ++round)
    {
      const double scale = std::array{1.0, 1.0, 1000.0, 1e-3}.at(draws.Next() % 4);
      const double far   = scale * std::array{0.0, 1.0, 10.0}.at(draws.Next() % 3);
      const Point  about = {.X = draws.Near(0.0, 1.0) * far, .Y = draws.Near(0.0, 1.0) * far};
      const auto   spot  = [&draws, about, scale]() -> Point
      {
        return {.X = about.X + (draws.Near(0.0, 1.0) * scale),
                .Y = about.Y + (draws.Near(0.0, 1.0) * scale)};
      };
      const auto slight = [&draws]() -> double
      {
        return (draws.Next() % 2 == 0 ? 1.0 : -1.0) * std::pow(10.0, -10.0 + draws.Near(0.0, 7.0));
      };
      const auto strung = [&draws, &spot, &slight]() -> Polyline
      {
        Polyline          polyline;
        const std::size_t count = 2 + (draws.Next() % 3);
        for (std::size_t index = 0; index < count; ++index)
        {
          const double bulge =
              draws.Next() % 5 == 0 ? std::array{0.0, 0.3, -0.6}.at(draws.Next() % 3) : slight();
          polyline.Vertices.push_back({.Position = spot(), .Bulge = bulge});
        }
        polyline.Vertices.front().Bulge = slight();
        return polyline;
      };
      const Polyline first = strung();
      const Point    start = first.Vertices[0].Position;
      const Point    end   = first.Vertices[1].Position;
      const Vector   chord = end - start;
      const double   sag   = first.Vertices[0].Bulge * Length(chord) / 2.0;
      const Vector   aside = Vector{.X = chord.Y, .Y = -chord.X} / Length(chord);
      Curve          other = Line{.Start = spot(), .End = spot()};
      switch (draws.Next() % 4)
      {
      case 0:
        other = Line{.Start = (start - (chord * 0.3)) + (aside * (sag * draws.Near(0.5, 1.0))),
                     .End   = (end + (chord * 0.3)) + (aside * (sag * draws.Near(0.5, 1.0)))};
        break;
      case 1:
        other = Arc{.Center = spot(),
                    .Radius = scale * (0.05 + std::abs(draws.Near(0.0, 0.9))),
                    .Start  = Degrees(0),
                    .Sweep  = Turn};
        break;
      case 2:
        other = strung();
        break;
      default:
        break;
      }
      const double    size      = std::max(std::abs(about.X), std::abs(about.Y)) + scale;
      const Tolerance tolerance = {.Length  = size * std::pow(10.0, -10.5 + draws.Near(0.0, 3.5)),
                                   .Radians = 1e-10};
      const Curve     bowed     = first;
      /**
       *  The length used is the tolerance's, or the grain of the numbers where that is coarser: an
       *  arc not remote counts by its circle, which is at most 4096 times the size.
       */
      const double               length = std::max(tolerance.Length, 0x1p-36 * size);
      const std::array           ready  = {Prepared{bowed}, Prepared{other}};
      std::array<std::size_t, 2> counts = {};
      for (const bool flip : {false, true})
      {
        const Curves         curves = flip ? Curves{other, bowed} : Curves{bowed, other};
        const Prepared::Pair pair = flip ? Prepared::Pair{std::cref(ready[1]), std::cref(ready[0])}
                                         : Prepared::Pair{std::cref(ready[0]), std::cref(ready[1])};
        const Intersections  met  = Intersect(curves[0], curves[1], tolerance);
        Intersections        readied;
        Intersect(pair, tolerance, readied);
        counts.at(flip ? 1 : 0)  = met.Hits.size() + (100 * met.Overlaps.size());
        hits                    += met.Hits.size();
        REQUIRE(readied.Hits.size()   == met.Hits.size());
        CHECK(readied.Overlaps.size() == met.Overlaps.size());
        for (std::size_t index = 0; index < met.Hits.size(); ++index)
        {
          const Hit &hit = met.Hits[index];
          CHECK(readied.Hits[index].At   == hit.At);
          CHECK(readied.Hits[index].Kind == hit.Kind);
          CHECK(Distance(hit.First.At, At(curves[0], hit.First.Along)) <= 64 * 0x1p-52 * size);
          CHECK(Distance(hit.Second.At, At(curves[1], hit.Second.Along)) <= 64 * 0x1p-52 * size);
          if (hit.Kind == Contact::Cross)
          {
            CHECK(Distance(hit.First.At, hit.Second.At) <= length);
            CHECK(Distance(hit.At, hit.First.At) <= length);
          }
        }
      }
      uneven += counts[0] == counts[1] ? 0U : 1U;
    }
    MESSAGE("hits: " << hits << ", pairs not met alike either way round: " << uneven);
    CHECK(hits > 100);
  }

  TEST_CASE("intersect: remote arcs whose chords do not tell them are left to the thorough way")
  {
    const Tolerance usual = {.Length = 1e-10, .Radians = 1e-10};

    /**
     *  Twenty arcs each standing 5e-7 off its chord, in a zigzag a line runs through: past sixteen
     *  chords to weigh, the thorough way meets them, and finds each arc's crossing.
     */
    Polyline zigzag;
    for (int index = 0; index <= 20; ++index)
    {
      zigzag.Vertices.push_back(
          {.Position = {.X = static_cast<double>(index), .Y = index % 2 == 0 ? -0.5 : 0.5},
           .Bulge    = 1e-6});
    }
    const Curve level = Line{.Start = {.X = -1, .Y = 0.0625}, .End = {.X = 21, .Y = 0.0625}};
    const Intersections across = Intersect(Curve{zigzag}, level, usual);
    CHECK(across.Hits.size() == 20);
    CHECK(across.Overlaps.empty());
    for (const Hit &hit : across.Hits)
    {
      CHECK(hit.Kind == Contact::Cross);
      CHECK(Distance(hit.First.At, hit.Second.At) <= 0x1p-48);
    }

    /**
     *  An arc coiled past the turns a meeting looks at leaves the chord's meeting unsure, and the
     *  arc to the thorough way, which says so too.
     */
    const Curve bow    = Polyline{.Vertices = {{.Position = {.X = 0, .Y = 0}, .Bulge = 1e-6},
                                               {.Position = {.X = 1, .Y = 0}}}};
    const Curve coiled = Arc{.Center = {.X = 0.5, .Y = 0},
                             .Radius = 0.25,
                             .Start  = Degrees(0),
                             .Sweep  = Degrees(360.0 * 3000)};
    CHECK(Intersect(bow, coiled, usual).Spent);

    /**
     *  Made ready and far apart, nothing.
     */
    const Curve      away  = Line{.Start = {.X = 5, .Y = 5}, .End = {.X = 6, .Y = 7}};
    const std::array ready = {Prepared{bow}, Prepared{away}};
    Intersections    none;
    Intersect(Prepared::Pair{std::cref(ready[0]), std::cref(ready[1])}, usual, none);
    CHECK(none.Hits.empty());
    CHECK(none.Overlaps.empty());
    Marking marks;
    Intersect(Prepared::Pair{std::cref(ready[0]), std::cref(ready[1])}, usual, marks);
    CHECK(marks.Marks.empty());

    /**
     *  A line through the middle of a chord of a thousand at a sine of 1e-5, under sixteen times
     *  the arc's bulge of 1e-6: no steeper than the arc turns, it may cross the arc twice, and the
     *  chord's one crossing does not tell. It crosses the arc once, at x = 450.49, where the arc
     *  stands 4.95e-4 off its chord.
     */
    const Curve wide    = Polyline{.Vertices = {{.Position = {.X = 0, .Y = 0}, .Bulge = 1e-6},
                                                {.Position = {.X = 1000, .Y = 0}}}};
    const Curve shallow = Line{.Start = {.X = 400, .Y = -1e-3}, .End = {.X = 600, .Y = 1e-3}};
    for (const Intersections &met :
         {Intersect(wide, shallow, usual), Intersect(shallow, wide, usual)})
    {
      REQUIRE(met.Hits.size() == 1);
      CHECK(met.Hits[0].Kind  == Contact::Cross);
      CHECK(met.Hits[0].At    == Point{.X = 450.49024320360712, .Y = -0.00049509756796392893});
    }

    /**
     *  Two legs of one polyline cross the chord a five thousandth apart, within the width that
     *  takes the arc in: the chord alone meets them in fewer hits than the chords do, and does not
     *  tell. Each leg crosses the arc, half a thousandth below the chord.
     */
    const Curve legs = Polyline{.Vertices = {{.Position = {.X = 499.9999, .Y = -1}},
                                             {.Position = {.X = 499.9999, .Y = 1}},
                                             {.Position = {.X = 500.0001, .Y = 1}},
                                             {.Position = {.X = 500.0001, .Y = -1}}}};
    for (const Intersections &met : {Intersect(wide, legs, usual), Intersect(legs, wide, usual)})
    {
      REQUIRE(met.Hits.size() == 2);
      for (const Hit &hit : met.Hits)
      {
        CHECK(hit.Kind == Contact::Cross);
        CHECK(std::abs(hit.At.Y + 0.0005) <= 1e-15);
        CHECK(Distance(hit.First.At, hit.Second.At) <= 0x1p-42);
      }
    }

    /**
     *  A spline of a weight no polish holds crosses the arc twice: the crossings found on the chord
     *  are not put on the arc by a polish, and the thorough way meets the arc itself.
     */
    const Curve heavy =
        Spline{.Degree  = 2,
               .Knots   = {0, 0, 0, 1, 1, 1},
               .Poles   = {{.X = 300, .Y = -1}, {.X = 500, .Y = 3}, {.X = 700, .Y = -1}},
               .Weights = {1, 0x1p300, 1}};
    for (const Intersections &met : {Intersect(wide, heavy, usual), Intersect(heavy, wide, usual)})
    {
      REQUIRE(met.Hits.size() == 2);
      CHECK(met.Hits[0].Kind  == Contact::Cross);
      CHECK(met.Hits[1].Kind  == Contact::Cross);
    }
  }
} // namespace Geoxy
