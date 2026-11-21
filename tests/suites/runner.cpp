/***************************************************************************************************
 *  @file      tests/suites/runner.cpp
 *  @brief     Suites of cases run against geoxy, its own made on an exact oracle and other engines'
 *             as references: each case two curves and the hits and shared stretches expected.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "cases.hpp"
#include "respelled.hpp"

#include "geoxy/at.hpp"
#include "geoxy/curve.hpp"
#include "geoxy/intersect.hpp"
#include "geoxy/point.hpp"

#include <doctest/doctest.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <filesystem>
#include <optional>
#include <span>
#include <sstream>
#include <string>
#include <string_view>
#include <vector>

namespace Geoxy::Suites
{
  namespace
  {
    /**
     *  How far apart two points may be for a case: its tolerance, scaled by their size where the
     *  case says so, and never under a few steps of a double there.
     *
     *  @param  one   The case.
     *  @param  point A point.
     *  @return       The distance.
     */
    auto Slack(const Case &one, Point point) -> double
    {
      const double size  = std::max({1.0, std::abs(point.X), std::abs(point.Y)});
      const double grain = 64.0 * size * 2.220446049250313e-16;
      return std::max(one.Relative ? one.Within * size : one.Within, grain);
    }

    /**
     *  A case's curves, read, and the tolerance to meet them by.
     */
    struct Posed
    {
      Curve  First  = {};  /**< The first curve. */
      Curve  Second = {};  /**< The second.      */
      double Within = 0.0; /**< The tolerance.   */
    };

    /**
     *  A case posed: its curves, and its tolerance, scaled by the case's size where it says so, or
     *  where it is nothing a millionth of a millionth of that size.
     *
     *  @param  one The case.
     *  @return     It posed; nothing where its curves do not read.
     */
    auto Pose(const Case &one) -> std::optional<Posed>
    {
      if (!one.First || !one.Second)
      {
        return std::nullopt;
      }
      Posed  posed = {.First = *one.First, .Second = *one.Second};
      double size  = 1.0;
      for (const Expected &point : one.Points)
      {
        size = std::max({size, std::abs(point.At.X), std::abs(point.At.Y)});
      }
      for (const Curve *curve : {&posed.First, &posed.Second})
      {
        const Point start = At(*curve, 0.0);
        size              = std::max({size, std::abs(start.X), std::abs(start.Y)});
      }
      const double given = one.Relative ? one.Within * size : one.Within;
      posed.Within       = given > 0.0 ? given : 1e-12 * size;
      return posed;
    }

    /**
     *  A contact's name, as cases write it.
     *
     *  @param  kind The contact.
     *  @return      Its name.
     */
    auto Named(Contact kind) -> std::string_view
    {
      switch (kind)
      {
      case Contact::Cross:
        return "Cross";
      case Contact::Touch:
        return "Touch";
      case Contact::End:
        return "End";
      }
      return "?";
    }

    /**
     *  The stretches shared, each by its ends on the first curve, those that run on joined, past a
     *  closed curve's seam too.
     *
     *  @param  curve   The first curve.
     *  @param  meeting What meeting it found.
     *  @param  within  The tolerance.
     *  @return         The stretches' ends.
     */
    auto Stretches(const Curve &curve, const Intersections &meeting, double within)
        -> std::vector<std::array<Point, 2>>
    {
      std::vector<std::array<Point, 2>> shared;
      for (const Overlap &overlap : meeting.Overlaps)
      {
        const std::array<Point, 2> ends = {At(curve, overlap.First.Start),
                                           At(curve, overlap.First.End)};
        if (!shared.empty() && Distance(shared.back()[1], ends[0]) <= within)
        {
          shared.back()[1] = ends[1];
        }
        else
        {
          shared.push_back(ends);
        }
      }
      if (shared.size() >= 2 && Distance(shared.back()[1], shared.front()[0]) <= within)
      {
        shared.front()[0] = shared.back()[0];
        shared.pop_back();
      }
      return shared;
    }

    /**
     *  Whether meeting a case's curves the other way round finds the same hits, of the same kinds:
     *  each within twice the slack of the hit expected nearest it, as both may be off by that much,
     *  or twice the tolerance where none is expected.
     *
     *  @param  one     The case.
     *  @param  posed   It posed.
     *  @param  meeting What meeting its curves found.
     *  @return         Whether it does.
     */
    auto Swapped(const Case &one, const Posed &posed, const Intersections &meeting) -> bool
    {
      const double        within = posed.Within;
      const Intersections back =
          Intersect(posed.Second, posed.First, {.Length = within, .Radians = 1e-9});
      const auto slack = [&one, within](Point point) -> double
      {
        const auto nearest =
            std::ranges::min_element(one.Points, {}, [point](const Expected &expected) -> double
                                     { return Distance(expected.At, point); });
        return 2 * std::max(within, nearest == one.Points.end() ? within : nearest->Slack);
      };
      return back.Hits.size() == meeting.Hits.size()
          && std::ranges::all_of(
                 back.Hits,
                 [&meeting, &slack](const Hit &hit) -> bool
                 {
                   const double reach = slack(hit.At);
                   return std::ranges::any_of(
                       meeting.Hits, [&hit, reach](const Hit &other) -> bool
                       { return other.Kind == hit.Kind && Distance(hit.At, other.At) <= reach; });
                 });
    }

    /**
     *  Whether a hit found is one a case expects: of its kind, if it says, and within a slack.
     *
     *  @param  expected The hit expected.
     *  @param  hits     The hits found.
     *  @param  slack    How far off it may be.
     *  @return          Whether one of them is.
     */
    auto Found(const Expected &expected, std::span<const Hit> hits, double slack) -> bool
    {
      return std::ranges::any_of(hits,
                                 [&expected, slack](const Hit &hit) -> bool
                                 {
                                   const bool kind =
                                       expected.Kind.empty() || expected.Kind == Named(hit.Kind);
                                   return kind && Distance(hit.At, expected.At) <= slack;
                                 });
    }

    /**
     *  Whether a stretch found has the ends a case expects, either way round, an end at infinity
     *  matched by one there.
     *
     *  @param  found The stretch found.
     *  @param  ends  The ends expected.
     *  @param  slack How far off an end may be.
     *  @return       Whether it does.
     */
    auto Ends(const std::array<Point, 2> &found, const std::array<Point, 2> &ends, double slack)
        -> bool
    {
      const auto same = [slack](Point lhs, Point rhs) -> bool
      {
        if (!std::isfinite(rhs.X) || !std::isfinite(rhs.Y))
        {
          return !std::isfinite(lhs.X) || !std::isfinite(lhs.Y);
        }
        return Distance(lhs, rhs) <= slack;
      };
      return (same(found[0], ends[0]) && same(found[1], ends[1]))
          || (same(found[0], ends[1]) && same(found[1], ends[0]));
    }

    /**
     *  Two points in the case files' text.
     *
     *  @param told Where to write.
     *  @param ends The points.
     */
    void Write(std::ostringstream &told, const std::array<Point, 2> &ends)
    {
      told
          << " ["
          << ends[0].X
          << ' '
          << ends[0].Y
          << " .. "
          << ends[1].X
          << ' '
          << ends[1].Y
          << ']';
    }

    /**
     *  What geoxy found for a case, and what the case expected.
     *
     *  @param  one     The case.
     *  @param  posed   It posed.
     *  @param  meeting What geoxy found.
     *  @return         The two, side by side.
     */
    auto Told(const Case &one, const Posed &posed, const Intersections &meeting) -> std::string
    {
      std::ostringstream told;
      told.precision(17);
      told << "found";
      for (const Hit &hit : meeting.Hits)
      {
        told << " (" << hit.At.X << ' ' << hit.At.Y << ' ' << Named(hit.Kind) << ')';
      }
      for (const auto &found : Stretches(posed.First, meeting, posed.Within))
      {
        Write(told, found);
      }
      told
          << (meeting.Spent ? " spent" : "")
          << (Swapped(one, posed, meeting) ? "" : " (not the same swapped)")
          << "; expected";
      for (const Expected &point : one.Points)
      {
        told << " (" << point.At.X << ' ' << point.At.Y << ' ' << point.Kind << ')';
      }
      for (const auto &ends : one.Shared)
      {
        Write(told, ends);
      }
      return told.str();
    }

    /**
     *  What is wrong with geoxy's answer to a case.
     *
     *  @param  one The case.
     *  @return     Nothing when it agrees; otherwise what it found and what was expected.
     */
    auto Judge(const Case &one) -> std::string
    {
      const std::optional<Posed> posed = Pose(one);
      if (!posed)
      {
        return "its curves do not read";
      }
      const double        within = posed->Within;
      const Intersections meeting =
          Intersect(posed->First, posed->Second, {.Length = within, .Radians = 1e-9});
      const auto shared = Stretches(posed->First, meeting, within);
      const bool hits =
          meeting.Hits.size() == one.Points.size()
          && std::ranges::all_of(one.Points,
                                 [&one, &meeting, within](const Expected &expected) -> bool
                                 {
                                   const double slack =
                                       std::max({within, expected.Slack, Slack(one, expected.At)});
                                   return Found(expected, meeting.Hits, slack);
                                 });
      const bool stretches =
          shared.size() == one.Shared.size()
          && std::ranges::all_of(one.Shared,
                                 [&one, &shared, within](const std::array<Point, 2> &ends) -> bool
                                 {
                                   const double slack = std::max(within, Slack(one, ends[0]));
                                   return std::ranges::any_of(
                                       shared,
                                       [&ends, slack](const std::array<Point, 2> &found) -> bool
                                       { return Ends(found, ends, slack); });
                                 });
      if (hits && stretches && Swapped(one, *posed, meeting))
      {
        return {};
      }
      return Told(one, *posed, meeting);
    }

    /**
     *  What is wrong with geoxy's answer to a case once its curves are spelled another way, the
     *  same geometry: a line as a Bézier curve, an arc or a polyline as a spline. The answer must
     *  not change with the spelling.
     *
     *  @param  one The case.
     *  @return     Nothing when it agrees, or where neither curve has another spelling.
     */
    auto Respelling(const Case &one) -> std::string
    {
      if (!one.First || !one.Second)
      {
        return {};
      }
      Case                       other  = one;
      const std::optional<Curve> first  = Respelled(*one.First);
      const std::optional<Curve> second = Respelled(*one.Second);
      if (!first && !second)
      {
        return {};
      }
      other.First           = first ? first : one.First;
      other.Second          = second ? second : one.Second;
      const std::string why = Judge(other);
      return why.empty() ? std::string{} : "; spelled as splines, " + why;
    }

    /**
     *  How the cases came out.
     */
    struct Tally
    {
      std::size_t Total    = 0; /**< The cases run.                                    */
      std::size_t Agreed   = 0; /**< Those geoxy agrees with.                          */
      std::size_t Disabled = 0; /**< Those their source skips, where geoxy differs.    */
      std::size_t Limited  = 0; /**< Those past what geoxy holds, where geoxy differs. */
    };

    /**
     *  A case run and counted: a failure unless geoxy agrees, its source skips it, or it is past
     *  what geoxy can hold, each of those told.
     *
     *  @param one   The case.
     *  @param tally The count so far.
     */
    void Weigh(const Case &one, Tally &tally)
    {
      ++tally.Total;
      const std::string wrong = Judge(one).append(Respelling(one));
      if (wrong.empty())
      {
        ++tally.Agreed;
      }
      else if (one.Disabled)
      {
        ++tally.Disabled;
        MESSAGE("disabled upstream, and geoxy differs: " << one.Name << ": " << wrong);
      }
      else if (!one.Limit.empty())
      {
        ++tally.Limited;
        MESSAGE("past what geoxy holds, and geoxy differs: "
                << one.Name
                << ": "
                << one.Limit
                << ": "
                << wrong);
      }
      else
      {
        CHECK_MESSAGE(wrong.empty(), one.Name << ": " << wrong);
      }
    }
  } // namespace

  TEST_CASE("suites: every case meets as its suite says")
  {
    Tally                                      tally;
    const std::array<std::filesystem::path, 1> roots = {std::filesystem::path{GEOXY_SUITES_TREE}};
    for (const std::filesystem::path &root : roots)
    {
      if (!std::filesystem::is_directory(root))
      {
        MESSAGE("no suites at " << root.string());
      }
    }
    for (const std::filesystem::path &file : Files(roots))
    {
      for (const Case &one : Read(file))
      {
        Weigh(one, tally);
      }
    }
    MESSAGE(tally.Agreed
            << " of "
            << tally.Total
            << " cases agree; "
            << tally.Disabled
            << " that their source skips differ; "
            << tally.Limited
            << " past what geoxy can hold differ");
  }
} // namespace Geoxy::Suites
