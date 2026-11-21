/***************************************************************************************************
 *  @file      tests/suites/score.cpp
 *  @brief     How near and how fast geoxy meets every case of the suites: for each crossing and
 *             each end, how many steps of a double its point is from the exact one; for each case,
 *             how long a meeting takes. With --check, untimed, it fails where a hit is past half a
 *             step, but for the cases known to be past it, each within its budget.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "cases.hpp"

#include "geoxy/at.hpp"
#include "geoxy/curve.hpp"
#include "geoxy/domain.hpp"
#include "geoxy/intersect.hpp"
#include "geoxy/point.hpp"
#include "geoxy/tolerance.hpp"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <exception>
#include <filesystem>
#include <limits>
#include <map>
#include <print>
#include <span>
#include <string>
#include <string_view>
#include <vector>

namespace
{
  using namespace Geoxy;
  using namespace Geoxy::Suites;

  /**
   *  What one group of cases came to.
   */
  struct Score
  {
    std::size_t         Cases   = 0;    /**< The cases.                              */
    std::vector<double> Errors  = {};   /**< Each hit's error, in steps of a double. */
    double              Seconds = 0;    /**< The time all its meetings took.         */
    bool                Timed   = true; /**< Whether its meetings are timed.         */
  };

  /**
   *  A hit's error, and the case it is of.
   */
  struct Worst
  {
    double      Error = 0;  /**< The error, in steps of a double. */
    std::string Name  = {}; /**< The case.                        */
  };

  /**
   *  The distance between neighbouring doubles at a size: one step of a double there.
   *
   *  @param  size The size.
   *  @return      The step.
   */
  auto Step(double size) -> double
  {
    const double magnitude = std::abs(size);
    return std::nextafter(magnitude, std::numeric_limits<double>::infinity()) - magnitude;
  }

  /**
   *  The size of a case's coordinates: the largest magnitude of its curves' starts and its points.
   *
   *  @param  first  One curve.
   *  @param  second The other.
   *  @param  one    The case.
   *  @return        The size, at least 1.
   */
  auto Size(const Curve &first, const Curve &second, const Case &one) -> double
  {
    double size = 1.0;
    for (const Curve *curve : {&first, &second})
    {
      const Point start = At(*curve, Domain(*curve).Start);
      if (std::isfinite(start.X) && std::isfinite(start.Y))
      {
        size = std::max({size, std::abs(start.X), std::abs(start.Y)});
      }
    }
    for (const Expected &point : one.Points)
    {
      size = std::max({size, std::abs(point.At.X), std::abs(point.At.Y)});
    }
    return size;
  }

  /**
   *  How long one meeting takes: the least, over three rounds of at least two milliseconds each, of
   *  a round's time over its meetings.
   *
   *  @param  first     One curve.
   *  @param  second    The other.
   *  @param  tolerance The tolerance.
   *  @return           The seconds.
   */
  auto Timed(const Curve &first, const Curve &second, const Tolerance &tolerance) -> double
  {
    using Clock  = std::chrono::steady_clock;
    double least = std::numeric_limits<double>::infinity();
    for (int round = 0; round < 3; ++round)
    {
      std::size_t count = 0;
      const auto  start = Clock::now();
      auto        now   = start;
      while (count == 0 || now - start < std::chrono::milliseconds(2))
      {
        count += Intersect(first, second, tolerance).Hits.size() + 1 > 0 ? 1U : 0U;
        now    = Clock::now();
      }
      least = std::min(least, std::chrono::duration<double>(now - start).count()
                                  / static_cast<double>(count));
    }
    return least;
  }

  /**
   *  A case's group: its source and family, the part of its name before the colon.
   *
   *  @param  one The case.
   *  @return     The group.
   */
  auto Group(const Case &one) -> std::string
  {
    return one.Name.substr(0, one.Name.find(':'));
  }

  /**
   *  The share of errors at most a bound.
   *
   *  @param  errors The errors.
   *  @param  bound  The bound.
   *  @return        The percentage.
   */
  auto Within(const std::vector<double> &errors, double bound) -> double
  {
    const auto count =
        std::ranges::count_if(errors, [bound](double error) -> bool { return error <= bound; });
    return errors.empty() ? 100.0
                          : 100.0 * static_cast<double>(count) / static_cast<double>(errors.size());
  }

  /**
   *  A case scored: its meeting's time, and each crossing's and end's error.
   *
   *  @param one   The case.
   *  @param score Its group's score so far.
   *  @param worst Every hit's error so far.
   */
  void Scored(const Case &one, Score &score, std::vector<Worst> &worst)
  {
    if (!one.First || !one.Second || !one.Limit.empty() || one.Disabled)
    {
      return;
    }
    const Curve        &first     = *one.First;
    const Curve        &second    = *one.Second;
    const double        size      = Size(first, second, one);
    const double        given     = one.Relative ? one.Within * size : one.Within;
    const Tolerance     tolerance = {.Length = given > 0.0 ? given : 1e-12 * size, .Radians = 1e-9};
    const Intersections meeting   = Intersect(first, second, tolerance);
    ++score.Cases;
    score.Seconds += score.Timed ? Timed(first, second, tolerance) : 0.0;
    for (const Expected &expected : one.Points)
    {
      if (expected.Kind == "Touch")
      {
        continue;
      }
      double nearest = std::numeric_limits<double>::infinity();
      for (const Hit &hit : meeting.Hits)
      {
        nearest = std::min(nearest, std::max(std::abs(hit.At.X - expected.At.X),
                                             std::abs(hit.At.Y - expected.At.Y)));
      }
      score.Errors.push_back(nearest / Step(size));
      worst.push_back({.Error = nearest / Step(size), .Name = one.Name});
    }
  }

  /**
   *  The scores printed, a group a row, and the worst hits.
   *
   *  @param scores Each group's score.
   *  @param worst  Every hit's error.
   */
  void Printed(std::map<std::string, Score> &scores, std::vector<Worst> &worst)
  {
    std::println("{:<22} {:>6} {:>6} {:>8} {:>8} {:>8} {:>10} {:>10} {:>10}", "group", "cases",
                 "hits", "exact %", "<=1 %", "<=16 %", "median", "worst", "us/meet");
    for (auto &[group, score] : scores)
    {
      std::ranges::sort(score.Errors);
      const double median = score.Errors.empty() ? 0 : score.Errors[score.Errors.size() / 2];
      const double most   = score.Errors.empty() ? 0 : score.Errors.back();
      std::println("{:<22} {:>6} {:>6} {:>8.1f} {:>8.1f} {:>8.1f} {:>10.3g} {:>10.3g} {:>10.2f}",
                   group, score.Cases, score.Errors.size(), Within(score.Errors, 0.5),
                   Within(score.Errors, 1), Within(score.Errors, 16), median, most,
                   1e6 * score.Seconds / static_cast<double>(score.Cases));
    }
    std::ranges::sort(worst, [](const Worst &lhs, const Worst &rhs) -> bool
                      { return lhs.Error > rhs.Error; });
    std::println("\nworst hits, in steps of a double at the case's size:");
    for (std::size_t index = 0; index < std::min<std::size_t>(worst.size(), 15); ++index)
    {
      std::println("{:>12.4g}  {}", worst[index].Error, worst[index].Name);
    }
  }

  /**
   *  Whether every hit is as near as allowed: the double nearest the exact point, within half a
   *  step of a double at the case's size; each that is not printed.
   *
   *  @param  worst Every hit's error.
   *  @return       Whether all are.
   */
  auto Checked(const std::vector<Worst> &worst) -> bool
  {
    bool kept = true;
    for (const Worst &hit : worst)
    {
      const double allowed = 0.5;
      if (!(hit.Error <= allowed))
      {
        std::println(stderr, "past its budget of {} steps: {} steps, {}", allowed, hit.Error,
                     hit.Name);
        kept = false;
      }
    }
    return kept;
  }
} // namespace

auto main(int argc, char **argv) -> int
{
  try
  {
    const std::span<char *> arguments{argv, static_cast<std::size_t>(argc)};
    const bool check = arguments.size() > 1 && std::string_view{arguments[1]} == "--check";
    const std::array<std::filesystem::path, 1> roots = {std::filesystem::path{GEOXY_SUITES_TREE}};
    std::map<std::string, Score>               scores;
    std::vector<Worst>                         worst;
    for (const std::filesystem::path &file : Files(roots))
    {
      for (const Case &one : Read(file))
      {
        Score &score = scores[Group(one)];
        score.Timed  = !check;
        Scored(one, score, worst);
      }
    }
    const bool kept = !check || Checked(worst);
    Printed(scores, worst);
    return kept ? 0 : 1;
  }
  catch (const std::exception &error)
  {
    std::fputs(error.what(), stderr);
    std::fputc('\n', stderr);
    return 1;
  }
  catch (...)
  {
    return 1;
  }
}
