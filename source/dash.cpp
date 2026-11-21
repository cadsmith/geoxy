/***************************************************************************************************
 *  @file      source/dash.cpp
 *  @brief     A linetype's pattern laid along a curve, walked piece by piece.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "geoxy/dash.hpp"

#include "lap.hpp"
#include "real.hpp"

#include "geoxy/accuracy.hpp"
#include "geoxy/at.hpp"
#include "geoxy/curve.hpp"
#include "geoxy/line.hpp"
#include "geoxy/part.hpp"
#include "geoxy/point.hpp"

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
     *  A pattern made ready to lay: where each element starts within it, its length, where the
     *  curve's start falls in it, and how it is laid.
     */
    struct Laying
    {
      const std::vector<double> *Elements = nullptr; /**< The elements.                        */
      std::vector<double>        Starts   = {};      /**< Where each starts in the pattern.    */
      double                     Period   = 0.0;     /**< The pattern's length.                */
      double                     Phase    = 0.0;     /**< The phase, within the period.        */
      Geoxy::Accuracy            Accuracy = {};      /**< How far a place may be off.          */
      bool                       Alone    = false;   /**< Whether one walk is the whole curve. */
    };

    /**
     *  A pattern made ready to lay.
     *
     *  @param  pattern  The pattern.
     *  @param  accuracy How far a place may be off.
     *  @return          Its starts, period and phase; a period of no length, or none, for a pattern
     *                   that lays nothing.
     */
    auto Readied(const Pattern &pattern, Accuracy accuracy) -> Laying
    {
      Laying laying = {.Elements = &pattern.Elements,
                       .Starts   = {},
                       .Period   = 0.0,
                       .Phase    = 0.0,
                       .Accuracy = accuracy,
                       .Alone    = false};
      for (const double element : pattern.Elements)
      {
        laying.Starts.push_back(laying.Period);
        laying.Period += std::abs(element);
      }
      const double phase  = std::isfinite(pattern.Phase) ? pattern.Phase : 0.0;
      laying.Phase        = std::fmod(phase, laying.Period);
      laying.Phase       += laying.Phase < 0.0 ? laying.Period : 0.0;
      return laying;
    }

    /**
     *  A stretch of a walk to cut: from one length to another, a dot where they are one.
     */
    using Interval = std::array<double, 2>;

    /**
     *  The runs of a pattern along a walk of a length: every dash clipped to it and every dot
     *  within it, in order.
     *
     *  @param  laying The pattern, ready.
     *  @param  length The walk's length, at least a period.
     *  @return        The runs; none where the pattern would be laid more than MostPatterns times.
     */
    auto Runs(const Laying &laying, double length) -> std::vector<Interval>
    {
      std::vector<Interval>      runs;
      const std::vector<double> &elements = *laying.Elements;
      if (!(length / laying.Period < static_cast<double>(MostPatterns)))
      {
        return runs;
      }
      for (std::size_t cycle = 0;; ++cycle)
      {
        for (std::size_t index = 0; index < elements.size(); ++index)
        {
          const double start =
              (static_cast<double>(cycle) * laying.Period) + laying.Starts[index] - laying.Phase;
          if (!(start <= length))
          {
            return runs;
          }
          const double element = elements[index];
          const double from    = std::max(start, 0.0);
          const double until   = std::min(start + element, length);
          if (element > 0.0 && from < until)
          {
            runs.push_back({from, until});
          }
          if (Zero(element) && !(start < 0.0))
          {
            runs.push_back({start, start});
          }
        }
      }
    }

    /**
     *  The dashes of a pattern along a walk.
     *
     *  @param  curve  The curve.
     *  @param  walk   The walk along it.
     *  @param  laying The pattern, ready.
     *  @return        The dashes; the walk whole where it is shorter than the pattern, or the
     *                 pattern too fine, the curve itself where the walk is all of it.
     */
    auto Laid(const Curve &curve, const Walk &walk, const Laying &laying) -> std::vector<Curve>
    {
      const double          length = Long(walk);
      Walker                walker = {.Laps = &walk, .Accuracy = laying.Accuracy, .Index = 0};
      std::vector<Interval> runs =
          length < laying.Period ? std::vector<Interval>{} : Runs(laying, length);
      std::vector<Curve> into;
      if (runs.empty() && laying.Alone)
      {
        return {curve};
      }
      if (runs.empty())
      {
        runs.push_back({0.0, length});
      }
      for (const auto &[from, until] : runs)
      {
        const double start = Trodden(walker, from);
        if (Exactly(from, until))
        {
          const Point point = At(curve, start);
          into.emplace_back(Line{.Start = point, .End = point});
          continue;
        }
        into.push_back(
            Part(curve, {.Start = start, .End = Trodden(walker, until)}).value_or(curve));
      }
      return into;
    }
  } // namespace

  auto Dashes(const Curve &curve, const Pattern &pattern, Accuracy accuracy) -> std::vector<Curve>
  {
    Laying laying = Readied(pattern, accuracy);
    if (!(laying.Period > 0.0) || !std::isfinite(laying.Period))
    {
      return {curve};
    }
    const std::vector<Walk> walks = Walks(curve, pattern.Across, accuracy);
    if (walks.empty())
    {
      return {curve};
    }
    laying.Alone = walks.size() == 1;
    std::vector<Curve> dashes;
    for (const Walk &walk : walks)
    {
      const std::vector<Curve> laid = Laid(curve, walk, laying);
      dashes.insert(dashes.end(), laid.begin(), laid.end());
    }
    return dashes;
  }
} // namespace Geoxy
