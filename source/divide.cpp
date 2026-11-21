/***************************************************************************************************
 *  @file      source/divide.cpp
 *  @brief     Places at equal lengths along a curve, walked once.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "geoxy/divide.hpp"

#include "lap.hpp"

#include "geoxy/accuracy.hpp"
#include "geoxy/at.hpp"
#include "geoxy/curve.hpp"
#include "geoxy/place.hpp"

#include <cmath>
#include <cstddef>
#include <vector>

namespace Geoxy
{
  namespace
  {
    /**
     *  The lengths along a curve to set places at, rising, and how far each may be off.
     */
    struct Marks
    {
      std::vector<double> Distances = {}; /**< The lengths from the start. */
      Geoxy::Accuracy     Accuracy  = {}; /**< How far each may be off.    */
    };

    /**
     *  The places along a curve's walk at lengths from its start.
     *
     *  @param  curve The curve.
     *  @param  walk  Its walk.
     *  @param  marks The lengths, rising, and the accuracy.
     *  @return       The places.
     */
    auto Set(const Curve &curve, const Walk &walk, const Marks &marks) -> std::vector<Place>
    {
      Walker             walker = {.Laps = &walk, .Accuracy = marks.Accuracy, .Index = 0};
      std::vector<Place> places;
      places.reserve(marks.Distances.size());
      for (const double distance : marks.Distances)
      {
        const double along = Trodden(walker, distance);
        places.push_back({.At = At(curve, along), .Along = along});
      }
      return places;
    }
  } // namespace

  auto Divide(const Curve &curve, std::size_t parts, Accuracy accuracy) -> std::vector<Place>
  {
    if (parts < 2 || parts > MostPlaces)
    {
      return {};
    }
    const std::vector<Walk> walks = Walks(curve, true, accuracy);
    if (walks.empty())
    {
      return {};
    }
    const double length = Long(walks.front());
    Marks        marks  = {.Distances = {}, .Accuracy = accuracy};
    marks.Distances.reserve(parts - 1);
    for (std::size_t index = 1; index < parts; ++index)
    {
      marks.Distances.push_back(length * static_cast<double>(index) / static_cast<double>(parts));
    }
    return Set(curve, walks.front(), marks);
  }

  auto Measure(const Curve &curve, double length, Accuracy accuracy) -> std::vector<Place>
  {
    const std::vector<Walk> walks = Walks(curve, true, accuracy);
    if (walks.empty() || !(length > 0.0) || !std::isfinite(length))
    {
      return {};
    }
    const double total = std::floor(Long(walks.front()) / length);
    if (total > static_cast<double>(MostPlaces))
    {
      return {};
    }
    const auto count = static_cast<std::size_t>(total);
    Marks      marks = {.Distances = {}, .Accuracy = accuracy};
    marks.Distances.reserve(count);
    for (std::size_t index = 1; index <= count; ++index)
    {
      marks.Distances.push_back(length * static_cast<double>(index));
    }
    return Set(curve, walks.front(), marks);
  }
} // namespace Geoxy
