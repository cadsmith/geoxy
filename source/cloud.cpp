/***************************************************************************************************
 *  @file      source/cloud.cpp
 *  @brief     Points kept in a tree: what a box holds and what is nearest a point, asked of the few
 *             points the tree leaves.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "geoxy/cloud.hpp"

#include "boughs.hpp"

#include "geoxy/box.hpp"
#include "geoxy/point.hpp"

#include <algorithm>
#include <cstddef>
#include <memory>
#include <span>
#include <vector>

namespace Geoxy
{
  Cloud::Cloud() : m_Boughs{std::make_unique<Boughs>()} {}

  Cloud::Cloud(std::span<const Point> points)
      : m_Points{points.begin(), points.end()}, m_Held(points.size(), 1), m_Count{points.size()}
  {
    std::vector<Box> boxes;
    boxes.reserve(points.size());
    for (const Point point : points)
    {
      boxes.push_back({.Min = point, .Max = point});
    }
    m_Boughs = std::make_unique<Boughs>(Boughs::Packed(boxes));
  }

  Cloud::Cloud(const Cloud &other)
      : m_Points{other.m_Points}, m_Held{other.m_Held}, m_Spare{other.m_Spare},
        m_Count{other.m_Count}, m_Boughs{std::make_unique<Boughs>(*other.m_Boughs)}
  {
  }

  Cloud::Cloud(Cloud &&other) noexcept = default;

  auto Cloud::operator=(const Cloud &other) -> Cloud &
  {
    *this = Cloud{other};
    return *this;
  }

  auto Cloud::operator=(Cloud &&other) noexcept -> Cloud & = default;

  Cloud::~Cloud() = default;

  auto Cloud::Insert(Point point) -> std::size_t
  {
    std::size_t index = m_Points.size();
    if (m_Spare.empty())
    {
      m_Points.push_back(point);
      m_Held.push_back(1);
    }
    else
    {
      index = m_Spare.back();
      m_Spare.pop_back();
      m_Points[index] = point;
      m_Held[index]   = 1;
    }
    ++m_Count;
    m_Boughs->Added(index, {.Min = point, .Max = point});
    return index;
  }

  void Cloud::Erase(std::size_t index)
  {
    if (!Holds(index))
    {
      return;
    }
    m_Boughs->Dropped(index);
    m_Held[index] = 0;
    m_Spare.push_back(index);
    --m_Count;
  }

  auto Cloud::Within(const Box &box) const -> std::vector<std::size_t>
  {
    std::vector<std::size_t> within;
    m_Boughs->Each([&box](const Box &bounds) -> bool { return Overlaps(bounds, box); },
                   [this, &box, &within](std::size_t index) -> void
                   {
                     if (Contains(box, m_Points[index]))
                     {
                       within.push_back(index);
                     }
                   });
    std::ranges::sort(within);
    return within;
  }

  auto Cloud::Nearest(Point point, std::size_t count, double within) const -> std::vector<Close>
  {
    const auto         found = m_Boughs->Closest(point, {.Count = count, .Within = within},
                                                 [this, point](std::size_t index) -> double
                                                 { return Distance(m_Points[index], point); });
    std::vector<Close> nearest;
    nearest.reserve(found.size());
    for (const Boughs::Nearer &one : found)
    {
      nearest.push_back({.Index = one.Mark, .Distance = one.Distance});
    }
    return nearest;
  }
} // namespace Geoxy
