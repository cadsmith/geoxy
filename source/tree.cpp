/***************************************************************************************************
 *  @file      source/tree.cpp
 *  @brief     Curves kept in a tree of their boxes: what is near a box, nearest a point and struck
 *             by a ray, asked of the few curves the tree leaves.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "geoxy/tree.hpp"

#include "boughs.hpp"
#include "real.hpp"
#include "thorough.hpp"

#include "geoxy/at.hpp"
#include "geoxy/box.hpp"
#include "geoxy/contain.hpp"
#include "geoxy/curve.hpp"
#include "geoxy/extents.hpp"
#include "geoxy/intersect.hpp"
#include "geoxy/nearest.hpp"
#include "geoxy/place.hpp"
#include "geoxy/point.hpp"
#include "geoxy/ray.hpp"
#include "geoxy/tolerance.hpp"
#include "geoxy/vector.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <memory>
#include <span>
#include <tuple>
#include <utility>
#include <vector>

namespace Geoxy
{
  namespace
  {
    /**
     *  How much of its size a product is taken to be off by at most, well past its rounding.
     */
    constexpr double Slack = 0x1p-48;

    /**
     *  Whether a ray may come within a distance of a box: not when the box is wholly to one side of
     *  the ray's line, nor when it is wholly behind the ray's start, each told only where rounding
     *  cannot have told it wrong.
     *
     *  @param  ray    The ray.
     *  @param  box    The box.
     *  @param  within The distance.
     *  @return        Whether it may.
     */
    auto Aimed(const Ray &ray, const Box &box, double within) -> bool
    {
      const std::array<Point, 4> corners = {
          Point{.X = box.Min.X - within, .Y = box.Min.Y - within},
          Point{.X = box.Max.X + within, .Y = box.Min.Y - within},
          Point{.X = box.Max.X + within, .Y = box.Max.Y + within},
          Point{.X = box.Min.X - within, .Y = box.Max.Y + within}};
      int left   = 0;
      int right  = 0;
      int behind = 0;
      for (const Point corner : corners)
      {
        const Vector reach = corner - ray.Start;
        const double side  = Cross(ray.Way, reach);
        const double along = Dot(ray.Way, reach);
        const double sway = Slack * (std::abs(ray.Way.X * reach.Y) + std::abs(ray.Way.Y * reach.X));
        const double give = Slack * (std::abs(ray.Way.X * reach.X) + std::abs(ray.Way.Y * reach.Y));
        left   += side > sway ? 1 : 0;
        right  += side < -sway ? 1 : 0;
        behind += along < -give ? 1 : 0;
      }
      return left < 4 && right < 4 && behind < 4;
    }
  } // namespace

  Tree::Tree() : m_Boughs{std::make_unique<Boughs>()} {}

  Tree::Tree(std::span<const Curve> curves)
      : m_Curves{curves.begin(), curves.end()}, m_Held(curves.size(), 1), m_Count{curves.size()}
  {
    std::vector<Box> boxes;
    boxes.reserve(curves.size());
    for (const Curve &curve : curves)
    {
      boxes.push_back(Extents(curve));
      m_Size = std::max(m_Size, Sized(curve));
    }
    m_Boughs = std::make_unique<Boughs>(Boughs::Packed(boxes));
  }

  Tree::Tree(const Tree &other)
      : m_Curves{other.m_Curves}, m_Held{other.m_Held}, m_Spare{other.m_Spare},
        m_Count{other.m_Count}, m_Size{other.m_Size},
        m_Boughs{std::make_unique<Boughs>(*other.m_Boughs)}
  {
  }

  Tree::Tree(Tree &&other) noexcept = default;

  auto Tree::operator=(const Tree &other) -> Tree &
  {
    *this = Tree{other};
    return *this;
  }

  auto Tree::operator=(Tree &&other) noexcept -> Tree & = default;

  Tree::~Tree() = default;

  auto Tree::Insert(Curve curve) -> std::size_t
  {
    const Box   box   = Extents(curve);
    std::size_t index = m_Curves.size();
    m_Size            = std::max(m_Size, Sized(curve));
    if (m_Spare.empty())
    {
      m_Curves.push_back(std::move(curve));
      m_Held.push_back(1);
    }
    else
    {
      index = m_Spare.back();
      m_Spare.pop_back();
      m_Curves[index] = std::move(curve);
      m_Held[index]   = 1;
    }
    ++m_Count;
    m_Boughs->Added(index, box);
    return index;
  }

  void Tree::Erase(std::size_t index)
  {
    if (!Holds(index))
    {
      return;
    }
    m_Boughs->Dropped(index);
    m_Curves[index] = Curve{};
    m_Held[index]   = 0;
    m_Spare.push_back(index);
    --m_Count;
  }

  auto Tree::Near(const Box &box) const -> std::vector<std::size_t>
  {
    std::vector<std::size_t> near;
    m_Boughs->Each([&box](const Box &bounds) -> bool { return Overlaps(bounds, box); },
                   [&near](std::size_t index) -> void { near.push_back(index); });
    std::ranges::sort(near);
    return near;
  }

  auto Tree::Window(const Box &box) const -> std::vector<std::size_t>
  {
    std::vector<std::size_t> near = Near(box);
    std::erase_if(near, [this, &box](std::size_t index) -> bool
                  { return !Inside(m_Curves[index], box); });
    return near;
  }

  auto Tree::Crossing(const Box &box) const -> std::vector<std::size_t>
  {
    std::vector<std::size_t> near = Near(box);
    std::erase_if(near, [this, &box](std::size_t index) -> bool
                  { return !Intersects(m_Curves[index], box); });
    return near;
  }

  auto Tree::Nearest(Point point, std::size_t count, double within) const -> std::vector<Close>
  {
    const auto found = m_Boughs->Closest(
        point, {.Count = count, .Within = within}, [this, point](std::size_t index) -> double
        { return Distance(Geoxy::Nearest(m_Curves[index], point).At, point); });
    std::vector<Close> nearest;
    nearest.reserve(found.size());
    for (const Boughs::Nearer &one : found)
    {
      nearest.push_back({.Index    = one.Mark,
                         .Nearest  = Geoxy::Nearest(m_Curves[one.Mark], point),
                         .Distance = one.Distance});
    }
    return nearest;
  }

  auto Tree::Cast(const Ray &ray, const Tolerance &tolerance) const -> std::vector<Strike>
  {
    /**
     *  Intersect takes a length under 16 steps of a double at the size of two curves' numbers as
     *  that, so a box is passed by only where it lies farther from the ray than the most that can
     *  be for any curve held.
     */
    std::vector<Strike> strikes;
    const Curve         aimed = ray;
    const double        size  = std::max(m_Size, Sized(aimed));
    const double        reach = std::max(tolerance.Length, 16.0 * Above(size));
    m_Boughs->Each([&ray, reach](const Box &bounds) -> bool { return Aimed(ray, bounds, reach); },
                   [this, &aimed, &tolerance, &strikes](std::size_t index) -> void
                   {
                     const Intersections met = Intersect(aimed, m_Curves[index], tolerance);
                     for (const Hit &hit : met.Hits)
                     {
                       strikes.push_back({.Index = index, .At = hit.At, .Along = hit.First.Along});
                     }
                     for (const Overlap &overlap : met.Overlaps)
                     {
                       for (const double along : {overlap.First.Start, overlap.First.End})
                       {
                         strikes.push_back(
                             {.Index = index, .At = Geoxy::At(aimed, along), .Along = along});
                       }
                     }
                   });
    std::ranges::sort(strikes, {}, [](const Strike &one) -> std::tuple<double, std::size_t>
                      { return {one.Along, one.Index}; });
    return strikes;
  }
} // namespace Geoxy
