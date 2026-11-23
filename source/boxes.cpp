/***************************************************************************************************
 *  @file      source/boxes.cpp
 *  @brief     Boxes kept in a tree by the caller's tags.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "geoxy/boxes.hpp"

#include "boughs.hpp"

#include "geoxy/box.hpp"
#include "geoxy/point.hpp"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <optional>
#include <span>
#include <type_traits>
#include <utility>
#include <vector>

namespace Geoxy
{
  Boxes::Boxes() : m_Boughs{std::make_unique<Boughs>()} {}

  Boxes::Boxes(std::span<const Entry> entries)
  {
    for (const Entry &entry : entries)
    {
      const auto [found, fresh] = m_Marks.try_emplace(entry.Tag, m_Tags.size());
      if (fresh)
      {
        m_Tags.push_back(entry.Tag);
        m_Kept.push_back(entry.Bounds);
      }
      else
      {
        m_Kept[found->second] = entry.Bounds;
      }
    }
    m_Boughs = std::make_unique<Boughs>(Boughs::Packed(m_Kept));
  }

  Boxes::Boxes(const Boxes &other)
      : m_Marks{other.m_Marks}, m_Tags{other.m_Tags}, m_Kept{other.m_Kept}, m_Spare{other.m_Spare},
        m_Boughs{std::make_unique<Boughs>(*other.m_Boughs)}
  {
  }

  Boxes::Boxes(Boxes &&other) noexcept = default;

  auto Boxes::operator=(const Boxes &other) -> Boxes &
  {
    *this = Boxes{other};
    return *this;
  }

  auto Boxes::operator=(Boxes &&other) noexcept -> Boxes & = default;

  Boxes::~Boxes() = default;

  auto Boxes::Bounds(std::uint64_t tag) const -> std::optional<Box>
  {
    const auto found = m_Marks.find(tag);
    if (found == m_Marks.end())
    {
      return std::nullopt;
    }
    return m_Kept[found->second];
  }

  void Boxes::Insert(std::uint64_t tag, const Box &box)
  {
    const auto found = m_Marks.find(tag);
    if (found != m_Marks.end())
    {
      m_Boughs->Dropped(found->second);
      m_Kept[found->second] = box;
      m_Boughs->Added(found->second, box);
      return;
    }
    std::size_t mark = m_Tags.size();
    if (m_Spare.empty())
    {
      m_Tags.push_back(tag);
      m_Kept.push_back(box);
    }
    else
    {
      mark = m_Spare.back();
      m_Spare.pop_back();
      m_Tags[mark] = tag;
      m_Kept[mark] = box;
    }
    m_Marks.emplace(tag, mark);
    m_Boughs->Added(mark, box);
  }

  void Boxes::Erase(std::uint64_t tag)
  {
    const auto found = m_Marks.find(tag);
    if (found == m_Marks.end())
    {
      return;
    }
    m_Boughs->Dropped(found->second);
    m_Spare.push_back(found->second);
    m_Marks.erase(found);
  }

  auto Boxes::Extents() const -> Box
  {
    return m_Boughs->Whole();
  }

  void Boxes::Overlapping(const Box &box, const Visit &visit) const
  {
    m_Boughs->Each([&box](const Box &bounds) -> bool { return Overlaps(bounds, box); },
                   [this, &box, &visit](std::size_t mark) -> bool
                   { return !Overlaps(m_Kept[mark], box) || visit(m_Tags[mark]); });
  }

  auto Boxes::Overlapping(const Box &box) const -> std::vector<std::uint64_t>
  {
    std::vector<std::uint64_t> found;
    Overlapping(box,
                [&found](std::uint64_t tag) -> bool
                {
                  found.push_back(tag);
                  return true;
                });
    std::ranges::sort(found);
    return found;
  }

  void Boxes::Within(const Box &box, const Visit &visit) const
  {
    m_Boughs->Each([&box](const Box &bounds) -> bool { return Overlaps(bounds, box); },
                   [this, &box, &visit](std::size_t mark) -> bool
                   {
                     const Box &kept = m_Kept[mark];
                     const bool inside =
                         !IsEmpty(kept) && Contains(box, kept.Min) && Contains(box, kept.Max);
                     return !inside || visit(m_Tags[mark]);
                   });
  }

  auto Boxes::Within(const Box &box) const -> std::vector<std::uint64_t>
  {
    std::vector<std::uint64_t> found;
    Within(box,
           [&found](std::uint64_t tag) -> bool
           {
             found.push_back(tag);
             return true;
           });
    std::ranges::sort(found);
    return found;
  }

  void Boxes::Nearest(Point point, const Measure &distance, const Walk &walk) const
  {
    const auto far = [this, &distance](std::size_t mark) -> double
    {
      return distance(m_Tags[mark]);
    };
    const auto tag = [this](std::size_t mark) -> std::uint64_t
    {
      return m_Tags[mark];
    };
    m_Boughs->Outward(
        point,
        Boughs::Ranking<std::remove_cvref_t<decltype(far)>, std::remove_cvref_t<decltype(tag)>>{
            .Far = far, .Key = tag},
        [&walk](const Boughs::Nearer &one) -> bool
        { return walk({.Tag = one.Order, .Distance = one.Distance}); });
  }

  auto Boxes::Nearest(Point point, const Seeking &seeking, const Measure &distance) const
      -> std::vector<Close>
  {
    std::vector<Close> nearest;
    if (seeking.Count == 0)
    {
      return nearest;
    }
    Nearest(point, distance,
            [&nearest, &seeking](const Close &close) -> bool
            {
              if (!(close.Distance <= seeking.Within))
              {
                return false;
              }
              nearest.push_back(close);
              return nearest.size() < seeking.Count;
            });
    return nearest;
  }
} // namespace Geoxy
