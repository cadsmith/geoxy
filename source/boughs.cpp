/***************************************************************************************************
 *  @file      source/boughs.cpp
 *  @brief     Boxes kept in a tree of boxes: packed all at once, one put in, one taken out.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "boughs.hpp"

#include "geoxy/box.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <optional>
#include <span>
#include <utility>
#include <vector>

namespace Geoxy
{
  namespace
  {
    /**
     *  How much of its size a box is kept larger by: four steps of a double.
     */
    constexpr double Pad = 0x1p-50;

    /**
     *  Whether a box can be kept in a limb: not empty, and finite.
     *
     *  @param  box The box.
     *  @return     Whether it can.
     */
    auto Bounded(const Box &box) -> bool
    {
      return !IsEmpty(box) && std::isfinite((box.Max.X - box.Min.X) + (box.Max.Y - box.Min.Y));
    }

  } // namespace

  auto Boughs::Padded(const Box &box) -> Box
  {
    const double size = std::max(
        {std::abs(box.Min.X), std::abs(box.Min.Y), std::abs(box.Max.X), std::abs(box.Max.Y)});
    const double pad = size * Pad;
    return {.Min = {.X = box.Min.X - pad, .Y = box.Min.Y - pad},
            .Max = {.X = box.Max.X + pad, .Y = box.Max.Y + pad}};
  }

  namespace
  {
    /**
     *  The middle of a box, across or upward.
     *
     *  @param  box    The box.
     *  @param  upward Whether upward.
     *  @return        Its middle that way.
     */
    auto Middle(const Box &box, bool upward) -> double
    {
      return upward ? (0.5 * box.Min.Y) + (0.5 * box.Max.Y) : (0.5 * box.Min.X) + (0.5 * box.Max.X);
    }

    /**
     *  The area of a box.
     *
     *  @param  box The box.
     *  @return     Its area.
     */
    auto Girth(const Box &box) -> double
    {
      return (box.Max.X - box.Min.X) * (box.Max.Y - box.Min.Y);
    }

    /**
     *  The box that holds all a limb holds.
     *
     *  @param  limb The limb.
     *  @return      Its box.
     */
    auto Bounds(const Boughs::Bough &limb) -> Box
    {
      Box all;
      for (std::size_t slot = 0; slot < limb.Count; ++slot)
      {
        all = Joined(all, limb.Boxes.at(slot));
      }
      return all;
    }

    /**
     *  The limb within a limb that a box goes in: the one whose box grows least to hold it, and of
     *  those the smallest.
     *
     *  @param  limb The limb, not a leaf.
     *  @param  box  The box.
     *  @return      Where in the limb.
     */
    auto Chosen(const Boughs::Bough &limb, const Box &box) -> std::size_t
    {
      std::size_t               best  = 0;
      std::pair<double, double> least = {std::numeric_limits<double>::infinity(),
                                         std::numeric_limits<double>::infinity()};
      for (std::size_t slot = 0; slot < limb.Count; ++slot)
      {
        const double                    girth = Girth(limb.Boxes.at(slot));
        const std::pair<double, double> cost  = {Girth(Joined(limb.Boxes.at(slot), box)) - girth,
                                                 girth};
        if (cost < least)
        {
          least = cost;
          best  = slot;
        }
      }
      return best;
    }
  } // namespace

  auto Boughs::Fresh(bool leaf) -> std::size_t
  {
    if (Spare.empty())
    {
      Limbs.push_back({.Leaf = leaf});
      return Limbs.size() - 1;
    }
    const std::size_t limb = Spare.back();
    Spare.pop_back();
    Limbs[limb] = Bough{.Leaf = leaf};
    return limb;
  }

  void Boughs::Housed(std::size_t limb, const Entry &entry)
  {
    Bough &into               = Limbs[limb];
    into.Boxes.at(into.Count) = entry.Bounds;
    into.Marks.at(into.Count) = entry.Mark;
    ++into.Count;
    if (into.Leaf)
    {
      Homes[entry.Mark] = limb;
    }
    else
    {
      Limbs[entry.Mark].Parent = limb;
    }
  }

  auto Boughs::Put(std::size_t limb, const Entry &entry) -> std::optional<Entry>
  {
    if (Limbs[limb].Count < Most)
    {
      Housed(limb, entry);
      return std::nullopt;
    }
    std::vector<Entry> all   = {entry};
    Box                whole = entry.Bounds;
    for (std::size_t slot = 0; slot < Most; ++slot)
    {
      all.push_back({.Bounds = Limbs[limb].Boxes.at(slot), .Mark = Limbs[limb].Marks.at(slot)});
      whole = Joined(whole, all.back().Bounds);
    }
    const bool upward = whole.Max.Y - whole.Min.Y > whole.Max.X - whole.Min.X;
    std::ranges::sort(all, {},
                      [upward](const Entry &one) -> double { return Middle(one.Bounds, upward); });
    const std::size_t twin = Fresh(Limbs[limb].Leaf);
    Limbs[limb].Count      = 0;
    for (std::size_t index = 0; index < all.size(); ++index)
    {
      Housed(index <= Most / 2 ? limb : twin, all[index]);
    }
    return Entry{.Bounds = Bounds(Limbs[twin]), .Mark = twin};
  }

  auto Boughs::Slot(std::size_t holder, std::size_t held) const -> std::size_t
  {
    const Bough &within = Limbs[holder];
    std::size_t  slot   = 0;
    while (within.Marks.at(slot) != held)
    {
      ++slot;
    }
    return slot;
  }

  void Boughs::Taken(std::size_t holder, std::size_t held)
  {
    const std::size_t slot = Slot(holder, held);
    Bough            &from = Limbs[holder];
    --from.Count;
    from.Boxes.at(slot) = from.Boxes.at(from.Count);
    from.Marks.at(slot) = from.Marks.at(from.Count);
  }

  auto Boughs::Tier(std::vector<Entry> row, bool leaf) -> std::vector<Entry>
  {
    const std::size_t limbs = (row.size() + Most - 1) / Most;
    const auto strips = static_cast<std::size_t>(std::ceil(std::sqrt(static_cast<double>(limbs))));
    std::ranges::sort(row, {},
                      [](const Entry &one) -> double { return Middle(one.Bounds, false); });
    std::vector<Entry>     above;
    const std::span<Entry> all = row;
    for (std::size_t start = 0; start < all.size(); start += strips * Most)
    {
      const std::span<Entry> strip =
          all.subspan(start, std::min(strips * Most, all.size() - start));
      std::ranges::sort(strip, {},
                        [](const Entry &one) -> double { return Middle(one.Bounds, true); });
      for (std::size_t index = 0; index < strip.size(); ++index)
      {
        if (index % Most == 0)
        {
          above.push_back({.Bounds = {}, .Mark = Fresh(leaf)});
        }
        Housed(above.back().Mark, strip[index]);
      }
    }
    for (Entry &one : above)
    {
      one.Bounds = Bounds(Limbs[one.Mark]);
    }
    return above;
  }

  auto Boughs::Packed(std::span<const Box> boxes) -> Boughs
  {
    Boughs made;
    made.Homes.assign(boxes.size(), None);
    std::vector<Entry> row;
    for (std::size_t mark = 0; mark < boxes.size(); ++mark)
    {
      if (Bounded(boxes[mark]))
      {
        row.push_back({.Bounds = boxes[mark], .Mark = mark});
      }
      else
      {
        made.Loose.push_back(mark);
        made.Homes[mark] = Adrift;
      }
    }
    if (row.size() <= Most)
    {
      for (const Entry &one : row)
      {
        made.Housed(made.Root, one);
      }
      return made;
    }
    made.Limbs.clear();
    row = made.Tier(std::move(row), true);
    while (row.size() > 1)
    {
      row = made.Tier(std::move(row), false);
    }
    made.Root = row.front().Mark;
    return made;
  }

  void Boughs::Added(std::size_t mark, const Box &box)
  {
    Homes.resize(std::max(Homes.size(), mark + 1), None);
    if (!Bounded(box))
    {
      Loose.insert(std::ranges::lower_bound(Loose, mark), mark);
      Homes[mark] = Adrift;
      return;
    }
    std::size_t limb = Root;
    while (!Limbs[limb].Leaf)
    {
      limb = Limbs[limb].Marks.at(Chosen(Limbs[limb], box));
    }
    std::optional<Entry> spill = Put(limb, {.Bounds = box, .Mark = mark});
    while (limb != Root)
    {
      const std::size_t parent                   = Limbs[limb].Parent;
      Limbs[parent].Boxes.at(Slot(parent, limb)) = Bounds(Limbs[limb]);
      spill                                      = spill ? Put(parent, *spill) : std::nullopt;
      limb                                       = parent;
    }
    if (spill)
    {
      const std::size_t top = Fresh(false);
      Housed(top, {.Bounds = Bounds(Limbs[Root]), .Mark = Root});
      Housed(top, *spill);
      Root = top;
    }
  }

  auto Boughs::Whole() const -> Box
  {
    return Bounds(Limbs[Root]);
  }

  void Boughs::Dropped(std::size_t mark)
  {
    const std::size_t home = Homes[mark];
    Homes[mark]            = None;
    if (home == Adrift)
    {
      Loose.erase(std::ranges::lower_bound(Loose, mark));
      return;
    }
    Taken(home, mark);
    std::size_t limb = home;
    while (limb != Root)
    {
      const std::size_t parent = Limbs[limb].Parent;
      if (Limbs[limb].Count == 0)
      {
        Taken(parent, limb);
        Spare.push_back(limb);
      }
      else
      {
        Limbs[parent].Boxes.at(Slot(parent, limb)) = Bounds(Limbs[limb]);
      }
      limb = parent;
    }
    while (!Limbs[Root].Leaf && Limbs[Root].Count == 1)
    {
      Spare.push_back(Root);
      Root               = Limbs[Root].Marks.at(0);
      Limbs[Root].Parent = None;
    }
    if (Limbs[Root].Count == 0)
    {
      Limbs[Root].Leaf = true;
    }
  }
} // namespace Geoxy
