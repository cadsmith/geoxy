/***************************************************************************************************
 *  @file      source/boughs.hpp
 *  @brief     Boxes kept in a tree of boxes, each limb's box holding its own: what a Tree keeps its
 *             curves in and a Cloud its points, to find those near a box or a point without asking
 *             every one.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#ifndef GEOXY_BOUGHS_HPP
#define GEOXY_BOUGHS_HPP

#include "geoxy/box.hpp"
#include "geoxy/point.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <limits>
#include <optional>
#include <queue>
#include <span>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

namespace Geoxy
{
  /**
   *  Boxes by their marks, in an R-tree: every leaf holds up to Most boxes with their marks, every
   *  other limb up to Most limbs with the boxes that hold them, and all leaves are as deep.
   *
   *  A box is kept a little larger than given, by a few steps of a double at its size, so that a
   *  point worked out on what it bounds is never outside it by rounding. A box that is empty or not
   *  finite is kept apart, adrift, and told with every answer: nothing can rule it out.
   */
  struct Boughs
  {
    /**
     *  How many boxes a limb holds at most.
     */
    static constexpr std::size_t Most = 16;

    /**
     *  No limb: the root's parent, and the home of a mark that is not kept.
     */
    static constexpr std::size_t None = std::numeric_limits<std::size_t>::max();

    /**
     *  The home of a mark kept adrift.
     */
    static constexpr std::size_t Adrift = None - 1;

    /**
     *  How far short of its rounded distance a box is told to be, so that no rounding puts it
     *  farther than what it holds.
     */
    static constexpr double Shy = 1.0 - 0x1p-50;

    /**
     *  A limb: its boxes, and what each bounds, a mark in a leaf and a limb elsewhere.
     */
    struct Bough
    {
      std::array<Box, Most>         Boxes  = {};   /**< Its boxes.                              */
      std::array<std::size_t, Most> Marks  = {};   /**< What each bounds: a mark, or a limb.    */
      std::size_t                   Count  = 0;    /**< How many it holds.                      */
      std::size_t                   Parent = None; /**< The limb that holds it; none, the root. */
      bool                          Leaf   = true; /**< Whether it holds marks.                 */
    };

    /**
     *  A box and what it bounds.
     */
    struct Entry
    {
      Box         Bounds = {}; /**< The box.                     */
      std::size_t Mark   = 0;  /**< What it bounds, as a limb's. */
    };

    /**
     *  A mark and how far it is.
     */
    struct Nearer
    {
      double        Distance = 0.0; /**< How far.                               */
      std::size_t   Mark     = 0;   /**< Which.                                 */
      std::uint64_t Order    = 0;   /**< Its place among those as far: its key. */

      /**
       *  What two are ordered by: the distance, then the key.
       *
       *  @param  one One.
       *  @return     Its distance and its key.
       */
      [[nodiscard]] static constexpr auto Rank(const Nearer &one)
          -> std::pair<double, std::uint64_t>
      {
        return {one.Distance, one.Order};
      }
    };

    /**
     *  How many marks are asked for at most, and how far at most.
     */
    struct Limit
    {
      std::size_t Count  = 0;   /**< How many at most. */
      double      Within = 0.0; /**< How far at most.  */
    };

    /**
     *  A limb or a mark yet to be looked at, and the least distance anything in it can be at.
     */
    /**
     *  How the marks nearest a point are ranked: how far each is, and its key among those as far.
     *
     *  @tparam FarT How far a mark is.
     *  @tparam KeyT A mark's key.
     */
    template <typename FarT, typename KeyT>
    struct Ranking
    {
      FarT Far; /**< How far a mark is. */
      KeyT Key; /**< A mark's key.      */
    };

    struct Ahead
    {
      double      Bound = 0.0;   /**< The least distance.    */
      std::size_t Mark  = 0;     /**< The limb, or the mark. */
      bool        Limb  = false; /**< Whether it is a limb.  */

      /**
       *  Whether one is to be looked at after another: the farther later.
       *
       *  @param  lhs One.
       *  @param  rhs The other.
       *  @return     Whether the first is farther.
       */
      [[nodiscard]] friend constexpr auto operator>(const Ahead &lhs, const Ahead &rhs) -> bool
      {
        return lhs.Bound > rhs.Bound;
      }
    };

    std::vector<Bough>       Limbs = {Bough{}}; /**< The limbs, the spare ones among them.   */
    std::vector<std::size_t> Spare = {};        /**< The limbs given back, to be used again. */
    std::vector<std::size_t> Homes = {};        /**< Each mark's leaf; Adrift, or None.      */
    std::vector<std::size_t> Loose = {};        /**< The marks adrift, in order.             */
    std::size_t              Root  = 0;         /**< The root.                               */

    /**
     *  Boxes packed all at once, each marked by its place among them: sorted across, cut in strips,
     *  each strip sorted upward and cut in leaves, and so again for the limbs above.
     *
     *  @param  boxes The boxes.
     *  @return       Them, kept.
     */
    [[nodiscard]] static auto Packed(std::span<const Box> boxes) -> Boughs;

    /**
     *  A box padded by 2^-50 of its largest coordinate's magnitude on every side, as every query
     *  sees the boxes kept, so that no rounding of what is asked about leaves one out.
     *
     *  @param  box The box, kept exact.
     *  @return     It padded.
     */
    [[nodiscard]] static auto Padded(const Box &box) -> Box;

    /**
     *  The exact box of every box kept with finite sides: the root's.
     *
     *  @return The box; empty where none is kept.
     */
    [[nodiscard]] auto Whole() const -> Box;

    /**
     *  A box more.
     *
     *  @pre        The mark is not kept.
     *  @param mark Its mark.
     *  @param box  The box.
     */
    void Added(std::size_t mark, const Box &box);

    /**
     *  A box no more.
     *
     *  @pre        The mark is kept.
     *  @param mark Its mark.
     */
    void Dropped(std::size_t mark);

    /**
     *  A limb to fill: a spare one, or one more.
     *
     *  @param  leaf Whether it is to hold marks.
     *  @return      The limb, empty.
     */
    [[nodiscard]] auto Fresh(bool leaf) -> std::size_t;

    /**
     *  An entry put last in a limb that has room, and told where it is.
     *
     *  @param limb  The limb.
     *  @param entry The entry.
     */
    void Housed(std::size_t limb, const Entry &entry);

    /**
     *  An entry put in a limb: last where it has room; where it has none, the limb's entries and
     *  this one sorted along the longer side of their box and shared between it and a new limb.
     *
     *  @param  limb  The limb.
     *  @param  entry The entry.
     *  @return       The new limb with its box, where one was made.
     */
    [[nodiscard]] auto Put(std::size_t limb, const Entry &entry) -> std::optional<Entry>;

    /**
     *  Where a limb keeps what it bounds.
     *
     *  @pre           The limb holds it.
     *  @param  holder The limb.
     *  @param  held   What it bounds.
     *  @return        Its place in the limb.
     */
    [[nodiscard]] auto Slot(std::size_t holder, std::size_t held) const -> std::size_t;

    /**
     *  An entry taken out of a limb, its last put in its place.
     *
     *  @pre          The limb holds it.
     *  @param holder The limb.
     *  @param held   What the entry bounds.
     */
    void Taken(std::size_t holder, std::size_t held);

    /**
     *  A row of entries packed in limbs, and the row of those limbs.
     *
     *  @param  row  The entries.
     *  @param  leaf Whether they are marks.
     *  @return      The limbs made, each with its box.
     */
    [[nodiscard]] auto Tier(std::vector<Entry> row, bool leaf) -> std::vector<Entry>;

    /**
     *  How far a box is from a point at least, told a little short.
     *
     *  @param  box   The box.
     *  @param  point The point.
     *  @return       No more than the distance from the point to anything in the box.
     */
    [[nodiscard]] static auto Apart(const Box &box, Point point) -> double
    {
      const double across = std::max({box.Min.X - point.X, 0.0, point.X - box.Max.X});
      const double upward = std::max({box.Min.Y - point.Y, 0.0, point.Y - box.Max.Y});
      return std::hypot(across, upward) * Shy;
    }

    /**
     *  Every mark whose box is reached, and every mark adrift, in no order.
     *
     *  @tparam ReachesT Whether a box is reached: a limb's, then those within it.
     *  @tparam SinkT    What takes each mark.
     *  @param  reaches  Tells whether a box is reached.
     *  @param  sink     Takes each mark.
     */
    template <typename ReachesT, typename SinkT>
    void Each(const ReachesT &reaches, const SinkT &sink) const
    {
      const auto told = [&sink](std::size_t mark) -> bool
      {
        if constexpr (std::is_same_v<decltype(sink(mark)), bool>)
        {
          return sink(mark);
        }
        else
        {
          sink(mark);
          return true;
        }
      };
      for (const std::size_t mark : Loose)
      {
        if (!told(mark))
        {
          return;
        }
      }
      std::vector<std::size_t> walk = {Root};
      while (!walk.empty())
      {
        const Bough &limb = Limbs[walk.back()];
        walk.pop_back();
        for (std::size_t slot = 0; slot < limb.Count; ++slot)
        {
          if (!reaches(Padded(limb.Boxes.at(slot))))
          {
            continue;
          }
          if (limb.Leaf)
          {
            if (!told(limb.Marks.at(slot)))
            {
              return;
            }
          }
          else
          {
            walk.push_back(limb.Marks.at(slot));
          }
        }
      }
    }

    /**
     *  Every mark in order of distance from a point, those as far in the order of their keys, told
     *  one at a time until the walk says stop: limbs taken nearest box first, and a limb no farther
     *  than a mark before it, so that the order is exact without looking at every mark.
     *
     *  @tparam FarT    How far a mark is from the point.
     *  @tparam KeyT    A mark's key.
     *  @tparam WalkT   What is told of each mark, in order: false to stop.
     *  @param  point   The point.
     *  @param  ranking Tells how far a mark is, never less than its box, a mark it tells no number
     *                  for left out, and its key.
     *  @param  walk    Told each mark and its distance, nearest first.
     */
    template <typename FarT, typename KeyT, typename WalkT>
    void Outward(Point point, const Ranking<FarT, KeyT> &ranking, const WalkT &walk) const
    {
      using Step = std::tuple<double, bool, std::uint64_t, std::size_t>;
      std::priority_queue<Step, std::vector<Step>, std::greater<>> ahead;
      const auto offer = [&ahead, &ranking](std::size_t mark) -> void
      {
        const double distance = ranking.Far(mark);
        if (distance >= 0.0)
        {
          ahead.emplace(distance, true, ranking.Key(mark), mark);
        }
      };
      for (const std::size_t mark : Loose)
      {
        offer(mark);
      }
      ahead.emplace(0.0, false, 0, Root);
      while (!ahead.empty())
      {
        const auto [bound, marked, key, which] = ahead.top();
        ahead.pop();
        if (marked)
        {
          if (!walk(Nearer{.Distance = bound, .Mark = which, .Order = key}))
          {
            return;
          }
          continue;
        }
        const Bough &limb = Limbs[which];
        for (std::size_t slot = 0; slot < limb.Count; ++slot)
        {
          if (limb.Leaf)
          {
            offer(limb.Marks.at(slot));
          }
          else
          {
            ahead.emplace(Apart(Padded(limb.Boxes.at(slot)), point), false, 0, limb.Marks.at(slot));
          }
        }
      }
    }

    /**
     *  The marks nearest a point, the nearest first: limbs and marks looked at nearest box first,
     *  until none left can be nearer than the farthest kept.
     *
     *  @tparam FarT  How far a mark is from the point.
     *  @param  point The point.
     *  @param  most  How many at most, and how far at most.
     *  @param  far   Tells how far a mark is; a mark it tells no number for is left out.
     *  @return       The marks, by distance, then by mark.
     */
    template <typename FarT>
    [[nodiscard]] auto Closest(Point point, const Limit &most, const FarT &far) const
        -> std::vector<Nearer>
    {
      const auto mark = [](std::size_t index) -> std::uint64_t
      {
        return index;
      };
      return Closest(point, most,
                     Ranking<FarT, std::remove_cvref_t<decltype(mark)>>{.Far = far, .Key = mark});
    }

    /**
     *  The marks nearest a point, the nearest first, those as near in the order of a key: as the
     *  other Closest, ranked by distance and then by the key.
     *
     *  @tparam FarT    How far a mark is from the point.
     *  @tparam KeyT    A mark's key.
     *  @param  point   The point.
     *  @param  most    How many at most, and how far at most.
     *  @param  ranking Tells how far a mark is, a mark it tells no number for left out, and its
     *                  key.
     *  @return         The marks, by distance, then by key.
     */
    template <typename FarT, typename KeyT>
    [[nodiscard]] auto Closest(Point point, const Limit &most,
                               const Ranking<FarT, KeyT> &ranking) const -> std::vector<Nearer>
    {
      const FarT         &far = ranking.Far;
      const KeyT         &key = ranking.Key;
      std::vector<Nearer> best;
      if (most.Count == 0)
      {
        return best;
      }
      const auto offer = [&best, &most, &far, &key](std::size_t mark) -> void
      {
        const Nearer one = {.Distance = far(mark), .Mark = mark, .Order = key(mark)};
        if (!(one.Distance <= most.Within))
        {
          return;
        }
        best.push_back(one);
        std::ranges::push_heap(best, {}, Nearer::Rank);
        if (best.size() > most.Count)
        {
          std::ranges::pop_heap(best, {}, Nearer::Rank);
          best.pop_back();
        }
      };
      for (const std::size_t mark : Loose)
      {
        offer(mark);
      }
      std::priority_queue<Ahead, std::vector<Ahead>, std::greater<>> ahead;
      ahead.push({.Bound = 0.0, .Mark = Root, .Limb = true});
      while (!ahead.empty())
      {
        const Ahead next = ahead.top();
        ahead.pop();
        if (next.Bound > most.Within
            || (best.size() == most.Count && next.Bound > best.front().Distance))
        {
          break;
        }
        if (!next.Limb)
        {
          offer(next.Mark);
          continue;
        }
        const Bough &limb = Limbs[next.Mark];
        for (std::size_t slot = 0; slot < limb.Count; ++slot)
        {
          ahead.push({.Bound = Apart(Padded(limb.Boxes.at(slot)), point),
                      .Mark  = limb.Marks.at(slot),
                      .Limb  = !limb.Leaf});
        }
      }
      std::ranges::sort(best, {}, Nearer::Rank);
      return best;
    }
  };
} // namespace Geoxy

#endif /* GEOXY_BOUGHS_HPP */
