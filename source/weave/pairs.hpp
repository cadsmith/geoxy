/***************************************************************************************************
 *  @file      source/weave/pairs.hpp
 *  @brief     The pairs of curves whose boxes come near, found without asking every pair.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#ifndef GEOXY_WEAVE_PAIRS_HPP
#define GEOXY_WEAVE_PAIRS_HPP

#include "real.hpp"
#include "scratch.hpp"
#include "weave/found.hpp"

#include "geoxy/angle.hpp"
#include "geoxy/arc.hpp"
#include "geoxy/box.hpp"
#include "geoxy/curve.hpp"
#include "geoxy/line.hpp"
#include "geoxy/point.hpp"
#include "geoxy/vector.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <iterator>
#include <numeric>
#include <ranges>
#include <span>
#include <variant>
#include <vector>

namespace Geoxy::Weave
{
  /**
   *  A line's or an arc's stance among pairs, as a drawing tells it: a line level or upright where
   *  its ends share one coordinate and not the other, any other line, or an arc.
   *
   *  @param  source The line or arc.
   *  @return        0 level, 1 upright, 2 any other line, 3 an arc.
   */
  inline auto Stance(const Curve &source) -> std::uint8_t
  {
    if (const auto *line = std::get_if<Line>(&source))
    {
      const bool level   = Exactly(line->Start.Y, line->End.Y);
      const bool upright = Exactly(line->Start.X, line->End.X);
      if (level == upright)
      {
        return 2;
      }
      return level ? 0 : 1;
    }
    return 3;
  }

  /**
   *  A box holding an arc, and some steps of a double more: its ends' as At gives them, with each
   *  point of its circle on an axis through its center that its sweep passes or comes within
   *  TurnSlack of, or every one of them where its start angle is past PlainTurn or its radius not
   *  above 0. Each piece of an arc between those points runs one way along each axis, so it lies
   *  between its ends.
   *
   *  @param  source The arc.
   *  @param  ends   Its ends, as At gives them.
   *  @return        The box.
   */
  inline auto Snug(const Arc &source, const std::array<Point, 2> &ends) -> Box
  {
    const Box    box    = Grown(Grown(Box{}, ends[0]), ends[1]);
    const auto  *arc    = &source;
    const double sweep  = arc->Sweep.InDegrees();
    const double from   = arc->Start.InDegrees() + std::min(sweep, 0.0);
    const bool   plain  = std::abs(from) <= PlainTurn && arc->Radius > 0.0;
    const double radius = std::abs(arc->Radius);
    Box          held   = box;
    for (std::size_t quarter = 0; quarter < 4; ++quarter)
    {
      const double toward  = 90.0 * static_cast<double>(quarter);
      double       turn    = std::fmod(toward - from + TurnSlack, 360.0);
      turn                += turn < 0.0 ? 360.0 : 0.0;
      if (!plain || turn <= std::abs(sweep) + (2.0 * TurnSlack))
      {
        const Vector way = {.X = quarter % 2 == 0 ? 1.0 - static_cast<double>(quarter) : 0.0,
                            .Y = quarter % 2 == 1 ? 2.0 - static_cast<double>(quarter) : 0.0};
        held             = Grown(held, arc->Center + (way * radius));
      }
    }
    return held;
  }

  /**
   *  A box as the sweep reads it: its least x, its span of y, and its curve.
   */
  struct Slat
  {
    double      Least  = 0.0; /**< Its least x. */
    double      Bottom = 0.0; /**< Its least y. */
    double      Top    = 0.0; /**< Its most y.  */
    std::size_t Index  = 0;   /**< Its curve.   */
  };

  /**
   *  The pairs of boxes that overlap, each pair once, the lower index first, in no order: the boxes
   *  laid in as many columns as there are boxes by their least x, by counting, side by side so that
   *  a column is read straight through, and each box asked against those whose least x lies within
   *  its own span of x and that come after it by their least x, then their index, so that every
   *  overlap is asked once.
   *
   *  @param  boxes The boxes, two or more.
   *  @return       The pairs.
   */
  inline auto Overlapping(std::span<const Box> boxes) -> List<Couple>
  {
    const std::size_t count  = boxes.size();
    double            lowest = boxes.front().Min.X;
    double            last   = lowest;
    for (const Box &box : boxes)
    {
      lowest = std::min(lowest, box.Min.X);
      last   = std::max(last, box.Min.X);
    }
    const double scale =
        static_cast<double>(count) / std::max(last - lowest, Above(std::abs(last)));
    const auto column = [lowest, scale, count](double x) -> std::size_t
    {
      const double place = std::min((x - lowest) * scale, static_cast<double>(count - 1));
      return static_cast<std::size_t>(std::max(place, 0.0));
    };
    List<std::size_t> starts(count + 1, 0);
    for (const Box &box : boxes)
    {
      ++starts[column(box.Min.X) + 1];
    }
    std::partial_sum(starts.begin(), starts.end(), starts.begin());
    List<Slat> slats(count);
    {
      List<std::size_t> filled(starts.begin(), starts.end() - 1);
      for (std::size_t index = 0; index < count; ++index)
      {
        const Box &box                     = boxes[index];
        slats[filled[column(box.Min.X)]++] = {
            .Least = box.Min.X, .Bottom = box.Min.Y, .Top = box.Max.Y, .Index = index};
      }
    }
    List<Couple> pairs;
    pairs.reserve(2 * count);
    for (std::size_t index = 0; index < count; ++index)
    {
      const Box        &box = boxes[index];
      const std::size_t end = starts[column(box.Max.X) + 1];
      for (std::size_t slot = starts[column(box.Min.X)]; slot < end; ++slot)
      {
        const Slat &near = slats[slot];
        if (near.Bottom <= box.Max.Y
            && box.Min.Y <= near.Top
            && near.Least <= box.Max.X
            && (box.Min.X < near.Least || (Exactly(box.Min.X, near.Least) && index < near.Index)))
        {
          pairs.push_back({std::min(index, near.Index), std::max(index, near.Index)});
        }
      }
    }
    return pairs;
  }

  /**
   *  How many boxes a ply holds, at the least on average, for plies to be worth laying; and up to
   *  how many of a ply's boxes are asked one by one rather than by halving.
   */
  inline constexpr std::size_t Plies = 4;

  /**
   *  A run of curves in order whose boxes each lie, along each axis, no lower than the one before,
   *  their least and their most both, or each no higher: so that the boxes of a ply that come near
   *  any box are a run of it, found by halving, and those after a box that come near it are the
   *  first after it.
   */
  struct Ply
  {
    std::size_t         First   = 0;  /**< Its first curve.   */
    std::size_t         End     = 0;  /**< One past its last. */
    Box                 Bounds  = {}; /**< Its boxes joined.  */
    std::array<bool, 2> Falling = {}; /**< Whether its boxes come lower, along x and
                                           along y; else they come no lower.      */
  };

  /**
   *  Boxes laid in plies, each as long as its boxes keep rising or falling along both axes.
   *
   *  @param  boxes The boxes, one or more.
   *  @return       The plies, in order, every box in one.
   */
  inline auto Plied(std::span<const Box> boxes) -> List<Ply>
  {
    List<Ply>           plies;
    std::size_t         first = 0;
    std::array<bool, 4> ways  = {true, true, true, true};
    Box                 held  = boxes.front();
    for (std::size_t index = 1; index <= boxes.size(); ++index)
    {
      if (index < boxes.size())
      {
        const Box                &last  = boxes[index - 1];
        const Box                &next  = boxes[index];
        const std::array<bool, 4> still = {
            ways[0] && last.Min.X <= next.Min.X && last.Max.X <= next.Max.X,
            ways[1] && next.Min.X <= last.Min.X && next.Max.X <= last.Max.X,
            ways[2] && last.Min.Y <= next.Min.Y && last.Max.Y <= next.Max.Y,
            ways[3] && next.Min.Y <= last.Min.Y && next.Max.Y <= last.Max.Y};
        if ((still[0] || still[1]) && (still[2] || still[3]))
        {
          ways = still;
          held = Joined(held, next);
          continue;
        }
      }
      plies.push_back(
          {.First = first, .End = index, .Bounds = held, .Falling = {!ways[0], !ways[2]}});
      if (index < boxes.size())
      {
        first = index;
        ways  = {true, true, true, true};
        held  = boxes[index];
      }
    }
    return plies;
  }

  /**
   *  The run of a ply's boxes that overlap a box: along each axis those of a ply that rises whose
   *  most is not below the box's least come after those whose most is, and those whose least is not
   *  above the box's most before those whose least is; and the other way round in a ply that falls.
   *  Each is found by halving, and the boxes in both runs of both axes are the ones.
   *
   *  @param  boxes The boxes.
   *  @param  ply   The ply, or a run of one.
   *  @param  box   The box.
   *  @return       The run, as a ply.
   */
  inline auto Covering(std::span<const Box> boxes, const Ply &ply, const Box &box) -> Ply
  {
    Ply run = ply;
    for (std::size_t axis = 0; axis < 2; ++axis)
    {
      const auto least = [axis](const Box &one) -> double
      {
        return axis == 0 ? one.Min.X : one.Min.Y;
      };
      const auto most = [axis](const Box &one) -> double
      {
        return axis == 0 ? one.Max.X : one.Max.Y;
      };
      const auto first = [&boxes, &run](const auto &before) -> std::size_t
      {
        std::size_t low  = run.First;
        std::size_t high = run.End;
        while (low < high)
        {
          const std::size_t middle = low + ((high - low) / 2);
          if (before(boxes[middle]))
          {
            low = middle + 1;
          }
          else
          {
            high = middle;
          }
        }
        return low;
      };
      const double floor   = least(box);
      const double ceiling = most(box);
      std::size_t  from    = 0;
      std::size_t  until   = 0;
      if (ply.Falling.at(axis))
      {
        from  = first([&least, ceiling](const Box &one) -> bool { return least(one) > ceiling; });
        until = first([&most, floor](const Box &one) -> bool { return most(one) >= floor; });
      }
      else
      {
        from  = first([&most, floor](const Box &one) -> bool { return most(one) < floor; });
        until = first([&least, ceiling](const Box &one) -> bool { return least(one) <= ceiling; });
      }
      run.First = std::max(run.First, from);
      run.End   = std::max(run.First, std::min(run.End, until));
    }
    return run;
  }

  /**
   *  The run of a ply's boxes that overlap a box, as Covering finds it, where the ply's last box
   *  overlaps it: those that overlap it are a run of the ply, which then holds its last, and so are
   *  the boxes back from its last while they overlap.
   *
   *  @param  boxes The boxes.
   *  @param  ply   The ply, its last box overlapping the box.
   *  @param  box   The box.
   *  @return       The run, as a ply.
   */
  inline auto Trailing(std::span<const Box> boxes, const Ply &ply, const Box &box) -> Ply
  {
    const std::span<const Box> before = boxes.subspan(ply.First, ply.End - 1 - ply.First);
    const auto                 stop =
        std::ranges::find_if_not(before | std::views::reverse,
                                 [&box](const Box &one) -> bool { return Overlaps(one, box); });
    return {.First =
                ply.First + static_cast<std::size_t>(std::distance(before.begin(), stop.base())),
            .End     = ply.End,
            .Bounds  = ply.Bounds,
            .Falling = ply.Falling};
  }

  /**
   *  The run of a ply's boxes that overlap a box, as Covering finds it, where the ply's first box
   *  overlaps it: the boxes on from its first while they overlap.
   *
   *  @param  boxes The boxes.
   *  @param  ply   The ply, its first box overlapping the box.
   *  @param  box   The box.
   *  @return       The run, as a ply.
   */
  inline auto Leading(std::span<const Box> boxes, const Ply &ply, const Box &box) -> Ply
  {
    const std::span<const Box> after = boxes.subspan(ply.First + 1, ply.End - 1 - ply.First);
    const auto                 stop = std::ranges::find_if_not(after, [&box](const Box &one) -> bool
                                                               { return Overlaps(one, box); });
    return {.First   = ply.First,
            .End     = ply.First + 1 + static_cast<std::size_t>(std::distance(after.begin(), stop)),
            .Bounds  = ply.Bounds,
            .Falling = ply.Falling};
  }

  /**
   *  The pairs of a ply's boxes that overlap, each box against those after it while they overlap
   *  it: those that do are the first after it, as along each axis the boxes past one that does not
   *  lie farther still.
   *
   *  @tparam SeenT The kind of what is done with a pair.
   *  @param  boxes The boxes.
   *  @param  ply   The ply.
   *  @param  seen  What is done with each pair, its two indices given in order.
   */
  template <typename SeenT>
  void Lengthwise(std::span<const Box> boxes, const Ply &ply, const SeenT &seen)
  {
    for (std::size_t index = ply.First; index < ply.End; ++index)
    {
      for (std::size_t other = index + 1; other < ply.End && Overlaps(boxes[index], boxes[other]);
           ++other)
      {
        seen(index, other);
      }
    }
  }

  /**
   *  The pairs of two plies' boxes that overlap: the boxes of the one that come near the other's
   *  bounds, each against the run of the other's that overlap it, or each of the other's where it
   *  has few.
   *
   *  @tparam SeenT The kind of what is done with a pair.
   *  @param  boxes The boxes.
   *  @param  plies The two plies.
   *  @param  seen  What is done with each pair, its two indices given in order.
   */
  template <typename SeenT>
  void Crosswise(std::span<const Box> boxes, const std::array<const Ply *, 2> &plies,
                 const SeenT &seen)
  {
    /**
     *  Two plies one after the other whose boxes where they join overlap: the last of the first
     *  overlaps the second's bounds, and the first of the second the first's.
     */
    const bool joined = plies[0]->End == plies[1]->First
                     && Overlaps(boxes[plies[0]->End - 1], boxes[plies[1]->First]);
    const Ply  one    = joined ? Trailing(boxes, *plies[0], plies[1]->Bounds)
                               : Covering(boxes, *plies[0], plies[1]->Bounds);
    const Ply  other  = joined ? Leading(boxes, *plies[1], plies[0]->Bounds)
                               : Covering(boxes, *plies[1], plies[0]->Bounds);
    const bool fewer  = one.End - one.First <= other.End - other.First;
    const Ply &outer  = fewer ? one : other;
    const Ply &inner  = fewer ? other : one;
    for (std::size_t index = outer.First; index < outer.End; ++index)
    {
      const Ply near =
          inner.End - inner.First > Plies ? Covering(boxes, inner, boxes[index]) : inner;
      for (std::size_t across = near.First; across < near.End; ++across)
      {
        if (Overlaps(boxes[index], boxes[across]))
        {
          seen(std::min(index, across), std::max(index, across));
        }
      }
    }
  }

  /**
   *  Every pair of boxes that overlap, each pair once, the lower index first, in no order: within
   *  each ply as Lengthwise finds them, and of each two plies whose bounds overlap, as Overlapping
   *  finds them, as Crosswise finds them. Where the plies are short, the boxes themselves are swept
   *  by Overlapping.
   *
   *  @tparam SeenT The kind of what is done with a pair.
   *  @param  boxes The boxes, two or more.
   *  @param  seen  What is done with each pair, its two indices given in order.
   */
  template <typename SeenT>
  void Woven(std::span<const Box> boxes, const SeenT &seen)
  {
    const List<Ply> plies = Plied(boxes);
    if (plies.size() * Plies > boxes.size())
    {
      for (const Couple &pair : Overlapping(boxes))
      {
        seen(pair[0], pair[1]);
      }
      return;
    }
    List<Box> bounds;
    bounds.reserve(plies.size());
    for (const Ply &ply : plies)
    {
      Lengthwise(boxes, ply, seen);
      bounds.push_back(ply.Bounds);
    }
    if (plies.size() < 2)
    {
      return;
    }
    for (const Couple &meeting : Overlapping(bounds))
    {
      Crosswise(boxes, {&plies[meeting[0]], &plies[meeting[1]]}, seen);
    }
  }

  /**
   *  Pairs of curves in the order IntersectAll meets the pairs that meet: by kind, then by the
   *  first index, then the second; sorted where they are few, else put in order by counting twice:
   *  each pair's first dealt by its second, then each second, in order, dealt by its first, so that
   *  each first's seconds come in order. A sort of each first's seconds would guess at every
   *  comparison.
   *
   *  @param  found   The pairs, the lower index first, in no order.
   *  @param  stances Each curve's stance.
   *  @return         The pairs of each kind, in order.
   */
  inline auto Coupled(std::span<const Couple> found, std::span<const std::uint8_t> stances)
      -> std::array<List<Couple>, Kinds>
  {
    std::array<List<Couple>, Kinds> kinds = {};
    if (found.size() * Plies < stances.size())
    {
      List<Couple> sorted(found.begin(), found.end());
      std::ranges::sort(sorted);
      for (const Couple &pair : sorted)
      {
        kinds.at(Couplings.at(stances[pair[0]]).at(stances[pair[1]])).push_back(pair);
      }
      return kinds;
    }
    const std::size_t count = stances.size();
    List<std::size_t> starts(count + 1, 0);
    List<std::size_t> later(count + 1, 0);
    for (const Couple &pair : found)
    {
      ++starts[pair[0] + 1];
      ++later[pair[1] + 1];
    }
    std::partial_sum(starts.begin(), starts.end(), starts.begin());
    std::partial_sum(later.begin(), later.end(), later.begin());
    List<std::size_t> firsts(found.size());
    for (const Couple &pair : found)
    {
      firsts[later[pair[1]]++] = pair[0];
    }
    List<std::size_t> seconds(found.size());
    {
      List<std::size_t> filled(starts.begin(), starts.end() - 1);
      std::size_t       from = 0;
      for (std::size_t other = 0; other < count; ++other)
      {
        for (; from < later[other]; ++from)
        {
          seconds[filled[firsts[from]]++] = other;
        }
      }
    }
    for (std::size_t index = 0; index < count; ++index)
    {
      const std::span<const std::size_t> near = std::span<const std::size_t>(seconds).subspan(
          starts[index], starts[index + 1] - starts[index]);
      for (const std::size_t other : near)
      {
        kinds.at(Couplings.at(stances[index]).at(stances[other])).push_back({index, other});
      }
    }
    return kinds;
  }
} // namespace Geoxy::Weave

#endif /* GEOXY_WEAVE_PAIRS_HPP */
