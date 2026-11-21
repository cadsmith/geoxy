/***************************************************************************************************
 *  @file      source/weave/order.cpp
 *  @brief     The pieces between a curve's notches, in order along it, and the pieces of one
 *             stretch made one edge.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "weave/order.hpp"

#include "arranged.hpp"
#include "real.hpp"
#include "scratch.hpp"
#include "unions.hpp"
#include "weave/found.hpp"
#include "weave/merge.hpp"

#include "geoxy/curve.hpp"
#include "geoxy/domain.hpp"
#include "geoxy/intersect.hpp"
#include "geoxy/line.hpp"
#include "geoxy/point.hpp"

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <map>
#include <numeric>
#include <ranges>
#include <span>
#include <utility>
#include <variant>
#include <vector>

namespace Geoxy::Weave
{
  namespace
  {
    /**
     *  Whether one cut comes before another along their curve: by place, then by vertex.
     *
     *  @param  lhs One cut.
     *  @param  rhs The other.
     *  @return     Whether it does.
     */
    auto Sooner(const Notch &lhs, const Notch &rhs) -> bool
    {
      return lhs.Along < rhs.Along || (Exactly(lhs.Along, rhs.Along) && lhs.Vertex < rhs.Vertex);
    }

    /**
     *  How many cuts of a curve, at the least, are dealt into bins rather than sorted whole; a bin
     *  of more than half as many is sorted whole too.
     */
    constexpr std::size_t Bins = 64;

    /**
     *  Room to put a curve's cuts in order: the cuts as they are dealt, and where each bin of them
     *  starts.
     */
    struct Dealt
    {
      List<Notch>       Cuts   = {}; /**< The cuts, dealt bin by bin. */
      List<std::size_t> Starts = {}; /**< Where each bin starts.      */
    };

    /**
     *  A curve's many cuts in the order Sooner puts them, as one sort of them all would leave them:
     *  dealt by counting into as many bins as there are cuts, by their place along the curve times
     *  their count, and each bin's few put in order one by one, its many sorted. A product by a
     *  positive number, rounded, and its floor never fall as the place rises, so a cut in an
     *  earlier bin comes sooner; two at one place share a bin. A curve's crossings lie spread along
     *  it, most often one or two to a bin, where a sort's every comparison is a guess the branches
     *  miss half the time.
     *
     *  @param run   The cuts.
     *  @param dealt Room to deal them in.
     */
    void Bucketed(std::span<Notch> run, Dealt &dealt)
    {
      const auto sooner = [](const Notch &lhs, const Notch &rhs) -> bool
      {
        return Sooner(lhs, rhs);
      };
      const std::size_t count = run.size();
      const auto        scale = static_cast<double>(count);
      const auto        bin   = [scale, count](double along) -> std::size_t
      {
        const double place = std::min(along * scale, static_cast<double>(count - 1));
        return static_cast<std::size_t>(std::max(place, 0.0));
      };
      List<std::size_t> &starts = dealt.Starts;
      starts.assign(count + 1, 0);
      for (const Notch &notch : run)
      {
        ++starts[bin(notch.Along) + 1];
      }
      std::partial_sum(starts.begin(), starts.end(), starts.begin());
      dealt.Cuts.resize(count);
      for (const Notch &notch : run)
      {
        dealt.Cuts[starts[bin(notch.Along)]++] = notch;
      }
      std::size_t from = 0;
      for (std::size_t index = 0; index < count; ++index)
      {
        const std::span<Notch> held =
            std::span<Notch>(dealt.Cuts).subspan(from, starts[index] - from);
        from = starts[index];
        if (held.size() > Bins / 2)
        {
          std::ranges::sort(held, sooner);
          continue;
        }
        for (std::size_t at = 1; at < held.size(); ++at)
        {
          const Notch cut  = held[at];
          std::size_t slot = at;
          for (; slot > 0 && Sooner(cut, held[slot - 1]); --slot)
          {
            held[slot] = held[slot - 1];
          }
          held[slot] = cut;
        }
      }
      std::ranges::copy(dealt.Cuts, run.begin());
    }

    /**
     *  Whether a cluster, its atoms each given the lines through it, has two lines that meet in no
     *  atom of it.
     *
     *  @param  lines Each atom's lines, sorted, the atoms one after another.
     *  @param  ends  Where each atom's lines end among them, in order.
     *  @return       Whether it has.
     */
    auto Askew(std::span<const std::size_t> lines, std::span<const std::size_t> ends) -> bool
    {
      List<std::size_t> all(lines.begin(), lines.end());
      std::ranges::sort(all);
      const auto [last, stop] = std::ranges::unique(all);
      all.erase(last, stop);
      for (std::size_t one = 0; one < all.size(); ++one)
      {
        for (std::size_t other = one + 1; other < all.size(); ++other)
        {
          bool        met  = false;
          std::size_t from = 0;
          for (const std::size_t end : ends)
          {
            const std::span<const std::size_t> atom = lines.subspan(from, end - from);
            met                                     = met
                                                   || (std::ranges::binary_search(atom, all[one])
                                                       && std::ranges::binary_search(atom, all[other]));
            from                                    = end;
          }
          if (!met)
          {
            return true;
          }
        }
      }
      return false;
    }

    /**
     *  The members of the clusters that a vertex merged into from another point: each cluster's
     *  first, then its vertices, in order.
     *
     *  @param  found   The findings.
     *  @param  sets    The merged vertices, every merge made.
     *  @param  strands The vertices merged elsewhere.
     *  @return         Each member's cluster, as its least vertex, and the member.
     */
    auto Clustered(const Findings &found, Unions &sets, const Tangle &strands) -> List<Couple>
    {
      const std::size_t  count = found.Points.size();
      List<std::uint8_t> clustered(count, 0);
      for (const std::size_t moved : strands.Moved)
      {
        clustered[sets.Root(moved)] = 1;
      }
      List<Couple> members;
      for (std::size_t index = 0; index < count; ++index)
      {
        const std::size_t root = sets.Root(index);
        if (clustered[root] != 0)
        {
          members.push_back({root, index});
        }
      }
      std::ranges::sort(members);
      return members;
    }

    /**
     *  The lines through the clusters' members: the line an end is of, the lines cut there.
     *
     *  @param  found   The findings.
     *  @param  members The members, by cluster.
     *  @param  strands The curves.
     *  @return         Each member and a line through it, in order.
     */
    auto Raveled(const Findings &found, std::span<const Couple> members, const Tangle &strands)
        -> List<Couple>
    {
      const std::size_t count = found.Points.size();
      const auto        line  = [&strands](std::size_t curve) -> bool
      {
        return std::holds_alternative<Line>(strands.Sources[curve]);
      };
      List<std::uint8_t> member(count, 0);
      List<Couple>       lines;
      const auto         firsts = std::span<const std::size_t>(found.Firsts);
      for (const Couple &pair : members)
      {
        member[pair[1]] = 1;
        if (pair[1] < firsts.back())
        {
          const auto        after = std::ranges::upper_bound(firsts, pair[1]);
          const std::size_t curve = static_cast<std::size_t>(after - firsts.begin()) - 1;
          if (line(curve))
          {
            lines.push_back({pair[1], curve});
          }
        }
      }
      for (const Notch &notch : found.Notches)
      {
        if (member[notch.Vertex] != 0 && line(notch.Source))
        {
          lines.push_back({notch.Vertex, notch.Source});
        }
      }
      std::ranges::sort(lines);
      return lines;
    }

    /**
     *  The clusters' members joined into atoms: a member with its twin, the vertices of two cuts at
     *  one place, and members of one cluster at one point.
     *
     *  @param  found   The findings.
     *  @param  members The members, by cluster.
     *  @param  strands The vertices of cuts at one place.
     *  @return         The atoms, over the members' places.
     */
    auto Atomized(const Findings &found, std::span<const Couple> members, const Tangle &strands)
        -> Unions
    {
      const std::size_t count = found.Points.size();
      List<std::size_t> place(count, count);
      for (std::size_t index = 0; index < members.size(); ++index)
      {
        place[members[index][1]] = index;
      }
      Unions atoms{members.size()};
      for (std::size_t index = 0; index < members.size(); ++index)
      {
        const std::size_t twin = found.Twins[members[index][1]];
        if (twin < count)
        {
          atoms.Join(place[twin], index);
        }
      }
      for (const Couple &join : strands.Joins)
      {
        if (place[join[0]] < count)
        {
          atoms.Join(place[join[0]], place[join[1]]);
        }
      }
      for (std::size_t one = 1; one < members.size(); ++one)
      {
        for (std::size_t other = one; other > 0 && members[other - 1][0] == members[one][0];
             --other)
        {
          if (Twinned(found.Points[members[one][1]], found.Points[members[other - 1][1]]))
          {
            atoms.Join(one, other - 1);
          }
        }
      }
      return atoms;
    }

    /**
     *  A cluster's members, their atoms and the lines through them.
     */
    struct Cluster
    {
      std::span<const Couple> Members = {};      /**< All the members, by cluster.    */
      std::size_t             From    = 0;       /**< The cluster's first place.      */
      std::size_t             Until   = 0;       /**< One past its last.              */
      Unions                 *Atoms   = nullptr; /**< The atoms, over all the places. */
    };

    /**
     *  Whether a cluster has two lines that meet in no atom of it: each atom's lines gathered, the
     *  atoms one after another, and asked of Askew.
     *
     *  @param  cluster The cluster.
     *  @param  lines   Each member and a line through it, in order.
     *  @return         Whether it has.
     */
    auto Twisted(const Cluster &cluster, std::span<const Couple> lines) -> bool
    {
      const std::size_t from  = cluster.From;
      const std::size_t until = cluster.Until;
      List<std::size_t> gathered;
      List<std::size_t> bounds;
      for (std::size_t atom = from; atom < until; ++atom)
      {
        if (cluster.Atoms->Root(atom) != atom)
        {
          continue;
        }
        const std::size_t start = gathered.size();
        for (std::size_t index = atom; index < until; ++index)
        {
          if (cluster.Atoms->Root(index) != atom)
          {
            continue;
          }
          const auto held = std::ranges::equal_range(
              lines, Couple{cluster.Members[index][1], 0},
              [](const Couple &lhs, const Couple &rhs) -> bool { return lhs[0] < rhs[0]; });
          for (const Couple &pair : held)
          {
            gathered.push_back(pair[1]);
          }
        }
        std::sort(gathered.begin() + static_cast<std::ptrdiff_t>(start), gathered.end());
        bounds.push_back(gathered.size());
      }
      return Askew(gathered, bounds);
    }
  } // namespace

  auto Sequenced(const Findings &found, Unions &sets, bool kept) -> Sequence
  {
    Sequence sequence = {.Starts = List<std::size_t>(found.Count + 1, 0)};
    if (found.Notches.empty())
    {
      return sequence;
    }
    List<std::size_t> &starts = sequence.Starts;
    for (const Notch &notch : found.Notches)
    {
      ++starts[notch.Source + 1];
    }
    std::size_t total = 0;
    for (std::size_t curve = 0; curve < found.Count; ++curve)
    {
      const std::size_t more = starts[curve + 1];
      starts[curve]          = total;
      total += more > 0 ? more + (found.Firsts[curve + 1] - found.Firsts[curve]) : 0;
    }
    starts[found.Count] = total;
    sequence.Cuts.resize(total);
    List<std::size_t> filled(starts.begin(), starts.end() - 1);
    for (std::size_t curve = 0; curve < found.Count; ++curve)
    {
      if (starts[curve] < starts[curve + 1])
      {
        const auto [ends, count] = Ended(found, curve);
        std::copy_n(ends.begin(), count,
                    sequence.Cuts.begin() + static_cast<std::ptrdiff_t>(filled[curve]));
        filled[curve] += count;
      }
    }
    for (const Notch &notch : found.Notches)
    {
      sequence.Cuts[filled[notch.Source]++] = notch;
    }
    Dealt dealt;
    for (std::size_t curve = 0; curve < found.Count; ++curve)
    {
      const std::span<Notch> run =
          std::span<Notch>(sequence.Cuts).subspan(starts[curve], starts[curve + 1] - starts[curve]);
      if (run.size() < Bins)
      {
        std::ranges::sort(run, [](const Notch &lhs, const Notch &rhs) -> bool
                          { return Sooner(lhs, rhs); });
      }
      else
      {
        Bucketed(run, dealt);
      }
      for (std::size_t index = 1; index < run.size(); ++index)
      {
        const Notch &before = run[index - 1];
        const Notch &after  = run[index];
        if (Exactly(before.Along, after.Along))
        {
          sequence.Joined = sequence.Joined || sets.Root(before.Vertex) != sets.Root(after.Vertex);
          sets.Join(before.Vertex, after.Vertex);
          if (kept)
          {
            sequence.Joins.push_back({before.Vertex, after.Vertex});
          }
        }
      }
    }
    return sequence;
  }

  auto Numbered(const Findings &found, Unions &sets, std::span<const std::size_t> leaders)
      -> List<std::size_t>
  {
    const std::size_t count = found.Points.size();
    List<std::size_t> ids(count);
    std::size_t       next = 0;
    for (std::size_t index = 0; index < count; ++index)
    {
      const std::size_t root =
          !leaders.empty() && leaders[index] < count ? leaders[index] : sets.Root(index);
      ids[index] = root == index ? next++ : ids[root];
    }
    return ids;
  }

  auto Untwisted(const Findings &found, Unions &sets, const Tangle &strands) -> Untwist
  {
    Untwist untwist;
    if (strands.Moved.empty())
    {
      return untwist;
    }
    const std::size_t  count   = found.Points.size();
    const List<Couple> members = Clustered(found, sets, strands);
    const List<Couple> lines   = Raveled(found, members, strands);
    Unions             atoms   = Atomized(found, members, strands);
    for (std::size_t from = 0; from < members.size();)
    {
      std::size_t until = from + 1;
      while (until < members.size() && members[until][0] == members[from][0])
      {
        ++until;
      }
      const Cluster cluster = {.Members = members, .From = from, .Until = until, .Atoms = &atoms};
      if (Twisted(cluster, lines))
      {
        untwist.Twisted = true;
        untwist.Leaders.resize(count, count);
        for (std::size_t index = from; index < until; ++index)
        {
          untwist.Leaders[members[index][1]] = members[atoms.Root(index)][1];
        }
      }
      from = until;
    }
    return untwist;
  }

  auto Coincide(const Findings &found, std::span<const Piece> pieces) -> Unions
  {
    Unions                                sets{pieces.size()};
    std::vector<std::vector<std::size_t>> owned(found.Count);
    for (std::size_t index = 0; index < pieces.size(); ++index)
    {
      owned[pieces[index].Source].push_back(index);
    }
    const auto inside = [&pieces](std::size_t piece, Range range) -> bool
    {
      const Range  stretch = pieces[piece].Stretch;
      const double middle  = (stretch.Start / 2.0) + (stretch.End / 2.0);
      return range.Start <= middle && middle <= range.End;
    };
    const auto key = [&pieces](std::size_t piece) -> std::pair<std::size_t, std::size_t>
    {
      const auto ends = pieces[piece].Ends;
      return std::minmax(ends[0], ends[1]);
    };
    for (const Common &common : found.Commons)
    {
      std::multimap<std::pair<std::size_t, std::size_t>, std::size_t> seconds;
      for (const std::size_t piece : owned[common.Second])
      {
        if (inside(piece, common.Shared.Second))
        {
          seconds.emplace(key(piece), piece);
        }
      }
      for (const std::size_t piece : owned[common.First])
      {
        if (!inside(piece, common.Shared.First))
        {
          continue;
        }
        const Range  stretch = pieces[piece].Stretch;
        const double middle = Reflected(common, (stretch.Start / 2.0) + (stretch.End / 2.0), false);
        const auto [first, last] = seconds.equal_range(key(piece));
        for (auto match = first; match != last; ++match)
        {
          const Range other = pieces[match->second].Stretch;
          if (other.Start <= middle && middle <= other.End)
          {
            sets.Join(piece, match->second);
          }
        }
      }
    }
    return sets;
  }

  auto Twistable(std::span<const Curve> sources) -> bool
  {
    return std::ranges::count_if(sources, [](const Curve &source) -> bool
                                 { return std::holds_alternative<Line>(source); })
         > 1;
  }

  void Lonesome(const Findings &found, std::span<const std::size_t> ids,
                std::pmr::vector<std::array<std::uint32_t, 2>> &lone)
  {
    const std::size_t                  points = found.Points.size();
    const std::size_t                  ends   = found.Firsts[found.Count];
    List<std::array<std::uint32_t, 2>> lines(points - ends, {Arranged::Crowded, Arranged::Crowded});
    for (const Notch &notch : found.Notches)
    {
      std::array<std::uint32_t, 2> &pair            = lines[notch.Vertex - ends];
      pair.at(pair[0] == Arranged::Crowded ? 0 : 1) = static_cast<std::uint32_t>(notch.Source);
    }
    List<std::uint8_t> crowds(lone.size(), 0);
    for (std::size_t index = 0; index < points; ++index)
    {
      const std::size_t vertex = ids[index];
      if (crowds[vertex] == 0)
      {
        lone[vertex] =
            index < ends ? std::array{Arranged::Joint, Arranged::Joint} : lines[index - ends];
      }
      else if (index >= ends)
      {
        lone[vertex] = {Arranged::Crowded, Arranged::Crowded};
      }
      crowds[vertex] = static_cast<std::uint8_t>(std::min(crowds[vertex] + 1, 3));
    }
    for (std::size_t vertex = 0; vertex < lone.size(); ++vertex)
    {
      const std::size_t wanted = lone[vertex][0] == Arranged::Joint ? 2 : 1;
      if (crowds[vertex] != wanted)
      {
        lone[vertex] = {Arranged::Crowded, Arranged::Crowded};
      }
    }
  }

  auto Solitary(const Inputs &inputs, const Findings &found, const Arranged &arranged) -> bool
  {
    return inputs.Chained
        && found.Commons.empty()
        && std::ranges::all_of(arranged.Sources, [](const Curve &source) -> bool
                               { return std::holds_alternative<Line>(source); });
  }
} // namespace Geoxy::Weave
