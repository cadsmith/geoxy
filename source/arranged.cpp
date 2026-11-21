/***************************************************************************************************
 *  @file      source/arranged.cpp
 *  @brief     Lines and arcs cut where they meet: every crossing and shared stretch found, cuts
 *             carried across shared stretches, vertices merged in a stated order, and the pieces of
 *             one stretch made one edge.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "arranged.hpp"

#include "sane.hpp"
#include "scratch.hpp"
#include "unions.hpp"
#include "weave/find.hpp"
#include "weave/found.hpp"
#include "weave/layout.hpp"
#include "weave/merge.hpp"
#include "weave/order.hpp"

#include "geoxy/angle.hpp"
#include "geoxy/arc.hpp"
#include "geoxy/box.hpp"
#include "geoxy/curve.hpp"
#include "geoxy/domain.hpp"
#include "geoxy/extents.hpp"
#include "geoxy/line.hpp"
#include "geoxy/tolerance.hpp"

#include <algorithm>
#include <array>
#include <cstddef>
#include <memory>
#include <optional>
#include <ranges>
#include <span>
#include <utility>
#include <variant>
#include <vector>

namespace Geoxy
{
  namespace
  {
    /**
     *  What an arrangement's vertices and edges are made from: the findings, the cuts in order and
     *  each vertex's number, and each curve's tally.
     */
    struct Making
    {
      const Weave::Findings *Found   = nullptr; /**< The findings.                       */
      Weave::Numbering       Numbers = {};      /**< The cuts and the vertices' numbers. */
      std::span<const Tally> Tallies = {};      /**< Each curve's tally.                 */
    };

    /**
     *  An arrangement's vertices, each the first point numbered so, and its edges: the pieces of
     *  the curves, those of one shared stretch made one, their tallies summed along the first.
     *
     *  @param making   What they are made from.
     *  @param arranged The arrangement, its curves laid.
     */
    void Built(const Making &making, Arranged &arranged)
    {
      const Weave::Findings   &found   = *making.Found;
      const List<std::size_t> &ids     = *making.Numbers.Ids;
      const auto               tallies = making.Tallies;
      arranged.Points.reserve(found.Points.size());
      for (std::size_t index = 0; index < found.Points.size(); ++index)
      {
        if (ids[index] == arranged.Points.size())
        {
          arranged.Points.push_back(found.Points[index]);
        }
      }
      arranged.Edges.reserve(found.Points.size());
      if (found.Commons.empty())
      {
        Weave::Hewn(found, making.Numbers,
                    [&arranged, &tallies](const Weave::Piece &piece) -> void
                    {
                      arranged.Edges.push_back({.Source  = piece.Source,
                                                .Stretch = piece.Stretch,
                                                .Ends    = piece.Ends,
                                                .Weighs  = tallies[piece.Source],
                                                .Lines   = piece.Lines});
                    });
        return;
      }
      List<Weave::Piece> pieces;
      Weave::Hewn(found, making.Numbers,
                  [&pieces](const Weave::Piece &piece) -> void { pieces.push_back(piece); });
      Unions                   same = Weave::Coincide(found, pieces);
      std::vector<std::size_t> edges(pieces.size(), pieces.size());
      for (std::size_t index = 0; index < pieces.size(); ++index)
      {
        const Weave::Piece &piece = pieces[index];
        const std::size_t   root  = same.Root(index);
        if (root == index)
        {
          edges[index] = arranged.Edges.size();
          arranged.Edges.push_back({.Source  = piece.Source,
                                    .Stretch = piece.Stretch,
                                    .Ends    = piece.Ends,
                                    .Lines   = piece.Lines});
        }
        Arranged::Edge &edge  = arranged.Edges[edges[root]];
        const int       along = piece.Ends == pieces[root].Ends ? 1 : -1;
        const Tally     tally = tallies[piece.Source];
        edge.Weighs = {edge.Weighs[0] + (along * tally[0]), edge.Weighs[1] + (along * tally[1])};
      }
    }
  } // namespace

  auto Carved(const Arranged &arranged, std::size_t edge) -> Curve
  {
    const Arranged::Edge &piece = arranged.Edges[edge];
    if (const auto *arc = std::get_if<Arc>(&arranged.Sources[piece.Source]))
    {
      return Arc{.Center = arc->Center,
                 .Radius = arc->Radius,
                 .Start  = arc->Start + (arc->Sweep * piece.Stretch.Start),
                 .Sweep  = arc->Sweep * (piece.Stretch.End - piece.Stretch.Start)};
    }
    return Line{.Start = arranged.Points[piece.Ends[0]], .End = arranged.Points[piece.Ends[1]]};
  }

  auto Arrange(const Inputs &inputs, const Tolerance &tolerance) -> Arranged
  {
    Weave::Layout layout{Sane(tolerance).Length, inputs.Memory};
    layout.Reserve(inputs);
    for (std::size_t run = 0; run < inputs.Runs.size(); ++run)
    {
      layout.Begin();
      for (const Curve &curve : inputs.Runs[run])
      {
        layout.Add(curve, inputs.Tallies[run]);
      }
    }
    const std::array<double, 2> lengths = layout.Lengths();
    auto [arranged, tallies]            = layout.Laid();
    const double length                 = lengths[0];
    arranged.Within                     = {.Length = length, .Radians = Sane(tolerance).Radians};
    arranged.Chained                    = inputs.Chained;
    if (arranged.Fault)
    {
      return std::move(arranged);
    }
    Weave::Findings found =
        Weave::Find(arranged.Sources, tolerance,
                    {.Reach = lengths[1], .Length = length, .Chained = inputs.Chained});
    arranged.Fault = Weave::Spending.at(std::min<std::size_t>(found.Spent, 1));
    Weave::Carry(found);
    Unions               vertices{found.Points.size()};
    const Weave::Merging merged = Weave::Merge(found, length, vertices);
    const bool kept = inputs.Untwisted || (inputs.Spare && Weave::Twistable(arranged.Sources));
    const Weave::Sequence sequence = Weave::Sequenced(found, vertices, kept);
    const Weave::Untwist  untwist =
        !kept ? Weave::Untwist{}
              : Weave::Untwisted(
                    found, vertices,
                    {.Sources = arranged.Sources, .Moved = merged.Moved, .Joins = sequence.Joins});
    const std::span<const std::size_t> leaders = inputs.Untwisted
                                                   ? std::span<const std::size_t>(untwist.Leaders)
                                                   : std::span<const std::size_t>{};
    arranged.Twisted                           = inputs.Untwisted && untwist.Twisted;
    const List<std::size_t> ids                = sequence.Joined || !leaders.empty()
                                                   ? Weave::Numbered(found, vertices, leaders)
                                                   : merged.Ids;
    const bool              solitary           = Weave::Solitary(inputs, found, arranged);
    if (inputs.Spare && untwist.Twisted)
    {
      /**
       *  The split arrangement, from the same cuts: each vertex numbered by its leader where its
       *  cluster is twisted, as Untwisted would number it.
       */
      auto split = std::make_shared<Arranged>(inputs.Memory);
      split->Sources.assign(arranged.Sources.begin(), arranged.Sources.end());
      split->Within                = arranged.Within;
      split->Chained               = arranged.Chained;
      split->Fault                 = arranged.Fault;
      split->Twisted               = true;
      const List<std::size_t> leds = Weave::Numbered(found, vertices, untwist.Leaders);
      Built({.Found = &found, .Numbers = {.Order = &sequence, .Ids = &leds}, .Tallies = tallies},
            *split);
      if (solitary)
      {
        split->Runs = arranged.Runs;
        split->Tallies.assign(tallies.begin(), tallies.end());
        split->Lone.resize(split->Points.size());
        Weave::Lonesome(found, leds, split->Lone);
      }
      arranged.Split = std::move(split);
    }
    Built({.Found = &found, .Numbers = {.Order = &sequence, .Ids = &ids}, .Tallies = tallies},
          arranged);
    if (solitary)
    {
      arranged.Tallies.assign(tallies.begin(), tallies.end());
      arranged.Lone.resize(arranged.Points.size());
      Weave::Lonesome(found, ids, arranged.Lone);
    }
    return std::move(arranged);
  }

  void Outfitted(Arranged &arranged)
  {
    arranged.Shapes.clear();
    arranged.Boxes.clear();
    arranged.Shapes.reserve(arranged.Edges.size());
    arranged.Boxes.reserve(arranged.Edges.size());
    for (std::size_t index = 0; index < arranged.Edges.size(); ++index)
    {
      arranged.Shapes.push_back(Carved(arranged, index));
      arranged.Boxes.push_back(Extents(arranged.Shapes.back()));
    }
    for (Arranged::Cycle &cycle : arranged.Cycles)
    {
      cycle.Bounds = {};
      Toured(arranged, cycle.First, [&arranged, &cycle](std::size_t half) -> void
             { cycle.Bounds = Joined(cycle.Bounds, arranged.Boxes[half / 2]); });
    }
  }

  auto Frayed(const Arranged &arranged) -> std::vector<std::size_t>
  {
    std::vector<std::size_t> degrees(arranged.Points.size(), 0);
    for (const Arranged::Edge &edge : arranged.Edges)
    {
      ++degrees[edge.Ends[0]];
      ++degrees[edge.Ends[1]];
    }
    std::vector<std::size_t> loose;
    for (std::size_t index = 0; index < degrees.size(); ++index)
    {
      if (degrees[index] == 1)
      {
        loose.push_back(index);
      }
    }
    return loose;
  }

  void Prune(Arranged &arranged)
  {
    std::vector<std::vector<std::size_t>> reaching(arranged.Points.size());
    for (std::size_t index = 0; index < arranged.Edges.size(); ++index)
    {
      reaching[arranged.Edges[index].Ends[0]].push_back(index);
      reaching[arranged.Edges[index].Ends[1]].push_back(index);
    }
    std::vector<bool>        gone(arranged.Edges.size(), false);
    std::vector<std::size_t> degrees(arranged.Points.size());
    std::vector<std::size_t> loose;
    for (std::size_t index = 0; index < reaching.size(); ++index)
    {
      degrees[index] = reaching[index].size();
      if (degrees[index] == 1)
      {
        loose.push_back(index);
      }
    }
    while (!loose.empty())
    {
      const std::size_t vertex = loose.back();
      loose.pop_back();
      for (const std::size_t index : reaching[vertex])
      {
        if (gone[index])
        {
          continue;
        }
        gone[index]                            = true;
        const std::array<std::size_t, 2> ends  = arranged.Edges[index].Ends;
        const std::size_t                other = ends[0] == vertex ? ends[1] : ends[0];
        --degrees[vertex];
        if (--degrees[other] == 1)
        {
          loose.push_back(other);
        }
      }
    }
    std::pmr::vector<Arranged::Edge> kept(arranged.Edges.get_allocator());
    for (std::size_t index = 0; index < arranged.Edges.size(); ++index)
    {
      if (!gone[index])
      {
        kept.push_back(arranged.Edges[index]);
      }
    }
    arranged.Edges = std::move(kept);
    arranged.Shapes.clear();
    arranged.Boxes.clear();
  }
} // namespace Geoxy
