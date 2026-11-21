/***************************************************************************************************
 *  @file      source/stitched.cpp
 *  @brief     What is read from an arrangement's faces: the loops that bound the faces kept, a
 *             cycle as a loop, the cycle about a point, and whether a point is on an edge.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "arranged.hpp"
#include "scratch.hpp"
#include "winding.hpp"

#include "geoxy/angle.hpp"
#include "geoxy/arc.hpp"
#include "geoxy/box.hpp"
#include "geoxy/curve.hpp"
#include "geoxy/domain.hpp"
#include "geoxy/line.hpp"
#include "geoxy/loop.hpp"
#include "geoxy/nearest.hpp"
#include "geoxy/point.hpp"
#include "geoxy/region.hpp"
#include "geoxy/reverse.hpp"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <span>
#include <variant>
#include <vector>

namespace Geoxy
{
  namespace
  {
    /**
     *  The vertex a half ends at.
     *
     *  @param  arranged The arrangement.
     *  @param  half     The half.
     *  @return          The vertex.
     */
    auto Head(const Arranged &arranged, std::size_t half) -> std::size_t
    {
      return arranged.Edges[half / 2].Ends.at(1 - (half % 2));
    }

    /**
     *  Whether one half runs on from another along the same curve. Two halves that follow each
     *  other in a cycle and are pieces of one curve run the same way, the one starting where the
     *  other ends, or where the pieces between them were dropped, their ends merged: a cycle turns
     *  back along an edge only where the edge has a loose end, and none is left.
     *
     *  @param  arranged The arrangement.
     *  @param  before   The half before.
     *  @param  after    The half after.
     *  @return          Whether it does.
     */
    auto Follows(const Arranged &arranged, std::size_t before, std::size_t after) -> bool
    {
      return arranged.Edges[before / 2].Source == arranged.Edges[after / 2].Source;
    }

    /**
     *  The curve of a run of halves along one curve, from the first's start to the last's end.
     *
     *  @param  arranged The arrangement.
     *  @param  first    The first half.
     *  @param  last     The last half.
     *  @return          A line between their vertices, or the arc's stretch, across its seam where
     *                   it runs on past it, turned about where the halves run back.
     */
    auto RunCurve(const Arranged &arranged, std::size_t first, std::size_t last) -> Curve
    {
      const Arranged::Edge &edge   = arranged.Edges[first / 2];
      const Curve          &source = arranged.Sources[edge.Source];
      const auto           *arc    = std::get_if<Arc>(&source);
      if (arc == nullptr)
      {
        return Line{.Start = arranged.Points[edge.Ends.at(first % 2)],
                    .End   = arranged.Points[Head(arranged, last)]};
      }
      const bool  back   = first % 2 == 1;
      const Range lower  = arranged.Edges[(back ? last : first) / 2].Stretch;
      const Range upper  = arranged.Edges[(back ? first : last) / 2].Stretch;
      double      end    = upper.End;
      end               += end > lower.Start ? 0.0 : 1.0;
      const Angle start  = arc->Start + (arc->Sweep * (back ? end : lower.Start));
      const Angle sweep  = arc->Sweep * (end - lower.Start);
      return Arc{.Center = arc->Center,
                 .Radius = arc->Radius,
                 .Start  = start,
                 .Sweep  = back ? -sweep : sweep};
    }

    /**
     *  The lines of a half along a run of lines added to a loop, in the order the half runs along
     *  them, each between the points it passes.
     *
     *  @param arranged The arrangement.
     *  @param half     The half.
     *  @param loop     The loop.
     */
    void Laced(const Arranged &arranged, std::size_t half, Loop &loop)
    {
      const std::size_t lines = arranged.Edges[half / 2].Lines;
      for (std::size_t line = 0; line < lines; ++line)
      {
        const bool        forward = half % 2 == 0;
        const std::size_t from    = forward ? line : lines - line;
        const std::size_t onto    = forward ? line + 1 : lines - line - 1;
        loop.Edges.emplace_back(Line{.Start = Beaded(arranged, half / 2, from),
                                     .End   = Beaded(arranged, half / 2, onto)});
      }
    }

    /**
     *  How many times a cycle winds about a point.
     *
     *  @param  arranged The arrangement.
     *  @param  cycle    The cycle.
     *  @param  point    The point.
     *  @return          The winding number.
     */
    auto Wound(const Arranged &arranged, std::size_t cycle, Point point) -> std::int64_t
    {
      std::int64_t count = 0;
      Toured(arranged, arranged.Cycles[cycle].First,
             [&](std::size_t half) -> void
             {
               const Arranged::Edge &edge    = arranged.Edges[half / 2];
               const std::int64_t    passes  = Ascents({.Shape = &arranged.Shapes[half / 2],
                                                        .From  = arranged.Points[edge.Ends[0]],
                                                        .To    = arranged.Points[edge.Ends[1]]},
                                                       point);
               count                        += half % 2 == 0 ? passes : -passes;
             });
      return count;
    }

    /**
     *  The half a boundary runs on into at a half's end: the first clockwise from its twin about
     *  that vertex that bounds what is kept. One is always found before the twin comes round again,
     *  as the face just clockwise of the twin is the half's own, kept, and the twin's is not.
     *
     *  @tparam BoundsT  The kind of the test.
     *  @param  arranged The traced arrangement.
     *  @param  half     The half.
     *  @param  bounds   Whether a half bounds what is kept.
     *  @return          The half it runs on into.
     */
    template <typename BoundsT>
    auto Onward(const Arranged &arranged, std::size_t half, const BoundsT &bounds) -> std::size_t
    {
      std::size_t next = Clockwise(arranged, half ^ 1U);
      while (!bounds(next))
      {
        next = Clockwise(arranged, next);
      }
      return next;
    }
  } // namespace

  auto Sewn(const Arranged &arranged, std::span<std::size_t> halves) -> Loop
  {
    const std::size_t count = halves.size();
    std::size_t       start = 0;
    while (start < count && Follows(arranged, halves[(start + count - 1) % count], halves[start]))
    {
      ++start;
    }
    if (start == count)
    {
      const Arranged::Edge &edge = arranged.Edges[halves.front() / 2];
      if (const auto *arc = std::get_if<Arc>(&arranged.Sources[edge.Source]))
      {
        if (halves.front() % 2 == 0)
        {
          return {.Edges = {*arc}};
        }
        return {.Edges = {Reversed(*arc)}};
      }
      /**
       *  A line's pieces that come round to where they began, the end of the last merged into the
       *  start of the first across a tolerance wider than they are apart: no run of them is one
       *  line, so each is its own.
       */
      Loop loop = {};
      loop.Edges.reserve(count);
      for (const std::size_t half : halves)
      {
        loop.Edges.push_back(RunCurve(arranged, half, half));
      }
      return loop;
    }
    std::ranges::rotate(halves, halves.begin() + static_cast<std::ptrdiff_t>(start));
    std::size_t runs = 0;
    for (std::size_t index = 0; index < count; ++index)
    {
      const Arranged::Edge &edge = arranged.Edges[halves[index] / 2];
      const bool            onto =
          index + 1 < count && edge.Source == arranged.Edges[halves[index + 1] / 2].Source;
      runs += static_cast<std::size_t>(!onto) + edge.Lines - 1;
    }
    Loop loop = {};
    loop.Edges.reserve(runs);
    std::size_t first = 0;
    for (std::size_t index = 0; index < count; ++index)
    {
      if (index + 1 < count && Follows(arranged, halves[index], halves[index + 1]))
      {
        continue;
      }
      if (arranged.Edges[halves[index] / 2].Lines > 1)
      {
        Laced(arranged, halves[index], loop);
      }
      else
      {
        loop.Edges.push_back(RunCurve(arranged, halves[first], halves[index]));
      }
      first = index + 1;
    }
    return loop;
  }

  auto Stitched(const Arranged &arranged, const std::vector<bool> &kept) -> Region
  {
    const std::size_t  halves = arranged.Next.size();
    List<std::uint8_t> keeps(kept.size());
    for (std::size_t cycle = 0; cycle < kept.size(); ++cycle)
    {
      keeps[cycle] = static_cast<std::uint8_t>(kept[cycle]);
    }
    /**
     *  Each half's state: 0 where it bounds nothing kept, 1 where it does, 2 once it is stitched.
     */
    List<std::uint8_t> states(halves);
    for (std::size_t half = 0; half < halves; ++half)
    {
      states[half] = static_cast<std::uint8_t>(keeps[arranged.Faces[half]]
                                               & (keeps[arranged.Faces[half ^ 1U]] ^ 1U));
    }
    const auto bounding = [&states](std::size_t half) -> bool
    {
      return states[half] != 0;
    };
    List<std::size_t> run;
    Region            region = {};
    for (std::size_t first = 0; first < halves; ++first)
    {
      if (states[first] != 1)
      {
        continue;
      }
      run.assign(1, first);
      states[first] = 2;
      for (std::size_t half = Onward(arranged, first, bounding); half != first;
           half             = Onward(arranged, half, bounding))
      {
        states[half] = 2;
        run.push_back(half);
      }
      region.Loops.push_back(Sewn(arranged, run));
    }
    return region;
  }

  auto Ringed(const Arranged &arranged, std::size_t cycle) -> Loop
  {
    List<std::size_t> halves;
    Toured(arranged, arranged.Cycles[cycle].First,
           [&halves](std::size_t half) -> void { halves.push_back(half); });
    return Sewn(arranged, halves);
  }

  auto Enclosing(const Arranged &arranged, Point point, std::size_t skipped)
      -> std::optional<std::size_t>
  {
    std::optional<std::size_t> best;
    for (std::size_t index = 0; index < arranged.Cycles.size(); ++index)
    {
      const Arranged::Cycle &cycle = arranged.Cycles[index];
      const bool             fit   = cycle.Area > 0.0
                                  && cycle.Component != skipped
                                  && Contains(cycle.Bounds, point)
                                  && (!best || cycle.Area < arranged.Cycles[*best].Area);
      if (fit && Wound(arranged, index, point) != 0)
      {
        best = index;
      }
    }
    return best;
  }

  auto Grazes(const Arranged &arranged, Point point) -> bool
  {
    const double length = arranged.Within.Length;
    for (std::size_t index = 0; index < arranged.Edges.size(); ++index)
    {
      if (Contains(Padded(arranged.Boxes[index], length), point)
          && Distance(Nearest(arranged.Shapes[index], point).At, point) <= length)
      {
        return true;
      }
    }
    return false;
  }
} // namespace Geoxy
