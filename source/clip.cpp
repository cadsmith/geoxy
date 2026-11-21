/***************************************************************************************************
 *  @file      source/clip.cpp
 *  @brief     Curves cut to a box: lines by their slabs, every other curve by where it crosses the
 *             box's edges, each stretch between kept by its middle.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "geoxy/clip.hpp"

#include "eased.hpp"
#include "elliptic.hpp"
#include "kernel/piece.hpp"
#include "real.hpp"
#include "scratch.hpp"
#include "straight.hpp"
#include "trig.hpp"

#include "geoxy/arc.hpp"
#include "geoxy/at.hpp"
#include "geoxy/bezier/cubic.hpp"
#include "geoxy/bezier/quadratic.hpp"
#include "geoxy/box.hpp"
#include "geoxy/curve.hpp"
#include "geoxy/domain.hpp"
#include "geoxy/ellipse.hpp"
#include "geoxy/line.hpp"
#include "geoxy/point.hpp"
#include "geoxy/polyline.hpp"
#include "geoxy/ray.hpp"
#include "geoxy/segment.hpp"
#include "geoxy/spline.hpp"
#include "geoxy/xline.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <optional>
#include <span>
#include <utility>
#include <variant>
#include <vector>

namespace Geoxy
{
  namespace
  {
    using Kernel::Axis;
    using Kernel::Coordinate;
    using Kernel::Level;

    /**
     *  The lines along a box's four edges.
     *
     *  @param  box The box.
     *  @return     Its least and greatest x, then its least and greatest y.
     */
    auto Edges(const Box &box) -> std::array<Level, 4>
    {
      return {
          Level{.Along = Axis::X, .Value = box.Min.X},
          Level{.Along = Axis::X, .Value = box.Max.X},
          Level{.Along = Axis::Y, .Value = box.Min.Y}, Level{.Along = Axis::Y, .Value = box.Max.Y}};
    }

    /**
     *  The stretch of a straight run within a box, by Liang and Barsky's slabs.
     *
     *  @param  run The run.
     *  @param  box The box.
     *  @return     At most one stretch, longer than nothing.
     */
    auto Clipped(const Run &run, const Box &box) -> std::vector<Range>
    {
      if (IsEmpty(box))
      {
        return {};
      }
      double enter = run.Along.Start;
      double leave = run.Along.End;
      for (const Axis axis : {Axis::X, Axis::Y})
      {
        const double start = Coordinate(run.Start, axis);
        const double way   = Coordinate(run.Way, axis);
        const double low   = Coordinate(box.Min, axis);
        const double high  = Coordinate(box.Max, axis);
        if (!Zero(way))
        {
          const double first  = (low - start) / way;
          const double second = (high - start) / way;
          enter               = std::max(enter, std::min(first, second));
          leave               = std::min(leave, std::max(first, second));
        }
        else if (!(low <= start) || !(start <= high))
        {
          return {};
        }
      }
      if (enter < leave)
      {
        return {Range{.Start = enter, .End = leave}};
      }
      return {};
    }

    /**
     *  Where an ellipse crosses the lines of a box's edges.
     *
     *  Along an axis the ellipse is a cos(angle) + b sin(angle) + c, which meets a value v at the
     *  phase of a,b, plus or less acos((v - c) / hypot(a, b)), once in each turn it sweeps.
     *
     *  @param  ellipse The ellipse.
     *  @param  box     The box.
     *  @return         The parameters, from 0 to 1, unordered; none where it only touches a line.
     */
    auto Cuts(const Ellipse &ellipse, const Box &box) -> std::vector<double>
    {
      std::vector<double> cuts;
      for (const Level edge : Edges(box))
      {
        const double cosine = Coordinate(ellipse.Map.X, edge.Along);
        const double sine   = Coordinate(ellipse.Map.Y, edge.Along);
        const double offset = edge.Value - Coordinate(ellipse.Map.Origin, edge.Along);
        const double reach  = Trig::Hypot(cosine, sine);
        if (!(std::abs(offset) < reach))
        {
          continue;
        }
        const double phase  = Trig::Atan2(sine, cosine);
        const double spread = Trig::Acos(offset / reach);
        for (const double angle : {phase - spread, phase + spread})
        {
          const List<double> passes = Passes(ellipse, angle);
          cuts.insert(cuts.end(), passes.begin(), passes.end());
        }
      }
      return cuts;
    }

    /**
     *  Where a piece crosses the lines of a box's edges.
     *
     *  @param  piece The piece.
     *  @param  box   The box.
     *  @return       The piece's parameters, from 0 to 1, unordered.
     */
    auto Cuts(const Kernel::Piece &piece, const Box &box) -> std::vector<double>
    {
      std::vector<double> cuts;
      for (const Level edge : Edges(box))
      {
        const List<double> roots = Kernel::Roots(piece, edge, Kernel::Budget).Roots;
        cuts.insert(cuts.end(), roots.begin(), roots.end());
      }
      return cuts;
    }

    /**
     *  Where a polyline's segment crosses the lines of a box's edges.
     *
     *  @param  segment The segment.
     *  @param  box     The box.
     *  @return         The segment's parameters, from 0 to 1, unordered.
     */
    auto Cuts(const Segment &segment, const Box &box) -> std::vector<double>
    {
      if (const Line *line = std::get_if<Line>(&segment))
      {
        return Cuts(Kernel::Straight(line->Start, line->End), box);
      }
      return Cuts(Elliptic(std::get<Arc>(segment)), box);
    }

    /**
     *  A curve's crossings within its domain, in order, its domain's ends among them.
     *
     *  @param  domain The curve's domain.
     *  @param  cuts   Its crossings, any of them, unordered; those outside its domain are let go.
     *  @return        The crossings, each once.
     */
    auto Sorted(Range domain, std::vector<double> cuts) -> std::vector<double>
    {
      std::erase_if(cuts, [domain](double cut) -> bool
                    { return !(domain.Start < cut) || !(cut < domain.End); });
      cuts.push_back(domain.Start);
      cuts.push_back(domain.End);
      Ordered(cuts);
      return cuts;
    }

    /**
     *  The stretches between crossings that are kept, those that meet joined.
     *
     *  @param  cuts The crossings, in order.
     *  @param  kept Whether each stretch between two of them is kept.
     *  @return      The stretches kept, in order.
     */
    auto Welded(const std::vector<double> &cuts, const std::vector<bool> &kept)
        -> std::vector<Range>
    {
      std::vector<Range> ranges;
      for (std::size_t index = 0; index < kept.size(); ++index)
      {
        if (!kept[index])
        {
          continue;
        }
        if (!ranges.empty() && Exactly(ranges.back().End, cuts[index]))
        {
          ranges.back().End = cuts[index + 1];
        }
        else
        {
          ranges.push_back({.Start = cuts[index], .End = cuts[index + 1]});
        }
      }
      return ranges;
    }

    /**
     *  The stretches of a curve within a box, between its crossings of the box's edges: each
     *  stretch kept when its middle is within, and kept stretches that meet joined.
     *
     *  @tparam CurveT The curve's kind.
     *  @param  curve  The curve.
     *  @param  cuts   Its crossings, any of them, unordered; those outside its domain are let go.
     *  @param  box    The box; an empty one holds no middle.
     *  @return        The stretches, in order.
     */
    template <typename CurveT>
    auto Kept(const CurveT &curve, std::vector<double> cuts, const Box &box) -> std::vector<Range>
    {
      const std::vector<double> sorted = Sorted(Domain(curve), std::move(cuts));
      std::vector<bool>         kept;
      for (std::size_t index = 0; index + 1 < sorted.size(); ++index)
      {
        const double middle = sorted[index] + ((sorted[index + 1] - sorted[index]) / 2.0);
        kept.push_back(Contains(box, At(curve, middle)));
      }
      return Welded(sorted, kept);
    }

    /**
     *  Where a spline's paving crosses the lines of a box's edges, and its two ends.
     *
     *  @param  pieces The paving's pieces.
     *  @param  box    The box.
     *  @return        The places, in order, each once: an end of a piece as the start of the next.
     */
    auto Cuts(std::span<const Kernel::Tread> pieces, const Box &box) -> std::vector<Kernel::Spot>
    {
      std::vector<Kernel::Spot> cuts = {Kernel::Spot{.Index = 0, .Along = 0.0},
                                        Kernel::Spot{.Index = pieces.size() - 1, .Along = 1.0}};
      for (const Level edge : Edges(box))
      {
        List<Kernel::Solution> solutions;
        solutions.reserve(pieces.size());
        for (const Kernel::Tread &piece : pieces)
        {
          solutions.push_back(Kernel::Roots(piece.Shape, edge, Kernel::Budget));
        }
        for (Kernel::Spot spot : Kernel::Chained(solutions, false))
        {
          if (Exactly(spot.Along, 1.0) && spot.Index + 1 < pieces.size())
          {
            spot = {.Index = spot.Index + 1, .Along = 0.0};
          }
          cuts.push_back(spot);
        }
      }
      const auto before = [](const Kernel::Spot &lhs, const Kernel::Spot &rhs) -> bool
      {
        return lhs.Index < rhs.Index || (lhs.Index == rhs.Index && lhs.Along < rhs.Along);
      };
      const auto same = [](const Kernel::Spot &lhs, const Kernel::Spot &rhs) -> bool
      {
        return lhs.Index == rhs.Index && Exactly(lhs.Along, rhs.Along);
      };
      std::ranges::sort(cuts, before);
      const auto repeated = std::ranges::unique(cuts, same);
      cuts.erase(repeated.begin(), repeated.end());
      return cuts;
    }

    /**
     *  A place of a paving between two of its places that are apart, as far from both as its pieces
     *  let it be told: a place a step of a double from a crossing is, to doubles, on the box's
     *  edge, and tells nothing of the stretch.
     *
     *  @param  start The first place.
     *  @param  end   The second, after the first.
     *  @return       Halfway between them where one piece has both; the middle of a piece between
     *                theirs where there is one; else halfway from the first to its piece's end, or
     *                from the second's piece's start to it, whichever is the longer way.
     */
    auto Midway(Kernel::Spot start, Kernel::Spot end) -> Kernel::Spot
    {
      if (start.Index == end.Index)
      {
        return {.Index = start.Index, .Along = start.Along + ((end.Along - start.Along) / 2.0)};
      }
      if (start.Index + 1 < end.Index)
      {
        return {.Index = start.Index + 1, .Along = 0.5};
      }
      if (1.0 - start.Along >= end.Along)
      {
        return {.Index = start.Index, .Along = start.Along + ((1.0 - start.Along) / 2.0)};
      }
      return {.Index = end.Index, .Along = end.Along / 2.0};
    }

    /**
     *  The stretches of a spline within a box, found on its paving: between the paving's crossings
     *  of the box's edges, each stretch kept when a place of the paving between them is within, and
     *  told in the spline's parameter.
     *
     *  Beside a heavy weight the spline's parameter tells places apart on the curve as one double:
     *  a stretch whose ends it tells as one is nothing it can tell and is left out, and stretches
     *  kept on either side of one left out so are joined.
     *
     *  @param  paving The paving.
     *  @param  box    The box; an empty one holds no place.
     *  @return        The stretches, in order.
     */
    auto OnPaving(const Paving &paving, const Box &box) -> std::vector<Range>
    {
      const std::vector<Kernel::Spot> cuts = Cuts(paving.Treads, box);
      std::vector<double> ends = {Homed(paving.Pieces[cuts.front().Index], cuts.front().Along)};
      std::vector<bool>   kept;
      for (std::size_t index = 0; index + 1 < cuts.size(); ++index)
      {
        const Kernel::Spot next = cuts[index + 1];
        const double       end  = Homed(paving.Pieces[next.Index], next.Along);
        if (!(ends.back() < end))
        {
          continue;
        }
        const Kernel::Spot middle = Midway(cuts[index], next);
        ends.push_back(end);
        kept.push_back(
            Contains(box, Kernel::Reached(paving.Treads[middle.Index].Shape, middle.Along)));
      }
      return Welded(ends, kept);
    }
  } // namespace

  auto Clip(const Line &line, const Box &box) -> std::vector<Range>
  {
    return Clipped(Running(line), box);
  }

  auto Clip(const Ray &ray, const Box &box) -> std::vector<Range>
  {
    return Clipped(Running(ray), box);
  }

  auto Clip(const Xline &xline, const Box &box) -> std::vector<Range>
  {
    return Clipped(Running(xline), box);
  }

  auto Clip(const Arc &arc, const Box &box) -> std::vector<Range>
  {
    return Kept(arc, Cuts(Elliptic(arc), box), box);
  }

  auto Clip(const Ellipse &ellipse, const Box &box) -> std::vector<Range>
  {
    return Kept(ellipse, Cuts(ellipse, box), box);
  }

  auto Clip(const Polyline &polyline, const Box &box) -> std::vector<Range>
  {
    const std::size_t count = Count(polyline);
    if (count == 0)
    {
      if (!polyline.Vertices.empty() && Contains(box, polyline.Vertices.front().Position))
      {
        return {Range{.Start = 0.0, .End = 0.0}};
      }
      return {};
    }
    std::vector<double> cuts;
    for (std::size_t index = 0; index < count; ++index)
    {
      for (const double cut : Cuts(Nth(polyline, index), box))
      {
        cuts.push_back(static_cast<double>(index) + cut);
      }
    }
    return Kept(polyline, cuts, box);
  }

  auto Clip(const Spline &spline, const Box &box) -> std::vector<Range>
  {
    if (!Valid(spline))
    {
      return {};
    }
    if (const std::optional<Paving> paving = Paved(spline); paving.has_value())
    {
      return OnPaving(*paving, box);
    }
    std::vector<double> cuts;
    for (const Level edge : Edges(box))
    {
      const List<double> alongs =
          Kernel::Solved(spline, [edge](const Kernel::Piece &piece) -> Kernel::Solution
                         { return Kernel::Roots(piece, edge, Kernel::Budget); });
      cuts.insert(cuts.end(), alongs.begin(), alongs.end());
    }
    return Kept(spline, cuts, box);
  }

  auto Clip(const Bezier::Quadratic &quadratic, const Box &box) -> std::vector<Range>
  {
    return Kept(quadratic, Cuts(Kernel::Polynomial(quadratic), box), box);
  }

  auto Clip(const Bezier::Cubic &cubic, const Box &box) -> std::vector<Range>
  {
    return Kept(cubic, Cuts(Kernel::Polynomial(cubic), box), box);
  }

  auto Clip(const Curve &curve, const Box &box) -> std::vector<Range>
  {
    return std::visit([&box](const auto &shape) -> std::vector<Range> { return Clip(shape, box); },
                      curve);
  }
} // namespace Geoxy
