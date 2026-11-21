/***************************************************************************************************
 *  @file      source/contain.cpp
 *  @brief     Whether a curve lies within a box, from its extents; and whether it meets one: a line
 *             by the sides of it the box's corners lie on, any other curve by its clip.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "geoxy/contain.hpp"

#include "sound.hpp"

#include "geoxy/at.hpp"
#include "geoxy/box.hpp"
#include "geoxy/clip.hpp"
#include "geoxy/curve.hpp"
#include "geoxy/domain.hpp"
#include "geoxy/extents.hpp"
#include "geoxy/intersect.hpp"
#include "geoxy/line.hpp"
#include "geoxy/loop.hpp"
#include "geoxy/orient.hpp"
#include "geoxy/point.hpp"
#include "geoxy/polyline.hpp"
#include "geoxy/region.hpp"
#include "geoxy/tolerance.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <variant>
#include <vector>

namespace Geoxy
{
  namespace
  {
    /**
     *  Whether a box holds another, edges included.
     *
     *  @param  outer The box that holds.
     *  @param  inner The box held, not empty.
     *  @return       Whether it does.
     */
    auto Holds(const Box &outer, const Box &inner) -> bool
    {
      return Contains(outer, inner.Min) && Contains(outer, inner.Max);
    }

    /**
     *  Whether a line meets a box.
     *
     *  A line and a box, both closed and convex, are apart only where an axis parts them: x or y,
     *  where the box of the line's ends does not overlap the box, or the line's normal, where every
     *  corner of the box lies on one side of the line, none on it. The box is first cut to the box
     *  of the ends, which it meets where the line does, so that its corners are finite where the
     *  line's are; a corner of numbers that are not finite is told on the line, which keeps any
     *  doubt on the side of meeting.
     *
     *  @param  line The line.
     *  @param  box  The box.
     *  @return      Whether they meet.
     */
    auto Crossed(const Line &line, const Box &box) -> bool
    {
      const Box ends = Grown(Grown(Box{}, line.Start), line.End);
      const Box cut  = {
          .Min = {.X = std::max(ends.Min.X, box.Min.X), .Y = std::max(ends.Min.Y, box.Min.Y)},
          .Max = {.X = std::min(ends.Max.X, box.Max.X), .Y = std::min(ends.Max.Y, box.Max.Y)},
      };
      if (IsEmpty(cut) || IsEmpty(box))
      {
        return false;
      }
      const std::array<Point, 4> corners = {
          cut.Min,
          Point{.X = cut.Max.X, .Y = cut.Min.Y},
          cut.Max,
          Point{.X = cut.Min.X, .Y = cut.Max.Y},
      };
      int sides = 0;
      for (const Point corner : corners)
      {
        sides += static_cast<int>(Orient(line.Start, line.End, corner));
      }
      return std::abs(sides) < 4;
    }

    /**
     *  The box of a region's edges.
     *
     *  @param  region The region.
     *  @return        The box of every edge's extents.
     */
    auto Bounds(const Region &region) -> Box
    {
      Box box = {};
      for (const Loop &loop : region.Loops)
      {
        for (const Curve &edge : loop.Edges)
        {
          box = Joined(box, Extents(edge));
        }
      }
      return box;
    }

    /**
     *  The places where a curve meets a region's edges.
     *
     *  @param  curve     The curve.
     *  @param  region    The region.
     *  @param  tolerance The tolerance, as Intersect takes it.
     *  @return           The curve's parameters at every hit, and at both ends of every stretch it
     *                    shares with an edge, in no order.
     */
    auto Cuts(const Curve &curve, const Region &region, const Tolerance &tolerance)
        -> std::vector<double>
    {
      std::vector<double> alongs;
      for (const Loop &loop : region.Loops)
      {
        for (const Curve &edge : loop.Edges)
        {
          const Intersections meeting = Intersect(curve, edge, tolerance);
          for (const Hit &hit : meeting.Hits)
          {
            alongs.push_back(hit.First.Along);
          }
          for (const Overlap &overlap : meeting.Overlaps)
          {
            alongs.push_back(overlap.First.Start);
            alongs.push_back(overlap.First.End);
          }
        }
      }
      return alongs;
    }

    /**
     *  A polygon's region: one loop, the polygon.
     *
     *  @param  polygon The polygon.
     *  @return         The region.
     */
    auto Filling(const Polyline &polygon) -> Region
    {
      return {.Loops = {Loop{.Edges = {polygon}}}};
    }
  } // namespace

  auto Inside(const Curve &curve, const Box &box) -> bool
  {
    const Box extents = Extents(curve);
    return !IsEmpty(extents) && Holds(box, extents);
  }

  auto Intersects(const Curve &curve, const Box &box) -> bool
  {
    if (const auto *line = std::get_if<Line>(&curve))
    {
      return Crossed(*line, box);
    }
    const Box extents = Extents(curve);
    if (!Overlaps(extents, box))
    {
      return false;
    }
    return Holds(box, extents) || !Clip(curve, box).empty();
  }

  auto Inside(const Curve &curve, const Region &region, const Tolerance &tolerance) -> bool
  {
    const Box extents = Extents(curve);
    const Box bounds  = Padded(Bounds(region), tolerance.Length);
    if (!Sound(curve)
        || IsEmpty(extents)
        || IsEmpty(bounds)
        || !Contains(bounds, extents.Min)
        || !Contains(bounds, extents.Max))
    {
      return false;
    }
    const Range         domain = Domain(curve);
    std::vector<double> alongs = Cuts(curve, region, tolerance);
    alongs.push_back(domain.Start);
    alongs.push_back(domain.End);
    std::ranges::sort(alongs);
    for (std::size_t index = 0; index + 1 < alongs.size(); ++index)
    {
      const double start = std::max(alongs[index], domain.Start);
      const double end   = std::min(alongs[index + 1], domain.End);
      if (end > start
          && Locate(region, At(curve, start + ((end - start) / 2.0)), tolerance)
                 == Location::Outside)
      {
        return false;
      }
    }
    return true;
  }

  auto Intersects(const Curve &curve, const Region &region, const Tolerance &tolerance) -> bool
  {
    const Box bounds = Padded(Bounds(region), tolerance.Length);
    if (!Sound(curve) || IsEmpty(bounds) || !Overlaps(Extents(curve), bounds))
    {
      return false;
    }
    if (!Cuts(curve, region, tolerance).empty())
    {
      return true;
    }
    const double start = Domain(curve).Start;
    const double along = std::isfinite(start) ? start : 0.0;
    return Locate(region, At(curve, along), tolerance) != Location::Outside;
  }

  auto Inside(const Curve &curve, const Polyline &polygon, const Tolerance &tolerance) -> bool
  {
    return Inside(curve, Filling(polygon), tolerance);
  }

  auto Intersects(const Curve &curve, const Polyline &polygon, const Tolerance &tolerance) -> bool
  {
    return Intersects(curve, Filling(polygon), tolerance);
  }
} // namespace Geoxy
