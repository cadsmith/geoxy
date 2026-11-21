/***************************************************************************************************
 *  @file      source/corner/rework.cpp
 *  @brief     A polyline's corners replaced in order: each segment's part kept narrowed by the
 *             corners at its ends, a corner left where its cut falls outside it, and the vertices
 *             laid again.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "corner/rework.hpp"

#include "corner/choose.hpp"
#include "corner/guide.hpp"
#include "real.hpp"
#include "trig.hpp"

#include "geoxy/arc.hpp"
#include "geoxy/curve.hpp"
#include "geoxy/fillet.hpp"
#include "geoxy/line.hpp"
#include "geoxy/point.hpp"
#include "geoxy/polyline.hpp"
#include "geoxy/segment.hpp"
#include "geoxy/tolerance.hpp"
#include "geoxy/vector.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <expected>
#include <optional>
#include <span>
#include <variant>
#include <vector>

namespace Geoxy::Corner
{
  namespace
  {
    /**
     *  A segment's length.
     *
     *  @param  piece The segment.
     *  @return       Its chord's length for a line, its arc's for an arc.
     */
    auto Long(const Segment &piece) -> double
    {
      if (const auto *line = std::get_if<Line>(&piece))
      {
        return Distance(line->Start, line->End);
      }
      const Arc &arc = std::get<Arc>(piece);
      return std::abs(arc.Radius * arc.Sweep.InRadians());
    }

    /**
     *  Whether two segments run on smoothly at a corner: the way the first ends and the second
     *  starts within the tolerance's angle of each other.
     *
     *  @param  guides    The two segments' guides.
     *  @param  corner    The corner.
     *  @param  tolerance The tolerance.
     *  @return           Whether they do.
     */
    auto Smooth(const Guides &guides, Point corner, const Tolerance &tolerance) -> bool
    {
      const Vector leaving  = Unit(Ahead(*guides[0], 1.0, corner)).value_or(Vector{});
      const Vector entering = Unit(Ahead(*guides[1], 0.0, corner)).value_or(Vector{});
      return Dot(leaving, entering)             > 0.0
          && std::abs(Cross(leaving, entering)) <= tolerance.Radians;
    }

    /**
     *  A segment as a curve.
     *
     *  @param  piece The segment.
     *  @return       The line or the arc.
     */
    auto Curved(const Segment &piece) -> Curve
    {
      if (const auto *line = std::get_if<Line>(&piece))
      {
        return *line;
      }
      return std::get<Arc>(piece);
    }

    /**
     *  The bulge of a part of a segment.
     *
     *  @param  bulge The segment's bulge.
     *  @param  kept  The part of its sweep kept.
     *  @return       The part's bulge; the segment's own where it is kept whole.
     */
    auto Narrowed(double bulge, double kept) -> double
    {
      return Exactly(kept, 1.0) ? bulge : Trig::Tan(Trig::Atan(bulge) * kept);
    }

    /**
     *  What each corner did to a polyline: the joints, and the part of each segment kept.
     */
    struct Narrowing
    {
      std::span<const std::optional<Joint>> Joints  = {}; /**< The joint at each vertex, if any. */
      std::span<const double>               Starts  = {}; /**< Where each segment's part starts. */
      std::span<const double>               Ends    = {}; /**< Where it ends.                    */
      std::span<const double>               Lengths = {}; /**< Each segment's length.            */
    };

    /**
     *  A polyline's vertices laid again: each segment's part kept, a joint's arc or line after it,
     *  a part shorter than the tolerance between two joints dropped.
     *
     *  @param  polyline  The polyline.
     *  @param  narrowing The joints and the parts kept.
     *  @param  tolerance Within which a part is of no length.
     *  @return           The polyline.
     */
    auto Relaid(const Polyline &polyline, const Narrowing &narrowing, const Tolerance &tolerance)
        -> Polyline
    {
      const std::vector<Vertex> &vertices = polyline.Vertices;
      const std::size_t          count    = narrowing.Starts.size();
      const std::size_t          size     = vertices.size();
      Polyline                   laid     = {.Vertices = {}, .Closed = polyline.Closed};
      std::vector<Vertex>       &out      = laid.Vertices;
      for (std::size_t index = 0; index < count; ++index)
      {
        const std::size_t           next  = (index + 1) % size;
        const std::optional<Joint> &here  = narrowing.Joints[index];
        const std::optional<Joint> &there = narrowing.Joints[next];
        const Point                 start = here ? here->Ends[1] : vertices[index].Position;
        const double                kept  = narrowing.Ends[index] - narrowing.Starts[index];
        if (here && there && kept * narrowing.Lengths[index] <= tolerance.Length)
        {
          out.push_back({.Position = start, .Bulge = there->Bulge});
          continue;
        }
        out.push_back({.Position = start, .Bulge = Narrowed(vertices[index].Bulge, kept)});
        if (there)
        {
          out.push_back({.Position = there->Ends[0], .Bulge = there->Bulge});
        }
      }
      if (!polyline.Closed && size > 0)
      {
        out.push_back({.Position = vertices.back().Position, .Bulge = 0.0});
      }
      return laid;
    }
  } // namespace

  auto Reworked(const Polyline &polyline, const Joiner &joiner, const Tolerance &tolerance)
      -> Reshaped
  {
    const std::vector<Vertex> &vertices = polyline.Vertices;
    const std::size_t          count    = Count(polyline);
    const std::size_t          size     = vertices.size();
    std::vector<Curve>         pieces;
    std::vector<double>        lengths;
    for (std::size_t index = 0; index < count; ++index)
    {
      const Segment piece = Nth(polyline, index);
      pieces.push_back(Curved(piece));
      lengths.push_back(Long(piece));
    }
    std::vector<double>               starts(count, 0.0);
    std::vector<double>               ends(count, 1.0);
    std::vector<std::optional<Joint>> joints(size);
    Reshaped                          reshaped;
    const std::size_t                 first = polyline.Closed && count > 1 ? 0 : 1;
    const std::size_t                 last =
        polyline.Closed && count > 1 ? size : size - std::min<std::size_t>(size, 1);
    for (std::size_t vertex = first; vertex < last; ++vertex)
    {
      const std::size_t                  before = (vertex + count - 1) % count;
      const std::size_t                  after  = vertex;
      const Point                        corner = vertices[vertex].Position;
      const std::expected<Guide, Misfit> leaving =
          Guided(pieces[before], vertices[before].Position);
      const std::expected<Guide, Misfit> entering =
          Guided(pieces[after], vertices[(after + 1) % size].Position);
      if (!leaving || !entering)
      {
        reshaped.Left.push_back({.Vertex = vertex, .Why = Misfit::Degenerate});
        continue;
      }
      const Guides guides = {&*leaving, &*entering};
      if (Smooth(guides, corner, tolerance))
      {
        continue;
      }
      const std::expected<Joint, Misfit> joint = joiner(guides, corner);
      if (!joint)
      {
        reshaped.Left.push_back({.Vertex = vertex, .Why = joint.error()});
        continue;
      }
      const double slackBefore = tolerance.Length / lengths[before];
      const double slackAfter  = tolerance.Length / lengths[after];
      if (joint->Alongs[0]    < starts[before] - slackBefore
          || joint->Alongs[1] > ends[after] + slackAfter)
      {
        reshaped.Left.push_back({.Vertex = vertex, .Why = Misfit::TooLarge});
        continue;
      }
      ends[before]   = std::max(joint->Alongs[0], starts[before]);
      starts[after]  = std::min(joint->Alongs[1], ends[after]);
      joints[vertex] = *joint;
    }
    reshaped.Shape =
        Relaid(polyline, {.Joints = joints, .Starts = starts, .Ends = ends, .Lengths = lengths},
               tolerance);
    return reshaped;
  }
} // namespace Geoxy::Corner
