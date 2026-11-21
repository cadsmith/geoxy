/***************************************************************************************************
 *  @file      source/outline.cpp
 *  @brief     The bands of a wide polyline, segment by segment, straight corners mitered.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "geoxy/outline.hpp"

#include "real.hpp"

#include "geoxy/accuracy.hpp"
#include "geoxy/angle.hpp"
#include "geoxy/arc.hpp"
#include "geoxy/flatten.hpp"
#include "geoxy/line.hpp"
#include "geoxy/point.hpp"
#include "geoxy/polyline.hpp"
#include "geoxy/segment.hpp"
#include "geoxy/vector.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <optional>
#include <span>
#include <type_traits>
#include <variant>
#include <vector>

namespace Geoxy
{
  namespace
  {
    /**
     *  A segment's band: its right edge and its left edge, each running the segment's way, and the
     *  way of a straight one.
     */
    struct Ribbon
    {
      std::vector<Vertex>   Right   = {}; /**< The right edge, bulges and all. */
      std::vector<Vertex>   Left    = {}; /**< The left edge, bulges and all.  */
      std::optional<Vector> Way     = {}; /**< A straight band's way, unit.    */
      std::size_t           Segment = 0;  /**< The segment it is the band of.  */
    };

    /**
     *  What a band is asked to be: its widths, and how near a spiral's chords must keep.
     */
    struct Brief
    {
      Geoxy::Width    Width    = {}; /**< Its widths.                      */
      Geoxy::Accuracy Accuracy = {}; /**< How near a spiral's chords keep. */
    };

    /**
     *  The band of a line.
     *
     *  @param  line The line, of some length.
     *  @param  half Its half widths at its start and its end.
     *  @return      The band.
     */
    auto Railed(const Line &line, Width half) -> Ribbon
    {
      const Vector way    = (line.End - line.Start) / Length(line.End - line.Start);
      const Vector across = Perpendicular(way);
      return {.Right   = {{.Position = line.Start - (across * half.Start), .Bulge = 0.0},
                          {.Position = line.End - (across * half.End), .Bulge = 0.0}},
              .Left    = {{.Position = line.Start + (across * half.Start), .Bulge = 0.0},
                          {.Position = line.End + (across * half.End), .Bulge = 0.0}},
              .Way     = way,
              .Segment = 0};
    }

    /**
     *  How many chords keep each edge of an arc's band within the accuracy, its width changing
     *  along it: the edge is r(t) u(t), whose second derivative by the angle is at most 2 |r'| +
     *  |r|, so chords over angles of sqrt(8 e / that) stray at most e.
     *
     *  @param  arc    The arc.
     *  @param  half   Its half widths.
     *  @param  within How near the chords keep.
     *  @return        The count, from 1 to MostChords.
     */
    auto Chording(const Arc &arc, Width half, Accuracy within) -> std::size_t
    {
      const double sweep = std::abs(arc.Sweep.InRadians());
      const double bend  = (2.0 * std::abs(half.End - half.Start) / sweep)
                         + arc.Radius
                         + std::max(half.Start, half.End);
      const auto   most  = static_cast<double>(MostChords);
      const double many =
          within.Error > 0.0 ? std::ceil(sweep / std::sqrt(8.0 * within.Error / bend)) : most;
      return static_cast<std::size_t>(std::clamp(many, 1.0, most));
    }

    /**
     *  The band of an arc: concentric arcs for one width, chords along spirals for a width that
     *  changes.
     *
     *  @param  arc   The arc.
     *  @param  ends  Its vertices, the first with its bulge.
     *  @param  asked Its half widths, and the accuracy.
     *  @return       The band.
     */
    auto Curled(const Arc &arc, std::span<const Vertex, 2> ends, Brief asked) -> Ribbon
    {
      const Width  half = asked.Width;
      const double side = arc.Sweep.InRadians() > 0.0 ? 1.0 : -1.0;
      Ribbon       band = {};
      if (Exactly(half.Start, half.End))
      {
        const double share = side * half.Start / arc.Radius;
        const Vector start = (ends[0].Position - arc.Center) * share;
        const Vector end   = (ends[1].Position - arc.Center) * share;
        band.Right         = {{.Position = ends[0].Position + start, .Bulge = ends[0].Bulge},
                              {.Position = ends[1].Position + end, .Bulge = 0.0}};
        band.Left          = {{.Position = ends[0].Position - start, .Bulge = ends[0].Bulge},
                              {.Position = ends[1].Position - end, .Bulge = 0.0}};
        return band;
      }
      const std::size_t chords = Chording(arc, half, asked.Accuracy);
      for (std::size_t index = 0; index <= chords; ++index)
      {
        const double fraction = static_cast<double>(index) / static_cast<double>(chords);
        const double width    = half.Start + ((half.End - half.Start) * fraction);
        const Vector way      = Direction(arc.Start + (arc.Sweep * fraction));
        band.Right.push_back(
            {.Position = arc.Center + (way * (arc.Radius + (side * width))), .Bulge = 0.0});
        band.Left.push_back(
            {.Position = arc.Center + (way * (arc.Radius - (side * width))), .Bulge = 0.0});
      }
      return band;
    }

    /**
     *  The band of a polyline's segment.
     *
     *  @param  polyline The polyline.
     *  @param  index    The segment.
     *  @param  asked    Its widths, and the accuracy.
     *  @return          The band; none for a segment of no width or no length.
     */
    auto Ribboned(const Polyline &polyline, std::size_t index, Brief asked) -> std::optional<Ribbon>
    {
      const std::array<Vertex, 2> ends = {
          polyline.Vertices[index], polyline.Vertices[(index + 1) % polyline.Vertices.size()]};
      const Width half = {.Start = std::abs(asked.Width.Start) / 2.0,
                          .End   = std::abs(asked.Width.End) / 2.0};
      if ((Zero(half.Start) && Zero(half.End)) || ends[0].Position == ends[1].Position)
      {
        return std::nullopt;
      }
      Ribbon band = std::visit(
          [&ends, half, &asked](const auto &segment) -> Ribbon
          {
            if constexpr (std::is_same_v<std::decay_t<decltype(segment)>, Line>)
            {
              return Railed(segment, half);
            }
            else
            {
              return Curled(segment, ends, {.Width = half, .Accuracy = asked.Accuracy});
            }
          },
          Nth(polyline, index));
      band.Segment = index;
      return band;
    }

    /**
     *  Two straight bands mitered where they meet: each side's edges run on to where they cross,
     *  and both end there, one point for both.
     *
     *  @param first  The band before the corner.
     *  @param second The band after it.
     */
    void Mitered(Ribbon &first, Ribbon &second)
    {
      const Vector before = first.Way.value_or(Vector{});
      const Vector after  = second.Way.value_or(Vector{});
      const double turn   = Cross(before, after);
      if (Zero(turn))
      {
        return;
      }
      for (const bool right : {true, false})
      {
        std::vector<Vertex> &ending   = right ? first.Right : first.Left;
        std::vector<Vertex> &starting = right ? second.Right : second.Left;
        const Point          from     = ending.front().Position;
        const double         reach    = Cross(starting.front().Position - from, after) / turn;
        const Point          corner   = from + (before * reach);
        ending.back().Position        = corner;
        starting.front().Position     = corner;
      }
    }

    /**
     *  A band as a closed polyline, counterclockwise: its right edge forward, its left edge back.
     *
     *  @param  band The band.
     *  @return      The polyline.
     */
    auto Outlined(const Ribbon &band) -> Polyline
    {
      Polyline outline = {.Vertices = band.Right, .Closed = true};
      for (std::size_t index = band.Left.size(); index-- > 0;)
      {
        const double bulge = index > 0 ? -band.Left[index - 1].Bulge : 0.0;
        outline.Vertices.push_back({.Position = band.Left[index].Position, .Bulge = bulge});
      }
      return outline;
    }
  } // namespace

  auto Outline(const Polyline &polyline, std::span<const Width> widths, Accuracy accuracy)
      -> std::vector<Polyline>
  {
    std::vector<Ribbon> bands;
    const std::size_t   count = widths.empty() ? 0 : Count(polyline);
    for (std::size_t index = 0; index < count; ++index)
    {
      const Width                 width = widths[std::min(index, widths.size() - 1)];
      const std::optional<Ribbon> band =
          Ribboned(polyline, index, {.Width = width, .Accuracy = accuracy});
      if (band)
      {
        bands.push_back(*band);
      }
    }
    for (std::size_t index = 0; index < bands.size(); ++index)
    {
      const std::size_t next    = (index + 1) % bands.size();
      Ribbon           &first   = bands[index];
      Ribbon           &second  = bands[next];
      const bool        follows = second.Segment == first.Segment + 1;
      const bool wraps  = polyline.Closed && second.Segment == 0 && first.Segment + 1 == count;
      const bool joined = follows || wraps;
      if (joined && first.Way && second.Way)
      {
        Mitered(first, second);
      }
    }
    std::vector<Polyline> outlines;
    outlines.reserve(bands.size());
    for (const Ribbon &band : bands)
    {
      outlines.push_back(Outlined(band));
    }
    return outlines;
  }
} // namespace Geoxy
