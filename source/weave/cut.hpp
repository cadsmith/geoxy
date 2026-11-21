/***************************************************************************************************
 *  @file      source/weave/cut.hpp
 *  @brief     A curve's notches gathered as they are found, each placed along it.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#ifndef GEOXY_WEAVE_CUT_HPP
#define GEOXY_WEAVE_CUT_HPP

#include "real.hpp"
#include "scratch.hpp"
#include "seek.hpp"
#include "weave/found.hpp"

#include "geoxy/angle.hpp"
#include "geoxy/arc.hpp"
#include "geoxy/at.hpp"
#include "geoxy/curve.hpp"
#include "geoxy/domain.hpp"
#include "geoxy/intersect.hpp"
#include "geoxy/line.hpp"
#include "geoxy/nearest.hpp"
#include "geoxy/place.hpp"
#include "geoxy/point.hpp"
#include "geoxy/prepared.hpp"

#include <array>
#include <cmath>
#include <cstddef>
#include <span>
#include <utility>
#include <variant>
#include <vector>

namespace Geoxy::Weave
{
  /**
   *  Where a point lies along a line or an arc, as Nearest places it: a foot on a line, an angle on
   *  an arc, asked of the curve's own kind.
   *
   *  @param  curve The line or arc.
   *  @param  point The point.
   *  @param  ends  Its start and end, as At gives them.
   *  @return       The curve's parameter there.
   */
  inline auto Footed(const Curve &curve, Point point, const std::array<Point, 2> &ends) -> double
  {
    if (const auto *arc = std::get_if<Arc>(&curve))
    {
      return Nearest(*arc, point, ends).Along;
    }
    return Nearest(std::get<Line>(curve), point).Along;
  }

  /**
   *  A vertex added, and a cut of a curve at it.
   */
  class Cutter
  {
  public:
    /**
     *  Cuts of some curves, each at first a curve of no cut, with room for their ends and a few
     *  meetings each.
     *
     *  @param count How many curves.
     */
    explicit Cutter(std::size_t count)
        : m_Found{.Firsts = List<std::size_t>(count + 1, 0), .Count = count}
    {
      m_Found.Points.reserve(3 * count);
      m_Found.Twins.reserve(3 * count);
      m_Found.Notches.reserve(count / 4);
    }

    /**
     *  What was found, taken.
     *
     *  @return The findings.
     */
    [[nodiscard]] auto Found() -> Findings { return std::move(m_Found); }

    /**
     *  A vertex added at a point, of no twin known.
     *
     *  @param  point The point.
     *  @return       Its index.
     */
    auto Vertex(Point point) -> std::size_t { return Twin(point, Unknown); }

    /**
     *  A vertex added at a point, and an earlier one at the same point, if any.
     *
     *  @param  point The point.
     *  @param  twin  The earlier vertex at the same point, or Unknown.
     *  @return       Its index.
     */
    auto Twin(Point point, std::size_t twin) -> std::size_t
    {
      m_Found.Points.push_back(point);
      m_Found.Twins.push_back(twin);
      return m_Found.Points.size() - 1;
    }

    /**
     *  A curve's ends made its first two vertices, but those Lonely leaves out, an arc's middle a
     *  third where it is past a half turn, so that no edge is a loop of its own; its start the twin
     *  of the last curve's end where they are the same point, as a loop's run of curves most often
     *  has it.
     *
     *  @param curve The curve's index.
     *  @param ends  Its start and its end, as At gives them.
     *  @param bare  Whether its start, and whether its end, is a lone joint: no vertex.
     */
    void Placed(std::size_t curve, const std::array<Point, 2> &ends,
                const std::array<bool, 2> &bare)
    {
      const std::size_t last = m_Last;
      const bool        twin = last != Unknown && Twinned(ends[0], m_Found.Points[last]);
      m_Found.Firsts[curve]  = m_Found.Points.size();
      m_Last                 = Unknown;
      if (!bare[0])
      {
        Twin(ends[0], twin ? last : Unknown);
      }
      if (!bare[1])
      {
        m_Last = Vertex(ends[1]);
      }
      m_Found.Firsts[curve + 1] = m_Found.Points.size();
    }

    /**
     *  An arc's middle made its third vertex, where it is past a half turn.
     *
     *  @param curve The arc's index.
     *  @param arc   The arc.
     */
    void Halved(std::size_t curve, const Arc &arc)
    {
      if (std::abs(arc.Sweep.InDegrees()) > 180.0)
      {
        Vertex(At(arc, 0.5));
        m_Found.Firsts[curve + 1] = m_Found.Points.size();
      }
    }

    /**
     *  A cut of a curve at a vertex.
     *
     *  @param curve The curve's index.
     *  @param notch Where along it, and the vertex.
     */
    void Cut(std::size_t curve, const Notch &notch)
    {
      m_Found.Notches.push_back({.Source = curve, .Along = notch.Along, .Vertex = notch.Vertex});
    }

    /**
     *  What meeting two curves found, cut into both: each mark a vertex of both, at the place of
     *  each nearest it, as Intersect places a hit, a foot on a line, an angle on an arc; each
     *  stretch's ends vertices of both, the stretch kept.
     *
     *  @param pair    The two curves.
     *  @param marking What they met in.
     *  @param sources The curves.
     */
    void Record(const Couple &pair, const Marking &marking, std::span<const Curve> sources)
    {
      const Curve &first   = sources[pair[0]];
      const Curve &second  = sources[pair[1]];
      m_Found.Spent       += static_cast<std::size_t>(marking.Spent);
      for (const Mark &mark : marking.Marks)
      {
        const std::array<double, 2> places = {Footed(first, mark.At, Terminals(pair[0])),
                                              Footed(second, mark.At, Terminals(pair[1]))};
        const std::array<bool, 2>   idle   = {Idle(pair[0], {.Along = places[0]}, mark.At),
                                              Idle(pair[1], {.Along = places[1]}, mark.At)};
        if (idle[0] && idle[1])
        {
          continue;
        }
        const std::size_t vertex = Twin(mark.At, Matched(pair, mark.At));
        for (std::size_t side = 0; side < 2; ++side)
        {
          if (!idle.at(side))
          {
            Cut(pair.at(side), {.Along = places.at(side), .Vertex = vertex});
          }
        }
      }
      for (const Overlap &overlap : marking.Overlaps)
      {
        const std::size_t start = Vertex(At(first, overlap.First.Start));
        const std::size_t end   = Vertex(At(first, overlap.First.End));
        const Range       other = overlap.Second;
        Cut(pair[0], {.Along = overlap.First.Start, .Vertex = start});
        Cut(pair[0], {.Along = overlap.First.End, .Vertex = end});
        Cut(pair[1], {.Along = overlap.Along ? other.Start : other.End, .Vertex = start});
        Cut(pair[1], {.Along = overlap.Along ? other.End : other.Start, .Vertex = end});
        m_Found.Commons.push_back({.First = pair[0], .Second = pair[1], .Shared = overlap});
      }
    }

  private:
    /**
     *  A curve's start and end, as its first two vertices hold them.
     *
     *  @param  curve The curve's index.
     *  @return       The two points.
     */
    [[nodiscard]] auto Terminals(std::size_t curve) const -> std::array<Point, 2>
    {
      const std::size_t first = m_Found.Firsts[curve];
      return {m_Found.Points[first], m_Found.Points[first + 1]};
    }

    /**
     *  The index that stands for no twin, past any vertex.
     */
    static constexpr std::size_t Unknown = static_cast<std::size_t>(-1);

    /**
     *  The end of either of two curves at a point, if any is there.
     *
     *  @param  pair  The two curves.
     *  @param  point The point.
     *  @return       The end's vertex, or Unknown.
     */
    [[nodiscard]] auto Matched(const Couple &pair, Point point) const -> std::size_t
    {
      for (const std::size_t curve : pair)
      {
        const std::size_t first = m_Found.Firsts[curve];
        for (const std::size_t end : {first, first + 1})
        {
          if (Twinned(m_Found.Points[end], point))
          {
            return end;
          }
        }
      }
      return Unknown;
    }

    /**
     *  Whether a cut at a mark would change nothing: the mark at a curve's end's very point, and
     *  its place there the end's very place, 0 or 1. Its vertex would join the end's, being at the
     *  same point after it; the cut would come beside the end's, at one place, joining their
     *  vertices again; and no piece would lie between them.
     *
     *  @param  curve The curve's index.
     *  @param  notch The cut's place along it.
     *  @param  point The mark's point.
     *  @return       Whether it would.
     */
    [[nodiscard]] auto Idle(std::size_t curve, const Notch &notch, Point point) const -> bool
    {
      const std::size_t start = m_Found.Firsts[curve];
      return (Exactly(notch.Along, 0.0) && Twinned(m_Found.Points[start], point))
          || (Exactly(notch.Along, 1.0) && Twinned(m_Found.Points[start + 1], point));
    }

    Findings    m_Found = {};      /**< The findings. */
    std::size_t m_Last  = Unknown; /**< The last curve's end's vertex, or
                                        Unknown where it has none.    */
  };
} // namespace Geoxy::Weave

#endif /* GEOXY_WEAVE_CUT_HPP */
