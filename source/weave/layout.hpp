/***************************************************************************************************
 *  @file      source/weave/layout.hpp
 *  @brief     The curves of an arrangement laid out: each run's curves cut into lines and arcs,
 *             with what each weighs.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#ifndef GEOXY_WEAVE_LAYOUT_HPP
#define GEOXY_WEAVE_LAYOUT_HPP

#include "arranged.hpp"
#include "real.hpp"
#include "sound.hpp"
#include "weave/found.hpp"

#include "geoxy/angle.hpp"
#include "geoxy/arc.hpp"
#include "geoxy/curve.hpp"
#include "geoxy/error.hpp"
#include "geoxy/line.hpp"
#include "geoxy/point.hpp"
#include "geoxy/polyline.hpp"
#include "geoxy/segment.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <memory_resource>
#include <optional>
#include <span>
#include <utility>
#include <variant>
#include <vector>

namespace Geoxy::Weave
{
  /**
   *  The curves an arrangement is cut from, laid out as lines and arcs, and the sizes of their
   *  numbers gathered as they are laid.
   */
  class Layout
  {
  public:
    /**
     *  A layout of no curve yet.
     *
     *  @param length The tolerance's length, made fit.
     *  @param memory What the arrangement's lists draw on.
     */
    explicit Layout(double length, std::pmr::memory_resource *memory)
        : m_Arranged(memory), m_Length(length)
    {
    }

    /**
     *  Room made for the lines and arcs of some runs of curves: a polyline's segments, one for any
     *  other curve.
     *
     *  @param inputs The runs.
     */
    void Reserve(const Inputs &inputs)
    {
      std::size_t count = 0;
      for (const std::span<const Curve> run : inputs.Runs)
      {
        for (const Curve &curve : run)
        {
          const auto *polyline  = std::get_if<Polyline>(&curve);
          count                += polyline != nullptr ? Count(*polyline) : 1;
        }
      }
      m_Arranged.Sources.reserve(count);
      m_Tallies.reserve(count);
    }

    /**
     *  A curve laid out: a line or an arc as it is, a polyline as its segments, each with the
     *  curve's tally; any other kind, or numbers that are not finite, a fault.
     *
     *  @param curve The curve.
     *  @param tally Its tally.
     */
    void Add(const Curve &curve, const Tally &tally)
    {
      if (!Sound(curve))
      {
        Fail(Error::NotFinite);
        return;
      }
      if (const auto *polyline = std::get_if<Polyline>(&curve))
      {
        /**
         *  Each segment as Nth makes it, a straight one made a line in place as Bent makes it.
         */
        const std::span<const Vertex> vertices = polyline->Vertices;
        const std::size_t             count    = Count(*polyline);
        for (std::size_t index = 0; index < count; ++index)
        {
          const Vertex &start = vertices[index];
          const Point   end   = vertices[index + 1 < vertices.size() ? index + 1 : 0].Position;
          if (!(std::abs(start.Bulge) >= Straight) || start.Position == end)
          {
            Put(Line{.Start = start.Position, .End = end}, tally);
            continue;
          }
          std::visit([this, &tally](const auto &shape) -> void { Put(shape, tally); },
                     Bent(start.Position, end, start.Bulge));
        }
        return;
      }
      const auto *arc   = std::get_if<Arc>(&curve);
      const bool  plain = std::holds_alternative<Line>(curve)
                       || (arc != nullptr && !(std::abs(arc->Sweep.InDegrees()) > 360.0));
      if (!plain)
      {
        Fail(Error::Degenerate);
        return;
      }
      Put(curve, tally);
    }

    /**
     *  The arrangement begun: its curves and their tallies, or its fault.
     *
     *  @return The arrangement, and the tallies.
     */
    [[nodiscard]] auto Laid() -> std::pair<Arranged, std::vector<Tally>>
    {
      m_Arranged.Runs.push_back(m_Arranged.Sources.size());
      return {std::move(m_Arranged), std::move(m_Tallies)};
    }

    /**
     *  The length within which two vertices are one: the tolerance's, made fit, and no less than 16
     *  steps of a double at the size of the curves' numbers, as Intersect takes it; and the length
     *  the curves are met within, at its most for any pair of them, at the size Intersect reads.
     *
     *  @return The two lengths, above 0.
     */
    [[nodiscard]] auto Lengths() const -> std::array<double, 2>
    {
      return {std::max(m_Length, 16.0 * Above(m_Size)), std::max(m_Length, 16.0 * Above(m_Reach))};
    }

    /**
     *  A run of curves begun: its first curve's index kept.
     */
    void Begin() { m_Arranged.Runs.push_back(m_Arranged.Sources.size()); }

  private:
    /**
     *  A line or an arc put down, made in place.
     *
     *  @tparam ShapeT The kind of the curve: a line, an arc, or a curve holding one.
     *  @param  shape  The curve.
     *  @param  tally  Its tally.
     */
    template <typename ShapeT>
    void Put(const ShapeT &shape, const Tally &tally)
    {
      m_Arranged.Sources.emplace_back(shape);
      m_Tallies.push_back(tally);
      const std::array<double, 2> sizes = Sized(m_Arranged.Sources.back(), m_Length);
      m_Size                            = std::max(m_Size, sizes[0]);
      m_Reach                           = std::max(m_Reach, sizes[1]);
    }

    /**
     *  A fault kept, the first one found.
     *
     *  @param error The fault.
     */
    void Fail(Error error)
    {
      if (!m_Arranged.Fault)
      {
        m_Arranged.Fault = error;
      }
    }

    Arranged           m_Arranged;      /**< The arrangement begun.                            */
    std::vector<Tally> m_Tallies = {};  /**< Each source's tally.                              */
    double             m_Length  = 0.0; /**< The tolerance's length, made fit.                 */
    double             m_Size    = 0.0; /**< The largest of the curves' numbers so far.        */
    double             m_Reach   = 0.0; /**< The largest each may reach as Intersect sizes it. */
  };
} // namespace Geoxy::Weave

#endif /* GEOXY_WEAVE_LAYOUT_HPP */
