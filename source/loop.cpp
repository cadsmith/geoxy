/***************************************************************************************************
 *  @file      source/loop.cpp
 *  @brief     Loops: how many times one winds about a point, by an exact ray cast over its edges
 *             and the lines across its gaps; its area and its perimeter.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "geoxy/loop.hpp"

#include "exact/two.hpp"
#include "finite.hpp"
#include "sane.hpp"
#include "sound.hpp"
#include "winding.hpp"

#include "geoxy/accuracy.hpp"
#include "geoxy/area.hpp"
#include "geoxy/at.hpp"
#include "geoxy/box.hpp"
#include "geoxy/curve.hpp"
#include "geoxy/domain.hpp"
#include "geoxy/extents.hpp"
#include "geoxy/length.hpp"
#include "geoxy/line.hpp"
#include "geoxy/nearest.hpp"
#include "geoxy/point.hpp"
#include "geoxy/polyline.hpp"
#include "geoxy/ray.hpp"
#include "geoxy/segment.hpp"
#include "geoxy/tolerance.hpp"
#include "geoxy/xline.hpp"

#include <cstddef>
#include <cstdint>
#include <deque>
#include <limits>
#include <optional>
#include <variant>
#include <vector>

namespace Geoxy
{
  namespace
  {
    /**
     *  Whether a curve runs without end.
     *
     *  @param  curve The curve.
     *  @return       Whether it is a ray or an endless line.
     */
    auto Endless(const Curve &curve) -> bool
    {
      return std::holds_alternative<Ray>(curve) || std::holds_alternative<Xline>(curve);
    }

    /**
     *  Whether a curve has ends and finite numbers.
     *
     *  @param  curve The curve.
     *  @return       Whether it has.
     */
    auto Finished(const Curve &curve) -> bool
    {
      return !Endless(curve) && Sound(curve);
    }

    /**
     *  Where an edge starts and where it ends.
     *
     *  @param  edge The edge.
     *  @return      Its points at the start and the end of its parameter.
     */
    auto Termini(const Curve &edge) -> Traverse
    {
      const Range domain = Domain(edge);
      return {.Shape = &edge, .From = At(edge, domain.Start), .To = At(edge, domain.End)};
    }

    /**
     *  A loop's courses: each edge from its start to its end, a polyline's segments one by one from
     *  vertex to vertex, and the line across each gap from an edge's end to the next one's start.
     */
    class Legs
    {
    public:
      /**
       *  The courses of a loop.
       *
       *  @param loop The loop, which must outlive this.
       */
      explicit Legs(const Loop &loop)
      {
        const std::size_t count = loop.Edges.size();
        for (std::size_t index = 0; index < count; ++index)
        {
          const Curve &edge   = loop.Edges[index];
          m_Finite            = m_Finite && Finished(edge);
          const Traverse ends = Termini(edge);
          Walk(edge, ends);
          const Point next = Termini(loop.Edges[(index + 1) % count]).From;
          if (next != ends.To)
          {
            m_Parts.emplace_back(Line{.Start = ends.To, .End = next});
            m_Courses.push_back({.Shape = &m_Parts.back(), .From = ends.To, .To = next});
          }
        }
      }

      /**
       *  Whether every edge has ends and finite numbers.
       *
       *  @return Whether they have.
       */
      [[nodiscard]] auto Finite() const -> bool { return m_Finite; }

      /**
       *  The courses, in order along the loop.
       *
       *  @return Them.
       */
      [[nodiscard]] auto Courses() const -> const std::vector<Traverse> & { return m_Courses; }

    private:
      /**
       *  An edge's courses: a polyline's segments, each from its vertex to the next, or the edge.
       *
       *  @param edge The edge.
       *  @param ends Where it starts and ends.
       */
      void Walk(const Curve &edge, const Traverse &ends)
      {
        const auto *polyline = std::get_if<Polyline>(&edge);
        if (polyline == nullptr)
        {
          m_Courses.push_back(ends);
          return;
        }
        const std::size_t segments = Count(*polyline);
        const std::size_t vertices = polyline->Vertices.size();
        for (std::size_t index = 0; index < segments; ++index)
        {
          m_Parts.push_back(
              std::visit([](const auto &shape) -> Curve { return shape; }, Nth(*polyline, index)));
          m_Courses.push_back({.Shape = &m_Parts.back(),
                               .From  = polyline->Vertices[index].Position,
                               .To    = polyline->Vertices[(index + 1) % vertices].Position});
        }
      }

      std::deque<Curve>     m_Parts   = {};   /**< The segments and gaps' lines, kept in place. */
      std::vector<Traverse> m_Courses = {};   /**< The courses.                                 */
      bool                  m_Finite  = true; /**< Whether every edge has ends, finite.         */
    };

    /**
     *  Whether a point is within a length of a traverse's curve.
     *
     *  @param  traverse The traverse.
     *  @param  point    The point.
     *  @param  length   The length.
     *  @return          Whether it is.
     */
    auto Near(const Traverse &traverse, Point point, double length) -> bool
    {
      if (!Contains(Padded(Extents(*traverse.Shape), length), point))
      {
        return false;
      }
      return Distance(Nearest(*traverse.Shape, point).At, point) <= length;
    }

    /**
     *  A sum of doubles carried with its rounding error.
     */
    class Total
    {
    public:
      /**
       *  A double added.
       *
       *  @param value The double.
       */
      void Add(double value)
      {
        const Exact::Two sum  = Exact::TwoSum(m_Big, value);
        m_Big                 = sum.Big;
        m_Small              += sum.Small;
      }

      /**
       *  The sum.
       *
       *  @return It, rounded once from the two parts.
       */
      [[nodiscard]] auto Sum() const -> double { return m_Big + m_Small; }

    private:
      double m_Big   = 0.0; /**< The sum, rounded.   */
      double m_Small = 0.0; /**< What rounding took. */
    };
  } // namespace

  auto Winding(const Loop &loop, Point point, const Tolerance &tolerance) -> std::optional<int>
  {
    const Legs legs{loop};
    if (!Finite(point) || !legs.Finite())
    {
      return std::nullopt;
    }
    const double length = Sane(tolerance).Length;
    std::int64_t count  = 0;
    for (const Traverse &traverse : legs.Courses())
    {
      if (Near(traverse, point, length))
      {
        return std::nullopt;
      }
      count += Ascents(traverse, point);
    }
    if (count < std::numeric_limits<int>::min() || count > std::numeric_limits<int>::max())
    {
      return std::nullopt;
    }
    return static_cast<int>(count);
  }

  auto Area(const Loop &loop, Accuracy accuracy) -> double
  {
    Total    total   = {};
    Polyline polygon = {.Closed = true};
    for (const Curve &edge : loop.Edges)
    {
      if (Endless(edge))
      {
        return 0.0;
      }
      total.Add(Area(edge, accuracy));
      const Traverse ends = Termini(edge);
      polygon.Vertices.push_back({.Position = ends.From});
      polygon.Vertices.push_back({.Position = ends.To});
    }
    total.Add(Area(polygon));
    return total.Sum();
  }

  auto Length(const Loop &loop, Accuracy accuracy) -> double
  {
    const std::size_t count = loop.Edges.size();
    Total             total = {};
    for (std::size_t index = 0; index < count; ++index)
    {
      const Curve &edge = loop.Edges[index];
      if (Endless(edge))
      {
        return std::numeric_limits<double>::infinity();
      }
      total.Add(Length(edge, accuracy));
      total.Add(Distance(Termini(edge).To, Termini(loop.Edges[(index + 1) % count]).From));
    }
    return total.Sum();
  }

  auto Closed(const Loop &loop, const Tolerance &tolerance) -> bool
  {
    const std::size_t count = loop.Edges.size();
    bool              shut  = count > 0;
    for (std::size_t index = 0; index < count; ++index)
    {
      const Curve &edge = loop.Edges[index];
      shut              = shut
                       && Finished(edge)
                       && Same(Termini(edge).To, Termini(loop.Edges[(index + 1) % count]).From, tolerance);
    }
    return shut;
  }
} // namespace Geoxy
