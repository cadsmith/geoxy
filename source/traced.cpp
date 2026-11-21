/***************************************************************************************************
 *  @file      source/traced.cpp
 *  @brief     An arrangement's faces: the halves of its edges ordered about each vertex by their
 *             exact tangents and their curvature, its cycles traced, and how many times each
 *             operand winds about each face.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "traced.hpp"

#include "arranged.hpp"
#include "exact/two.hpp"
#include "scratch.hpp"
#include "unions.hpp"
#include "winding.hpp"

#include "geoxy/angle.hpp"
#include "geoxy/arc.hpp"
#include "geoxy/area.hpp"
#include "geoxy/at.hpp"
#include "geoxy/box.hpp"
#include "geoxy/curve.hpp"
#include "geoxy/derivative.hpp"
#include "geoxy/extents.hpp"
#include "geoxy/line.hpp"
#include "geoxy/point.hpp"
#include "geoxy/vector.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <compare>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <numeric>
#include <span>
#include <utility>
#include <variant>
#include <vector>

namespace Geoxy
{
  namespace
  {
    /**
     *  The area a closed run of points encloses, as Area takes a polygon's: each side's cross
     *  product about the first point halved, its two products exact and summed with what rounding
     *  took carried apart, so that a cycle's polygon needs no polyline made of it.
     *
     *  Or plainly, from the same two differences, its two products and every sum rounded once, with
     *  a slack it lies within of the exact way's: of n sides and M the sum of the products'
     *  magnitudes, halved, the plain sum lies within (n + 3) u M of the exact one, the exact way's
     *  within u of it and 2 n u^2 M, u the unit roundoff, 2^-53, a product below the normal doubles
     *  within 2^-1074 of its own; the slack is twice all of that and more, so that its own
     *  roundings leave it a bound.
     */
    class Shoelace
    {
    public:
      /**
       *  A run begun at a point, the one the areas are taken about.
       *
       *  @param origin The point.
       *  @param plain  Whether the area is worked out plainly, with its slack.
       */
      Shoelace(Point origin, bool plain) : m_Origin(origin), m_Last(origin), m_Plain(plain) {}

      /**
       *  The run taken on to a point, its side from the last point added.
       *
       *  @param point The point.
       */
      void Add(Point point)
      {
        Side(m_Last, point);
        m_Last = point;
      }

      /**
       *  The area, the run closed back to its first point: a cycle has two points or more, and a
       *  run of one would be 0 as a polygon of one point is. With it its slack, 0 where exact.
       *
       *  @return The area and its slack.
       */
      [[nodiscard]] auto Closed() -> std::array<double, 2>
      {
        Side(m_Last, m_Origin);
        if (!m_Plain)
        {
          return {m_Area + m_Error, 0.0};
        }
        return {m_Area / 2.0,
                (0x1p-52 * ((m_Sides + 8.0) * (m_Size / 2.0))) + (m_Sides * 0x1p-1070)};
      }

    private:
      /**
       *  A side's share added.
       *
       *  @param start Where it starts.
       *  @param end   Where it ends.
       */
      void Side(Point start, Point end)
      {
        const Vector p = start - m_Origin;
        const Vector q = end - m_Origin;
        if (m_Plain)
        {
          const double first   = p.X * q.Y;
          const double second  = p.Y * q.X;
          m_Area              += first - second;
          m_Size              += std::abs(first) + std::abs(second);
          m_Sides             += 1.0;
          return;
        }
        const Exact::Two first  = Exact::TwoProduct(p.X, q.Y);
        const Exact::Two second = Exact::TwoProduct(p.Y, q.X);
        Carry(first.Big / 2.0);
        Carry(-second.Big / 2.0);
        m_Error += (first.Small - second.Small) / 2.0;
      }

      /**
       *  A value added to the area, what rounding took kept apart.
       *
       *  @param value The value.
       */
      void Carry(double value)
      {
        const Exact::Two sum  = Exact::TwoSum(m_Area, value);
        m_Area                = sum.Big;
        m_Error              += sum.Small;
      }

      Point  m_Origin = {};    /**< The point areas are taken about.          */
      Point  m_Last   = {};    /**< The last point added.                     */
      double m_Area   = 0.0;   /**< The area so far, rounded; twice it plain. */
      double m_Error  = 0.0;   /**< The rounding errors so far.               */
      double m_Size   = 0.0;   /**< The products' magnitudes so far, summed.  */
      double m_Sides  = 0.0;   /**< How many sides so far.                    */
      bool   m_Plain  = false; /**< Whether worked out plainly.               */
    };

    /**
     *  The vertex a half starts from.
     *
     *  @param  arranged The arrangement.
     *  @param  half     The half.
     *  @return          The vertex.
     */
    auto Tail(const Arranged &arranged, std::size_t half) -> std::size_t
    {
      return arranged.Edges[half / 2].Ends.at(half % 2);
    }

    /**
     *  How a half leaves its vertex: a line's piece along the line's own way, an arc's along its
     *  tangent there and turning by the arc's curvature, either the other way back.
     *
     *  @param  arranged The arrangement.
     *  @param  half     The half.
     *  @return          Its spoke.
     */
    auto Spoked(const Arranged &arranged, std::size_t half) -> Spoke
    {
      const Arranged::Edge &edge   = arranged.Edges[half / 2];
      const bool            back   = half % 2 == 1;
      const Curve          &source = arranged.Sources[edge.Source];
      Spoke                 spoke  = {.Half = half};
      if (const auto *arc = std::get_if<Arc>(&source))
      {
        spoke.Way  = Derivative(*arc, back ? edge.Stretch.End : edge.Stretch.Start);
        spoke.Curl = (arc->Sweep.InDegrees() > 0.0 ? 1.0 : -1.0) / std::abs(arc->Radius);
      }
      else
      {
        /**
         *  A run of lines leaves its end back along its last.
         */
        const std::size_t last = back ? edge.Lines - 1 : 0;
        const Line       &line = std::get<Line>(arranged.Sources[edge.Source + last]);
        spoke.Way              = line.End - line.Start;
      }
      if (back)
      {
        spoke.Way  = -spoke.Way;
        spoke.Curl = -spoke.Curl;
      }
      return spoke;
    }

    /**
     *  Whether one spoke turns less to the left than another; by the halves where they turn alike.
     *
     *  @param  lhs One spoke.
     *  @param  rhs The other.
     *  @return     Whether it does.
     */
    auto Straighter(const Spoke &lhs, const Spoke &rhs) -> bool
    {
      if (lhs.Curl < rhs.Curl || rhs.Curl < lhs.Curl)
      {
        return lhs.Curl < rhs.Curl;
      }
      return lhs.Half < rhs.Half;
    }

    /**
     *  The spokes of a vertex counterclockwise: by the exact heading of their tangents, and where
     *  two tangents are within an angle of each other, by curvature, the one that turns more to the
     *  left after the other, as it lies left of it once they leave.
     *
     *  Where banded, two tangents are as near also within the root of the band times the difference
     *  of their curvatures: two curves that leave at an angle t and bend apart by k stay within t^2
     *  / 2k of each other until they cross at 2t / k, so within the band, twice the length, where
     *  t^2 is within twice the band times k; Intersect takes such a band for one hit, and they part
     *  on the sides their curvature, not their tangents, puts them on.
     *
     *  @param spokes The spokes, three or more.
     *  @param angle  The angle.
     *  @param band   Four lengths where banded, else 0.
     */
    void Fanned(std::span<Spoke> spokes, double angle, double band)
    {
      const std::size_t count = spokes.size();
      const auto        apart = [&spokes, count, angle, band](std::size_t index) -> bool
      {
        const Spoke &spoke = spokes[index];
        const Spoke &after = spokes[index + 1 < count ? index + 1 : 0];
        const double wide = band > 0.0
                              ? std::max(angle, std::sqrt(band * std::abs(spoke.Curl - after.Curl)))
                              : angle;
        return !(Dot(spoke.Way, after.Way) > 0.0)
            || std::abs(Cross(spoke.Way, after.Way)) > wide * Length(spoke.Way) * Length(after.Way);
      };
      std::size_t start = 0;
      while (start < count && !apart(start))
      {
        ++start;
      }
      const std::size_t after = start + 1 < count ? start + 1 : start + 1 - count;
      std::ranges::rotate(spokes, spokes.begin() + static_cast<std::ptrdiff_t>(after));
      std::size_t first = 0;
      for (std::size_t index = 0; index < count; ++index)
      {
        if (index + 1 < count && !apart(index))
        {
          continue;
        }
        if (index > first)
        {
          std::sort(spokes.begin() + static_cast<std::ptrdiff_t>(first),
                    spokes.begin() + static_cast<std::ptrdiff_t>(index + 1), Straighter);
        }
        first = index + 1;
      }
    }

    /**
     *  The half of the plane a way points into, counterclockwise from the x-axis's negative end, as
     *  Order counts them: below the axis, along it forward, above it, along it back.
     *
     *  @param  way The way.
     *  @return     0, 1, 2 or 3.
     */
    auto Halfplane(Vector way) -> int
    {
      if (way.Y < 0.0)
      {
        return 0;
      }
      if (way.Y > 0.0)
      {
        return 2;
      }
      return way.X < 0.0 ? 3 : 1;
    }

    /**
     *  The area each edge adds to the cycle of its first half with its chord: an arc's piece's, as
     *  Area gives it; a line's none.
     *
     *  @param  arranged The arrangement.
     *  @return          Each edge's area, along it.
     */
    auto Bulges(const Arranged &arranged) -> List<double>
    {
      List<double> bulges(arranged.Edges.size(), 0.0);
      for (std::size_t index = 0; index < arranged.Edges.size(); ++index)
      {
        if (std::holds_alternative<Arc>(arranged.Sources[arranged.Edges[index].Source]))
        {
          bulges[index] = Area(std::get<Arc>(Carved(arranged, index)));
        }
      }
      return bulges;
    }

    /**
     *  A cycle of halves as Cycled tours it: the cycle's number, what each edge adds as an arc's
     *  piece with its chord, the arcs' sum so far, and the shoelace of its points.
     */
    struct Circuit
    {
      std::size_t             Cycle  = 0;  /**< The cycle's number. */
      std::span<const double> Bulges = {}; /**< Each edge's area with its chord,
                                                or none where no curve is an arc. */
      double   Arcs = 0.0;                 /**< The arcs' sum so far.       */
      Shoelace Polygon;                    /**< The shoelace of its points. */
    };

    /**
     *  A cycle toured from its first half: each half given the cycle, its arc's area added, and the
     *  points the halves leave after the first added to the shoelace, which begins where the first
     *  leaves; a run of lines adds each of its joints, as its lines' halves would. The first half
     *  is never a run's: at a run's start, a vertex of two edges, the curve before it ends, and its
     *  last piece, of a lower number, is in every cycle the run is.
     *
     *  @param arranged The arrangement, its next halves made.
     *  @param first    The cycle's first half.
     *  @param circuit  The cycle, its shoelace begun.
     */
    void Encircled(Arranged &arranged, std::size_t first, Circuit &circuit)
    {
      Toured(arranged, first,
             [&arranged, first, &circuit](std::size_t half) -> void
             {
               arranged.Faces[half] = circuit.Cycle;
               if (!circuit.Bulges.empty())
               {
                 const double bulge  = circuit.Bulges[half / 2];
                 circuit.Arcs       += half % 2 == 1 ? -bulge : bulge;
               }
               if (half == first)
               {
                 return;
               }
               const Arranged::Edge &edge = arranged.Edges[half / 2];
               if (edge.Lines == 1)
               {
                 circuit.Polygon.Add(arranged.Points[edge.Ends.at(half % 2)]);
                 return;
               }
               for (std::size_t line = 0; line < edge.Lines; ++line)
               {
                 const std::size_t joint = half % 2 == 0 ? line : edge.Lines - line;
                 circuit.Polygon.Add(Beaded(arranged, half / 2, joint));
               }
             });
    }

    /**
     *  The cycles of a traced arrangement's halves, each half followed by its next, each with its
     *  signed area: its polygon's, through the vertices its halves leave, and its arcs' with their
     *  chords. A line's piece adds 0 or -0 to the second, as Area gives a line none: the sum, begun
     *  at 0, is never -0, so adding either changes no bit, and where every curve is a line nothing
     *  is added. A chained arrangement's are worked out plainly, each with its slack: the
     *  shoelace's, and twice the rounding of the last sum each way, the arcs' sum being the same.
     *
     *  @param arranged The arrangement.
     *  @param bulges   Each edge's area with its chord, or none where no curve is an arc.
     */
    void Cycled(Arranged &arranged, std::span<const double> bulges)
    {
      const std::size_t halves = 2 * arranged.Edges.size();
      arranged.Faces.assign(halves, halves);
      for (std::size_t first = 0; first < halves; ++first)
      {
        if (arranged.Faces[first] != halves)
        {
          continue;
        }
        Circuit circuit = {.Cycle  = arranged.Cycles.size(),
                           .Bulges = bulges,
                           .Polygon =
                               Shoelace{arranged.Points[Tail(arranged, first)], arranged.Chained}};
        Encircled(arranged, first, circuit);
        const auto [area, slack] = circuit.Polygon.Closed();
        const double sum         = area + circuit.Arcs;
        arranged.Cycles.push_back({.First = first,
                                   .Area  = sum,
                                   .Slack = arranged.Chained
                                              ? slack + (0x1p-51 * (std::abs(sum) + std::abs(area)))
                                              : 0.0});
      }
    }

    /**
     *  Each component's outer cycle, its least by area, made certain where the areas were worked
     *  out plainly: one whose area and slack lie below every other's area less its slack is least
     *  by the exact areas too, by these bounds, and the only one; in a component where any other
     *  comes nearer, every cycle's area is worked out exactly, as Cycled works it out where not
     *  chained, and the least taken from them, the first of the least, as Parted takes it.
     *
     *  @param arranged The arrangement, its outer cycles taken by the plain areas.
     *  @param bulges   Each edge's area with its chord, or none where no curve is an arc.
     */
    void Settled(Arranged &arranged, std::span<const double> bulges)
    {
      const std::size_t  cycles = arranged.Cycles.size();
      List<std::uint8_t> unsure(arranged.Outer.size(), 0);
      bool               any = false;
      for (std::size_t index = 0; index < cycles; ++index)
      {
        const std::size_t component = arranged.Cycles[index].Component;
        const std::size_t outer     = arranged.Outer[component];
        const bool        clear     = index == outer
                                   || arranged.Cycles[index].Area - arranged.Cycles[index].Slack
                                          > arranged.Cycles[outer].Area + arranged.Cycles[outer].Slack;
        if (!clear)
        {
          unsure[component] = 1;
          any               = true;
        }
      }
      if (!any)
      {
        return;
      }
      for (std::size_t index = 0; index < cycles; ++index)
      {
        Arranged::Cycle &cycle = arranged.Cycles[index];
        if (unsure[cycle.Component] != 0)
        {
          Circuit circuit = {.Cycle  = index,
                             .Bulges = bulges,
                             .Polygon =
                                 Shoelace{arranged.Points[Tail(arranged, cycle.First)], false}};
          Encircled(arranged, cycle.First, circuit);
          cycle.Area                      = circuit.Polygon.Closed()[0] + circuit.Arcs;
          cycle.Slack                     = 0.0;
          arranged.Outer[cycle.Component] = cycles;
        }
      }
      for (std::size_t index = 0; index < cycles; ++index)
      {
        const Arranged::Cycle &cycle = arranged.Cycles[index];
        std::size_t           &outer = arranged.Outer[cycle.Component];
        if (unsure[cycle.Component] != 0
            && (outer == cycles || cycle.Area < arranged.Cycles[outer].Area))
        {
          outer = index;
        }
      }
    }

    /**
     *  The components of a traced arrangement: each vertex an edge reaches in the component of its
     *  least, each cycle in its vertices', and each component's outer cycle the one of least area.
     *  The cycles of an edge's two halves are of one component, as are those of every half about a
     *  vertex, each following the one just clockwise of its twin's: so the cycles' sets, joined
     *  across each edge, are the components, and a vertex's is that of the first half leaving it.
     *  An edge between the two cycles the last one joined joins nothing more.
     *
     *  @param arranged The arrangement.
     *  @param firsts   The first two halves leaving each vertex, each past every half where there
     *                  is none.
     */
    void Parted(Arranged &arranged, std::span<const std::size_t> firsts)
    {
      const std::size_t          cycles = arranged.Cycles.size();
      Unions                     unions{cycles};
      std::array<std::size_t, 2> last = {cycles, cycles};
      for (std::size_t edge = 0; edge < arranged.Edges.size(); ++edge)
      {
        const std::array<std::size_t, 2> faces = {arranged.Faces[2 * edge],
                                                  arranged.Faces[(2 * edge) + 1]};
        if (faces != last)
        {
          unions.Join(faces[0], faces[1]);
          last = faces;
        }
      }
      const std::size_t halves = arranged.Faces.size();
      List<std::size_t> components(cycles, cycles);
      for (std::size_t vertex = 0; vertex < firsts.size() / 2; ++vertex)
      {
        const std::size_t leaving = firsts[2 * vertex];
        if (leaving == halves)
        {
          continue;
        }
        const std::size_t root = unions.Root(arranged.Faces[leaving]);
        if (components[root] == cycles)
        {
          components[root] = arranged.Samples.size();
          arranged.Samples.push_back(vertex);
        }
      }
      arranged.Outer.assign(arranged.Samples.size(), cycles);
      for (std::size_t index = 0; index < cycles; ++index)
      {
        Arranged::Cycle &cycle = arranged.Cycles[index];
        cycle.Component        = components[unions.Root(index)];
        std::size_t &outer     = arranged.Outer[cycle.Component];
        if (outer == cycles || cycle.Area < arranged.Cycles[outer].Area)
        {
          outer = index;
        }
      }
    }

    /**
     *  Each piece of an edge as a curve of its own, and its ends: the edge's piece, as Carved makes
     *  it, or each line of a run of lines in turn, between the points it passes.
     *
     *  @tparam VisitT   The kind of what is done with a piece.
     *  @param  arranged The arrangement.
     *  @param  edge     The edge.
     *  @param  visit    What is done with each piece and its two ends.
     */
    template <typename VisitT>
    void Pieces(const Arranged &arranged, std::size_t edge, const VisitT &visit)
    {
      const Arranged::Edge &piece = arranged.Edges[edge];
      if (piece.Lines == 1)
      {
        visit(Carved(arranged, edge),
              std::array{arranged.Points[piece.Ends[0]], arranged.Points[piece.Ends[1]]});
        return;
      }
      for (std::size_t line = 0; line < piece.Lines; ++line)
      {
        const std::array<Point, 2> ends = {Beaded(arranged, edge, line),
                                           Beaded(arranged, edge, line + 1)};
        visit(Curve{Line{.Start = ends[0], .End = ends[1]}}, ends);
      }
    }

    /**
     *  How many times each operand's loops wind about each component's least vertex, those of its
     *  own component left out: every other component's edges' ascents about it, each component
     *  passed over whose outer cycle's box does not hold the point.
     *
     *  @param  arranged The traced arrangement, of two components or more.
     *  @return          Each component's windings.
     */
    auto Based(const Arranged &arranged) -> List<Tally>
    {
      const std::size_t components = arranged.Outer.size();
      List<std::size_t> starts(components + 1, 0);
      for (std::size_t index = 0; index < arranged.Edges.size(); ++index)
      {
        ++starts[arranged.Cycles[arranged.Faces[2 * index]].Component + 1];
      }
      std::partial_sum(starts.begin(), starts.end(), starts.begin());
      List<std::size_t> edges(arranged.Edges.size());
      {
        List<std::size_t> filled(starts.begin(), starts.end() - 1);
        for (std::size_t index = 0; index < arranged.Edges.size(); ++index)
        {
          edges[filled[arranged.Cycles[arranged.Faces[2 * index]].Component]++] = index;
        }
      }
      List<Box> boxes(components);
      for (std::size_t component = 0; component < components; ++component)
      {
        Toured(arranged, arranged.Cycles[arranged.Outer[component]].First,
               [&arranged, &boxes, component](std::size_t half) -> void
               {
                 Pieces(
                     arranged, half / 2,
                     [&boxes, component](const Curve &shape, const std::array<Point, 2> &) -> void
                     { boxes[component] = Joined(boxes[component], Extents(shape)); });
               });
      }
      List<Tally> bases(components, Tally{});
      for (std::size_t component = 0; component < components; ++component)
      {
        const Point point = arranged.Points[arranged.Samples[component]];
        Tally       base  = {};
        for (std::size_t other = 0; other < components; ++other)
        {
          if (other == component || !Contains(boxes[other], point))
          {
            continue;
          }
          for (std::size_t slot = starts[other]; slot < starts[other + 1]; ++slot)
          {
            const Arranged::Edge &edge   = arranged.Edges[edges[slot]];
            int                   passes = 0;
            Pieces(arranged, edges[slot],
                   [&passes, point](const Curve &shape, const std::array<Point, 2> &ends) -> void
                   {
                     passes += static_cast<int>(
                         Ascents({.Shape = &shape, .From = ends[0], .To = ends[1]}, point));
                   });
            base = {base[0] + (passes * edge.Weighs[0]), base[1] + (passes * edge.Weighs[1])};
          }
        }
        bases[component] = base;
      }
      return bases;
    }

    /**
     *  How many halves leave each vertex, and the first two that do, in the order of the halves.
     *
     *  @param arranged The arrangement.
     *  @param degrees  How many halves leave each vertex, counted.
     *  @param firsts   The first two halves leaving each vertex, filled where there are any.
     */
    void Leaving(const Arranged &arranged, std::span<std::size_t> degrees,
                 std::span<std::size_t> firsts)
    {
      for (std::size_t edge = 0; edge < arranged.Edges.size(); ++edge)
      {
        const std::array<std::size_t, 2> &ends = arranged.Edges[edge].Ends;
        for (std::size_t side = 0; side < 2; ++side)
        {
          const std::size_t vertex = ends.at(side);
          const std::size_t degree = degrees[vertex]++;
          if (degree < 2)
          {
            firsts[(2 * vertex) + degree] = (2 * edge) + side;
          }
        }
      }
    }

    /**
     *  The next of each half about a vertex of one or two halves: about such a vertex the half just
     *  clockwise of each is the other, or itself, whichever way they lie.
     *
     *  @param  arranged The arrangement, its next halves sized.
     *  @param  degrees  How many halves leave each vertex.
     *  @param  firsts   The first two halves leaving each vertex.
     *  @return          How many halves leave the vertices of three or more, left to Thronged.
     */
    auto Join(Arranged &arranged, std::span<const std::size_t> degrees,
              std::span<const std::size_t> firsts) -> std::size_t
    {
      std::size_t busy = 0;
      for (std::size_t vertex = 0; vertex < degrees.size(); ++vertex)
      {
        const std::size_t degree = degrees[vertex];
        const std::size_t one    = firsts[2 * vertex];
        const std::size_t other  = firsts[(2 * vertex) + 1];
        if (degree == 1)
        {
          arranged.Next[one ^ 1U] = one;
        }
        if (degree == 2)
        {
          arranged.Next[one ^ 1U]   = other;
          arranged.Next[other ^ 1U] = one;
        }
        busy += degree > 2 ? degree : 0;
      }
      return busy;
    }

    /**
     *  The next of each half about a vertex of three halves or more: the halves of such vertices
     *  gathered vertex by vertex, by counting, each vertex's put in order as Ordered puts them,
     *  each half's next the one before its twin about the vertex, the last for the first.
     *
     *  @param arranged The arrangement, its next halves sized.
     *  @param degrees  How many halves leave each vertex.
     *  @param busy     How many halves leave the vertices of three or more.
     */
    void Thronged(Arranged &arranged, std::span<const std::size_t> degrees, std::size_t busy)
    {
      const std::size_t count = degrees.size();
      List<std::size_t> fans(count + 1, 0);
      for (std::size_t vertex = 0; vertex < count; ++vertex)
      {
        fans[vertex + 1] = fans[vertex] + (degrees[vertex] > 2 ? degrees[vertex] : 0);
      }
      List<std::size_t> around(busy);
      List<std::size_t> filled(fans.begin(), fans.end() - 1);
      for (std::size_t half = 0; half < 2 * arranged.Edges.size(); ++half)
      {
        const std::size_t vertex = Tail(arranged, half);
        if (degrees[vertex] > 2)
        {
          around[filled[vertex]++] = half;
        }
      }
      List<Spoke> spokes;
      for (std::size_t vertex = 0; vertex < count; ++vertex)
      {
        const std::size_t first = fans[vertex];
        const std::size_t end   = fans[vertex + 1];
        if (first == end)
        {
          continue;
        }
        Ordered(arranged, std::span<std::size_t>(around).subspan(first, end - first), spokes);
        for (std::size_t place = first; place < end; ++place)
        {
          arranged.Next[around[place] ^ 1U] = around[(place > first ? place : end) - 1];
        }
      }
    }

    /**
     *  The point halfway along an edge's piece: an arc's at the middle of its stretch, a line's
     *  between its vertices, as Carved draws it.
     *
     *  @param  arranged The arrangement.
     *  @param  edge     The edge, of one line or arc.
     *  @return          The point.
     */
    auto Pivot(const Arranged &arranged, std::size_t edge) -> Point
    {
      const Arranged::Edge &piece = arranged.Edges[edge];
      if (const auto *arc = std::get_if<Arc>(&arranged.Sources[piece.Source]))
      {
        return At(*arc, (piece.Stretch.Start / 2.0) + (piece.Stretch.End / 2.0));
      }
      return Lerp(arranged.Points[piece.Ends[0]], arranged.Points[piece.Ends[1]], 0.5);
    }

    /**
     *  Whether two edges between the same two vertices bound a face of their own: a half of the
     *  first followed by a half of the second, which is followed by it.
     *
     *  @param  arranged The traced arrangement.
     *  @param  pair     The two edges.
     *  @return          Whether they do.
     */
    auto Lensed(const Arranged &arranged, const std::array<std::size_t, 2> &pair) -> bool
    {
      return std::ranges::any_of(std::array{2 * pair[0], (2 * pair[0]) + 1},
                                 [&arranged, &pair](std::size_t half) -> bool
                                 {
                                   const std::size_t next = arranged.Next[half];
                                   return next / 2 == pair[1] && arranged.Next[next] == half;
                                 });
    }

    /**
     *  Whether two edges between the same two vertices are one, as Melded tells them: two pieces of
     *  lines; or two of which one is of an arc, halfway along within four lengths of each other,
     *  that bound no face of their own.
     *
     *  Neither is a run of lines: a run's lines are in no pair found but with their neighbours, and
     *  a curve with a piece at both its ends would be in a pair with its first and its last.
     *
     *  @param  arranged The traced arrangement.
     *  @param  pair     The two edges.
     *  @return          Whether they are.
     */
    auto Fusible(const Arranged &arranged, const std::array<std::size_t, 2> &pair) -> bool
    {
      if (std::holds_alternative<Line>(arranged.Sources[arranged.Edges[pair[0]].Source])
          && std::holds_alternative<Line>(arranged.Sources[arranged.Edges[pair[1]].Source]))
      {
        return true;
      }
      return Distance(Pivot(arranged, pair[0]), Pivot(arranged, pair[1]))
              <= 4.0 * arranged.Within.Length
          && !Lensed(arranged, pair);
    }

    /**
     *  The edge an edge is made one with, followed to the end: itself where it is made one with
     *  none.
     *
     *  @param  into Each edge's, or past the last edge for none.
     *  @param  edge The edge.
     *  @return      The edge at the end.
     */
    auto Fused(std::span<const std::size_t> into, std::size_t edge) -> std::size_t
    {
      while (into[edge] < into.size())
      {
        edge = into[edge];
      }
      return edge;
    }
  } // namespace

  auto Preceding(Vector lhs, Vector rhs) -> int
  {
    const int one   = Halfplane(lhs);
    const int other = Halfplane(rhs);
    if (one != other)
    {
      return one - other;
    }
    const std::strong_ordering order = Order(lhs, rhs);
    return static_cast<int>(order > 0) - static_cast<int>(order < 0);
  }

  void Ordered(const Arranged &arranged, std::span<std::size_t> around, List<Spoke> &spokes)
  {
    const double angle = std::max(arranged.Within.Radians, LeastAngle);
    spokes.clear();
    for (const std::size_t half : around)
    {
      spokes.push_back(Spoked(arranged, half));
    }
    std::ranges::sort(spokes,
                      [](const Spoke &lhs, const Spoke &rhs) -> bool
                      {
                        const int order = Preceding(lhs.Way, rhs.Way);
                        return order == 0 ? Straighter(lhs, rhs) : order < 0;
                      });
    Fanned(spokes, angle, arranged.Banded ? 4.0 * arranged.Within.Length : 0.0);
    for (std::size_t place = 0; place < around.size(); ++place)
    {
      around[place] = spokes[place].Half;
    }
  }

  /**
   *  About a vertex of one or two halves the half just clockwise of each is the other, or itself;
   *  only those of three or more are put in order.
   */
  void Trace(Arranged &arranged)
  {
    const std::size_t halves = 2 * arranged.Edges.size();
    const std::size_t count  = arranged.Points.size();
    arranged.Cycles.clear();
    arranged.Samples.clear();
    List<std::size_t> degrees(count, 0);
    List<std::size_t> firsts(2 * count, halves);
    Leaving(arranged, degrees, firsts);
    arranged.Next.resize(halves);
    const std::size_t busy = Join(arranged, degrees, firsts);
    if (busy > 0)
    {
      Thronged(arranged, degrees, busy);
    }
    const bool         arced = std::ranges::any_of(arranged.Sources, [](const Curve &source) -> bool
                                                   { return std::holds_alternative<Arc>(source); });
    const List<double> bulges = arced ? Bulges(arranged) : List<double>{};
    Cycled(arranged, bulges);
    Parted(arranged, firsts);
    if (arranged.Chained)
    {
      Settled(arranged, bulges);
    }
  }

  /**
   *  The halves of a connected graph, each followed by the one just clockwise of its twin about its
   *  end, however they are ordered about each vertex, trace the faces of a closed surface, of 2 -
   *  2g for g the surface's handles, 0 or more: so the sum over the components is twice their
   *  number just where every component's is 2, the vertices reached by an edge counted once.
   */
  auto Planar(const Arranged &arranged) -> bool
  {
    List<std::uint8_t> reached(arranged.Points.size(), 0);
    for (const Arranged::Edge &edge : arranged.Edges)
    {
      reached[edge.Ends[0]] = 1;
      reached[edge.Ends[1]] = 1;
    }
    const auto vertices = std::ranges::count(reached, std::uint8_t{1});
    const auto edges    = static_cast<std::ptrdiff_t>(arranged.Edges.size());
    const auto faces    = static_cast<std::ptrdiff_t>(arranged.Cycles.size());
    return vertices - edges + faces == 2 * static_cast<std::ptrdiff_t>(arranged.Outer.size());
  }

  /**
   *  The edges sorted by their two vertices, the lesser first, then by index, so that those between
   *  the same two lie in a run; each edge of a run made one with the first before it in the run
   *  that it is one with, and so with the first of them all that edge is made one with.
   */
  auto Melded(Arranged &arranged) -> bool
  {
    const std::size_t count = arranged.Edges.size();
    const auto        key   = [&arranged](std::size_t edge) -> std::array<std::size_t, 3>
    {
      const std::array<std::size_t, 2> &ends = arranged.Edges[edge].Ends;
      return {std::min(ends[0], ends[1]), std::max(ends[0], ends[1]), edge};
    };
    List<std::size_t> order(count);
    std::ranges::iota(order, std::size_t{0});
    std::ranges::sort(order, {}, key);
    List<std::size_t> into(count, count);
    bool              any   = false;
    std::size_t       first = 0;
    for (std::size_t place = 1; place < count; ++place)
    {
      const std::size_t edge = order[place];
      first = key(edge)[1] == key(order[first])[1] && key(edge)[0] == key(order[first])[0] ? first
                                                                                           : place;
      const auto before = std::span<const std::size_t>(order).subspan(first, place - first);
      const auto match  = std::ranges::find_if(before, [&arranged, edge](std::size_t other) -> bool
                                               { return Fusible(arranged, {other, edge}); });
      if (match != before.end())
      {
        into[edge] = *match;
        any        = true;
      }
    }
    if (!any)
    {
      return false;
    }
    for (std::size_t edge = 0; edge < count; ++edge)
    {
      if (into[edge] < count)
      {
        Arranged::Edge &kept  = arranged.Edges[Fused(into, edge)];
        const Tally     tally = arranged.Edges[edge].Weighs;
        const int       along = arranged.Edges[edge].Ends == kept.Ends ? 1 : -1;
        kept.Weighs = {kept.Weighs[0] + (along * tally[0]), kept.Weighs[1] + (along * tally[1])};
      }
    }
    std::size_t placed = 0;
    for (std::size_t edge = 0; edge < count; ++edge)
    {
      if (into[edge] == count)
      {
        arranged.Edges[placed++] = arranged.Edges[edge];
      }
    }
    arranged.Edges.resize(placed);
    Trace(arranged);
    return true;
  }

  auto Mended(Arranged arranged, bool pruned) -> Arranged
  {
    if (Planar(arranged))
    {
      return arranged;
    }
    if (arranged.Split)
    {
      if (pruned)
      {
        Prune(*arranged.Split);
      }
      Trace(*arranged.Split);
      if (Planar(*arranged.Split))
      {
        return std::move(*arranged.Split);
      }
    }
    for (const bool banded : {false, true})
    {
      Arranged mending = arranged;
      mending.Banded   = banded;
      Trace(mending);
      if (Planar(mending) || (Melded(mending) && Planar(mending)))
      {
        return mending;
      }
    }
    return arranged;
  }

  auto Windings(const Arranged &arranged) -> std::vector<Tally>
  {
    const std::size_t  components = arranged.Outer.size();
    const List<Tally>  bases      = components > 1 ? Based(arranged) : List<Tally>(components);
    std::vector<Tally> windings(arranged.Cycles.size());
    List<std::uint8_t> seen(arranged.Cycles.size(), 0);
    List<std::size_t>  queue;
    queue.reserve(arranged.Cycles.size());
    for (std::size_t component = 0; component < components; ++component)
    {
      std::size_t head = queue.size();
      queue.push_back(arranged.Outer[component]);
      windings[queue.back()] = bases[component];
      seen[queue.back()]     = 1;
      while (head < queue.size())
      {
        const std::size_t face = queue[head++];
        Toured(arranged, arranged.Cycles[face].First,
               [&](std::size_t half) -> void
               {
                 const std::size_t     across = arranged.Faces[half ^ 1U];
                 const Arranged::Edge &edge   = arranged.Edges[half / 2];
                 const int             way    = half % 2 == 0 ? 1 : -1;
                 if (seen[across] == 0)
                 {
                   seen[across]     = 1;
                   windings[across] = {windings[face][0] - (way * edge.Weighs[0]),
                                       windings[face][1] - (way * edge.Weighs[1])};
                   queue.push_back(across);
                 }
               });
      }
    }
    return windings;
  }
} // namespace Geoxy
