/***************************************************************************************************
 *  @file      source/weave/found.hpp
 *  @brief     What weaving an arrangement finds: each curve's notches, the stretches two curves
 *             share, and the points they are told by.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#ifndef GEOXY_WEAVE_FOUND_HPP
#define GEOXY_WEAVE_FOUND_HPP

#include "real.hpp"
#include "scratch.hpp"

#include "geoxy/angle.hpp"
#include "geoxy/arc.hpp"
#include "geoxy/box.hpp"
#include "geoxy/curve.hpp"
#include "geoxy/error.hpp"
#include "geoxy/extents.hpp"
#include "geoxy/intersect.hpp"
#include "geoxy/line.hpp"
#include "geoxy/point.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <utility>
#include <variant>
#include <vector>

namespace Geoxy::Weave
{
  /**
   *  The fault of an arrangement whose meetings ran out of budget, by how many did: none, or some.
   *  Lines and arcs never run out; the kinds still to come may.
   */
  inline constexpr std::array<std::optional<Error>, 2> Spending = {std::nullopt, Error::Spent};

  /**
   *  The start angle, in degrees, up to which an arc's box is known to lie within twice the sum of
   *  its center's largest coordinate and its radius, Extents's pads being far smaller there.
   */
  inline constexpr double PlainTurn = 1e6;

  /**
   *  How far, in degrees, an arc's sweep is taken to pass a quarter turn it only nears: far more
   *  than its angles in degrees are off by under PlainTurn.
   */
  inline constexpr double TurnSlack = 1e-6;

  /**
   *  Two curves, by their indices, the first first.
   */
  using Couple = std::array<std::size_t, 2>;

  /**
   *  The kinds of pairs of lines and arcs, in the order IntersectAll meets them: a level line and
   *  an upright one, two other lines, a line and an arc, two arcs.
   */
  inline constexpr std::size_t Kinds = 4;

  /**
   *  The kind of a pair of each two stances: level, upright, other line, arc.
   */
  inline constexpr std::array<std::array<std::uint8_t, Kinds>, Kinds> Couplings = {{
      {1, 0, 1, 2},
      {0, 1, 1, 2},
      {1, 1, 1, 2},
      {2, 2, 2, 3},
  }};

  /**
   *  A cut of a curve: where along it, and the vertex there.
   */
  struct Notch
  {
    std::size_t Source = 0;   /**< The curve.                   */
    double      Along  = 0.0; /**< The curve's parameter there. */
    std::size_t Vertex = 0;   /**< The vertex.                  */
  };

  /**
   *  A stretch two curves share, as Intersect gave it.
   */
  struct Common
  {
    std::size_t First  = 0;  /**< The first curve.              */
    std::size_t Second = 0;  /**< The second.                   */
    Overlap     Shared = {}; /**< The stretch, on each of them. */
  };

  /**
   *  What meeting the curves found: the vertices, each curve's cuts, and the stretches shared, the
   *  vertices and cuts in the thread's scratch memory. Each curve's ends are its first two
   *  vertices, and an arc's middle past a half turn its third, all of them before any other vertex;
   *  the cuts there, at 0, 1 and 0.5 along the curve, are held by its first vertex alone, and the
   *  cuts kept are the rest. A lone joint is no vertex: a line next to one has only its other
   *  end's, or none.
   */
  struct Findings
  {
    List<Point>       Points = {};    /**< The vertices, in order. */
    List<std::size_t> Twins  = {};    /**< Each vertex's earlier one at the same
                                           point, or one past the last vertex. */
    List<std::size_t> Firsts = {};    /**< Each curve's first vertex, and past
                                           the last, the first that is no end. */
    List<Notch>         Notches = {}; /**< The curves' other cuts.             */
    std::vector<Common> Commons = {}; /**< The stretches shared.               */
    std::size_t         Spent   = 0;  /**< How many meetings' budgets ran out. */
    std::size_t         Count   = 0;  /**< How many curves.                    */
    List<std::uint8_t>  Strung  = {}; /**< Whether each curve's end and the next's
                                           start are a lone joint, as Lonely tells
                                           it, of which neither is a vertex; none
                                           where not chained.                    */
  };

  /**
   *  A curve's cuts at its ends, and at its middle where it has one, in order along it: at 0, 0.5
   *  and 1.
   *
   *  @param  found The findings.
   *  @param  curve The curve's index.
   *  @return       The cuts, two or three, and how many.
   */
  inline auto Ended(const Findings &found, std::size_t curve)
      -> std::pair<std::array<Notch, 3>, std::size_t>
  {
    const std::size_t first = found.Firsts[curve];
    if (found.Firsts[curve + 1] - first == 2)
    {
      return {{Notch{.Source = curve, .Along = 0.0, .Vertex = first},
               Notch{.Source = curve, .Along = 1.0, .Vertex = first + 1}},
              2};
    }
    return {{Notch{.Source = curve, .Along = 0.0, .Vertex = first},
             Notch{.Source = curve, .Along = 0.5, .Vertex = first + 2},
             Notch{.Source = curve, .Along = 1.0, .Vertex = first + 1}},
            3};
  }

  /**
   *  Whether two points are the same, coordinate by coordinate.
   *
   *  @param  first  One point.
   *  @param  second The other.
   *  @return        Whether they are.
   */
  inline auto Twinned(Point first, Point second) -> bool
  {
    return Exactly(first.X, second.X) && Exactly(first.Y, second.Y);
  }

  /**
   *  Two magnitudes of a line or an arc: the largest of its box, as Extents gives the box, a line's
   *  ends'; or none, 0, for an arc whose box lies within a bound so small that 16 steps of a double
   *  there are within a length: within twice its center's largest coordinate and its radius
   *  together, as the box of an arc of a start angle under PlainTurn and of up to a turn is. Then
   *  the largest it may reach as Intersect sizes it: a line's ends', an arc's center's largest
   *  coordinate and its radius together.
   *
   *  @param  source The line or arc.
   *  @param  length The length.
   *  @return        The two magnitudes.
   */
  inline auto Sized(const Curve &source, double length) -> std::array<double, 2>
  {
    if (const auto *line = std::get_if<Line>(&source))
    {
      const double size = std::max({std::abs(line->Start.X), std::abs(line->Start.Y),
                                    std::abs(line->End.X), std::abs(line->End.Y)});
      return {size, size};
    }
    const Arc   &arc = std::get<Arc>(source);
    const double reach =
        std::max(std::abs(arc.Center.X), std::abs(arc.Center.Y)) + std::abs(arc.Radius);
    if (std::abs(arc.Start.InDegrees()) <= PlainTurn && 16.0 * Above(2.0 * reach) <= length)
    {
      return {0.0, reach};
    }
    const Box box = Extents(arc);
    return {std::max({std::abs(box.Min.X), std::abs(box.Min.Y), std::abs(box.Max.X),
                      std::abs(box.Max.Y)}),
            reach};
  }
} // namespace Geoxy::Weave

#endif /* GEOXY_WEAVE_FOUND_HPP */
