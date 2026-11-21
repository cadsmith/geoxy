/***************************************************************************************************
 *  @file      source/met/slight.cpp
 *  @brief     A polyline's remote arcs, too slight to be worked out through their centers, met
 *             through their chords and told again on the arcs themselves.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "met/slight.hpp"

#include "bow.hpp"
#include "crossing.hpp"
#include "met/judge.hpp"
#include "met/parts.hpp"
#include "real.hpp"

#include "geoxy/at.hpp"
#include "geoxy/box.hpp"
#include "geoxy/curve.hpp"
#include "geoxy/derivative.hpp"
#include "geoxy/extents.hpp"
#include "geoxy/intersect.hpp"
#include "geoxy/line.hpp"
#include "geoxy/point.hpp"
#include "geoxy/polyline.hpp"
#include "geoxy/segment.hpp"
#include "geoxy/tolerance.hpp"
#include "geoxy/vector.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <optional>
#include <utility>
#include <variant>

namespace Geoxy::Met
{
  namespace
  {
    /**
     *  A polyline's segment at a place, where it is a remote arc.
     */
    struct Bowing
    {
      double      Bulge = 0.0; /**< Its bulge.                     */
      double      Chord = 0.0; /**< Its chord's length.            */
      double      Part  = 0.0; /**< How far along it the place is. */
      std::size_t Index = 0;   /**< Which segment it is.           */

      /**
       *  How far the arc stands off its chord at the most: half its bulge times the chord's length,
       *  to within a hundred millionth of that.
       *
       *  @return How far.
       */
      [[nodiscard]] auto Sagitta() const -> double { return std::abs(Bulge) * Chord / 2.0; }
    };

    /**
     *  A polyline's segment at a place, where the segment is a remote arc.
     *
     *  @pre             The polyline has a segment.
     *  @param  polyline The polyline.
     *  @param  along    The place.
     *  @return          The arc; nothing where the segment there is none.
     */
    auto Bowed(const Polyline &polyline, double along) -> std::optional<Bowing>
    {
      const auto    last   = static_cast<double>(Count(polyline) - 1);
      const double  whole  = std::clamp(std::floor(along), 0.0, last);
      const auto    index  = static_cast<std::size_t>(whole);
      const Vertex &vertex = polyline.Vertices[index];
      const Point   next   = polyline.Vertices[(index + 1) % polyline.Vertices.size()].Position;
      if (!Remote(vertex.Position, next, vertex.Bulge))
      {
        return std::nullopt;
      }
      return Bowing{.Bulge = vertex.Bulge,
                    .Chord = Distance(vertex.Position, next),
                    .Part  = along - whole,
                    .Index = index};
    }

    /**
     *  A polyline with its remote arcs within a length of their chords made chords.
     *
     *  @param  curve  The curve.
     *  @param  length The length; past every length for every remote arc.
     *  @return        What stands for it; nothing where the curve is no polyline or has no such
     *                 arc.
     */
    auto Chorded(const Curve &curve, double length) -> std::optional<Standing>
    {
      const auto             *polyline = std::get_if<Polyline>(&curve);
      const std::size_t       count    = polyline != nullptr ? Count(*polyline) : 0;
      std::optional<Standing> out;
      for (std::size_t index = 0; index < count; ++index)
      {
        const Vertex &vertex = polyline->Vertices[index];
        const Point   next   = polyline->Vertices[(index + 1) % polyline->Vertices.size()].Position;
        if (!Remote(vertex.Position, next, vertex.Bulge))
        {
          continue;
        }
        const double sagitta = std::abs(vertex.Bulge) * Distance(vertex.Position, next) / 2.0;
        if (sagitta <= length)
        {
          out = out ? std::move(out) : std::optional<Standing>{Standing{.Shape = *polyline}};
          std::get<Polyline>(out->Shape).Vertices[index].Bulge = 0.0;
          out->Deepest                                         = std::max(out->Deepest, sagitta);
        }
      }
      return out;
    }

    /**
     *  Two curves, what stands for each, and the meeting being told.
     */
    struct Telling
    {
      std::array<const Curve *, 2>    Given   = {};  /**< The curves themselves.         */
      std::array<const Standing *, 2> Stand   = {};  /**< What stands for each, or none. */
      double                          Length  = 0.0; /**< The length used.               */
      double                          Radians = 0.0; /**< The tolerance's angle.         */
    };

    /**
     *  The curves a meeting is asked of: what stands for each, or the curve itself.
     *
     *  @param  telling The curves and what stands for each.
     *  @return         Them.
     */
    auto Shapes(const Telling &telling) -> Curves
    {
      return {telling.Stand[0] != nullptr ? telling.Stand[0]->Shape : *telling.Given[0],
              telling.Stand[1] != nullptr ? telling.Stand[1]->Shape : *telling.Given[1]};
    }

    /**
     *  A hit's place on a curve, where it lies on a chord that stands for a remote arc.
     *
     *  @param  telling The curves.
     *  @param  side    Which curve.
     *  @param  along   The place.
     *  @return         The arc; nothing where the curve has no stand-in or the place is on no
     *                  chord.
     */
    auto Chordal(const Telling &telling, std::size_t side, double along) -> std::optional<Bowing>
    {
      if (telling.Stand.at(side) == nullptr)
      {
        return std::nullopt;
      }
      const std::optional<Bowing> bow   = Bowed(std::get<Polyline>(*telling.Given.at(side)), along);
      const auto                 &stood = std::get<Polyline>(telling.Stand.at(side)->Shape);
      return bow && Zero(stood.Vertices[bow->Index].Bulge) ? bow : std::nullopt;
    }

    /**
     *  Whether a remote arc stands off its chord by more than the length used.
     *
     *  @param  telling The curves.
     *  @param  bow     The arc, or none.
     *  @return         Whether it does; never for none.
     */
    auto Deep(const Telling &telling, const std::optional<Bowing> &bow) -> bool
    {
      return bow && bow->Sagitta() > telling.Length;
    }

    /**
     *  Whether a crossing on a chord is one crossing of the arc the chord stands for, and nothing
     *  else of the arc is near: the curves cross there more than sixteen times steeper than the
     *  arcs turn, and the crossing lies farther from the arcs' ends than the arcs' standing off and
     *  the width looked within may move it along them.
     *
     *  @param  telling The curves.
     *  @param  hit     The crossing, as met on the chords.
     *  @param  width   The width looked within: the widened tolerance's length.
     *  @return         Whether it is.
     */
    auto Plain(const Telling &telling, const Hit &hit, double width) -> bool
    {
      const std::array<std::optional<Bowing>, 2> bows = {Chordal(telling, 0, hit.First.Along),
                                                         Chordal(telling, 1, hit.Second.Along)};
      const Vector way   = Derivative(*telling.Given[0], hit.First.Along);
      const Vector back  = Derivative(*telling.Given[1], hit.Second.Along);
      const double speed = Length(way) * Length(back);
      const double sine  = speed > 0.0 ? std::abs(Cross(way, back)) / speed : 0.0;
      double       turn  = 0.0;
      double       stray = width;
      for (const std::optional<Bowing> &bow : bows)
      {
        turn  += bow ? std::abs(bow->Bulge) : 0.0;
        stray += bow ? bow->Sagitta() : 0.0;
      }
      if (hit.Kind != Contact::Cross || !(sine > 16.0 * turn))
      {
        return false;
      }
      const double reach = 2.0 * stray / sine;
      return std::ranges::all_of(bows,
                                 [reach](const std::optional<Bowing> &bow) -> bool
                                 {
                                   return !bow
                                       || (bow->Part * bow->Chord > reach
                                           && (1.0 - bow->Part) * bow->Chord > reach);
                                 });
    }

    /**
     *  One curve's remote arc that stands off its chord by more than the length used, and what its
     *  chord is to be weighed against.
     */
    struct Weighing
    {
      std::size_t Side  = 0;   /**< Which curve the arc is of. */
      std::size_t Index = 0;   /**< Which segment of it.       */
      double      Width = 0.0; /**< The width looked within.   */
    };

    /**
     *  Whether what the chords meet in tells what one arc meets in, where it stands off its chord
     *  by more than the length used: its chord alone, met with what stands for the other curve
     *  within a width that takes the arc in, shares no stretch with it and is hit only at
     *  crossings, each plain. The chord is met alone so that nothing beside it is told in its
     *  place: next to another segment, a long shallow approach is told as one hit at the joint.
     *  Where the chord alone is so met, the chords met together are hit on that segment at those
     *  crossings and nowhere else, each of which a polish then puts on the arc.
     *
     *  @param  telling  The curves.
     *  @param  weighing The arc, what the chords meet in, and the width.
     *  @return          Whether it does.
     */
    auto Certain(const Telling &telling, const Weighing &weighing) -> bool
    {
      const std::size_t side  = weighing.Side;
      const auto       &given = std::get<Polyline>(*telling.Given.at(side));
      const Point       start = given.Vertices[weighing.Index].Position;
      const Point       end = given.Vertices[(weighing.Index + 1) % given.Vertices.size()].Position;
      const Curve       chord  = Line{.Start = start, .End = end};
      const Curves      shapes = Shapes(telling);
      thread_local Intersections alone;
      Intersect(Curves{chord, shapes.at(1 - side)},
                {.Length = weighing.Width, .Radians = telling.Radians}, alone);
      const auto placed = [&weighing, side](Hit hit) -> Hit
      {
        hit.First.Along += static_cast<double>(weighing.Index);
        if (side == 1)
        {
          std::swap(hit.First, hit.Second);
        }
        return hit;
      };
      return alone.Overlaps.empty()
          && !alone.Spent
          && std::ranges::all_of(alone.Hits, [&](const Hit &hit) -> bool
                                 { return Plain(telling, placed(hit), weighing.Width); });
    }

    /**
     *  Whether what the chords meet in tells what the arcs meet in: every arc that stands off its
     *  chord by more than the length used, and whose chord's box, grown by the width that takes the
     *  arcs in, reaches the other curve's, is told so, as Certain tells one.
     *
     *  @param  telling The curves.
     *  @param  width   The width.
     *  @return         Whether it does.
     */
    auto Certain(const Telling &telling, double width) -> bool
    {
      int weighed = 0;
      for (const std::size_t side : {0UZ, 1UZ})
      {
        const auto *polyline = std::get_if<Polyline>(telling.Given.at(side));
        if (polyline == nullptr || telling.Stand.at(side) == nullptr)
        {
          continue;
        }
        const Box         other = Extents(Shapes(telling).at(1 - side));
        const std::size_t count = Count(*polyline);
        for (std::size_t index = 0; index < count; ++index)
        {
          const std::optional<Bowing> bow =
              Chordal(telling, side, static_cast<double>(index) + 0.5);
          if (!Deep(telling, bow))
          {
            continue;
          }
          const Point start = polyline->Vertices[index].Position;
          const Point end   = polyline->Vertices[(index + 1) % polyline->Vertices.size()].Position;
          const bool  apart = std::min(start.X, end.X) - width > other.Max.X
                           || std::max(start.X, end.X) + width < other.Min.X
                           || std::min(start.Y, end.Y) - width > other.Max.Y
                           || std::max(start.Y, end.Y) + width < other.Min.Y;
          /**
           *  No more than sixteen chords are weighed: past that the thorough way is the faster.
           */
          weighed += apart ? 0 : 1;
          if (!apart
              && (weighed > 16
                  || !Certain(telling, {.Side = side, .Index = index, .Width = width})))
          {
            return false;
          }
        }
      }
      return true;
    }

    /**
     *  Where the curves themselves cross near where their chords do, polished in double-double from
     *  places brought near by Newton's method in doubles, as a candidate's are: the chords'
     *  crossing lies as far from the arcs' as the arcs stand off over the crossing's sine, farther
     *  than a polish takes for its own.
     *
     *  @param  telling The curves.
     *  @param  places  Where the chords cross, on each curve.
     *  @return         The crossing; nothing where it is not found so.
     */
    auto Polished(const Telling &telling, std::array<double, 2> places)
        -> std::optional<Crossings::Crossed>
    {
      const Local near =
          Neared(*telling.Given[0], *telling.Given[1], {.First = places[0], .Second = places[1]});
      places = {near.First, near.Second};
      return Crossings::Polished(telling.Given, places, At(*telling.Given[0], places[0]));
    }

    /**
     *  Each hit on a chord told on the arc the chord stands for: a crossing put where the curves
     *  themselves cross, polished in double-double; any other hit's point on the chord made the
     *  arc's own at its place.
     *
     *  @param telling The curves.
     *  @param into    The meeting, as met on the chords.
     */
    void Retold(const Telling &telling, Intersections &into)
    {
      for (Hit &hit : into.Hits)
      {
        const std::array<bool, 2>   stand  = {Chordal(telling, 0, hit.First.Along).has_value(),
                                              Chordal(telling, 1, hit.Second.Along).has_value()};
        const std::array<double, 2> places = {hit.First.Along, hit.Second.Along};
        const std::optional<Crossings::Crossed> exact =
            hit.Kind == Contact::Cross && (stand[0] || stand[1]) ? Polished(telling, places)
                                                                 : std::nullopt;
        if (exact)
        {
          hit = {.At     = exact->At,
                 .First  = {.At = exact->On[0], .Along = exact->Along[0]},
                 .Second = {.At = exact->On[1], .Along = exact->Along[1]},
                 .Kind   = hit.Kind};
          continue;
        }
        hit.First.At  = stand[0] ? At(*telling.Given[0], hit.First.Along) : hit.First.At;
        hit.Second.At = stand[1] ? At(*telling.Given[1], hit.Second.Along) : hit.Second.At;
      }
    }

    /**
     *  What a meeting under way knows of its curves' remote arcs.
     */
    struct Knowing
    {
      bool Declined = false; /**< Whether the thorough way is meeting curves with remote arcs, asked
                                  by Slightly: it is then left to meet them as they are.         */
      bool Vetted = false;   /**< Whether the curves being met are known to have none. */
    };

    /**
     *  What the meeting under way on this thread knows.
     *
     *  @return It.
     */
    auto Known() -> Knowing &
    {
      thread_local Knowing known;
      return known;
    }
  } // namespace

  auto Slightish(const Polyline &polyline) -> bool
  {
    const std::size_t count = Count(polyline);
    for (std::size_t index = 0; index < count; ++index)
    {
      const Vertex &vertex = polyline.Vertices[index];
      const Point   next   = polyline.Vertices[(index + 1) % polyline.Vertices.size()].Position;
      if (Remote(vertex.Position, next, vertex.Bulge))
      {
        return true;
      }
    }
    return false;
  }

  auto Slightish(const Curve &curve) -> bool
  {
    const auto *polyline = std::get_if<Polyline>(&curve);
    return polyline != nullptr && Slightish(*polyline);
  }

  auto Remotely(const Curve &first, const Curve &second) -> bool
  {
    const bool strung =
        std::holds_alternative<Polyline>(first) || std::holds_alternative<Polyline>(second);
    return strung && !Known().Vetted && (Slightish(first) || Slightish(second));
  }

  Vetting::Vetting() : m_Before(Known().Vetted)
  {
    Known().Vetted = !Known().Declined;
  }

  Vetting::~Vetting()
  {
    Known().Vetted = m_Before;
  }

  auto Stood(const Curve &curve) -> std::optional<Standing>
  {
    return Chorded(curve, std::numeric_limits<double>::infinity());
  }

  auto Slightly(const Curves &curves, const Chording &asked, Intersections &into) -> bool
  {
    if (Known().Declined)
    {
      return false;
    }
    const std::array<std::optional<Standing>, 2> made = {
        asked.Kept[0] != nullptr ? std::nullopt : Stood(curves[0]),
        asked.Kept[1] != nullptr ? std::nullopt : Stood(curves[1])};
    if (asked.Kept[0] == nullptr && asked.Kept[1] == nullptr && !made[0] && !made[1])
    {
      return false;
    }
    const auto held = [](const std::optional<Standing> &stood) -> const Standing *
    {
      return stood ? &*stood : nullptr;
    };
    Telling      telling = {.Given   = {&curves[0].get(), &curves[1].get()},
                            .Stand   = {asked.Kept[0] != nullptr ? asked.Kept[0] : held(made[0]),
                                        asked.Kept[1] != nullptr ? asked.Kept[1] : held(made[1])},
                            .Length  = asked.Length,
                            .Radians = asked.Sane.Radians};
    const double deepest = std::max(telling.Stand[0] != nullptr ? telling.Stand[0]->Deepest : 0.0,
                                    telling.Stand[1] != nullptr ? telling.Stand[1]->Deepest : 0.0);
    asked.Intersect(asked.Context, Shapes(telling), asked.Sane, into);
    const bool sure = !(deepest > asked.Length) || Certain(telling, asked.Length + deepest);
    if (sure)
    {
      Retold(telling, into);
      return true;
    }
    /**
     *  The thorough way, each arc that stands off its chord by more than the length used a piece,
     *  the others still their chords.
     */
    const std::array<std::optional<Standing>, 2> shallow = {Chorded(curves[0], asked.Length),
                                                            Chorded(curves[1], asked.Length)};
    telling.Stand                                        = {held(shallow[0]), held(shallow[1])};
    Known().Declined                                     = true;
    asked.Deep(Shapes(telling), asked.Sane, into);
    Known().Declined = false;
    Retold(telling, into);
    return true;
  }
} // namespace Geoxy::Met
