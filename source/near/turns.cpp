/***************************************************************************************************
 *  @file      source/near/turns.cpp
 *  @brief     The place of an arc nearest a point, told by its angles where they tell it, and by
 *             its feet and ends where not.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "near/turns.hpp"

#include "elliptic.hpp"
#include "near/feet.hpp"
#include "near/plumb.hpp"
#include "real.hpp"
#include "straight.hpp"

#include "geoxy/angle.hpp"
#include "geoxy/arc.hpp"
#include "geoxy/at.hpp"
#include "geoxy/ellipse.hpp"
#include "geoxy/line.hpp"
#include "geoxy/place.hpp"
#include "geoxy/point.hpp"
#include "geoxy/polyline.hpp"
#include "geoxy/segment.hpp"
#include "geoxy/vector.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <initializer_list>
#include <optional>
#include <span>
#include <utility>
#include <variant>

namespace Geoxy::Near
{
  namespace
  {
    /**
     *  The parameter at which an arc swept round no more than once passes an angle of its circle,
     *  worked out as Passes works it out, to the bit, and how many times it passes it.
     */
    struct Pass
    {
      double Along = 0.0; /**< The parameter of the first pass. */
      double Count = 0.0; /**< How many passes, 0 to 2.         */
    };

    /**
     *  A whole turn in radians, as Turn.InRadians() gives it: the double nearest 2 pi, twice the
     *  double nearest pi.
     */
    constexpr double FullTurn = 2.0 * Pi;

    /**
     *  An arc's start and sweep in radians, as InRadians gives them.
     */
    struct Swing
    {
      double Start = 0.0; /**< The start. */
      double Sweep = 0.0; /**< The sweep. */
    };

    /**
     *  Where an arc passes an angle of its circle, as Passes finds it but with no list: the first
     *  turn's parameter, to the bit, and how many turns pass it.
     *
     *  @param  swing The arc's start and sweep in radians.
     *  @param  angle The angle, in radians.
     *  @return       The pass.
     */
    auto PassedOnce(Swing swing, double angle) -> Pass
    {
      const double turned = angle - swing.Start;
      const double first  = std::ceil((std::min(0.0, swing.Sweep) - turned) / FullTurn);
      const double last   = std::floor((std::max(0.0, swing.Sweep) - turned) / FullTurn);
      return {.Along = std::clamp((turned + (first * FullTurn)) / swing.Sweep, 0.0, 1.0),
              .Count = std::max(0.0, std::min(last - first + 1.0, MostTurns))};
    }

    /**
     *  An angle less the whole turns nearest it.
     *
     *  @param  radians The angle.
     *  @return         What is left, within a half turn and a hair of nothing.
     */
    auto TurnsOff(double radians) -> double
    {
      return radians - (FullTurn * std::rint(radians / FullTurn));
    }

    /**
     *  The place of an arc swept round once at most nearest a point, where the angles alone tell
     *  it, with no foot listed, one heading and one point worked out.
     *
     *  The places weighed are the arc's start, its end and its feet, the near one where the way
     *  from the center to the point, or for a negative radius the way back, crosses it, and the far
     *  one a half turn round; each is as far from the point as its angle from that way, a, is wide,
     *  the farther the wider (CloserBy). The angles are worked out in doubles, within 2^-40 (4 +
     *  |start| + |sweep|) radians of the angles At turns to, some 2^9 times what their roundings
     *  come to; At's points and their distances are within 2^-48 (|center| + d + 2 r) of the
     *  circle's, some 2^5 times what theirs come to. Where the nearest and the next differ by more
     *  than twice that, a foot not by an end, the nearest is the only place whose distance is least
     *  as Distance works it out, and so the one Nearest weighs out of all of them. A foot so far
     *  from the ends is passed once, where Passes finds it on the correctly rounded heading, which
     *  the angles are worked from; its point is the circle's toward the point (Dropped), as
     *  Perpendiculars gives it, the double nearest the exact foot.
     *
     *  @param  arc   The arc, swept round once at most.
     *  @param  point The point.
     *  @return       The place; nothing where the angles do not tell it, a foot by an end, the
     *                point at the center, the arc of no sweep or radius or swept round more than
     *                once.
     */
    auto ByAngles(const Arc &arc, Point point) -> std::optional<Place>
    {
      const Swing  swing  = {.Start = arc.Start.InRadians(), .Sweep = arc.Sweep.InRadians()};
      const Vector away   = point - arc.Center;
      const double reach  = Length(away);
      const double radius = std::abs(arc.Radius);
      if (!(reach > 0.0)
          || !(radius > 0.0)
          || Zero(swing.Sweep)
          || !(std::abs(swing.Sweep) <= FullTurn))
      {
        return std::nullopt;
      }
      const double toward = Heading(away).InRadians();
      const double way    = arc.Radius < 0.0 ? toward + Pi : toward;
      const double slack  = 0x1p-40 * (4.0 + std::abs(swing.Start) + std::abs(swing.Sweep));
      const double error =
          0x1p-48 * (std::abs(arc.Center.X) + std::abs(arc.Center.Y) + reach + (2.0 * radius));
      const double behind = TurnsOff(swing.Start - way);
      const double first  = std::abs(behind);
      const double last   = std::abs(TurnsOff(behind + swing.Sweep));
      if (!(std::min(first, last) > slack))
      {
        return std::nullopt;
      }
      /**
       *  The near foot is on the arc where the way lies within its sweep, at least the slack from
       *  either end.
       */
      double ahead = -behind;
      if (swing.Sweep > 0.0 && ahead < 0.0)
      {
        ahead += FullTurn;
      }
      if (swing.Sweep < 0.0 && ahead > 0.0)
      {
        ahead -= FullTurn;
      }
      const bool   foot    = std::abs(ahead) < std::abs(swing.Sweep);
      const double nearest = foot ? 0.0 : std::min(first, last);
      const double next    = foot ? std::min(first, last) : std::min(std::max(first, last), Pi);
      if (!(CloserBy({.X = reach, .Y = radius}, nearest + slack, next - slack)
            > 2.0 * error * (1.0 + 0x1p-20)))
      {
        return std::nullopt;
      }
      if (foot)
      {
        return Place{.At    = Dropped(arc.Center, radius, point),
                     .Along = PassedOnce(swing, way).Along};
      }
      const double along = first < last ? 0.0 : 1.0;
      return Place{.At = At(arc, along), .Along = along};
    }

    /**
     *  An arc's start or end as At gives it, or as given where it is.
     *
     *  @param  arc   The arc.
     *  @param  along 0 for its start, 1 for its end.
     *  @param  ends  Its start and end, as At gives them, or none.
     *  @return       The point.
     */
    auto Pointing(const Arc &arc, double along, const std::array<Point, 2> *ends) -> Point
    {
      if (ends != nullptr)
      {
        return (*ends).at(along > 0.0 ? 1 : 0);
      }
      return At(arc, along);
    }

    /**
     *  The feet of an arc swept no more than once round, as Perpendiculars gives them, with no
     *  list: Passes's parameters for the way from its center to the point and the way back, by the
     *  same steps, each way passed at most twice as the sweep is no more than a turn, each at the
     *  circle's point that way (Dropped); none for the center or a sweep of nothing.
     *
     *  @param  arc   The arc, swept no more than once round.
     *  @param  point The point.
     *  @return       The feet, unordered, and how many.
     */
    auto Footfalls(const Arc &arc, Point point) -> std::pair<std::array<Place, 4>, std::size_t>
    {
      std::pair<std::array<Place, 4>, std::size_t> falls = {{}, 0};
      if (point == arc.Center)
      {
        return falls;
      }
      const double toward = Heading(point - arc.Center).InRadians();
      const double sweep  = arc.Sweep.InRadians();
      if (Zero(sweep))
      {
        return falls;
      }
      const double whole = Turn.InRadians();
      for (const double side : {1.0, -1.0})
      {
        const Point  foot   = Dropped(arc.Center, side * arc.Radius, point);
        const double turned = (side > 0.0 ? toward : toward + Pi) - arc.Start.InRadians();
        const double first  = std::ceil((std::min(0.0, sweep) - turned) / whole);
        const double last   = std::floor((std::max(0.0, sweep) - turned) / whole);
        const auto   count =
            static_cast<std::size_t>(std::max(0.0, std::min(last - first + 1.0, MostTurns)));
        for (std::size_t turn = 0; turn < count; ++turn)
        {
          const double times             = first + static_cast<double>(turn);
          falls.first.at(falls.second++) = {
              .At = foot, .Along = std::clamp((turned + (times * whole)) / sweep, 0.0, 1.0)};
        }
      }
      return falls;
    }
  } // namespace

  auto Bypassed(const Vertex &from, Point end, const Asking &asking) -> bool
  {
    const Point  start    = from.Position;
    const Point  point    = asking.From;
    const double bulge    = from.Bulge;
    const Vector chord    = end - start;
    const double size     = asking.Size;
    const double limit    = asking.Least * (1.0 + 0x1p-40);
    const bool   moderate = size > 0x1p-400 && size < 0x1p400;
    if (!moderate)
    {
      return false;
    }
    if (!(std::abs(bulge) >= Straight) || start == end)
    {
      const Vector away   = point - start;
      const double square = Dot(chord, chord);
      const double along  = Dot(away, chord);
      const double margin =
          0x1p-40 * (std::abs(away.X) + std::abs(away.Y)) * (std::abs(chord.X) + std::abs(chord.Y));
      if (along < -margin || along > square + margin)
      {
        return true;
      }
      const double cross = Cross(chord, away);
      const double reach = limit + (0x1p-44 * size);
      return along         > margin
          && along         < square - margin
          && cross * cross > reach * reach * square * (1.0 + 0x1p-30);
    }
    const double reach0 = std::abs(chord.X) + std::abs(chord.Y);
    if (std::abs(bulge) <= 1.0)
    {
      /**
       *  An arc of a half turn at most lies within its chord's box grown by its height, the bulge
       *  times half the chord.
       */
      const double grow   = (std::abs(bulge) * reach0 * (0.5 + 0x1p-20)) + (0x1p-36 * size);
      const double across = std::max({std::min(start.X, end.X) - grow - point.X,
                                      point.X - std::max(start.X, end.X) - grow, 0.0});
      const double rising = std::max({std::min(start.Y, end.Y) - grow - point.Y,
                                      point.Y - std::max(start.Y, end.Y) - grow, 0.0});
      if ((across * across) + (rising * rising) > limit * limit * (1.0 + 0x1p-30))
      {
        return true;
      }
    }
    const double square = bulge * bulge;
    const Point  center =
        Midpoint(start, end) + (Perpendicular(chord) * ((1.0 - square) / (4.0 * bulge)));
    const double radius = Length(chord) * ((1.0 + square) / (4.0 * std::abs(bulge)));
    const Vector away   = point - center;
    const double reach  = Length(away);
    if (!(reach > 0.0))
    {
      return false;
    }
    const double sign   = bulge > 0.0 ? 1.0 : -1.0;
    const double base   = sign * Cross(chord, center - start);
    const double turn   = sign * (radius / reach) * Cross(chord, away);
    const double margin = 0x1p-36 * Largest(chord) * (size + radius) * (1.0 + (radius / reach));
    const double error  = 0x1p-44 * (size + (2.0 * radius) + reach);
    const auto   clear  = [margin, limit, error](double side, double distance) -> bool
    {
      return side > margin || (side < -margin && distance > limit + error);
    };
    return clear(base + turn, std::abs(reach - radius)) && clear(base - turn, reach + radius);
  }

  auto Braced(const Segment &segment, Point point) -> Brace
  {
    Brace brace;
    if (const Line *line = std::get_if<Line>(&segment))
    {
      if (const std::optional<double> along = Footing(Running(*line), point))
      {
        brace.Alongs.at(brace.Count++) = *along;
      }
      return brace;
    }
    const Arc &arc = std::get<Arc>(segment);
    if (point == arc.Center)
    {
      return brace;
    }
    const double toward = Heading(point - arc.Center).InRadians();
    const Swing  swing  = {.Start = arc.Start.InRadians(), .Sweep = arc.Sweep.InRadians()};
    for (const double side : {1.0, -1.0})
    {
      const Pass pass = PassedOnce(swing, side > 0.0 ? toward : toward + Pi);
      if (pass.Count >= 1.0)
      {
        brace.Sides.at(brace.Count)    = side;
        brace.Alongs.at(brace.Count++) = pass.Along;
      }
    }
    return brace;
  }

  auto Endpoint(const Arc &arc, Point point, const std::array<Point, 2> *ends) -> Place
  {
    /**
     *  ByAngles answers as the weighing does, so it is only how fast: a point by an end it never
     *  answers for, and with the ends given, one within a few steps of either is weighed at once.
     */
    const bool ended = ends != nullptr
                    && std::min(Largest((*ends)[0] - point), Largest((*ends)[1] - point))
                           <= 0x1p-30 * std::abs(arc.Radius);
    if (const std::optional<Place> place = ended ? std::nullopt : ByAngles(arc, point))
    {
      return *place;
    }
    Place      best    = {.At = Pointing(arc, 0.0, ends), .Along = 0.0};
    double     least   = Distance(point, best.At);
    const auto weighed = [point, &best, &least](const Place &place) -> void
    {
      const double away = Distance(point, place.At);
      if (away < least)
      {
        best  = place;
        least = away;
      }
    };
    const Place   end   = {.At = Pointing(arc, 1.0, ends), .Along = 1.0};
    const Ellipse round = Elliptic(arc);
    if (Coiled(round) && !(point == arc.Center))
    {
      /**
       *  Swept round more than once it passes each foot on every turn, all as near: the first pass
       *  of each is the first of the nearest, and the rest need no weighing, of which there may be
       *  65,536 a foot.
       */
      const double toward = Heading(point - arc.Center).InRadians();
      const Place  near   = {.At    = Dropped(arc.Center, arc.Radius, point),
                             .Along = Earliest(round, toward).value_or(0.0)};
      const Place  far    = {.At    = Dropped(arc.Center, -arc.Radius, point),
                             .Along = Earliest(round, toward + Pi).value_or(0.0)};
      weighed(near.Along < far.Along ? near : far);
      weighed(near.Along < far.Along ? far : near);
      weighed(end);
      return best;
    }
    auto [falls, count]         = Footfalls(arc, point);
    const std::span<Place> feet = std::span<Place>(falls).first(count);
    Lineup(feet);
    std::ranges::for_each(feet, weighed);
    weighed(end);
    return best;
  }
} // namespace Geoxy::Near
