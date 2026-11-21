/***************************************************************************************************
 *  @file      source/met/judge.cpp
 *  @brief     The candidates every pair of parts found, gathered with the ends of parts, each
 *             polished, judged by how the curves pass there, and merged into hits.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "met/judge.hpp"

#include "met/floor.hpp"
#include "met/parts.hpp"
#include "met/solve.hpp"
#include "real.hpp"
#include "scratch.hpp"
#include "trig.hpp"

#include "geoxy/angle.hpp"
#include "geoxy/arc.hpp"
#include "geoxy/at.hpp"
#include "geoxy/curve.hpp"
#include "geoxy/derivative.hpp"
#include "geoxy/domain.hpp"
#include "geoxy/ellipse.hpp"
#include "geoxy/intersect.hpp"
#include "geoxy/line.hpp"
#include "geoxy/matrix.hpp"
#include "geoxy/nearest.hpp"
#include "geoxy/place.hpp"
#include "geoxy/point.hpp"
#include "geoxy/polyline.hpp"
#include "geoxy/ray.hpp"
#include "geoxy/segment.hpp"
#include "geoxy/spline.hpp"
#include "geoxy/tolerance.hpp"
#include "geoxy/vector.hpp"
#include "geoxy/xline.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <iterator>
#include <limits>
#include <numeric>
#include <optional>
#include <span>
#include <utility>
#include <variant>
#include <vector>

namespace Geoxy::Met
{
  namespace
  {
    /**
     *  A curve's point and derivative at a place, as At and Derivative give them: an arc's or an
     *  ellipse's from one direction, worked out once for both.
     */
    struct Traced
    {
      Point  At  = {}; /**< The point.      */
      Vector Way = {}; /**< The derivative. */

      /**
       *  A curve's point and derivative at a place.
       *
       *  @param curve The curve.
       *  @param along The place.
       */
      Traced(const Curve &curve, double along)
      {
        if (const auto *arc = std::get_if<Arc>(&curve))
        {
          const Vector toward = Direction(arc->Start + (arc->Sweep * along));
          At                  = arc->Center + (toward * arc->Radius);
          Way                 = Perpendicular(toward) * (arc->Radius * arc->Sweep.InRadians());
          return;
        }
        if (const auto *ellipse = std::get_if<Ellipse>(&curve))
        {
          const Vector toward = Direction(ellipse->Start + (ellipse->Sweep * along));
          At                  = ellipse->Map * (Point{} + toward);
          Way                 = ellipse->Map * (Perpendicular(toward) * ellipse->Sweep.InRadians());
          return;
        }
        At  = Geoxy::At(curve, along);
        Way = Derivative(curve, along);
      }
    };

  } // namespace

  /**
   *  A candidate put in the curves' domains and polished by Newton's method on the gap between
   *  them, kept in their domains; a step that does not narrow the gap is not taken, so a tangency
   *  stays where it was found.
   *
   *  @param  first  The first curve.
   *  @param  second The second.
   *  @param  local  The candidate, in the curves' parameters.
   *  @return        The candidate polished.
   */
  auto Neared(const Curve &first, const Curve &second, Local local) -> Local
  {
    const Range one   = Domain(first);
    const Range other = Domain(second);
    local.First       = std::clamp(local.First, one.Start, one.End);
    local.Second      = std::clamp(local.Second, other.Start, other.End);
    std::array here   = {Traced(first, local.First), Traced(second, local.Second)};
    double     gap    = Distance(here[0].At, here[1].At);
    for (int step = 0; step < 16 && gap > 0.0; ++step)
    {
      const Vector apart = here[0].At - here[1].At;
      const Vector way   = here[0].Way;
      const Vector back  = here[1].Way;
      const double turn  = Cross(back, way);
      if (!(std::abs(turn) > 0x1p-48 * Length(way) * Length(back)))
      {
        break;
      }
      const Local next = {
          .First  = std::clamp(local.First + (Cross(apart, back) / turn), one.Start, one.End),
          .Second = std::clamp(local.Second + (Cross(apart, way) / turn), other.Start, other.End)};
      const std::array there    = {Traced(first, next.First), Traced(second, next.Second)};
      const double     narrower = Distance(there[0].At, there[1].At);
      if (!(narrower < gap))
      {
        break;
      }
      local = next;
      here  = there;
      gap   = narrower;
    }
    return local;
  }

  namespace
  {

    /**
     *  Whether a parameter reaches another, to within a few steps of a double: the end of one part
     *  and the start of the next, worked out apart.
     *
     *  @param  from  Where something starts.
     *  @param  until Where what it may run on from ends.
     *  @return       Whether it starts no later.
     */
    auto Reaches(double from, double until) -> bool
    {
      return from <= until
                         + (64.0
                            * std::numeric_limits<double>::epsilon()
                            * std::max({1.0, std::abs(from), std::abs(until)}));
    }

  } // namespace

  auto Continues(const Overlap &before, const Overlap &after) -> bool
  {
    const double enter  = before.Along ? after.Second.Start : after.Second.End;
    const bool   second = Reaches(before.Second.Start, enter) && Reaches(enter, before.Second.End);
    return before.Along == after.Along && Reaches(after.First.Start, before.First.End) && second;
  }

  auto Gathered(const List<Part> &parts, const List<Part> &others, const Scene &scene) -> Found
  {
    Found      gathered;
    List<bool> spent;
    for (const Part &part : parts)
    {
      for (const Part &other : others)
      {
        const Found found = Matched(part, other, scene.Tolerance);
        spent.push_back(found.Spent);
        for (const Local &local : found.Locals)
        {
          gathered.Locals.push_back(
              {.First = Placed(part, local.First), .Second = Placed(other, local.Second)});
        }
        for (const Overlap &shared : found.Shared)
        {
          gathered.Shared.push_back({.First  = {.Start = Placed(part, shared.First.Start),
                                                .End   = Placed(part, shared.First.End)},
                                     .Second = {.Start = Placed(other, shared.Second.Start),
                                                .End   = Placed(other, shared.Second.End)},
                                     .Along  = shared.Along});
        }
      }
    }
    gathered.Spent = std::ranges::contains(spent, true);
    std::ranges::sort(gathered.Shared, {},
                      [](const Overlap &overlap) -> double { return overlap.First.Start; });
    return gathered;
  }

  namespace
  {
    /**
     *  How far from the other curve, in tolerances, a part's end may lie and still be met with
     *  whatever else is found: one farther is left out where the hits the candidates found without
     *  it are settled.
     */
    constexpr double Far = 4.0;

    /**
     *  A curve's end, and the place nearest it on the other curve.
     */
    struct Outpost
    {
      Point  At      = {};  /**< The end.                           */
      double Nearest = 0.0; /**< The place nearest it on the other. */
    };

    /**
     *  The places other than the nearest where a polyline passes within the tolerance of the other
     *  curve's end. Nearest gives one place, the nearest of all. A curve's end within the tolerance
     *  of two passes of a polyline, as where the polyline doubles back past it or passes it again,
     *  is an end at each; and where the end's part crosses only the farther pass's carrier, short
     *  of the pass, no other candidate lies there. Each segment but the nearest place's, whose
     *  bounds do not put it farther, gives its own place nearest the end, where that lies within
     *  the tolerance and is not the nearest place. Only a part's end at one of the curve's ends, as
     *  a vertex put where the curve starts is, is so met.
     *
     *  @param  scene   The curves being met.
     *  @param  side    Whose the polyline is: 0 for the first curve, 1 for the second.
     *  @param  outpost A part's end of the other curve, and the place nearest it on the polyline.
     *  @return         The places, in order along the polyline; none for any other curve.
     */
    auto Repasses(const Scene &scene, std::size_t side, Outpost outpost) -> List<double>
    {
      const Curve       &curve    = side == 0 ? *scene.First : *scene.Second;
      const List<Floor> &floors   = side == 0 ? scene.FirstFloors : scene.SecondFloors;
      const List<Point> &ends     = side == 0 ? scene.SecondEnds : scene.FirstEnds;
      const auto        *polyline = std::get_if<Polyline>(&curve);
      const double       length   = scene.Tolerance.Length;
      const bool         tip      = std::ranges::any_of(
          ends, [&outpost](Point end) -> bool
          { return Exactly(end.X, outpost.At.X) && Exactly(end.Y, outpost.At.Y); });
      List<double> places;
      if (!tip || polyline == nullptr || floors.size() < 2)
      {
        return places;
      }
      const auto seat = static_cast<std::size_t>(
          std::clamp(std::floor(outpost.Nearest), 0.0, static_cast<double>(floors.size() - 1)));
      for (std::size_t index = 0; index < floors.size(); ++index)
      {
        if (index == seat || Remote(std::span{floors}.subspan(index, 1), outpost.At, length))
        {
          continue;
        }
        const Place foot =
            std::visit([&outpost](const auto &segment) -> Place
                       { return Nearest(segment, outpost.At); }, Nth(*polyline, index));
        const double place = static_cast<double>(index) + foot.Along;
        if (!Exactly(place, outpost.Nearest) && Distance(foot.At, outpost.At) <= length)
        {
          places.push_back(place);
        }
      }
      return places;
    }

  } // namespace

  auto Ended(const Scene &scene, std::size_t side, bool all) -> Tips
  {
    const bool         first   = side == 0;
    const Curve       &curve   = first ? *scene.First : *scene.Second;
    const Curve       &against = first ? *scene.Second : *scene.First;
    const List<Part>  &parts   = first ? *scene.FirstParts : *scene.SecondParts;
    const List<Floor> &floors  = first ? scene.SecondFloors : scene.FirstFloors;
    Tips               tips;
    const auto         add = [&tips, first](double along, double other) -> void
    {
      tips.Locals.push_back(first ? Local{.First = along, .Second = other}
                                  : Local{.First = other, .Second = along});
    };
    for (const Part &part : parts)
    {
      const Range bounds = Bounds(part);
      for (const double end : {bounds.Start, bounds.End})
      {
        const double along = Placed(part, end);
        if (!std::isfinite(along))
        {
          continue;
        }
        const Point point = At(curve, along);
        if (!all && Remote(floors, point, Far * scene.Tolerance.Length))
        {
          tips.Skipped = true;
          continue;
        }
        const double other = Nearest(against, point).Along;
        add(along, other);
        for (const double pass : Repasses(scene, 1 - side, {.At = point, .Nearest = other}))
        {
          add(along, pass);
        }
      }
    }
    return tips;
  }

  namespace
  {
    /**
     *  Whether a point lies within the tolerance of one of a run of ends.
     *
     *  @param  ends      The ends.
     *  @param  point     The point.
     *  @param  tolerance The tolerance.
     *  @return           Whether it does.
     */
    auto Ending(const List<Point> &ends, Point point, const Tolerance &tolerance) -> bool
    {
      return std::ranges::any_of(ends, [&tolerance, point](Point end) -> bool
                                 { return Distance(end, point) <= tolerance.Length; });
    }

    /**
     *  A curve's ways into a place and out of it: its derivative just before and just after, so
     *  that at a corner, a polyline's vertex or a spline's knot, the two are its two parts'; one
     *  derivative for the kinds that have no corners.
     *
     *  @param  curve The curve.
     *  @param  along Its parameter at the place.
     *  @return       The way in and the way out.
     */
    auto Ways(const Curve &curve, double along) -> std::array<Vector, 2>
    {
      if (!std::holds_alternative<Polyline>(curve) && !std::holds_alternative<Spline>(curve))
      {
        const Vector way = Derivative(curve, along);
        return {way, way};
      }
      const Range  domain = Domain(curve);
      const double before = along > domain.Start ? std::nextafter(along, domain.Start) : along;
      const double after  = along < domain.End ? std::nextafter(along, domain.End) : along;
      return {Derivative(curve, before), Derivative(curve, after)};
    }

  } // namespace

  auto Passing(const std::array<Vector, 2> &ways, Vector other, double radians) -> int
  {
    std::array<int, 2> sides = {};
    for (std::size_t index = 0; index < ways.size(); ++index)
    {
      const double turn  = Cross(other, ways.at(index));
      const double apart = Trig::Sin(radians) * Length(other) * Length(ways.at(index));
      sides.at(index)    = static_cast<int>(turn > apart) - static_cast<int>(turn < -apart);
    }
    return sides[0] == sides[1] ? sides[0] : 0;
  }

  namespace
  {
    /**
     *  How far one step of a double of a curve's parameter moves the curve at a place: its speed
     *  there times the step.
     *
     *  @param  curve The curve.
     *  @param  along The place.
     *  @return       The distance.
     */
    auto Stride(const Curve &curve, double along) -> double
    {
      return Length(Derivative(curve, along)) * Above(std::abs(along));
    }

    /**
     *  A place on a curve and a stretch of it.
     */
    struct Stretched
    {
      Range  Stretch = {};  /**< The stretch. */
      double Along   = 0.0; /**< The place.   */
    };

    /**
     *  Whether a place of a curve lies in a stretch of it, or its point within the tolerance of
     *  either end's.
     *
     *  @param  curve The curve.
     *  @param  place The stretch and the place.
     *  @param  point The curve's point at the place, and the tolerance's length.
     *  @return       Whether it does.
     */
    auto Sheltered(const Curve &curve, Stretched place, std::pair<Point, double> point) -> bool
    {
      const auto [spot, length] = point;
      return (place.Stretch.Start <= place.Along && place.Along <= place.Stretch.End)
          || Distance(spot, At(curve, place.Stretch.Start)) <= length
          || Distance(spot, At(curve, place.Stretch.End)) <= length;
    }

  } // namespace

  auto Judged(const Scene &scene, const List<Overlap> &overlaps, Local candidate)
      -> std::optional<Seen>
  {
    const Curve    &first     = *scene.First;
    const Curve    &second    = *scene.Second;
    const Tolerance tolerance = scene.Tolerance;
    const Range     one       = Domain(first);
    const Range     two       = Domain(second);
    const Local     held      = {.First  = std::clamp(candidate.First, one.Start, one.End),
                                 .Second = std::clamp(candidate.Second, two.Start, two.End)};
    const Point     start     = At(first, held.First);
    const Point     other     = At(second, held.Second);
    const bool      ending    = Distance(start, other) <= tolerance.Length
                             && (Ending(scene.FirstEnds, start, tolerance)
                                 || Ending(scene.SecondEnds, other, tolerance));
    const Local     polished  = ending ? held : Neared(first, second, candidate);
    const Point     here      = At(first, polished.First);
    const Point     there     = At(second, polished.Second);
    /**
     *  A place is told by a double, and the curve's point there by where that double puts it:
     *  beside a heavy weight of a spline one step of the parameter moves the curve by more than the
     *  tolerance, and no place at all would be found within it of a crossing that is there. Where
     *  either curve is a spline, whose crossings are worked out through its parameter, the points
     *  are held to the tolerance or to twice what a step of each place moves its curve, whichever
     *  is more.
     */
    const double coarse =
        scene.Splined ? 2.0 * (Stride(first, polished.First) + Stride(second, polished.Second))
                      : 0.0;
    const double within = std::max(tolerance.Length, coarse);
    if (!(Distance(here, there) <= within))
    {
      return std::nullopt;
    }
    /**
     *  A candidate in a shared stretch is the stretch's, on both curves: one whose place on either
     *  lies outside it is another pass of that curve through the stretch, a crossing of its own, as
     *  it is met the other way round.
     */
    const bool inside = std::ranges::any_of(
        overlaps,
        [&scene, &polished, &here, &there](const Overlap &overlap) -> bool
        {
          return Sheltered(*scene.First, {.Stretch = overlap.First, .Along = polished.First},
                           {here, scene.Tolerance.Length})
              && Sheltered(*scene.Second, {.Stretch = overlap.Second, .Along = polished.Second},
                           {there, scene.Tolerance.Length});
        });
    if (inside)
    {
      return std::nullopt;
    }
    const auto ways   = Ways(first, polished.First);
    const auto others = Ways(second, polished.Second);
    const int  side   = Passing(ways, others[0], tolerance.Radians);
    const bool across = side != 0 && Passing(others, ways[0], tolerance.Radians) != 0;
    Contact    kind   = across ? Contact::Cross : Contact::Touch;
    if (Ending(scene.FirstEnds, here, tolerance) || Ending(scene.SecondEnds, there, tolerance))
    {
      kind = Contact::End;
    }
    return Seen{.Found  = {.At     = Midpoint(here, there),
                           .First  = {.At = here, .Along = polished.First},
                           .Second = {.At = there, .Along = polished.Second},
                           .Kind   = kind},
                .Side   = across ? side : 0,
                .Within = within};
  }

  namespace
  {
    /**
     *  The step from one of a curve's parameters to another, the short way round a closed curve,
     *  across its seam where that is shorter.
     *
     *  @param  curve  The curve.
     *  @param  closed Whether it is closed.
     *  @param  span   The two parameters, from and to.
     *  @return        The step.
     */
    auto Shortcut(const Curve &curve, bool closed, Range span) -> double
    {
      const Range  domain = Domain(curve);
      const double whole  = domain.End - domain.Start;
      const double step   = span.End - span.Start;
      return closed && std::abs(step) > whole / 2 ? step - std::copysign(whole, step) : step;
    }

    /**
     *  A parameter part of the way between two of a curve's, the short way round a closed curve.
     *
     *  @param  curve  The curve.
     *  @param  closed Whether it is closed.
     *  @param  span   The two parameters, from and to, and the part of the way between them.
     *  @return        The parameter there.
     */
    auto Partway(const Curve &curve, bool closed, std::array<double, 3> span) -> double
    {
      const auto [from, until, part] = span;
      const Range  domain            = Domain(curve);
      const double whole             = domain.End - domain.Start;
      const double along = from + (Shortcut(curve, closed, {.Start = from, .End = until}) * part);
      if (along > domain.End)
      {
        return along - whole;
      }
      return along < domain.Start ? along + whole : along;
    }

    /**
     *  At most how far a curve runs between two of its parameters near each other: the step times
     *  the faster of its speeds at the two.
     *
     *  @param  curve  The curve.
     *  @param  closed Whether it is closed.
     *  @param  span   The two parameters.
     *  @return        The distance.
     */
    auto Travel(const Curve &curve, bool closed, Range span) -> double
    {
      const double speed =
          std::max(Length(Derivative(curve, span.Start)), Length(Derivative(curve, span.End)));
      return speed * std::abs(Shortcut(curve, closed, span));
    }

    /**
     *  Whether two places of a curve lie half a turn or more apart on an arc or an ellipse swept
     *  round more than once. Such a curve comes back to each of its points on every turn, so the
     *  places between two of its passes of one point, a whole number of turns apart, are that point
     *  again wherever the count of turns divides evenly: no telling them apart by looking between.
     *
     *  @param  curve The curve.
     *  @param  span  The two places.
     *  @return       Whether they do.
     */
    auto Lapped(const Curve &curve, Range span) -> bool
    {
      const double sweep = Wound(curve);
      return std::abs(sweep) > Turn.InRadians() && std::abs(sweep * (span.End - span.Start)) >= Pi;
    }

    /**
     *  A curve's stretch between two places where it lies on one line or one arc of a circle: the
     *  line, ray, endless line or arc itself, or the one segment of a polyline it lies in.
     */
    struct Stint
    {
      Curve Shape = Line{}; /**< The line, ray, endless line or arc. */
      Range Span  = {};     /**< The stretch, from one place to the other, in its own
                                 parameter, the short way round a closed curve.       */
    };

    /**
     *  A curve's stretch between two places as one line or arc, where it is one.
     *
     *  @param  curve  The curve.
     *  @param  closed Whether it is closed.
     *  @param  span   The two places.
     *  @return        The stretch; nothing where the curve is of another kind or the stretch runs
     *                 over a vertex.
     */
    auto Stinted(const Curve &curve, bool closed, Range span) -> std::optional<Stint>
    {
      const double reach = span.Start + Shortcut(curve, closed, span);
      if (const auto *polyline = std::get_if<Polyline>(&curve))
      {
        const double index = std::min(std::floor(std::min(span.Start, reach)),
                                      static_cast<double>(Count(*polyline)) - 1.0);
        if (!(index >= 0.0) || !(std::max(span.Start, reach) <= index + 1.0))
        {
          return std::nullopt;
        }
        return std::visit(
            [&span, reach, index](const auto &segment) -> Stint
            {
              return {.Shape = Curve{segment},
                      .Span  = {.Start = span.Start - index, .End = reach - index}};
            },
            Nth(*polyline, static_cast<std::size_t>(index)));
      }
      if (std::holds_alternative<Line>(curve)
          || std::holds_alternative<Ray>(curve)
          || std::holds_alternative<Xline>(curve)
          || std::holds_alternative<Arc>(curve))
      {
        return Stint{.Shape = curve, .Span = {.Start = span.Start, .End = reach}};
      }
      return std::nullopt;
    }

    /**
     *  A line, a ray, an endless line or an arc as a stretch is measured from: its circle, if an
     *  arc, its ends, and its way, if straight.
     */
    struct Yard
    {
      const Arc           *Round = nullptr; /**< The arc, if it is one. */
      std::array<Point, 2> Ends  = {};      /**< Its ends.              */
      std::size_t          Many  = 0;       /**< How many: none for a whole circle or an
                                                 endless line, one for a ray.               */
      Vector Way = {};                      /**< Its way, where it is straight. */
    };

    /**
     *  A line, a ray, an endless line or an arc as a stretch is measured from.
     *
     *  @param  other The line, ray, endless line or arc.
     *  @return       Its circle, ends and way.
     */
    auto Yarded(const Curve &other) -> Yard
    {
      if (const auto *round = std::get_if<Arc>(&other))
      {
        return {.Round = round,
                .Ends  = {At(other, 0.0), At(other, 1.0)},
                .Many  = std::abs(round->Sweep.InRadians()) < 2.0 * Pi ? 2UZ : 0UZ,
                .Way   = {}};
      }
      if (const auto *line = std::get_if<Line>(&other))
      {
        return {.Round = nullptr,
                .Ends  = {line->Start, line->End},
                .Many  = 2,
                .Way   = line->End - line->Start};
      }
      if (const auto *ray = std::get_if<Ray>(&other))
      {
        return {.Round = nullptr, .Ends = {ray->Start, ray->Start}, .Many = 1, .Way = ray->Way};
      }
      return {.Round = nullptr, .Ends = {}, .Many = 0, .Way = std::get<Xline>(other).Way};
    }

    /**
     *  The places a stretch may lie farthest from another curve at, as points.
     */
    struct Peaks
    {
      std::array<Point, 6> At    = {}; /**< The points.         */
      std::size_t          Count = 0;  /**< How many there are. */
    };

    /**
     *  Where a stretch of a circle may lie farthest from a line, a ray, an endless line or an arc:
     *  where the circle's tangent is parallel to a straight other; on the line of centers with an
     *  arc; opposite each end of the other; and where the bisector of an arc's ends cuts it. Each
     *  is the circle's point a way from its center, kept where the stretch passes it: its angle
     *  past the stretch's start, the way the stretch turns, within the stretch's turn. A way of
     *  nothing has no point, and its angle, not a number, is within no turn.
     *
     *  @param  arc   The circle's arc, and the stretch of it.
     *  @param  gauge The other.
     *  @return       The places.
     */
    auto Circlewise(std::pair<const Arc *, Range> arc, const Yard &gauge) -> Peaks
    {
      const auto &[round, span] = arc;
      const Point  center       = round->Center;
      const double from = round->Start.InRadians() + (round->Sweep.InRadians() * span.Start);
      const double turn = round->Sweep.InRadians() * (span.End - span.Start);
      Peaks        spots;
      const auto   toward = [&](Vector direction) -> void
      {
        const double none = std::numeric_limits<double>::quiet_NaN();
        const Vector unit = Unit(direction).value_or(Vector{.X = none, .Y = none});
        const double past = std::copysign(1.0, turn) * (std::atan2(unit.Y, unit.X) - from);
        if (std::fmod(std::fmod(past, 2.0 * Pi) + (2.0 * Pi), 2.0 * Pi) <= std::abs(turn))
        {
          spots.At.at(spots.Count) = center + (unit * round->Radius);
          ++spots.Count;
        }
      };
      const Vector way = Perpendicular(gauge.Way);
      toward(way);
      toward(-way);
      for (const Point tip : std::span{gauge.Ends}.first(gauge.Many))
      {
        toward(center - tip);
      }
      if (gauge.Round != nullptr)
      {
        toward(gauge.Round->Center - center);
        toward(center - gauge.Round->Center);
        /**
         *  The bisector of the arc's ends, middle + across k, k a root of |across|^2 k^2 + 2
         *  (across . middle) k + |middle|^2 - r^2. A whole circle has no ends to bisect, and two
         *  ends at one point no bisector: no root, or roots not numbers.
         */
        const Vector across =
            Perpendicular(gauge.Ends[1] - gauge.Ends[0]) * (static_cast<double>(gauge.Many) / 2.0);
        const Vector middle = Midpoint(gauge.Ends[0], gauge.Ends[1]) - center;
        const double square = Dot(across, across);
        const double lean   = Dot(across, middle);
        const double left =
            (lean * lean) - (square * (Dot(middle, middle) - (round->Radius * round->Radius)));
        if (left >= 0.0 && square > 0.0)
        {
          toward(middle + (across * ((-lean - std::sqrt(left)) / square)));
          toward(middle + (across * ((-lean + std::sqrt(left)) / square)));
        }
      }
      return spots;
    }

    /**
     *  Where a straight stretch may lie farthest from an arc: the foot of the arc's center, and
     *  where the bisector of the arc's ends crosses it, each kept where it lies strictly between
     *  the stretch's ends, the part of the way along it a ratio turned to a positive denominator.
     *
     *  @param  start The stretch's start.
     *  @param  run   The way from its start to its end.
     *  @param  gauge The arc.
     *  @return       The places.
     */
    auto Lining(Point start, Vector run, const Yard &gauge) -> Peaks
    {
      Peaks      spots;
      const auto keep = [&](double above, double below) -> void
      {
        if (0.0 < above && above < below)
        {
          spots.At.at(spots.Count) = start + (run * (above / below));
          ++spots.Count;
        }
      };
      keep(Dot(gauge.Round->Center - start, run), Dot(run, run));
      const Vector chord =
          (gauge.Ends[1] - gauge.Ends[0]) * (static_cast<double>(gauge.Many) / 2.0);
      const double sign = std::copysign(1.0, Dot(run, chord));
      keep(sign * Dot(Midpoint(gauge.Ends[0], gauge.Ends[1]) - start, chord),
           sign * Dot(run, chord));
      return spots;
    }

    /**
     *  How far a plain stretch lies from a line, a ray, an endless line or an arc at its farthest,
     *  where that is more than at its ends, exactly: where it is tangent to a curve all of whose
     *  points lie as far from the other. Those curves are lines parallel to a straight other,
     *  circles about an arc's center, and circles about an end of the other; and where the other is
     *  an arc, the bisector of its ends, across which the end nearest changes. Along a line, the
     *  distance from a line or a point only falls to its least and rises again, and from a circle
     *  it rises to its most at the foot of the center inside it; along a circle, the distance from
     *  a line is most where the circle's tangent is parallel to it, from another circle on the line
     *  of centers, and from a point at the circle's point opposite it. Where the nearest point of
     *  the other moves from inside it to an end, the two distances are one and so are their slopes,
     *  and the farthest is one of the two's own. The stretch's ends, each a hit's point on its
     *  curve, lie within the tolerance of the other, and are not looked at.
     *
     *  @param  stretch The stretch.
     *  @param  other   The line, ray, endless line or arc.
     *  @return         The distance; 0 where no place inside the stretch is such a place.
     */
    auto Utmost(const Stint &stretch, const Curve &other) -> double
    {
      const auto &[shape, span] = stretch;
      const Yard gauge          = Yarded(other);
      Peaks      spots;
      if (const auto *arc = std::get_if<Arc>(&shape))
      {
        spots = Circlewise({arc, span}, gauge);
      }
      else if (gauge.Round != nullptr)
      {
        const Point start = At(shape, span.Start);
        spots             = Lining(start, At(shape, span.End) - start, gauge);
      }
      double farthest = 0.0;
      for (const Point spot : std::span{spots.At}.first(spots.Count))
      {
        farthest = std::max(farthest, Distance(spot, Nearest(other, spot).At));
      }
      return farthest;
    }

    /**
     *  How far apart two curves come between two hits, exactly, where each runs between them on one
     *  line or one arc: the farther of the farthest each stretch lies from the other's line or arc,
     *  the one the other's stretch lies on. A polyline's other segments are not looked at: a pass
     *  of one near the stretch, as at a loop's double point, makes no two hits one.
     *
     *  @param  scene  The curves being met.
     *  @param  pair   The two hits.
     *  @param  closed Whether each curve is closed.
     *  @return        The distance; nothing where a curve between them is not one line or arc.
     */
    auto Sundering(const Scene &scene, const std::array<Hit, 2> &pair, std::array<bool, 2> closed)
        -> std::optional<double>
    {
      const std::optional<Stint> one = Stinted(
          *scene.First, closed[0], {.Start = pair[0].First.Along, .End = pair[1].First.Along});
      const std::optional<Stint> other = Stinted(
          *scene.Second, closed[1], {.Start = pair[0].Second.Along, .End = pair[1].Second.Along});
      if (!one || !other)
      {
        return std::nullopt;
      }
      return std::max(Utmost(*one, other->Shape), Utmost(*other, one->Shape));
    }

    /**
     *  Halfway along each segment a polyline runs along between two of its places, as parts of the
     *  way from the first to the second, the short way round a closed one: none where it stays on
     *  one segment or is no polyline; the first eight and the last eight where there are more.
     *
     *  @param  curve  The curve.
     *  @param  closed Whether it is closed.
     *  @param  span   The two places.
     *  @return        The parts.
     */
    auto Middles(const Curve &curve, bool closed, Range span) -> List<double>
    {
      List<double> parts;
      const double step = Shortcut(curve, closed, span);
      if (!std::holds_alternative<Polyline>(curve) || !(std::abs(step) > 0.0))
      {
        return parts;
      }
      const double low    = std::min(span.Start, span.Start + step);
      const double high   = std::max(span.Start, span.Start + step);
      const double first  = std::floor(low) + 1.0;
      const auto   joints = first < high ? static_cast<std::size_t>(std::ceil(high - first)) : 0UZ;
      double       before = low;
      for (std::size_t index = 0; index < joints; ++index)
      {
        const double joint = first + static_cast<double>(index);
        parts.push_back((((before + joint) / 2.0) - span.Start) / step);
        before = joint;
      }
      if (parts.empty())
      {
        return parts;
      }
      parts.push_back((((before + high) / 2.0) - span.Start) / step);
      if (parts.size() > 16)
      {
        parts.erase(parts.begin() + 8, parts.end() - 8);
      }
      return parts;
    }

    /**
     *  Whether two hits are one: the curves within twice the tolerance all the way between them,
     *  unless both run less than the tolerance between them. Where each curve runs between them on
     *  one line or one arc, as Sundering tells it, exactly; otherwise at a quarter, a half and
     *  three quarters of the way along the first and along the second, halfway first, so that two
     *  places of one curve that pass through one point, as at a loop's double point, stay two hits,
     *  though halfway between them the curve passes the other again. Twice, as candidates at the
     *  edge of a touch's band are each within the tolerance, and the curves between two of them as
     *  near it as rounding leaves. Never two passes, turns apart, of an arc swept round more than
     *  once: four turns apart, every place looked at between them is the same point again.
     *
     *  @param  scene The curves being met.
     *  @param  seen  The two hits.
     *  @return       Whether they are one.
     */
    auto Together(const Scene &scene, const std::array<Seen, 2> &seen) -> bool
    {
      constexpr int            Places = 4;
      const std::array<Hit, 2> pair   = {seen[0].Found, seen[1].Found};
      if (Lapped(*scene.First, {.Start = pair[0].First.Along, .End = pair[1].First.Along})
          || Lapped(*scene.Second, {.Start = pair[0].Second.Along, .End = pair[1].Second.Along}))
      {
        return false;
      }
      const auto closed = [](const Curve &curve, const List<Point> &ends) -> bool
      {
        return ends.empty() && std::isfinite(Domain(curve).Start);
      };
      const bool   firstClosed  = closed(*scene.First, scene.FirstEnds);
      const bool   secondClosed = closed(*scene.Second, scene.SecondEnds);
      const double length       = std::max(seen[0].Within, seen[1].Within);
      const bool   close        = Distance(pair[0].At, pair[1].At) <= length
                               && Travel(*scene.First, firstClosed,
                                         {.Start = pair[0].First.Along, .End = pair[1].First.Along})
                                      <= length
                               && Travel(*scene.Second, secondClosed,
                                         {.Start = pair[0].Second.Along, .End = pair[1].Second.Along})
                                      <= length;
      if (close)
      {
        return true;
      }
      if (const std::optional<double> apart = Sundering(scene, pair, {firstClosed, secondClosed}))
      {
        return *apart <= 2.0 * length;
      }
      constexpr std::array<int, 3> Order   = {2, 1, 3};
      const bool                   sampled = std::ranges::all_of(
          Order,
          [&scene, &pair, firstClosed, secondClosed, length](int place) -> bool
          {
            const double part = static_cast<double>(place) / Places;
            const Point  first =
                At(*scene.First, Partway(*scene.First, firstClosed,
                                         {pair[0].First.Along, pair[1].First.Along, part}));
            const Point second =
                At(*scene.Second, Partway(*scene.Second, secondClosed,
                                          {pair[0].Second.Along, pair[1].Second.Along, part}));
            const double within = 2.0 * length;
            return !Remote(scene.SecondFloors, first, within)
                && !Remote(scene.FirstFloors, second, within)
                && Distance(first, Nearest(*scene.Second, first).At) <= within
                && Distance(second, Nearest(*scene.First, second).At) <= within;
          });
      if (!sampled)
      {
        return false;
      }
      /**
       *  A polyline that passes a vertex between the two hits is looked at on every segment it runs
       *  along between them, halfway along each: a quarter, a half and three quarters of the way
       *  may all fall on a segment that lies along the other curve, and the last stretch, on the
       *  next segment, leave it by thousands of tolerances.
       */
      const double within = 2.0 * length;
      const auto   beside = [within](const Curve &other, const List<Floor> &floors,
                                     Point point) -> bool
      {
        return !Remote(floors, point, within)
            && Distance(point, Nearest(other, point).At) <= within;
      };
      const auto along = [](const Curve &curve, bool shut, std::array<double, 2> span,
                            double part) -> Point
      {
        return At(curve, Partway(curve, shut, {span[0], span[1], part}));
      };
      return std::ranges::all_of(
                 Middles(*scene.First, firstClosed,
                         {.Start = pair[0].First.Along, .End = pair[1].First.Along}),
                 [&](double part) -> bool
                 {
                   return beside(*scene.Second, scene.SecondFloors,
                                 along(*scene.First, firstClosed,
                                       {pair[0].First.Along, pair[1].First.Along}, part));
                 })
          && std::ranges::all_of(
                 Middles(*scene.Second, secondClosed,
                         {.Start = pair[0].Second.Along, .End = pair[1].Second.Along}),
                 [&](double part) -> bool
                 {
                   return beside(*scene.First, scene.FirstFloors,
                                 along(*scene.Second, secondClosed,
                                       {pair[0].Second.Along, pair[1].Second.Along}, part));
                 });
    }

    /**
     *  Whether two hits in a row along the first curve are one: the curves together between them,
     *  but never across an end of a shared stretch along the first: one inside it, as another pass
     *  of the second curve through it, and one outside are together along the stretch only as the
     *  stretch is.
     *
     *  @param  scene    The curves being met.
     *  @param  overlaps The stretches they share.
     *  @param  pair     The hit before and the hit after.
     *  @return          Whether they are one.
     */
    auto Fused(const Scene &scene, const List<Overlap> &overlaps, const std::array<Seen, 2> &pair)
        -> bool
    {
      const bool across =
          std::ranges::any_of(overlaps,
                              [&pair](const Overlap &overlap) -> bool
                              {
                                const auto [low, high] = std::minmax(pair[0].Found.First.Along,
                                                                     pair[1].Found.First.Along);
                                const auto between     = [low, high](double end) -> bool
                                {
                                  return low < end && end < high;
                                };
                                return between(overlap.First.Start) || between(overlap.First.End);
                              });
      return !across && Together(scene, pair);
    }

    /**
     *  How a run's hit is chosen among its members: one within the tolerance of an end where the
     *  curves cross, its two points a few steps of a double apart at their size, before any, as an
     *  end hit is put where they cross within the tolerance of the end; then one at a curve's end
     *  itself, as the end's own candidate is; then one within the tolerance of an end polished
     *  short of a crossing, as along the band about a touch; then any.
     *
     *  @param  scene The curves being met.
     *  @param  hit   A member's hit.
     *  @return       3 where they cross near an end, 2 at an end, 1 near one, 0 otherwise.
     */
    auto Standing(const Scene &scene, const Hit &hit) -> int
    {
      const auto tipped = [](const List<Point> &ends, Point point) -> bool
      {
        return std::ranges::any_of(ends, [point](Point end) -> bool
                                   { return Exactly(end.X, point.X) && Exactly(end.Y, point.Y); });
      };
      const Point  here  = hit.First.At;
      const Point  there = hit.Second.At;
      const double size =
          std::max({std::abs(here.X), std::abs(here.Y), std::abs(there.X), std::abs(there.Y)});
      if (hit.Kind                 == Contact::End
          && Distance(here, there) <= 8.0 * std::numeric_limits<double>::epsilon() * size)
      {
        return 3;
      }
      if (tipped(scene.FirstEnds, here) || tipped(scene.SecondEnds, there))
      {
        return 2;
      }
      return hit.Kind == Contact::End ? 1 : 0;
    }

    /**
     *  One hit for a run of hits that are one: an end of a curve if one is, where the curves cross
     *  within the tolerance of it, else the end itself, as an end within the tolerance of the other
     *  is a hit where it stands, before a candidate polished near it along the band about a touch,
     *  which the order of the curves and where it started from move; else the one where the curves
     *  come nearest, the first of those as near; an end if any is, a touch if any is or if the
     *  first passes to both sides, going in and coming back out, else a crossing. Where along a
     *  shallow crossing its candidates were polished to never makes it a touch.
     *
     *  @param  scene The curves being met.
     *  @param  run   The run, in order along the first curve.
     *  @return       The hit.
     */
    auto Merged(const Scene &scene, std::span<const Seen> run) -> Hit
    {
      std::size_t pick  = 0;
      bool        left  = false;
      bool        right = false;
      Contact     kind  = Contact::Cross;
      for (std::size_t member = 0; member < run.size(); ++member)
      {
        const Hit &hit  = run[member].Found;
        const Hit &best = run[pick].Found;
        left            = left || run[member].Side > 0;
        right           = right || run[member].Side < 0;
        const int  rank = Standing(scene, hit);
        const int  over = Standing(scene, best);
        const bool closer =
            Distance(hit.First.At, hit.Second.At) < Distance(best.First.At, best.Second.At);
        if (rank > over || (rank == over && closer))
        {
          pick = member;
        }
        if (hit.Kind == Contact::End || (hit.Kind == Contact::Touch && kind == Contact::Cross))
        {
          kind = hit.Kind;
        }
      }
      Hit merged  = run[pick].Found;
      merged.Kind = left && right && kind == Contact::Cross ? Contact::Touch : kind;
      return merged;
    }

  } // namespace

  auto Hits(const Scene &scene, const List<Overlap> &overlaps, const List<Seen> &seen) -> List<Hit>
  {
    /**
     *  The sightings in order along the first curve, those at one place in the order they came:
     *  their indices sorted by place and then by index, the order a stable sort gives, with no room
     *  taken but the indices' own.
     */
    List<std::size_t> order(seen.size());
    std::ranges::iota(order, std::size_t{0});
    std::ranges::sort(order,
                      [&seen](std::size_t lhs, std::size_t rhs) -> bool
                      {
                        return std::pair{seen[lhs].Found.First.Along, lhs}
                             < std::pair{seen[rhs].Found.First.Along, rhs};
                      });
    const bool       shut = scene.FirstEnds.empty() && std::isfinite(Domain(*scene.First).Start);
    List<List<Seen>> runs;
    for (const std::size_t index : order)
    {
      const Seen &one   = seen[index];
      const auto  joins = [&scene, &overlaps, &one](const List<Seen> &run) -> bool
      {
        return Fused(scene, overlaps, {run.back(), one});
      };
      /**
       *  A sighting is one with the run before it, or with any earlier run the curves stay together
       *  from: the run of a long shallow approach goes on past a crossing of another segment that
       *  falls between two of its sightings along the first curve, which made the answer hang on
       *  which curve came first. An earlier run is weighed only where, halfway to it along the
       *  first curve, the second is near, as it is wherever the curves stay together.
       */
      const auto near = [&scene, &one, shut](const List<Seen> &run) -> bool
      {
        const double within = 2.0 * std::max(run.back().Within, one.Within);
        const Point  between =
            At(*scene.First, Partway(*scene.First, shut,
                                     {run.back().Found.First.Along, one.Found.First.Along, 0.5}));
        return !Remote(scene.SecondFloors, between, within);
      };
      const auto earlier = std::span{runs}.first(runs.empty() ? 0 : runs.size() - 1);
      auto       home    = runs.end();
      if (!runs.empty() && joins(runs.back()))
      {
        home = std::prev(runs.end());
      }
      else if (const auto found =
                   std::ranges::find_if(earlier, [&near, &joins](const List<Seen> &run) -> bool
                                        { return near(run) && joins(run); });
               found != earlier.end())
      {
        home = runs.begin() + (found - earlier.begin());
      }
      if (home == runs.end())
      {
        runs.push_back({one});
      }
      else
      {
        home->push_back(one);
      }
    }
    List<Hit> hits;
    hits.reserve(runs.size());
    for (const List<Seen> &run : runs)
    {
      hits.push_back(Merged(scene, run));
    }
    return hits;
  }

  namespace
  {
    /**
     *  One of the two curves being met, and the other's floors about it.
     */
    struct Beside
    {
      const Curve       *Own    = nullptr; /**< The curve the places are of. */
      const Curve       *Other  = nullptr; /**< The other curve.             */
      const List<Floor> *Floors = nullptr; /**< The other's floors.          */
      double             Length = 0.0;     /**< The tolerance's length.      */
    };

    /**
     *  Whether the curves stay within the tolerance halfway between two places of one of them.
     *
     *  @param  beside The curve, the other and the tolerance.
     *  @param  from   One parameter of the curve.
     *  @param  until  The other.
     *  @return        Whether they do.
     */
    auto Stays(const Beside &beside, double from, double until) -> bool
    {
      const Point between = At(*beside.Own, from + ((until - from) / 2));
      return !Remote(*beside.Floors, between, beside.Length)
          && Distance(between, Nearest(*beside.Other, between).At) <= beside.Length;
    }

    /**
     *  Whether a place of a curve lies just past a stretch of it, the curves within the tolerance
     *  from the stretch's end to the place.
     *
     *  @param  beside  The curve, the other and the tolerance.
     *  @param  stretch The stretch.
     *  @param  along   The place.
     *  @return         Whether it does.
     */
    auto Trailing(const Beside &beside, Range stretch, double along) -> bool
    {
      return (along > stretch.End && Stays(beside, stretch.End, along))
          || (along < stretch.Start && Stays(beside, along, stretch.Start));
    }

  } // namespace

  void Absorb(const Scene &scene, Intersections &meeting)
  {
    /**
     *  A hit just past a shared stretch's end along either curve is the stretch's: told along the
     *  first alone, a hit past the end along the second but within the stretch along the first, as
     *  where the first runs on a step past a joint of the second, was kept one way round and not
     *  the other.
     */
    const Beside first  = {.Own    = scene.First,
                           .Other  = scene.Second,
                           .Floors = &scene.SecondFloors,
                           .Length = scene.Tolerance.Length};
    const Beside second = {.Own    = scene.Second,
                           .Other  = scene.First,
                           .Floors = &scene.FirstFloors,
                           .Length = scene.Tolerance.Length};
    for (const Overlap &overlap : meeting.Overlaps)
    {
      std::erase_if(meeting.Hits,
                    [&first, &second, &overlap](const Hit &hit) -> bool
                    {
                      return Trailing(first, overlap.First, hit.First.Along)
                          || Trailing(second, overlap.Second, hit.Second.Along);
                    });
    }
  }
} // namespace Geoxy::Met
