/***************************************************************************************************
 *  @file      source/met/ellipse.cpp
 *  @brief     Two parts met where one is an elliptic arc: with a piece, with a run, with another
 *             elliptic arc, of one ellipse or of two.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "met/ellipse.hpp"

#include "elliptic.hpp"
#include "kernel/piece.hpp"
#include "met/floor.hpp"
#include "met/parts.hpp"
#include "met/solve.hpp"
#include "real.hpp"
#include "scratch.hpp"
#include "straight.hpp"
#include "trig.hpp"

#include "geoxy/angle.hpp"
#include "geoxy/at.hpp"
#include "geoxy/domain.hpp"
#include "geoxy/ellipse.hpp"
#include "geoxy/intersect.hpp"
#include "geoxy/matrix.hpp"
#include "geoxy/nearest.hpp"
#include "geoxy/place.hpp"
#include "geoxy/point.hpp"
#include "geoxy/tolerance.hpp"
#include "geoxy/vector.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <optional>
#include <span>
#include <utility>

namespace Geoxy::Met
{
  namespace
  {
    /**
     *  The candidate meetings of a piece and an ellipse, solved on the ellipse's quarters as
     *  pieces, their parameters turned to the ellipse's angles and passed in every turn.
     *
     *  @param  piece     The piece.
     *  @param  ellipse   The ellipse.
     *  @param  tolerance How near the parts must come.
     *  @return           What they found, the piece first.
     */
    auto Quartered(const Piece &piece, const Ellipse &ellipse, const Tolerance &tolerance) -> Found
    {
      Found               found;
      std::array<bool, 4> spent = {};
      for (std::size_t index = 0; index < 4; ++index)
      {
        const Found quarter = Against(piece, Kernel::Quarter(ellipse.Map, index), tolerance);
        spent.at(index)     = quarter.Spent;
        for (const Local &local : quarter.Locals)
        {
          for (const double along : Laps(ellipse, Kernel::Quartered(index, local.Second)))
          {
            found.Locals.push_back({.First = local.First, .Second = along});
          }
        }
      }
      found.Spent = std::ranges::contains(spent, true);
      return found;
    }

  } // namespace

  auto Overlapped(const Ellipse &first, const Ellipse &second, double offset) -> Found
  {
    Found        found;
    const double whole       = Turn.InRadians();
    const double firstStart  = first.Start.InRadians();
    const double firstSweep  = first.Sweep.InRadians();
    const double secondStart = second.Start.InRadians() + offset;
    const double secondSweep = second.Sweep.InRadians();
    if (Zero(firstSweep) || Zero(secondSweep))
    {
      return found;
    }
    const double lowFirst   = std::min(firstStart, firstStart + firstSweep);
    const double highFirst  = std::max(firstStart, firstStart + firstSweep);
    const double lowSecond  = std::min(secondStart, secondStart + secondSweep);
    const double highSecond = std::max(secondStart, secondStart + secondSweep);
    const double earliest   = std::ceil((lowFirst - highSecond) / whole);
    const double latest     = std::floor((highFirst - lowSecond) / whole);
    const double laps       = latest - earliest + 1.0;
    /**
     *  An angle that is not finite makes the count of turns no number, and no turn is tried.
     */
    const auto count =
        static_cast<std::size_t>(laps >= 0.0 ? std::min(laps, static_cast<double>(MostLaps)) : 0.0);
    for (std::size_t turn = 0; turn < count; ++turn)
    {
      const double shift   = (earliest + static_cast<double>(turn)) * whole;
      const double low     = std::max(lowFirst, lowSecond + shift);
      const double high    = std::min(highFirst, highSecond + shift);
      const auto   onFirst = [firstStart, firstSweep](double angle) -> double
      {
        return (angle - firstStart) / firstSweep;
      };
      const auto onSecond = [secondStart, secondSweep, shift](double angle) -> double
      {
        return (angle - shift - secondStart) / secondSweep;
      };
      if (low < high)
      {
        const double a = onFirst(low);
        const double b = onFirst(high);
        const double c = onSecond(low);
        const double d = onSecond(high);
        found.Shared.push_back({.First  = {.Start = std::min(a, b), .End = std::max(a, b)},
                                .Second = {.Start = std::min(c, d), .End = std::max(c, d)},
                                .Along  = (firstSweep > 0.0) == (secondSweep > 0.0)});
      }
      else
      {
        found.Locals.push_back({.First = onFirst(low), .Second = onSecond(low)});
      }
    }
    return found;
  }

  namespace
  {
    /**
     *  The stretches a piece lying on an ellipse shares with it: the arc of the ellipse the piece
     *  runs over, from its start's angle to its end's the way its middle lies, met with the
     *  ellipse's own sweep, its stretches turned back to the piece's parameter.
     *
     *  @param  piece   The piece, on the ellipse.
     *  @param  ellipse The ellipse.
     *  @param  unit    The map that takes the ellipse to the unit circle.
     *  @return         The stretches shared, the piece first.
     */
    auto Lying(const Piece &piece, const Ellipse &ellipse, const Matrix &unit) -> Found
    {
      const double whole  = Turn.InRadians();
      const auto   turned = [whole](double angle) -> double
      {
        return angle - (whole * std::floor(angle / whole));
      };
      const double  start  = Angled(unit, Kernel::Start(piece));
      const double  ahead  = turned(Angled(unit, Kernel::End(piece)) - start);
      const double  middle = turned(Angled(unit, Evaluated(piece, 0.5)) - start);
      const double  sweep  = middle <= ahead ? ahead : ahead - whole;
      const Ellipse arc    = {.Map = ellipse.Map, .Start = Radians(start), .Sweep = Radians(sweep)};
      Found         found;
      for (const Overlap &overlap : Overlapped(arc, ellipse, 0.0).Shared)
      {
        const auto back = [&piece, &arc](double along) -> double
        {
          return Kernel::Closest(piece, At(arc, along));
        };
        const double from  = back(overlap.First.Start);
        const double until = back(overlap.First.End);
        found.Shared.push_back(
            {.First  = {.Start = std::min(from, until), .End = std::max(from, until)},
             .Second = overlap.Second,
             .Along  = overlap.Along});
      }
      return found;
    }

    /**
     *  Stretches shared no longer than four tolerances made points: where two arcs of one circle
     *  meet end to end, their ends a rounding apart, the hit at its middle, as two curves are found
     *  to lie along each other only over more.
     *
     *  @param  found  What was found, along the first curve's parameter.
     *  @param  scale  How long the first curve is per unit of its parameter.
     *  @param  length The tolerance's length.
     *  @return        What was found, short stretches made points.
     */
    auto Pointed(Found found, double scale, double length) -> Found
    {
      std::erase_if(
          found.Shared,
          [&found, scale, length](const Overlap &overlap) -> bool
          {
            if ((overlap.First.End - overlap.First.Start) * scale > 4.0 * length)
            {
              return false;
            }
            found.Locals.push_back(
                {.First = overlap.First.Start + ((overlap.First.End - overlap.First.Start) / 2),
                 .Second =
                     overlap.Second.Start + ((overlap.Second.End - overlap.Second.Start) / 2)});
            return true;
          });
      return found;
    }

    /**
     *  Whether a piece lies along an ellipse's carrier within a length, as Along asks of two
     *  pieces: at nine places of it, the carrier's nearest point within the length. A place is
     *  settled by the carrier's point on the ray from its center, no nearer than the nearest; one
     *  that point leaves open, by the nearest point, unless the least stretch of the map times how
     *  far the place is off the unit circle in its frame, no farther than the nearest, already puts
     *  it out. The middle is tried first, as most pieces near an ellipse are nowhere along it.
     *
     *  @pre            The ellipse's map is invertible.
     *  @param  piece   The piece.
     *  @param  ellipse The ellipse.
     *  @param  length  The length.
     *  @return         Whether it does.
     */
    auto Hugging(const Piece &piece, const Ellipse &ellipse, double length) -> bool
    {
      constexpr int                Places = 8;
      constexpr std::array<int, 9> Order  = {4, 0, 8, 2, 6, 1, 3, 5, 7};
      const Matrix                 unit   = Inverse(ellipse.Map).value_or(Matrix{});
      const double                 least  = Stretches(ellipse.Map)[0];
      const Ellipse                whole  = {.Map = ellipse.Map, .Start = {}, .Sweep = Turn};
      return std::ranges::all_of(Order,
                                 [&piece, &ellipse, &unit, &whole, least, length](int place) -> bool
                                 {
                                   const Point point =
                                       Evaluated(piece, static_cast<double>(place) / Places);
                                   const Vector seen   = (unit * point) - Point{};
                                   const double radius = Length(seen);
                                   if (!(radius > 0.0) || least * std::abs(radius - 1.0) > length)
                                   {
                                     return false;
                                   }
                                   const Point foot = ellipse.Map * (Point{} + (seen / radius));
                                   return Distance(point, foot)                     <= length
                                       || Distance(point, Nearest(whole, point).At) <= length;
                                 });
    }

  } // namespace

  namespace
  {
    /**
     *  Whether an ellipse's map is a circle's: its axes square and as long.
     *
     *  @param  ellipse The ellipse.
     *  @return         Whether it is a circle, exactly.
     */
    auto Circular(const Ellipse &ellipse) -> bool
    {
      return Zero(Dot(ellipse.Map.X, ellipse.Map.Y))
          && Exactly(Dot(ellipse.Map.X, ellipse.Map.X), Dot(ellipse.Map.Y, ellipse.Map.Y))
          && Cross(ellipse.Map.X, ellipse.Map.Y) > 0.0;
    }

    /**
     *  A run, or an arc of a circle, as Hugged reads it: the run; or the arc, the map that takes it
     *  to the unit circle, its center and its radius.
     */
    struct Girdle
    {
      const Ellipse *Arc    = nullptr; /**< The arc, or nothing for a run.                 */
      const Run     *Line   = nullptr; /**< The run, or nothing for an arc.                */
      Matrix         Unit   = {};      /**< The map that takes the arc to the unit circle. */
      Point          Center = {};      /**< The arc's center.                              */
      double         Radius = 0.0;     /**< The arc's radius.                              */
    };

    /**
     *  An arc of a circle as Hugged reads it.
     *
     *  @param  arc The arc, of a circle.
     *  @return     It.
     */
    auto Girded(const Ellipse &arc) -> Girdle
    {
      return {.Arc    = &arc,
              .Line   = nullptr,
              .Unit   = Inverse(arc.Map).value_or(Matrix{}),
              .Center = arc.Map.Origin,
              .Radius = Length(arc.Map.X)};
    }

    /**
     *  A point of a run or an arc.
     *
     *  @param  girdle The run or arc.
     *  @param  along  Its parameter.
     *  @return        The point.
     */
    auto Plot(const Girdle &girdle, double along) -> Point
    {
      if (girdle.Arc != nullptr)
      {
        return At(*girdle.Arc, along);
      }
      return girdle.Line->Start + (girdle.Line->Way * along);
    }

    /**
     *  How far a point lies from a run's line or an arc's circle.
     *
     *  @param  girdle The run or arc.
     *  @param  point  The point.
     *  @return        The distance.
     */
    auto Fringe(const Girdle &girdle, Point point) -> double
    {
      if (girdle.Arc != nullptr)
      {
        return std::abs(Distance(point, girdle.Center) - girdle.Radius);
      }
      const Run &run = *girdle.Line;
      return std::abs(Cross(point - run.Start, run.Way)) / Length(run.Way);
    }

    /**
     *  Whether an angle of an arc's unit circle, any turn of it, lies within a stretch of its
     *  sweep.
     *
     *  @param  girdle  The arc.
     *  @param  stretch The stretch of its parameter.
     *  @param  angle   The angle, in radians.
     *  @return         Whether it does.
     */
    auto Brims(const Girdle &girdle, Range stretch, double angle) -> bool
    {
      const double start   = girdle.Arc->Start.InRadians();
      const double sweep   = girdle.Arc->Sweep.InRadians();
      const double from    = start + (sweep * stretch.Start);
      const double until   = start + (sweep * stretch.End);
      const double low     = std::min(from, until);
      const double shifted = angle + (Revolution * std::ceil((low - angle) / Revolution));
      return shifted <= std::max(from, until);
    }

    /**
     *  Where on a run or an arc a point within a length of its line or circle stands abreast: the
     *  foot on a run, the place of its angle on an arc, the first pass where a whole circle passes
     *  it twice.
     *
     *  @param  girdle The run or arc.
     *  @param  point  The point.
     *  @param  length The length.
     *  @return        The place; nothing where the point lies farther from the line or circle, or
     *                 beside no place of the run or arc.
     */
    auto Fringed(const Girdle &girdle, Point point, double length) -> std::optional<double>
    {
      if (!(Fringe(girdle, point) <= length))
      {
        return std::nullopt;
      }
      if (girdle.Arc != nullptr)
      {
        return Earliest(*girdle.Arc, Angled(girdle.Unit, point));
      }
      return Footing(*girdle.Line, point);
    }

    /**
     *  How far a stretch of a run or an arc strays from the other's line or circle at most: at its
     *  ends, or where between them it comes nearest or goes farthest, as the distance to a circle's
     *  center grows with the angle from it, a run comes nearest a center at its foot, and an arc
     *  goes farthest to either side of a line where its tangent is the line's.
     *
     *  @param  girdle  The run or arc.
     *  @param  stretch The stretch of its parameter.
     *  @param  against The other, whose line or circle.
     *  @return         The distance.
     */
    auto Roves(const Girdle &girdle, Range stretch, const Girdle &against) -> double
    {
      double most = std::max(Fringe(against, Plot(girdle, stretch.Start)),
                             Fringe(against, Plot(girdle, stretch.End)));
      if (girdle.Arc == nullptr)
      {
        const std::optional<double> foot = Footing(*girdle.Line, against.Center);
        if (foot && stretch.Start < *foot && *foot < stretch.End)
        {
          most = std::max(most, Fringe(against, Plot(girdle, *foot)));
        }
        return most;
      }
      const Point  across = against.Arc != nullptr
                              ? against.Center
                              : girdle.Center + Perpendicular(against.Line->Way);
      const double toward = Angled(girdle.Unit, across);
      for (const double angle : {toward, toward + Pi})
      {
        if (Brims(girdle, stretch, angle))
        {
          const Point point = girdle.Arc->Map * (Point{} + Direction(Radians(angle)));
          most              = std::max(most, Fringe(against, point));
        }
      }
      return most;
    }

    /**
     *  How long a stretch of a run or an arc is.
     *
     *  @param  girdle  The run or arc.
     *  @param  stretch The stretch of its parameter.
     *  @return         Its length.
     */
    auto Furlong(const Girdle &girdle, Range stretch) -> double
    {
      const double part = stretch.End - stretch.Start;
      if (girdle.Arc != nullptr)
      {
        return part * girdle.Radius * std::abs(girdle.Arc->Sweep.InRadians());
      }
      return part * Length(girdle.Line->Way);
    }

    /**
     *  The ends of a run that it has, those that are finite, or of an arc, as places of it.
     *
     *  @param  girdle The run or arc.
     *  @return        Their places.
     */
    auto Tipmost(const Girdle &girdle) -> List<double>
    {
      if (girdle.Arc != nullptr)
      {
        return {0.0, 1.0};
      }
      List<double> tips;
      for (const double end : {girdle.Line->Along.Start, girdle.Line->Along.End})
      {
        if (std::isfinite(end))
        {
          tips.push_back(end);
        }
      }
      return tips;
    }

    /**
     *  An end of either of two curves that lies within a length of the other, as a place of each.
     */
    struct Brink
    {
      double First  = 0.0;   /**< Its place on the first.                             */
      double Second = 0.0;   /**< Its place on the second.                            */
      bool   Seam   = false; /**< Whether it is a whole second circle's seam, 0 or 1. */
    };

    /**
     *  The ends of either of a run or an arc and an arc of another circle that lie within a length
     *  of the other, in order along the first: the first's start, the second's, the first's end,
     *  the second's, before sorting, so that ends at one place keep that order.
     *
     *  @param  first  The run or arc.
     *  @param  second The arc.
     *  @param  length The length.
     *  @return        The ends.
     */
    auto Brinks(const Girdle &first, const Girdle &second, double length) -> List<Brink>
    {
      const bool         round = !(std::abs(second.Arc->Sweep.InRadians()) < Revolution);
      const List<double> tips  = Tipmost(first);
      List<Brink>        ends;
      for (std::size_t index = 0; index < 2; ++index)
      {
        if (index < tips.size())
        {
          const double tip = tips[index];
          if (const std::optional<double> across = Fringed(second, Plot(first, tip), length))
          {
            ends.push_back({.First = tip, .Second = *across, .Seam = false});
          }
        }
        const double end = index == 0 ? 0.0 : 1.0;
        if (const std::optional<double> along = Fringed(first, Plot(second, end), length))
        {
          ends.push_back({.First = *along, .Second = end, .Seam = round});
        }
      }
      std::ranges::stable_sort(ends, {}, [](const Brink &end) -> double { return end.First; });
      return ends;
    }

    /**
     *  The stretches a run or an arc of a circle shares with an arc of another circle that it parts
     *  from: between two ends, of either, each within the length of the other, where both stay
     *  within the length of the other's line or circle all the way, longer than four lengths, as
     *  two arcs of one circle are found to share one. Such a stretch lies about where the two come
     *  nearest or cross shallowly, each of its ends inside the band where they stay within the
     *  length, and where one of them ends; one that runs out of the band ends where they part
     *  instead, and is the band about a touch or a shallow crossing, a hit, which the ends and the
     *  crossings find. Each end is placed on an arc by its angle and on a run by its foot, the same
     *  whichever is first; a whole circle's seam bounds a stretch, its place 0 or 1 as the stretch
     *  runs. Arcs swept round more than once are left to their passes.
     *
     *  @param  first  One run, or arc of a circle.
     *  @param  second The other, an arc of another circle.
     *  @param  length How near they must stay.
     *  @return        The stretches, along the first.
     */
    auto Hugged(const Girdle &first, const Girdle &second, double length) -> List<Overlap>
    {
      List<Overlap> shared;
      const auto    coiled = [](const Girdle &girdle) -> bool
      {
        return girdle.Arc != nullptr && Coiled(*girdle.Arc);
      };
      if (coiled(first) || coiled(second))
      {
        return shared;
      }
      const List<Brink> ends    = Brinks(first, second, length);
      const auto        choices = [](const Brink &end) -> std::array<double, 2>
      {
        return end.Seam ? std::array{0.0, 1.0} : std::array{end.Second, end.Second};
      };
      for (std::size_t index = 0; index + 1 < ends.size(); ++index)
      {
        const Brink                &from  = ends[index];
        const Brink                &until = ends[index + 1];
        const Range                 lap   = {.Start = from.First, .End = until.First};
        const std::optional<double> middle =
            Fringed(second, Plot(first, lap.Start + ((lap.End - lap.Start) / 2)), length);
        if (!middle)
        {
          continue;
        }
        /**
         *  Both stretches are weighed alike, the shorter and the farther of the two, so that the
         *  same stretches are found whichever is first.
         */
        for (const double start : choices(from))
        {
          for (const double stop : choices(until))
          {
            const Range  back = {.Start = std::min(start, stop), .End = std::max(start, stop)};
            const double span = std::min(Furlong(first, lap), Furlong(second, back));
            const double off  = std::max(Roves(first, lap, second), Roves(second, back, first));
            if (back.Start <= *middle
                && *middle <= back.End
                && span    > 4.0 * length
                && off     <= length)
            {
              shared.push_back({.First = lap, .Second = back, .Along = start < stop});
            }
          }
        }
      }
      return shared;
    }

  } // namespace

  auto Against(const Piece &piece, const Ellipse &ellipse, const Tolerance &tolerance) -> Found
  {
    const std::optional<Matrix> unit = Inverse(ellipse.Map);
    if (!unit)
    {
      return Quartered(piece, ellipse, tolerance);
    }
    const Kernel::Zeros zeros = Kernel::Zeroed(piece, *unit, Kernel::Budget);
    if (zeros.On)
    {
      return Lying(piece, ellipse, *unit);
    }
    if (Hugging(piece, ellipse, tolerance.Length))
    {
      /**
       *  A piece within the tolerance of the ellipse at nine places lies along it as two pieces do,
       *  though its implicit equation there is not nothing: the stretches are found as for one on
       *  it, through the angles of its ends and its middle. One no longer than four tolerances, as
       *  where the piece starts a rounding before a whole ellipse's seam, is a point, as for two
       *  circles within the tolerance: measured along the ellipse, no faster than its largest
       *  stretch times its sweep.
       */
      const double speed = Stretches(ellipse.Map)[1] * std::abs(ellipse.Sweep.InRadians());
      return Swapped(Pointed(Swapped(Lying(piece, ellipse, *unit)), speed, tolerance.Length));
    }
    Found      found;
    const auto add = [&found, &piece, &ellipse, &unit](double along) -> void
    {
      for (const double pass : Laps(ellipse, Angled(*unit, Evaluated(piece, along))))
      {
        found.Locals.push_back({.First = along, .Second = pass});
      }
    };
    for (const double along : zeros.Crossing.Roots)
    {
      add(along);
    }
    for (const double along : zeros.Turns)
    {
      add(along);
    }
    return found;
  }

  auto Against(const Run &run, const Ellipse &ellipse, const Tolerance &tolerance) -> Found
  {
    Found        found;
    const double square = Dot(run.Way, run.Way);
    if (Zero(square))
    {
      return found;
    }
    const std::optional<Matrix> unit = Inverse(ellipse.Map);
    if (!unit)
    {
      for (std::size_t index = 0; index < 4; ++index)
      {
        const Piece quarter = Kernel::Quarter(ellipse.Map, index);
        const Found met     = Swapped(Against(run, quarter, tolerance));
        for (const Local &local : met.Locals)
        {
          for (const double along : Laps(ellipse, Kernel::Quartered(index, local.First)))
          {
            found.Locals.push_back({.First = local.Second, .Second = along});
          }
        }
      }
      return found;
    }
    if (Circular(ellipse))
    {
      const Girdle along = {.Arc = nullptr, .Line = &run, .Unit = {}, .Center = {}, .Radius = 0.0};
      found.Shared       = Hugged(along, Girded(ellipse), tolerance.Length);
    }
    const Vector start = (*unit * run.Start) - Point{};
    const Vector way   = *unit * run.Way;
    const double a     = Dot(way, way);
    const double b     = 2.0 * Dot(start, way);
    const double c     = Dot(start, start) - 1.0;
    /**
     *  A run so short in the circle's frame that its length squared is nothing, as beside a radius
     *  near the largest double, is a run of no length there: it has no crossing to work out, and
     *  its ends are weighed as any curve's are.
     */
    if (!(a > 0.0))
    {
      return found;
    }
    const auto add = [&found, &run, &ellipse, &unit](double along) -> void
    {
      const Point point = run.Start + (run.Way * along);
      for (const double pass : Laps(ellipse, Angled(*unit, point)))
      {
        found.Locals.push_back({.First = along, .Second = pass});
      }
    };
    const double discriminant = (b * b) - (4.0 * a * c);
    if (discriminant >= 0.0)
    {
      const double q = -0.5 * (b + std::copysign(std::sqrt(discriminant), b));
      add(q / a);
      if (!Zero(q))
      {
        add(c / q);
      }
    }
    /**
     *  The closest approach is a candidate where the run comes near the ellipse without crossing
     *  it: where the carrier's point nearest the center, in the unit circle's frame, lies off the
     *  circle by more than four tolerances, as the least stretch tells, with the rounding Remote
     *  allows for, none is sought. Outside, every point of the run lies as far from the ellipse,
     *  and no candidate there is a hit; inside, the run cuts deep through the ellipse, and its
     *  crossings are the roots above.
     */
    const auto [least, most] = Stretches(ellipse.Map);
    const double outside     = std::abs(Cross(start, way)) / std::sqrt(a);
    const double clear = (least * std::abs(outside - 1.0)) - (Absolute * most * (outside + 1.0));
    const double slack = Absolute * (Magnitude(run.Start) + Magnitude(ellipse.Map.Origin) + most);
    if (clear * (1.0 - Relative) > (4.0 * tolerance.Length) + slack)
    {
      return found;
    }
    const double nearest = -b / (2.0 * a);
    const Place  foot    = Nearest(ellipse, run.Start + (run.Way * nearest));
    add(Dot(foot.At - run.Start, run.Way) / square);
    return found;
  }

  namespace
  {
    /**
     *  The candidate meetings at points found on both of two ellipses, every turn of each: no more
     *  than MostPairs of them at a point, what was found marked spent where there were more.
     *
     *  @param found  Where to add them.
     *  @param passes The first ellipse's parameters there, and the second's.
     */
    void Paired(Found &found, const std::array<List<double>, 2> &passes)
    {
      std::size_t made = 0;
      for (const double pass : passes[0])
      {
        for (const double meet : passes[1])
        {
          if (made == MostPairs)
          {
            found.Spent = true;
            return;
          }
          found.Locals.push_back({.First = pass, .Second = meet});
          ++made;
        }
      }
    }

    /**
     *  Where two circles meet, in closed form: two within a length of each other all the way round,
     *  their centers and radii apart by that much together at most, as two rounded from one circle
     *  are, share stretches; otherwise the stretches Hugged finds between their arcs' ends, and the
     *  chord of their crossing, or, where they do not reach each other, the points on the line of
     *  their centers.
     *
     *  @pre           Both are circles, of a radius above 0 and so invertible.
     *  @param  first  One circle.
     *  @param  second The other.
     *  @param  length How near they must stay to be one.
     *  @return        What they found.
     */
    auto Circles(const Ellipse &first, const Ellipse &second, double length) -> Found
    {
      const Point  center = first.Map.Origin;
      const Point  other  = second.Map.Origin;
      const double radius = Length(first.Map.X);
      const double outer  = Length(second.Map.X);
      if (Distance(center, other) + std::abs(radius - outer) <= length)
      {
        Found shared = Overlapped(
            first, second, Heading(second.Map.X).InRadians() - Heading(first.Map.X).InRadians());
        return Pointed(std::move(shared), radius * std::abs(first.Sweep.InRadians()), length);
      }
      Found found;
      found.Shared       = Hugged(Girded(first), Girded(second), length);
      const double apart = Distance(center, other);
      if (Zero(apart))
      {
        return found;
      }
      /**
       *  The chord is worked from the smaller circle where, from the larger, its rounding may reach
       *  a sixteenth of the length: a chord much shorter than the larger radius is the difference
       *  of two squares all but equal, off by many tolerances, and a crossing near an end of an arc
       *  was taken for one past it.
       */
      const double grain      = std::numeric_limits<double>::epsilon() * radius * radius;
      const bool   smaller    = outer < radius && grain > outer * length / 16.0;
      const double near       = smaller ? outer : radius;
      const double far        = smaller ? radius : outer;
      const Vector toward     = (other - center) / apart;
      const double along      = ((apart * apart) + (near * near) - (far * far)) / (2 * apart);
      const double squared    = (near * near) - (along * along);
      const Matrix firstUnit  = Inverse(first.Map).value_or(Matrix{});
      const Matrix secondUnit = Inverse(second.Map).value_or(Matrix{});
      const auto   add = [&found, &first, &second, &firstUnit, &secondUnit](Point point) -> void
      {
        Paired(found,
               {Laps(first, Angled(firstUnit, point)), Laps(second, Angled(secondUnit, point))});
      };
      if (squared >= 0.0)
      {
        const double height = std::sqrt(squared);
        const Point  base   = smaller ? other - (toward * along) : center + (toward * along);
        add(base + (Perpendicular(toward) * height));
        add(base - (Perpendicular(toward) * height));
      }
      else
      {
        add(center + (toward * radius));
        add(center - (toward * radius));
      }
      return found;
    }

    /**
     *  Where an ellipse meets another by the other's implicit equation on its quarters, chained
     *  round, every turn of each.
     *
     *  @param  cut   The ellipse cut into quarters.
     *  @param  whole The other.
     *  @param  unit  The map that takes the other to the unit circle.
     *  @return       What they found.
     */
    auto Implicit(const Ellipse &cut, const Ellipse &whole, const Matrix &unit) -> Found
    {
      std::array<Kernel::Solution, 4> quarters;
      std::array<List<double>, 4>     turns;
      for (std::size_t index = 0; index < quarters.size(); ++index)
      {
        const Kernel::Zeros zeros =
            Kernel::Zeroed(Kernel::Quarter(cut.Map, index), unit, Kernel::Budget);
        quarters.at(index) = zeros.Crossing;
        turns.at(index)    = zeros.Turns;
      }
      Found      found;
      const auto add = [&found, &cut, &whole, &unit](std::size_t index, double along) -> void
      {
        const double angle = Kernel::Quartered(index, along);
        const Point  point = cut.Map * (Point{} + Direction(Radians(angle)));
        Paired(found, {Laps(cut, angle), Laps(whole, Angled(unit, point))});
      };
      for (const Kernel::Spot spot : Kernel::Chained(quarters, true))
      {
        add(spot.Index, spot.Along);
      }
      for (std::size_t index = 0; index < turns.size(); ++index)
      {
        for (const double along : turns.at(index))
        {
          add(index, along);
        }
      }
      return found;
    }

    /**
     *  Whether two ellipses lie within a length of each other all the way round, point for point:
     *  the first's point at an angle and the second's at that angle turned by one offset. The two
     *  points differ by the gap of the centers and by the first's map less the second's turned,
     *  applied to the angle's direction; that is no longer than the gap and the map's largest
     *  stretch together, and so no longer than the gap and the root of the sum of its four entries'
     *  squares. The turn is the one nearest the second's map undone after the first's. Where the
     *  two maps turn opposite ways, the second is spelled the other way round first: its axis y
     *  turned back and its angles negated, which is the same curve, place for place.
     *
     *  @param  first  One ellipse.
     *  @param  second The other.
     *  @param  length How near they must stay to be one.
     *  @return        What to add to the second's angles to have the first's, and the second as
     *                 spelled to have them; nothing where they part by more, or either flattens.
     */
    auto Coinciding(const Ellipse &first, const Ellipse &second, double length)
        -> std::optional<std::pair<double, Ellipse>>
    {
      const bool    opposite = Determinant(first.Map) * Determinant(second.Map) < 0.0;
      const Ellipse turned =
          opposite
              ? Ellipse{.Map = {.X = second.Map.X, .Y = -second.Map.Y, .Origin = second.Map.Origin},
                        .Start = -second.Start,
                        .Sweep = -second.Sweep}
              : second;
      const std::optional<Matrix> unit = Inverse(turned.Map);
      if (!unit || !Inverse(first.Map))
      {
        return std::nullopt;
      }
      const Vector x      = *unit * first.Map.X;
      const Vector y      = *unit * first.Map.Y;
      const double angle  = Trig::Atan2(x.Y - y.X, x.X + y.Y);
      const Vector way    = Trig::Way(angle);
      const Vector across = first.Map.X - (turned.Map * way);
      const Vector upward = first.Map.Y - (turned.Map * Perpendicular(way));
      const double apart  = Distance(first.Map.Origin, turned.Map.Origin)
                          + std::sqrt(Dot(across, across) + Dot(upward, upward));
      if (!(apart <= length))
      {
        return std::nullopt;
      }
      return std::pair{-angle, turned};
    }

  } // namespace

  auto Against(const Ellipse &first, const Ellipse &second, const Tolerance &tolerance) -> Found
  {
    if (first.Map == second.Map)
    {
      return Overlapped(first, second, 0.0);
    }
    if (Circular(first) && Circular(second))
    {
      return Circles(first, second, tolerance.Length);
    }
    /**
     *  Two ellipses within the tolerance of each other all the way round share stretches as two on
     *  one carrier do, though their maps differ in a last bit: their implicit equations are nothing
     *  nowhere on each other, and found only where they cross.
     */
    if (const auto coinciding = Coinciding(first, second, tolerance.Length))
    {
      /**
       *  The angles turned by the offset round to places a step of a double past an arc's end: each
       *  is held to its curve's domain.
       */
      const double speed = Stretches(first.Map)[1] * std::abs(first.Sweep.InRadians());
      Found        found = Pointed(Overlapped(first, coinciding->second, coinciding->first), speed,
                                   tolerance.Length);
      const auto   held  = [](Range range) -> Range
      {
        return {.Start = std::clamp(range.Start, 0.0, 1.0), .End = std::clamp(range.End, 0.0, 1.0)};
      };
      for (Overlap &overlap : found.Shared)
      {
        overlap = {
            .First = held(overlap.First), .Second = held(overlap.Second), .Along = overlap.Along};
      }
      return found;
    }
    if (const std::optional<Matrix> unit = Inverse(second.Map))
    {
      return Implicit(first, second, *unit);
    }
    if (const std::optional<Matrix> unit = Inverse(first.Map))
    {
      return Swapped(Implicit(second, first, *unit));
    }
    return {.Locals = {}, .Shared = {}, .Spent = true};
  }
} // namespace Geoxy::Met
