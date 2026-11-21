/***************************************************************************************************
 *  @file      source/abut.cpp
 *  @brief     Two ready curves where an end of one lies on the other, told from what they keep as
 *             Intersect tells them, to the bit, or left to Intersect.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "abut.hpp"

#include "axis.hpp"
#include "crossing.hpp"
#include "lean.hpp"
#include "met/linear.hpp"
#include "real.hpp"
#include "straight.hpp"
#include "thorough.hpp"
#include "trig.hpp"
#include "twofold.hpp"

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
#include "geoxy/tolerance.hpp"
#include "geoxy/vector.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <optional>
#include <span>
#include <variant>

namespace Geoxy
{
  namespace
  {
    /**
     *  How far a number must lie past what Intersect compares it with, as a factor: more than any
     *  few roundings either side can close.
     */
    constexpr double Margin = 1.0 + 0x1p-40;

    /**
     *  A length, and a slack either side of it: what is nearer than the length less the slack is
     *  near, what is farther than the length and the slack is far, and the rest is not told.
     */
    struct Scope
    {
      double Length = 0.0; /**< The length. */
      double Slack  = 0.0; /**< The slack.  */
    };

    /**
     *  The square of a distance, rounded: within 4 u of it, each difference, product and sum
     *  rounded once.
     *
     *  @param  from  One point.
     *  @param  until The other.
     *  @return       The square.
     */
    auto Squared(Point from, Point until) -> double
    {
      const Vector apart = until - from;
      return Dot(apart, apart);
    }

    /**
     *  How a distance squared, rounded, stands to a reach: 1 near, -1 far, 0 not told; a square's
     *  roundings either side taken by Tight and Loose.
     *
     *  @param  square The distance squared.
     *  @param  reach  The reach.
     *  @return        Which.
     */
    auto Gauged(double square, Scope reach) -> int
    {
      const double near = reach.Length - reach.Slack;
      const double far  = reach.Length + reach.Slack;
      if (square <= near * near * Lean::Tight)
      {
        return 1;
      }
      return square > far * far * Lean::Loose ? -1 : 0;
    }

    /**
     *  Whether a point lies within a reach of one of some, told: 1 where one is near, -1 where
     *  every one is far, 0 where neither is told.
     *
     *  @param  ends  The points.
     *  @param  point The point.
     *  @param  reach The reach.
     *  @return       Which.
     */
    auto Nigh(std::span<const Point> ends, Point point, Scope reach) -> int
    {
      int told = -1;
      for (const Point end : ends)
      {
        told = std::max(told, Gauged(Squared(end, point), reach));
      }
      return told;
    }

    /**
     *  A coordinate of a point.
     *
     *  @param  point The point.
     *  @param  axis  Which: 0 for x, 1 for y.
     *  @return       It.
     */
    auto Coordinate(Point point, std::size_t axis) -> double
    {
      return axis == 0 ? point.X : point.Y;
    }

    /**
     *  A level or upright line and an arc as Rimmed reads them, and what it tells them within.
     */
    struct Rim
    {
      std::array<Point, 2> Ends   = {};      /**< The line's ends.                            */
      std::array<Point, 2> Cuts   = {};      /**< The arc's ends, as it keeps them.           */
      const Arc           *Round  = nullptr; /**< The arc.                                    */
      double               Sweep  = 0.0;     /**< Its sweep, in radians, as it keeps it.      */
      std::size_t          Fixed  = 0;       /**< The axis the line holds fixed: 1 for level. */
      double               Length = 0.0;     /**< The length they are met within.             */
      double               Slack  = 0.0;     /**< How far a number may be from Intersect's.   */
    };

    /**
     *  The crossings of a carrier and a circle, and whether each is the double nearest its own.
     */
    struct Spots
    {
      std::array<Point, 2> At = {};  /**< The crossings, in the order along the line Roundly
                                          finds them.                                      */
      std::array<bool, 2> Sure = {}; /**< Whether each is certain. */
    };

    /**
     *  The ends of a line and an arc at a crossing, as Roundly keeps them.
     */
    struct Held
    {
      std::array<Point, 4> Ends  = {};    /**< The ends.                          */
      std::size_t          Count = 0;     /**< How many there are.                */
      bool                 Line  = false; /**< Whether one of them is the line's. */

      /**
       *  One more end.
       *
       *  @param end The end.
       */
      void Add(Point end)
      {
        Ends.at(Count) = end;
        ++Count;
      }
    };

    /**
     *  The crossings of a level or upright line's carrier and an arc's circle, c -+ sqrt(r^2 - g^2)
     *  along the line, c the center's coordinate there and g the carrier's gap from the center,
     *  each one exact sum certified as Levelled certifies it, in the order along the line Roundly
     *  finds them.
     *
     *  @param  rim     The line and the arc.
     *  @param  squared The arc's radius squared, exactly.
     *  @return         The crossings; nothing where the carrier cuts the circle at a sine of 1 /
     *                  sqrt(32) or less.
     */
    auto Spotted(const Rim &rim, Lean::Double squared) -> std::optional<Spots>
    {
      const std::size_t  fixed   = rim.Fixed;
      const std::size_t  free    = 1 - fixed;
      const double       across  = Coordinate(rim.Ends[0], fixed);
      const double       centre  = Coordinate(rim.Round->Center, free);
      const Twofold::Two apart   = Twofold::TwoSum(across, -Coordinate(rim.Round->Center, fixed));
      const Lean::Double gap     = {.High = apart.High, .Low = apart.Low};
      const Lean::Double rest    = Lean::Plus(squared, Lean::Minus(Lean::Times(gap, gap)));
      const double       offRest = Lean::Grain * (squared.High + (gap.High * gap.High));
      if (!((rest.High * Lean::Tight) - offRest > squared.High * Lean::Loose / 32.0))
      {
        return std::nullopt;
      }
      const Lean::Rooted root = Lean::Root(rest);
      const double       offRoot =
          (offRest * root.Inverse * Lean::Loose) + (Lean::Grain * root.Value.High);
      const bool   rising = Coordinate(rim.Ends[1], free) > Coordinate(rim.Ends[0], free);
      const Vector way    = rim.Ends[1] - rim.Ends[0];
      Spots        spots;
      for (std::size_t side = 0; side < 2; ++side)
      {
        const double       sign  = (side == 0) == rising ? -1.0 : 1.0;
        const Twofold::Two sum   = Twofold::TwoSum(centre, sign * root.Value.High);
        const double       part  = sum.Low + (sign * root.Value.Low);
        const double       value = sum.High + part;
        const double bound = (2.0 * offRoot) + (Lean::Grain * (std::abs(centre) + root.Value.High));
        const bool   sure  = Lean::Certain(value, (sum.High - value) + part, bound);
        /**
         *  A crossing on an axis, the center's coordinate and the root's share cancelling exactly,
         *  Intersect rounds to +0, as Chord's Zeroed tells.
         */
        const int             asked = static_cast<int>(!sure)
                                    * static_cast<int>(!(std::abs(value) > bound))
                                    * static_cast<int>(centre * (sign * root.Value.High) < 0.0);
        std::array<double, 2> spot  = {};
        spot.at(fixed)              = across;
        const bool zero =
            asked > 0
            && static_cast<int>(Axis::Fits(Dot(way, way) * centre))
                       * static_cast<int>(Axis::On(*rim.Round, {.X = spot[0], .Y = spot[1]}))
                   > 0;
        spot.at(free)       = zero ? 0.0 : value;
        spots.At.at(side)   = {.X = spot[0], .Y = spot[1]};
        spots.Sure.at(side) = sure || zero;
      }
      return spots;
    }

    /**
     *  The point of a level or upright line nearest a point.
     *
     *  @param  rim   The line, and the arc.
     *  @param  point The point.
     *  @return       The nearest point.
     */
    auto Plummet(const Rim &rim, Point point) -> Point
    {
      const std::size_t     free  = 1 - rim.Fixed;
      const double          from  = Coordinate(rim.Ends[0], free);
      const double          until = Coordinate(rim.Ends[1], free);
      std::array<double, 2> foot  = {};
      foot.at(rim.Fixed)          = Coordinate(rim.Ends[0], rim.Fixed);
      foot.at(free) =
          std::clamp(Coordinate(point, free), std::min(from, until), std::max(from, until));
      return {.X = foot[0], .Y = foot[1]};
    }

    /**
     *  Every end at a crossing, kept, or far from the other curve, the line's from the circle and
     *  the arc's from the line, as Roundly settles them: each told, or nothing.
     *
     *  @param  rim   The line and the arc.
     *  @param  spots Their crossings.
     *  @return       The ends kept; nothing where an end is neither at a crossing nor far, or not
     *                told so.
     */
    auto Settled(const Rim &rim, const Spots &spots) -> std::optional<Held>
    {
      const double radius = rim.Round->Radius;
      const double clear  = 4.0 * rim.Length;
      const Scope  half   = {.Length = rim.Length / 2.0, .Slack = rim.Slack};
      const Scope  deep   = {.Length = radius - clear, .Slack = rim.Slack};
      const Scope  wide   = {.Length = radius + clear, .Slack = rim.Slack};
      const Scope  far    = {.Length = clear, .Slack = rim.Slack};
      Held         held;
      for (const Point end : rim.Ends)
      {
        const double square = Squared(end, rim.Round->Center);
        if (Nigh(spots.At, end, half) > 0)
        {
          held.Add(end);
          held.Line = true;
        }
        else if (Gauged(square, deep) < 1 && Gauged(square, wide) > -1)
        {
          return std::nullopt;
        }
      }
      for (const Point end : rim.Cuts)
      {
        if (Nigh(spots.At, end, half) > 0)
        {
          held.Add(end);
        }
        else if (Gauged(Squared(end, Plummet(rim, end)), far) > -1)
        {
          return std::nullopt;
        }
      }
      return held;
    }

    /**
     *  Whether a point lies within an arc's sweep, told: 1 within, -1 off, 0 where Intersect's test
     *  of a point near this one may answer either way. The test is Intersect's, on the turns from
     *  the arc's ends as it keeps them, each told by its sign where it lies past the slack times
     *  the radius: both turns must be told to tell a point within a sweep up to a half turn, one is
     *  enough to tell it off it, and the other way round for a longer sweep.
     *
     *  @param  rim   The line, and the arc, swept less than a whole turn.
     *  @param  point The point.
     *  @return       Which.
     */
    auto Swept(const Rim &rim, Point point) -> int
    {
      const Arc   &arc    = *rim.Round;
      const auto   back   = static_cast<std::size_t>(!(rim.Sweep > 0.0));
      const Vector start  = rim.Cuts.at(back) - arc.Center;
      const Vector end    = rim.Cuts.at(1 - back) - arc.Center;
      const Vector toward = point - arc.Center;
      const double edge   = rim.Slack * arc.Radius;
      const auto   sign   = [edge](double turn) -> int
      {
        return static_cast<int>(turn > edge) - static_cast<int>(turn < -edge);
      };
      const int left  = sign(Cross(start, toward));
      const int right = sign(Cross(toward, end));
      return std::abs(rim.Sweep) <= Pi ? std::min(left, right) : std::max(left, right);
    }

    /**
     *  Whether a place along a level or upright line lies within it, told: 1 within, -1 off, 0
     *  within the slack of an end.
     *
     *  @param  rim   The line, and the arc.
     *  @param  place The place, its free coordinate.
     *  @return       Which.
     */
    auto Housed(const Rim &rim, double place) -> int
    {
      const std::size_t free  = 1 - rim.Fixed;
      const double      from  = Coordinate(rim.Ends[0], free);
      const double      until = Coordinate(rim.Ends[1], free);
      const double      low   = std::min(from, until);
      const double      high  = std::max(from, until);
      if (place - low > rim.Slack && high - place > rim.Slack)
      {
        return 1;
      }
      return place < low - rim.Slack || place > high + rim.Slack ? -1 : 0;
    }

    /**
     *  Whether a crossing is a hit, told: on the line, within it or at an end of it, and on the
     *  arc, within its sweep or at an end of it, as Roundly asks.
     *
     *  @param  rim  The line and the arc.
     *  @param  spot The crossing.
     *  @return      1 where it is, -1 where it is not, 0 where it is not told.
     */
    auto Struck(const Rim &rim, Point spot) -> int
    {
      const Scope half = {.Length = rim.Length / 2.0, .Slack = rim.Slack};
      const int   onLine =
          std::max(Nigh(rim.Ends, spot, half), Housed(rim, Coordinate(spot, 1 - rim.Fixed)));
      const int onArc = std::max(Nigh(rim.Cuts, spot, half), Swept(rim, spot));
      return std::min(onLine, onArc);
    }

    /**
     *  How far past a line's ends a point on its carrier may lie and be its, as Intersect places a
     *  crossing: a few steps of a double.
     */
    constexpr double Overstep = 16.0 * std::numeric_limits<double>::epsilon();

    /**
     *  How far a point is from a run, as Intersect's Away works it out, by the same steps: to the
     *  foot on it, or to its nearer end.
     *
     *  @param  point The point.
     *  @param  run   The run.
     *  @return       The distance.
     */
    auto Away(Point point, const Run &run) -> double
    {
      const double square = Dot(run.Way, run.Way);
      const double along  = square > 0.0 ? std::clamp(Dot(point - run.Start, run.Way) / square,
                                                      run.Along.Start, run.Along.End)
                                         : run.Along.Start;
      return Distance(point, run.Start + (run.Way * along));
    }

    /**
     *  Whether a line lies wholly to one side of a run's carrier, farther from it than a distance,
     *  as Intersect's Beside tells it, by the same steps. Intersect's own asks too that the cross
     *  products be finite, which they are for every pair asked here: a length within Coarse, past
     *  sixteen steps of a double at the lines' size, keeps their numbers under 2^98.
     *
     *  @param  line  The line.
     *  @param  run   The run whose carrier.
     *  @param  least The distance times the length of the run's way.
     *  @return       Whether it does.
     */
    auto Beside(const Line &line, const Run &run, double least) -> bool
    {
      const double start  = Cross(line.Start - run.Start, run.Way);
      const double end    = Cross(line.End - run.Start, run.Way);
      const bool   sided  = (start > 0.0) == (end > 0.0);
      const double nearer = std::min(std::abs(start), std::abs(end));
      return sided && nearer > least;
    }

    /**
     *  How far a point is from a run's carrier, as Intersect's Off works it out, by the same steps.
     *
     *  @param  run   The run.
     *  @param  point The point.
     *  @return       The distance.
     */
    auto Off(const Run &run, Point point) -> double
    {
      return std::abs(Cross(point - run.Start, run.Way)) / Length(run.Way);
    }

    /**
     *  The stretch of one line along another Intersect's Aligned finds, by the same steps: the
     *  stretch of the first the second's shadow covers, longer than the tolerance, its ends within
     *  the tolerance of the second's carrier. Two lines' stretches are never endless.
     *
     *  @param  first  The first line's run.
     *  @param  second The second's.
     *  @param  length The tolerance.
     *  @return        The stretch of the first, from and to; nothing where it finds none.
     */
    auto Aligned(const Run &first, const Run &second, double length) -> std::optional<Range>
    {
      const Vector way    = first.Way;
      const Vector other  = second.Way;
      const double square = Dot(way, way);
      if (!(square > 0.0))
      {
        return std::nullopt;
      }
      const double base  = Dot(second.Start - first.Start, way) / square;
      const double rate  = Dot(other, way) / square;
      const double from  = base + (rate * second.Along.Start);
      const double until = base + (rate * second.Along.End);
      const Range  span  = {.Start = std::max(first.Along.Start, std::min(from, until)),
                            .End   = std::min(first.Along.End, std::max(from, until))};
      if (!((span.End - span.Start) * Length(way) > length))
      {
        return std::nullopt;
      }
      const bool near = Off(second, first.Start + (way * span.Start)) <= length
                     && Off(second, first.Start + (way * span.End))   <= length;
      return near ? std::optional<Range>{span} : std::nullopt;
    }

    /**
     *  Whether every end of a line lies farther than a distance from a run, as Intersect asks it of
     *  its Away.
     *
     *  @param  line  The line.
     *  @param  run   The run.
     *  @param  clear The distance.
     *  @return       Whether they do.
     */
    auto Far(const Line &line, const Run &run, double clear) -> bool
    {
      return Away(line.Start, run) > clear && Away(line.End, run) > clear;
    }

    /**
     *  Whether two points lie farther apart than a distance by Distance, Intersect's test of two
     *  ends: first by their largest coordinate's difference, which Distance is never below, as no
     *  square root of a rounded square is below the number squared, and no rounded sum below a term
     *  of it.
     *
     *  @param  first  One point.
     *  @param  second The other.
     *  @param  apart  The distance.
     *  @return        Whether Distance takes them as farther apart.
     */
    auto Asunder(Point first, Point second, double apart) -> bool
    {
      return Largest(second - first) > apart || !(Distance(first, second) <= apart);
    }

    /**
     *  Whether a point lies farther from a run than a distance by Away: first by its gap from the
     *  box of the run's start and its start moved along its whole way, which every point Away
     *  places on it lies within, as a product with a place within 0 and 1 is no larger than the
     *  way, and a rounded sum keeps its order; Away is no nearer, as Asunder tells.
     *
     *  @param  point The point.
     *  @param  run   The run, over 0 to 1.
     *  @param  clear The distance.
     *  @return       Whether Away takes it as farther.
     */
    auto Remote(Point point, const Run &run, double clear) -> bool
    {
      const Point last = run.Start + run.Way;
      const auto  gap  = [](double value, double one, double other) -> double
      {
        const double low  = std::min(one, other);
        const double high = std::max(one, other);
        if (value < low)
        {
          return low - value;
        }
        return value > high ? value - high : 0.0;
      };
      return std::max(gap(point.X, run.Start.X, last.X), gap(point.Y, run.Start.Y, last.Y)) > clear
          || Away(point, run)                                                               > clear;
    }

    /**
     *  Whether a point on a line's carrier lies on the line, as Intersect places a crossing on it:
     *  its place along the way within a few steps of a double of 0 and 1.
     *
     *  @param  line  The line.
     *  @param  point The point.
     *  @return       Whether it does.
     */
    auto Perched(const Line &line, Point point) -> bool
    {
      const Vector way    = line.End - line.Start;
      const double square = Dot(way, way);
      if (!(square > 0.0))
      {
        return false;
      }
      const double along = Dot(point - line.Start, way) / square;
      return -Overstep <= along && along <= 1.0 + Overstep;
    }

    /**
     *  Two lines as Intersect's direct way takes them: each as a curve and as a run, the tolerance,
     *  and the numbers its tests are made on, worked out as it works them out.
     */
    struct Straights
    {
      std::array<const Curve *, 2> Shapes = {};  /**< Each as a curve.               */
      std::array<const Line *, 2>  Lines  = {};  /**< Each as a line.                */
      std::array<Run, 2>           Runs   = {};  /**< Each as a run.                 */
      std::array<double, 2>        Long   = {};  /**< The length of each one's way.  */
      double                       Length = 0.0; /**< The tolerance.                 */
      double                       Cross  = 0.0; /**< The cross product of the ways. */
    };

    /**
     *  What a pass on two lines tells: a hit, that Intersect goes on to its next way, or that
     *  Intersect must decide.
     */
    enum class Outcome : std::uint8_t
    {
      Hit,    /**< One hit, at the point given. */
      Passed, /**< Intersect goes on.           */
      Unsure  /**< Intersect decides.           */
    };

    /**
     *  A pass's hit, where it tells one.
     */
    struct Answer
    {
      Outcome Kind = Outcome::Unsure; /**< What it tells. */
      Mark    Hit  = {};              /**< The hit.       */
    };

    /**
     *  A difference of two doubles, exactly, in two doubles.
     *
     *  @param  head The double.
     *  @param  tail The double taken from it.
     *  @return      The difference.
     */
    auto Sundered(double head, double tail) -> Twofold::Double
    {
      const Twofold::Two difference = Twofold::TwoSum(head, -tail);
      return {.High = difference.High, .Low = difference.Low};
    }

    /**
     *  Where two lines' carriers cross as Intersect's quick crossing works it out, each coordinate
     *  on its own: the double nearest it where its bound settles it.
     */
    struct Quick
    {
      std::array<std::optional<double>, 2> At   = {};  /**< Each coordinate, where it settles. */
      double                               Turn = 0.0; /**< The ways' cross product, roughly, its
                                                            sign certain wherever a pair is asked
                                                            here.                              */
    };

    /**
     *  Where two lines' carriers cross, as Intersect's quick crossing works it out, by the same
     *  steps in double-double: each way and the gap exact in two doubles, the cross products and
     *  their quotient each within 2^-98 of their terms' sizes, each coordinate the double nearest
     *  it where its bound leaves no halfway point near.
     *
     *  @param  first  One line.
     *  @param  second The other.
     *  @return        Each coordinate where it settles, and the ways' cross product.
     */
    auto Quickly(const Line &first, const Line &second) -> Quick
    {
      constexpr double                     Grain = 0x1p-98;
      const std::array<Twofold::Double, 2> way   = {Sundered(first.End.X, first.Start.X),
                                                    Sundered(first.End.Y, first.Start.Y)};
      const std::array<Twofold::Double, 2> other = {Sundered(second.End.X, second.Start.X),
                                                    Sundered(second.End.Y, second.Start.Y)};
      const std::array<Twofold::Double, 2> gap   = {Sundered(second.Start.X, first.Start.X),
                                                    Sundered(second.Start.Y, first.Start.Y)};
      const Twofold::Double turn  = Twofold::Difference(Twofold::Multiply(way[0], other[1]),
                                                        Twofold::Multiply(way[1], other[0]));
      const Twofold::Double reach = Twofold::Difference(Twofold::Multiply(gap[0], other[1]),
                                                        Twofold::Multiply(gap[1], other[0]));
      const double          turnBound =
          Grain * (std::abs(way[0].High * other[1].High) + std::abs(way[1].High * other[0].High));
      const double reachBound =
          Grain * (std::abs(gap[0].High * other[1].High) + std::abs(gap[1].High * other[0].High));
      /**
       *  Where Intersect's gives up at once, its turn not clear of its bound, nothing settles here:
       *  worked out all the same, never branched on, as no pair asked here crosses so shallowly.
       */
      const int                   clear = static_cast<int>(std::abs(turn.High) > 4.0 * turnBound);
      const Twofold::Double       along = Twofold::Divide(reach, turn);
      const double                size  = std::abs(along.High);
      const double                alongBound = (2.0
                                                * (reachBound + (2.0 * size * turnBound))
                                                / (std::abs(turn.High) - (2.0 * turnBound)))
                                             + (Grain * size);
      Quick                       quick      = {.At = {}, .Turn = turn.High};
      const std::array<double, 2> start      = {first.Start.X, first.Start.Y};
      for (std::size_t axis = 0; axis < 2; ++axis)
      {
        const Twofold::Double step = way.at(axis);
        const Twofold::Double value =
            Twofold::Add(Twofold::Lifted(start.at(axis)), Twofold::Multiply(step, along));
        const double bound =
            (2.0 * std::abs(step.High) * alongBound)
            + (Grain * (std::abs(start.at(axis)) + std::abs(step.High * along.High)));
        const int settled =
            clear
            * static_cast<int>(Spacing(value.High) - (2.0 * std::abs(value.Low)) > 2.0 * bound);
        if (settled > 0)
        {
          quick.At.at(axis) = value.High;
        }
      }
      return quick;
    }

    /**
     *  A coordinate of a point.
     *
     *  @param  point The point.
     *  @param  axis  Which: 0 for x, 1 for y.
     *  @return       It.
     */
    auto Along(Point point, std::size_t axis) -> double
    {
      return axis == 0 ? point.X : point.Y;
    }

    /**
     *  The end two lines share exactly, where they have one: an end of the first that is an end of
     *  the second, bit for bit.
     *
     *  @param  first  One line.
     *  @param  second The other.
     *  @return        The end; nothing where they share none.
     */
    auto Commonly(const Line &first, const Line &second) -> std::optional<Point>
    {
      const auto same = [](Point one, Point other) -> bool
      {
        return Exactly(one.X, other.X) && Exactly(one.Y, other.Y);
      };
      for (const Point end : {first.Start, first.End})
      {
        if (same(end, second.Start) || same(end, second.End))
        {
          return end;
        }
      }
      return std::nullopt;
    }

    /**
     *  Where two lines' carriers cross, as Intersect's Exactly finds it: its quick crossing where
     *  that settles both coordinates; else its exact way's, each coordinate the quotient's double
     *  nearest, which is the quick one's where that settles it, a line's own where the line holds
     *  it, an end's the lines share exactly where it is not zero, and zero where the carriers cross
     *  on the axis, as Axis tells it exactly. A zero there is the exact quotient's of a numerator
     *  of zero over the ways' cross product, of that sign.
     *
     *  @param  pair The two lines.
     *  @return      The crossing; nothing where Intersect's quick crossing gives up at once, or a
     *               coordinate does not settle and is not zero.
     */
    auto Crossed(const Straights &pair) -> std::optional<Point>
    {
      const Line &first  = *pair.Lines[0];
      const Line &second = *pair.Lines[1];
      const Quick quick  = Quickly(first, second);
      if (quick.At[0] && quick.At[1])
      {
        return Point{.X = quick.At[0].value_or(0.0), .Y = quick.At[1].value_or(0.0)};
      }
      const double               zero  = quick.Turn > 0.0 ? 0.0 : -0.0;
      std::array<double, 2>      point = {};
      std::array<int, 2>         asked = {};
      const std::optional<Point> joint = Commonly(first, second);
      for (std::size_t axis = 0; axis < 2; ++axis)
      {
        std::optional<double> held;
        for (const Line *line : pair.Lines)
        {
          if (Exactly(Along(line->Start, axis), Along(line->End, axis)))
          {
            held = Along(line->Start, axis);
          }
        }
        /**
         *  An end the two lines share exactly is where their carriers cross, as the exact way finds
         *  it to the bit; a zero there is left to the zero below, of the turn's sign.
         */
        if (joint && !Zero(Along(*joint, axis)))
        {
          held = Along(*joint, axis);
        }
        /**
         *  A coordinate that settles at zero, or does not settle, may be zero, asked below: a zero
         *  of Intersect's exact way's has the turn's sign.
         */
        const double value = held.value_or(quick.At.at(axis).value_or(0.0));
        point.at(axis)     = Zero(value) ? zero : value;
        asked.at(axis)     = static_cast<int>(!held) * static_cast<int>(Zero(value));
      }
      /**
       *  A coordinate asked, zero, is the crossing's where the point lies on both carriers exactly.
       */
      const Point spot = {.X = point[0], .Y = point[1]};
      const bool  told =
          asked[0] + asked[1] == 0 || (Axis::On(first, spot) && Axis::On(second, spot));
      return told ? std::optional<Point>{spot} : std::nullopt;
    }

    /**
     *  The carriers' crossing of two lines as Intersect's Exactly finds it, a hit of a kind where
     *  Intersect places it on both lines.
     *
     *  @param  pair The two lines.
     *  @param  kind The hit's kind.
     *  @return      The hit; Passed where Intersect places it off either line, so that its Exactly
     *               finds nothing; Unsure where its crossing is not told.
     */
    auto Exacted(const Straights &pair, Contact kind) -> Answer
    {
      const std::optional<Point> exact = Crossed(pair);
      if (!exact)
      {
        return {};
      }
      if (!Perched(*pair.Lines[0], *exact) || !Perched(*pair.Lines[1], *exact))
      {
        return {.Kind = Outcome::Passed};
      }
      return {.Kind = Outcome::Hit, .Hit = {.At = *exact, .Kind = kind}};
    }

    /**
     *  Two lines met as Intersect's Endly meets them where an end of one lies near the other: at a
     *  sine past Steep, every end farther than four tolerances from the other line or within half
     *  the tolerance of the rough crossing, a hit at an end at the exact one.
     *
     *  @param  pair The two lines.
     *  @return      The hit; Passed where Endly gives them up.
     */
    auto Endly(const Straights &pair) -> Answer
    {
      const Run   &one     = pair.Runs[0];
      const Run   &other   = pair.Runs[1];
      const double lengths = pair.Long[0] * pair.Long[1];
      if (!(lengths > 0.0))
      {
        return {.Kind = Outcome::Passed};
      }
      const double sine   = pair.Cross / lengths;
      const double length = pair.Length;
      if (!(std::abs(sine) > Met::Steep))
      {
        return {.Kind = Outcome::Passed};
      }
      const Vector gap   = other.Start - one.Start;
      const Point  point = one.Start + (one.Way * (Cross(gap, other.Way) / pair.Cross));
      const auto   ended = [&point, length](const Line &line, const Run &run) -> bool
      {
        const auto end = [&point, &run, length](Point tip) -> bool
        {
          return Away(tip, run) > 4.0 * length || Distance(tip, point) <= length / 2.0;
        };
        return end(line.Start) && end(line.End);
      };
      if (!ended(*pair.Lines[0], other) || !ended(*pair.Lines[1], one))
      {
        return {.Kind = Outcome::Passed};
      }
      return Exacted(pair, Contact::End);
    }

    /**
     *  Two lines met as Intersect's Butt meets them, where its direct ways leave them to their
     *  bands: one pair of ends, one of each, within half the tolerance of each other, every other
     *  end farther than four tolerances from the other line; a hit at an end halfway between the
     *  two.
     *
     *  @param  pair The two lines.
     *  @return      The hit; Unsure where Butt finds none, for Intersect's thorough way.
     */
    auto Butt(const Straights &pair) -> Answer
    {
      const double clear  = 4.0 * pair.Length;
      const Line  &first  = *pair.Lines[0];
      const Line  &second = *pair.Lines[1];
      if (!(Dot(pair.Runs[0].Way, pair.Runs[0].Way) > 0.0)
          || !(Dot(pair.Runs[1].Way, pair.Runs[1].Way) > 0.0))
      {
        return {};
      }
      const std::array<Point, 2> ones   = {first.Start, first.End};
      const std::array<Point, 2> others = {second.Start, second.End};
      std::array<bool, 2>        paired = {};
      std::size_t                made   = 0;
      Answer                     told;
      for (std::size_t one = 0; one < 2; ++one)
      {
        bool shared = false;
        for (std::size_t two = 0; two < 2; ++two)
        {
          if (Distance(ones.at(one), others.at(two)) <= pair.Length / 2)
          {
            if (made > 0)
            {
              return {};
            }
            const Point here  = Lerp(first.Start, first.End, static_cast<double>(one));
            const Point there = Lerp(second.Start, second.End, static_cast<double>(two));
            told              = {.Kind = Outcome::Hit,
                                 .Hit  = {.At = Midpoint(here, there), .Kind = Contact::End}};
            ++made;
            shared         = true;
            paired.at(two) = true;
          }
        }
        if (!shared && !(Away(ones.at(one), pair.Runs[1]) > clear))
        {
          return {};
        }
      }
      for (std::size_t two = 0; two < 2; ++two)
      {
        if (!paired.at(two) && !(Away(others.at(two), pair.Runs[0]) > clear))
        {
          return {};
        }
      }
      return told;
    }

    /**
     *  Two lines met as Intersect's Collinearly meets them, where Aligned finds a stretch of the
     *  first along the second, by the same steps: the stretch longer than four tolerances, each end
     *  of the first an end of the stretch or farther than four tolerances from the second, each of
     *  the second at an end of the stretch within half the tolerance of the first, or as far. There
     *  Intersect's answer is the stretch Between finds, its places on the second worked back from
     *  the first's, and no hit.
     *
     *  @param  pair The two lines.
     *  @param  span The stretch of the first along the second, as Aligned finds it.
     *  @return      What they meet in; nothing where Collinearly gives them up.
     */
    auto Lying(const Straights &pair, Range span) -> std::optional<Abutment>
    {
      const Run   &one    = pair.Runs[0];
      const Run   &other  = pair.Runs[1];
      const double length = pair.Length;
      const double clear  = 4.0 * length;
      if (!((span.End - span.Start) * pair.Long[0] > clear))
      {
        return std::nullopt;
      }
      const double square = Dot(one.Way, one.Way);
      const double base   = Dot(other.Start - one.Start, one.Way) / square;
      const double rate   = Dot(other.Way, one.Way) / square;
      const auto   bounds = [&span](double along) -> bool
      {
        return Exactly(along, span.Start) || Exactly(along, span.End);
      };
      const std::array<Point, 2> ones   = {pair.Lines[0]->Start, pair.Lines[0]->End};
      const std::array<Point, 2> others = {pair.Lines[1]->Start, pair.Lines[1]->End};
      for (std::size_t index = 0; index < 2; ++index)
      {
        const auto place = static_cast<double>(index);
        if (!bounds(place) && !(Away(ones.at(index), other) > clear))
        {
          return std::nullopt;
        }
      }
      for (std::size_t index = 0; index < 2; ++index)
      {
        const double along = base + (rate * static_cast<double>(index));
        const double away  = Away(others.at(index), one);
        const bool   bound = bounds(along) && away <= length / 2;
        if (!bound && !(away > clear))
        {
          return std::nullopt;
        }
      }
      /**
       *  Between: the stretch Aligned finds again, the same, worked back onto the second, whose
       *  way's square, which Between asks to be above zero, is: its shadow is past four tolerances.
       */
      Abutment   lying;
      const auto back = [base, rate](double along) -> double
      {
        return (along - base) / rate;
      };
      const double start = back(span.Start);
      const double end   = back(span.End);
      lying.Stretch      = {.First  = span,
                            .Second = {.Start = std::min(start, end), .End = std::max(start, end)},
                            .Along  = rate > 0.0};
      lying.Stretches    = 1;
      return lying;
    }

    /**
     *  Two lines as Squarely reads them: each one's ends, each one's way's length, and the largest
     *  magnitude of their coordinates.
     */
    struct Rails
    {
      std::array<std::array<Point, 2>, 2> Ends  = {};  /**< Each line's start and end. */
      std::array<double, 2>               Long  = {};  /**< Each one's way's length.   */
      double                              Reach = 0.0; /**< Their coordinates' largest
                                                            magnitude.                     */
    };

    /**
     *  A level line and an upright one, either first, met as Squarely meets them.
     *
     *  @param  rails  The two lines.
     *  @param  length The length they are met within, as Reached gives it.
     *  @return        What they meet in; nothing where Intersect must decide.
     */
    auto Railed(const Rails &rails, double length) -> std::optional<Abutment>
    {
      const bool                  flat    = Exactly(rails.Ends[0][0].Y, rails.Ends[0][1].Y);
      const std::array<Point, 2> &level   = rails.Ends.at(flat ? 0 : 1);
      const std::array<Point, 2> &upright = rails.Ends.at(flat ? 1 : 0);
      const double                clear   = 4.0 * length * Margin;
      const bool                  square  = Exactly(level[0].Y, level[1].Y)
                                         && Exactly(upright[0].X, upright[1].X)
                                         && rails.Long.at(flat ? 0 : 1) > clear
                                         && rails.Long.at(flat ? 1 : 0) > clear
                                         && length > 0x1p-47 * rails.Reach;
      if (!square)
      {
        return std::nullopt;
      }
      const double                x      = upright[0].X;
      const double                y      = level[0].Y;
      const std::array<double, 4> gaps   = {level[0].X - x, level[1].X - x, upright[0].Y - y,
                                            upright[1].Y - y};
      const auto                  beside = [clear](double start, double end) -> bool
      {
        const bool sided = (start > 0.0 && end > 0.0) || (start < 0.0 && end < 0.0);
        return sided && std::min(std::abs(start), std::abs(end)) > clear;
      };
      if (beside(gaps[0], gaps[1]) || beside(gaps[2], gaps[3]))
      {
        return Abutment{};
      }
      bool ended = false;
      for (const double gap : gaps)
      {
        if (Zero(gap))
        {
          ended = true;
          continue;
        }
        if (!(std::abs(gap) > clear))
        {
          return std::nullopt;
        }
      }
      if (Zero(x) || Zero(y))
      {
        return std::nullopt;
      }
      Abutment told;
      told.Marks[0] = {.At = {.X = x, .Y = y}, .Kind = ended ? Contact::End : Contact::Cross};
      told.Count    = 1;
      return told;
    }

    /**
     *  Whether two points are the same, coordinate by coordinate.
     *
     *  @param  one   One point.
     *  @param  other The other.
     *  @return       Whether they are.
     */
    auto Selfsame(Point one, Point other) -> bool
    {
      return Exactly(one.X, other.X) && Exactly(one.Y, other.Y);
    }

    /**
     *  The length Intersect made ready meets two lines within, where it asks Abutting of them: as
     *  Reached gives it from their boxes, the ends', and the largest of their coordinates, within
     *  Coarse and Fine, at an angle within half of Steep; each line plain, its box finite.
     *
     *  @param  limbs The two lines as limbs.
     *  @param  sane  The tolerance, as Sane makes it.
     *  @return       The length; nothing where Intersect asks Abutting nothing.
     */
    auto Joinable(const std::array<Limb, 2> &limbs, const Tolerance &sane) -> std::optional<double>
    {
      const double reach  = std::max(limbs[0].Reach, limbs[1].Reach);
      const double length = std::max(sane.Length, 16.0 * Above(reach));
      const bool   plain  = limbs[0].Finite && limbs[1].Finite;
      if (!plain || length > Coarse || length < Fine || !(sane.Radians <= Met::Steep / 2.0))
      {
        return std::nullopt;
      }
      return length;
    }

    /**
     *  Whether two lines joined end to end are met as Squarely meets a level line and an upright
     *  one, Met's first pass, and so in the joint alone: Squarely's mark is at the upright's x and
     *  the level's y, the joint, and it finds neither beside the other, the end at the joint of
     *  each lying on the other's carrier, its gap zero.
     *
     *  @param  pair The two lines, their lengths not yet worked out.
     *  @return      Whether Squarely tells them.
     */
    auto Squarer(Straights &pair) -> bool
    {
      const Line &first  = *pair.Lines[0];
      const Line &second = *pair.Lines[1];
      if (!Exactly(first.Start.Y, first.End.Y) && !Exactly(second.Start.Y, second.End.Y))
      {
        return false;
      }
      pair.Long = {Length(pair.Runs[0].Way), Length(pair.Runs[1].Way)};
      const double reach =
          std::max({std::abs(first.Start.X), std::abs(first.Start.Y), std::abs(first.End.X),
                    std::abs(first.End.Y), std::abs(second.Start.X), std::abs(second.Start.Y),
                    std::abs(second.End.X), std::abs(second.End.Y)});
      return Railed({.Ends  = {std::array<Point, 2>{first.Start, first.End},
                               std::array<Point, 2>{second.Start, second.End}},
                     .Long  = pair.Long,
                     .Reach = reach},
                    pair.Length)
          .has_value();
    }

    /**
     *  Two lines joined end to end met as Butt meets them where Endly passes them at once: its hit
     *  pairs the joint with itself, halfway between the lines' ends there, which Midpoint takes to
     *  the joint itself, where every other two ends lie farther apart than half the tolerance and
     *  each line's other end farther than four tolerances from the other line; else it finds none.
     *  Each test is Butt's, asked first of what bounds it from below, as Asunder and Remote ask it.
     *  Butt's first, that each way's square is above zero, holds where the ends lie farther apart
     *  than half the tolerance, no finer than Fine: the square is no less than its largest
     *  coordinate's, rounded.
     *
     *  @param  pair  The two lines.
     *  @param  ahead Whether the first's end is the joint, not its start.
     *  @return       Whether Butt finds its hit, at the joint's very point.
     */
    auto Propped(const Straights &pair, bool ahead) -> bool
    {
      const Line  &first  = *pair.Lines[0];
      const Line  &second = *pair.Lines[1];
      const Run   &one    = pair.Runs[0];
      const Run   &other  = pair.Runs[1];
      const Point  joint  = ahead ? first.End : first.Start;
      const Point  near   = ahead ? first.Start : first.End;
      const Point  far    = ahead ? second.End : second.Start;
      const double half   = pair.Length / 2;
      const double clear  = 4.0 * pair.Length;
      return Asunder(near, joint, half)
          && Asunder(near, far, half)
          && Asunder(joint, far, half)
          && Remote(near, other, clear)
          && Remote(far, one, clear);
    }
  } // namespace

  /**
   *  Intersect's direct way meets two lines by Straightly. The level line's way is (dx, 0) and the
   *  upright's (0, dy), each difference rounded and each zero exact; their lengths are |dx| and
   *  |dy| exactly, the cross product of the ways is dx dy rounded, and their dot product is zero. S
   *  is the largest magnitude of the lines' coordinates, and L, at least 2^-47 S, the length. Each
   *  end's gap from the other's carrier is g, its own coordinate less the carrier's.
   *
   *  - Beside works the cross product of each end's gap from the other's start with the other's
   *    way: g times dx or dy, rounded, of g's sign. It holds where both of a line's are of one sign
   *    and past 4 L times the other's length. Gaps past 4 L (1 + 2^-40) pass it however they round,
   *    as a few roundings cannot close a factor of 1 + 2^-40; a zero gap, or two of opposite signs,
   *    never does. Every product is a normal double: each line and gap is past 4 L, 2^-108 at the
   *    least, and below 2^100, as L is no coarser than Coarse.
   *  - Aligned finds no stretch shared: one way runs along the other at a rate of zero, so what
   *    each has of the other is a point, of no length.
   *  - An end past 4 L (1 + 2^-40) from the other's carrier, each line crossing the other's, is
   *    farther from the other line than 4 L by Away, which is at least its gap rounded less two
   *    roundings of it.
   *  - An end on the other's carrier lies at the carriers' crossing (x, y). Endly works that point
   *    out as the first line's start moved along its way by a quotient of rounded products, within
   *    11.1 u S of (x, y), u being 2^-53: under L / 2, so the end is at it; and Away finds the end
   *    within as much of the other line, so not far from it.
   *
   *  So where an end lies on the other's carrier, Straightly goes to Endly, whose sine is 1 and
   *  which takes every end as at the crossing or far; else every end is far, and Straightly crosses
   *  the lines where their carriers cross, each place along its line within 0 and 1 by 2^-46 at the
   *  least, past its roundings, the angle never shallow at a tolerance's angle within Steep / 2.
   *  Either way the hit is at the double nearest the carriers' crossing, (x, y) itself, placed on
   *  each line by its gap from the line's start over its way: a quotient within 0 and 1, which no
   *  rounding takes past them. A hit at an end is the end's kind, a crossing a crossing.
   */
  auto Abutting::Squarely(const Prepared::Pair &lines, double length) -> std::optional<Abutment>
  {
    const Prepared &one   = lines[0];
    const Prepared &other = lines[1];
    return Railed({.Ends  = {one.m_Ends, other.m_Ends},
                   .Long  = {one.m_Long, other.m_Long},
                   .Reach = std::max(one.m_Reach, other.m_Reach)},
                  length);
  }

  /**
   *  Intersect's direct way meets a straight curve and an arc by Roundly, then rounds each crossing
   *  it keeps to the double nearest the carriers' crossing and places it on both curves. Take a
   *  level line at y = c, from a0 to a1 along x, its way (dx, 0) rounded, and an arc of center (p,
   *  q) and radius r; an upright line likewise, x and y swapped. S is the largest magnitude of the
   *  line's coordinates, |p| + r and |q| + r; u is 2^-53; L, the length, is at least 2^-42 S, r is
   *  past 1024 L, and r^2 - g^2, for g = c - q, past r^2 / 32.
   *
   *  - Roundly's foot of the center on the carrier is within 11.1 u S of (p, c), so the carrier's
   *    distance from the center is within 17.1 u S of |g|, under 2^-16 r as r is past 2^-32 S: 1 -
   *    (off / r)^2 is past 1/32 less 2^-14, and r - off past r / 64 less 2^-16 r, so the carrier
   *    cuts the circle twice at a sine past Steep and deeper than four tolerances, as Roundly asks.
   *  - Its crossings are p -+ sqrt(r^2 - off^2) along the line, at c exactly across it: the
   *    radicand within 37.3 u S r of r^2 - g^2, the root within 108.3 u S of the true one, and each
   *    crossing within 128 u S of the true one, its place along the line likewise.
   *  - Each number Roundly compares with half a tolerance, four tolerances, the line's ends and the
   *    arc's ends' rays is then within the slack here, 2^-45 S, of what is worked out here from the
   *    true crossings: an end's distance from a crossing, a line's end's from the circle and an
   *    arc's end's from the line, a crossing's place along the line, and its turn from an arc's
   *    end, a cross product, over r.
   *
   *  So where every test here is told, Roundly keeps the same ends as at a crossing, finds every
   *  other end far, cuts no corner as Aloof would, and keeps the same crossings as hits. Each is
   *  then the double nearest the true crossing, as here, placed on both curves by Perched, as here,
   *  and a hit at an end where an end kept lies within the tolerance of it, on the same doubles.
   *  Two hits come in order along the first curve by their places, as Intersect sorts them.
   */
  auto Abutting::Rimmed(const Prepared::Pair &curves, double length) -> std::optional<Abutment>
  {
    const bool      first = std::holds_alternative<Line>(curves[0].get().Shape());
    const Prepared &line  = curves.at(first ? 0 : 1);
    const Prepared &round = curves.at(first ? 1 : 0);
    const Arc      &arc   = std::get<Arc>(round.Shape());
    const double    size  = std::max(
        {line.m_Reach, std::abs(arc.Center.X) + arc.Radius, std::abs(arc.Center.Y) + arc.Radius});
    const Rim    rim    = {.Ends   = line.m_Ends,
                           .Cuts   = round.m_Ends,
                           .Round  = &arc,
                           .Sweep  = round.m_Sweep,
                           .Fixed  = Exactly(line.m_Ends[0].Y, line.m_Ends[1].Y) ? 1U : 0U,
                           .Length = length,
                           .Slack  = 0x1p-45 * size};
    const double across = Coordinate(rim.Ends[0], rim.Fixed);
    const bool   plain  = Exactly(across, Coordinate(rim.Ends[1], rim.Fixed))
                       && line.m_Long > 4.0 * length
                       && round.m_Round
                       && std::abs(round.m_Sweep) < 2.0 * Pi
                       && arc.Radius > 1024.0 * length
                       && length >= 0x1p-42 * size
                       && !(Zero(across) && std::signbit(across));
    const std::optional<Spots> spots =
        plain ? Spotted(rim, {.High = round.m_Square, .Low = round.m_Under}) : std::nullopt;
    const std::optional<Held> held = spots ? Settled(rim, *spots) : std::nullopt;
    if (!held)
    {
      return std::nullopt;
    }
    /**
     *  A crossing within the line, or the line's end at one, keeps Aloof from telling them apart.
     */
    std::array<int, 2> hits   = {};
    bool               inside = held->Line;
    for (std::size_t side = 0; side < 2; ++side)
    {
      const Point spot = spots->At.at(side);
      hits.at(side)    = Struck(rim, spot);
      if (hits.at(side) == 0 || (hits.at(side) > 0 && !spots->Sure.at(side)))
      {
        return std::nullopt;
      }
      inside = inside || Housed(rim, Coordinate(spot, 1 - rim.Fixed)) > 0;
    }
    if (!inside)
    {
      return std::nullopt;
    }
    const std::array<const Curve *, 2> shapes = {&curves[0].get().Shape(),
                                                 &curves[1].get().Shape()};
    const std::span<const Point>       ends   = std::span{held->Ends}.first(held->Count);
    std::array<double, 2>              places = {};
    Abutment                           abutment;
    for (std::size_t side = 0; side < 2; ++side)
    {
      if (hits.at(side) < 0)
      {
        continue;
      }
      const Point                             spot   = spots->At.at(side);
      const std::optional<Crossings::Crossed> placed = Crossings::Placed(shapes, {}, spot);
      if (!placed)
      {
        return std::nullopt;
      }
      const bool ending         = std::ranges::any_of(ends, [spot, length](Point end) -> bool
                                                      { return Distance(end, spot) <= length; });
      places.at(abutment.Count) = placed->Along[0];
      abutment.Marks.at(abutment.Count) = {.At   = spot,
                                           .Kind = ending ? Contact::End : Contact::Cross};
      ++abutment.Count;
    }
    if (abutment.Count == 2 && places[1] < places[0])
    {
      std::swap(abutment.Marks[0], abutment.Marks[1]);
    }
    return abutment;
  }

  /**
   *  Intersect's direct way meets two lines by Straightly, every test made here on the same numbers
   *  by the same steps, in the same order:
   *  - a line wholly beside the other's carrier: nothing;
   *  - a stretch along the other found by Aligned: left to Intersect, which weighs what they share;
   *  - an end within four tolerances of the other line: Endly, its hit at the exact crossing, an
   *    end's kind; where Endly gives them up, Shallowly at an angle within twice the tolerance's,
   *    left to Intersect, else on through Intersect's direct ways, which ask nothing more of two
   *    lines, to Butt's band: a hit at two ends that meet;
   *  - every end farther: at an angle within twice the tolerance's, or a stretch the other way
   *    round, left to Intersect; else a crossing where the carriers' crossing lies within both,
   *    rough places worked out as Intersect works them, at the exact crossing.
   *
   *  The exact crossing is Intersect's quick crossing, the same function on the same lines, placed
   *  on each line as Intersect places it by its own steps; where it does not settle, Intersect's
   *  exact way decides. A Butt's hit is Intersect's, halfway between the ends as it works it out.
   */
  auto Abutting::Slanted(const Prepared::Pair &lines, const Tolerance &usable)
      -> std::optional<Abutment>
  {
    const Curve &firstCurve  = lines[0].get().Shape();
    const Curve &secondCurve = lines[1].get().Shape();
    const Line  &first       = std::get<Line>(firstCurve);
    const Line  &second      = std::get<Line>(secondCurve);
    Straights    pair        = {.Lines  = {&first, &second},
                                .Runs   = {Running(first), Running(second)},
                                .Long   = {},
                                .Length = usable.Length,
                                .Cross  = 0.0};
    const Run   &one         = pair.Runs[0];
    const Run   &other       = pair.Runs[1];
    pair.Long                = {Length(one.Way), Length(other.Way)};
    if (!(pair.Long[0] > 0.0) || !(pair.Long[1] > 0.0))
    {
      return std::nullopt;
    }
    const double clear = 4.0 * usable.Length;
    if (Beside(second, one, clear * pair.Long[0]) || Beside(first, other, clear * pair.Long[1]))
    {
      return Abutment{};
    }
    const std::optional<Range> span = Aligned(one, other, usable.Length);
    if (span)
    {
      const std::optional<Abutment> lying = Lying(pair, *span);
      if (lying)
      {
        return lying;
      }
    }
    pair.Cross = Cross(one.Way, other.Way);
    const bool shallow =
        Zero(pair.Cross)
        || !(std::abs(pair.Cross) > 2.0 * Trig::Sin(usable.Radians) * pair.Long[0] * pair.Long[1]);
    Answer told;
    if (!Far(first, other, clear) || !Far(second, one, clear))
    {
      /**
       *  Shallowly finds nothing for lines parallel exactly, or with an end within twice the
       *  tolerance of the other: Intersect goes on to Butt.
       */
      const bool still = !shallow
                      || span.has_value()
                      || Zero(pair.Cross)
                      || !Far(first, other, 2.0 * usable.Length)
                      || !Far(second, one, 2.0 * usable.Length);
      told             = Endly(pair);
      if (told.Kind == Outcome::Passed)
      {
        told = still ? Butt(pair) : Answer{};
      }
    }
    else
    {
      if (shallow || span || Aligned(other, one, usable.Length))
      {
        return std::nullopt;
      }
      const Vector gap    = other.Start - one.Start;
      const double place  = Cross(gap, other.Way) / pair.Cross;
      const double beyond = Cross(gap, one.Way) / pair.Cross;
      const bool   within = one.Along.Start   <= place
                         && place             <= one.Along.End
                         && other.Along.Start <= beyond
                         && beyond            <= other.Along.End;
      if (!within)
      {
        return Abutment{};
      }
      told = Exacted(pair, Contact::Cross);
    }
    if (told.Kind != Outcome::Hit)
    {
      return std::nullopt;
    }
    Abutment abutment;
    abutment.Marks[0] = told.Hit;
    abutment.Count    = 1;
    return abutment;
  }

  /**
   *  Take the first line from a to b and the second from c to d, joined at p, the first's end b and
   *  the second's start c, or the first's start a and the second's end d, the very doubles. Their
   *  ways are b - a and d - c rounded, as Running makes them, and rounding is odd: x - y is -(y -
   *  x) to the bit, as is a product of the negated, and a sum of two products either way round.
   *
   *  - Intersect made ready reads the lines' boxes and their ways' lengths as worked out here: a
   *    line's box is its ends', the largest magnitude of its edges the largest of its coordinates,
   *    and the length of its way the length of b - a rounded, as Lean's Apart rounds it. So the
   *    length it meets them within is the one here, and its passes are asked where they are here.
   *  - Its passes before Abutting's tell them nothing. Judged's gap from the first's start to the
   *    second's is the first's way itself (p joining b and c), or the second's way negated (p
   *    joining a and d): its cross product with the second's way is the ways' own, or zero, so the
   *    joint lies at place 1 or 0 along the first, no place Judged takes as clear of an end.
   *    Straight's place along the first is that cross product over the ways', 1 or 0 exactly, which
   *    Seated never takes as well within. Glancing asks nothing of two lines.
   *  - Squarely is asked as Met asks it, and answers where it answers: its mark is at the upright's
   *    x and the level's y, the joint itself.
   *  - Beside finds neither line beside the other: of each line, the end at the joint has a cross
   *    product with the other's way of a way with itself, or of a zero gap, zero exactly.
   *  - Aligned finds no stretch where the ways' dot product is not below zero: the second's shadow
   *    on the first runs from 1 to 1 + rate, or from -rate to 0, rate that dot product over the
   *    first's square, and is cut to the first's [0, 1] at a point. Below zero it is asked by the
   *    same steps; a stretch it finds is left to Intersect.
   *  - The end at the joint of one line or the other lies at a run's start, where Away's place
   *    along it is a zero gap's, 0, and its point the joint's: not farther than any multiple of the
   *    tolerance. So the near branch is taken, and Shallowly never: Endly, then Butt where it
   *    passes, as Slanted asks them, on the same lines, runs, lengths and cross product.
   *  - Where Endly passes at once, at its sine, Butt pairs the joint with itself, and its hit is
   *    halfway between the lines' ends there, where every other two ends lie farther apart than
   *    half the tolerance, and each line's other end farther than four tolerances from the other
   *    line; else it finds none. Each test is Butt's, asked first of what bounds it from below, as
   *    Asunder and Remote ask it.
   *
   *  Where that hit is at the joint's very point, the arrangement's cut at it is idle on both
   *  lines: Nearest places the joint at 1 on the line it ends, the foot's quotient a way's dot
   *  product with itself over itself, and at 0 on the line it starts, nothing nearer than its
   *  start.
   */
  auto Abutting::Jointed(const std::array<const Line *, 2> &lines, const Tolerance &sane) -> bool
  {
    return Jointed(std::array<Limb, 2>{Limbed(*lines[0]), Limbed(*lines[1])}, sane);
  }

  /**
   *  A limb's way is Running's, the same difference rounded, and its largest magnitudes are those
   *  the steps below take of it and of the two lines' ends together, a largest of largests being
   *  the largest of all.
   */
  auto Abutting::Jointed(const std::array<Limb, 2> &limbs, const Tolerance &sane) -> bool
  {
    const Line                 &first  = *limbs[0].Shape;
    const Line                 &second = *limbs[1].Shape;
    const bool                  ahead  = Selfsame(first.End, second.Start);
    const std::optional<double> length =
        ahead != Selfsame(first.Start, second.End) ? Joinable(limbs, sane) : std::nullopt;
    if (!length)
    {
      return false;
    }
    Straights pair = {
        .Lines  = {&first, &second},
        .Runs   = {Run{.Start = first.Start, .Way = limbs[0].Way, .Along = Domain(first)},
                   Run{.Start = second.Start, .Way = limbs[1].Way, .Along = Domain(second)}},
        .Long   = {},
        .Length = *length,
        .Cross  = 0.0};
    if (Squarer(pair))
    {
      return true;
    }
    const Run &one   = pair.Runs[0];
    const Run &other = pair.Runs[1];
    if (!(limbs[0].Largest > 0.0)
        || !(limbs[1].Largest > 0.0)
        || (!(Dot(one.Way, other.Way) >= 0.0) && Aligned(one, other, *length)))
    {
      return false;
    }
    pair.Cross = Cross(one.Way, other.Way);
    /**
     *  Endly passes at once where its sine is not past Steep, as it is where the cross product is
     *  no more than an eighth of the product of the ways' largest coordinates, normal, which their
     *  lengths and their product, rounded, are never below.
     */
    const double least = limbs[0].Largest * limbs[1].Largest;
    if (!(least > 0x1p-1000) || !(std::abs(pair.Cross) <= Met::Steep * least))
    {
      pair.Long            = {Length(one.Way), Length(other.Way)};
      const double lengths = pair.Long[0] * pair.Long[1];
      if (lengths > 0.0 && std::abs(pair.Cross / lengths) > Met::Steep)
      {
        return Endly(pair).Kind == Outcome::Hit;
      }
    }
    return Propped(pair, ahead);
  }

  auto Abutting::Crossing(const Prepared::Pair &lines) -> std::optional<Point>
  {
    const Straights pair = {
        .Lines = {&std::get<Line>(lines[0].get().Shape()), &std::get<Line>(lines[1].get().Shape())},
        .Runs  = {},
        .Long  = {},
        .Length = 0.0,
        .Cross  = 0.0};
    return Crossed(pair);
  }

  auto Abutting::Crossing(const Prepared &one, const Prepared &other) -> std::optional<Point>
  {
    return Crossing({one, other});
  }

  auto Abutting::Met(const Prepared::Pair &curves, const Tolerance &usable)
      -> std::optional<Abutment>
  {
    const Prepared &one    = curves[0];
    const Prepared &other  = curves[1];
    const auto     *line   = std::get_if<Line>(&one.Shape());
    const auto     *across = std::get_if<Line>(&other.Shape());
    if (!one.m_Plain || !other.m_Plain)
    {
      return std::nullopt;
    }
    if (line != nullptr && across != nullptr)
    {
      if (!(usable.Radians <= Met::Steep / 2.0))
      {
        return std::nullopt;
      }
      const std::optional<Abutment> square = Squarely(curves, usable.Length);
      return square ? square : Slanted(curves, usable);
    }
    /**
     *  Both plain, so each a line or an arc: one line, and the other an arc.
     */
    return line != nullptr || across != nullptr ? Rimmed(curves, usable.Length) : std::nullopt;
  }

  auto Abutting::Dotted(const Prepared::Pair &lines, const Tolerance &sane)
      -> std::optional<Abutment>
  {
    const Prepared &one    = lines[0];
    const Prepared &other  = lines[1];
    const auto     *first  = std::get_if<Line>(&one.Shape());
    const auto     *second = std::get_if<Line>(&other.Shape());
    if (first == nullptr
        || second == nullptr
        || !one.m_Plain
        || !other.m_Plain
        || !(sane.Radians <= Met::Steep / 2.0))
    {
      return std::nullopt;
    }
    const bool dot    = first->Start == first->End;
    const bool across = second->Start == second->End;
    if (dot == across)
    {
      return std::nullopt;
    }
    const Point     spot = dot ? first->Start : second->Start;
    const Prepared &line = dot ? other : one;
    const double    used = Usable({one.Shape(), other.Shape()}, sane).Length;
    if (!(line.m_Long > 64.0 * used))
    {
      return std::nullopt;
    }
    const Point  foot  = At(line.Shape(), Nearest(line.Shape(), spot).Along);
    const double apart = Distance(spot, foot);
    if (apart > used * Margin)
    {
      return Abutment{};
    }
    if (!(apart * Margin < used))
    {
      return std::nullopt;
    }
    Abutment told;
    told.Marks[0] = {.At = Midpoint(spot, foot), .Kind = Contact::End};
    told.Count    = 1;
    return told;
  }
} // namespace Geoxy
