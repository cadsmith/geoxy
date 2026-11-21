/***************************************************************************************************
 *  @file      source/prepared.cpp
 *  @brief     Curves made ready to be met many times: far apart told from boxes kept, two lines or
 *             arcs that cross plainly told with no place worked out, all else met as Intersect
 *             does.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "geoxy/prepared.hpp"

#include "abut.hpp"
#include "axis.hpp"
#include "chord.hpp"
#include "crossing.hpp"
#include "graze.hpp"
#include "lean.hpp"
#include "met/linear.hpp"
#include "met/round.hpp"
#include "met/slight.hpp"
#include "piecewise.hpp"
#include "plain.hpp"
#include "real.hpp"
#include "sane.hpp"
#include "thorough.hpp"
#include "trig.hpp"
#include "twofold.hpp"

#include "geoxy/angle.hpp"
#include "geoxy/arc.hpp"
#include "geoxy/box.hpp"
#include "geoxy/curve.hpp"
#include "geoxy/intersect.hpp"
#include "geoxy/line.hpp"
#include "geoxy/point.hpp"
#include "geoxy/ray.hpp"
#include "geoxy/tolerance.hpp"
#include "geoxy/vector.hpp"
#include "geoxy/xline.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <memory>
#include <optional>
#include <span>
#include <type_traits>
#include <utility>
#include <variant>

namespace Geoxy
{
  namespace
  {
    /**
     *  How many tolerances a crossing must lie from every end of both curves to be told here, by
     *  its point alone.
     *
     *  Two curves that cross above Steep part by an eighth of the way gone from their crossing, so
     *  an end this far along one lies eight tolerances from the other: no end is a hit, and the
     *  crossing is one, a plain one, whichever way Intersect goes about it.
     */
    constexpr double Wide = 64.0;

    /**
     *  Crossings told, two at the most, held in place: roughly by the rough passes, or exactly by
     *  the exact ones, each the double nearest its crossing and in the order Intersect gives them.
     */
    struct Sightings
    {
      std::array<Point, 2> At      = {};    /**< The points.         */
      std::size_t          Count   = 0;     /**< How many there are. */
      bool                 Further = false; /**< Whether a line's one crossing with a circle
                                                 is the further along it, of the carrier's two. */

      /**
       *  One more crossing.
       *
       *  @param point Its point.
       */
      void Add(Point point)
      {
        At.at(Count) = point;
        ++Count;
      }
    };

    /**
     *  Whether a place along a line lies well within it, 1, well off it, -1, or near an end, 0.
     *
     *  @param  along The place, 0 at the line's start and 1 at its end.
     *  @param  reach The line's length.
     *  @param  wide  How far from each end is well away.
     *  @return       Which.
     */
    auto Seated(double along, double reach, double wide) -> int
    {
      if (!(std::abs(along) * reach > wide) || !(std::abs(along - 1.0) * reach > wide))
      {
        return 0;
      }
      return along > 0.0 && along < 1.0 ? 1 : -1;
    }
  } // namespace

  /**
   *  What is told of ready curves from what they keep.
   */
  struct Standby
  {
    Met::Standing Chords; /**< The polyline with its remote arcs made chords. */
    Prepared      Kept;   /**< That polyline made ready.                      */

    /**
     *  What stands for a polyline, made ready.
     *
     *  @param chords The polyline with its remote arcs made chords.
     */
    explicit Standby(Met::Standing chords) : Chords(std::move(chords)), Kept(Chords.Shape) {}
  };

  struct Preparing
  {
    /**
     *  The length two ready curves are met within, where they may meet at all, as Reachable gives
     *  it for their curves, to the bit, from the boxes and sizes kept.
     *
     *  @param  curves    The two ready curves.
     *  @param  tolerance The tolerance asked for.
     *  @return           The length used; nothing where they lie too far apart to meet.
     */
    static auto Reached(const Prepared::Pair &curves, const Tolerance &tolerance)
        -> std::optional<double>
    {
      const Prepared &one    = curves[0];
      const Prepared &other  = curves[1];
      const auto      remote = [&one, &other](double clear) -> bool
      {
        const Box &box  = one.m_Box;
        const Box &next = other.m_Box;
        return box.Min.X - next.Max.X > clear
            || next.Min.X - box.Max.X > clear
            || box.Min.Y - next.Max.Y > clear
            || next.Min.Y - box.Max.Y > clear;
      };
      /**
       *  Two lines or arcs first, the commonest pair: their grain is their boxes', and they are far
       *  at four tolerances.
       */
      if (one.m_Plain && other.m_Plain)
      {
        const double size   = std::max(one.m_Reach, other.m_Reach);
        const double length = std::max(tolerance.Length, 16.0 * Above(size));
        return remote(4.0 * length) ? std::nullopt : std::optional<double>{length};
      }
      const double size   = std::max(one.m_Size, other.m_Size);
      const double length = std::max(tolerance.Length, 16.0 * Above(size));
      const bool   boxed  = one.m_Boxed && other.m_Boxed;
      return boxed && remote(16.0 * length) ? std::nullopt : std::optional<double>{length};
    }

    /**
     *  Whether a point of a ready arc's circle lies well within its sweep, 1, well off it, -1, or
     *  near an end of it, 0: by the turns from its ends seen from the center, no angle taken. A
     *  whole circle's seam is where the thorough way cuts it, and a point near it is not told here.
     *
     *  @param  round The ready arc.
     *  @param  point The point.
     *  @param  wide  How far from each end's ray is well away.
     *  @return       Which.
     */
    static auto Swept(const Prepared &round, Point point, double wide) -> int
    {
      const Arc &arc = std::get<Arc>(round.Shape());
      if (std::abs(round.m_Sweep) >= 2.0 * Pi)
      {
        return Distance(point, round.m_Ends[0]) > wide ? 1 : 0;
      }
      const bool   forward = round.m_Sweep > 0.0;
      const Vector start   = round.m_Ends.at(forward ? 0 : 1) - arc.Center;
      const Vector end     = round.m_Ends.at(forward ? 1 : 0) - arc.Center;
      const Vector toward  = point - arc.Center;
      const double left    = Cross(start, toward);
      const double right   = Cross(toward, end);
      const double edge    = wide * arc.Radius;
      if (!(std::abs(left) > edge) || !(std::abs(right) > edge))
      {
        return 0;
      }
      if (std::abs(round.m_Sweep) <= Pi)
      {
        return left > 0.0 && right > 0.0 ? 1 : -1;
      }
      return left < 0.0 && right < 0.0 ? -1 : 1;
    }

    /**
     *  How far round a ready arc a point of it lies, from its start the way it sweeps: a number
     *  that grows with the place along it, no angle taken.
     *
     *  @param  round The ready arc.
     *  @param  point The point.
     *  @return       From 0 at the start to 4 a whole turn on.
     */
    static auto Turned(const Prepared &round, Point point) -> double
    {
      const Arc   &arc    = std::get<Arc>(round.Shape());
      const Vector start  = round.m_Ends[0] - arc.Center;
      const Vector toward = point - arc.Center;
      const double side   = round.m_Sweep > 0.0 ? Cross(start, toward) : Cross(toward, start);
      const double level  = Dot(start, toward);
      const double share  = level / (std::abs(level) + std::abs(side));
      return side >= 0.0 ? 1.0 - share : 3.0 + share;
    }

    /**
     *  A ready line's way from its start to its end, rounded, as plain arithmetic gives it.
     *
     *  @param  line The ready line.
     *  @return      Its way.
     */
    static auto Way(const Prepared &line) -> Vector { return line.m_Ends[1] - line.m_Ends[0]; }

    /**
     *  What a ready curve keeps to be met piece by piece.
     *
     *  @param  curve The ready curve.
     *  @return       A polyline's legwork; nothing for any other curve.
     */
    static auto Worked(const Prepared &curve) -> const Legwork * { return curve.m_Legwork.get(); }

    /**
     *  What stands for a ready curve.
     *
     *  @param  curve The ready curve.
     *  @return       It; none for a curve with no remote arc.
     */
    static auto Stood(const Prepared &curve) -> const Standby * { return curve.m_Standby.get(); }

    /**
     *  Two ready curves met where one has remote arcs: as Slightly meets them, what stands for each
     *  kept ready already.
     *
     *  @param curves The two ready curves.
     *  @param sane   The tolerance, as Sane makes it.
     *  @param into   Where what they meet in goes, empty on entry.
     */
    static void Slightly(const Prepared::Pair &curves, const Tolerance &sane, Intersections &into)
    {
      const std::array<const Standby *, 2>  stood = {Stood(curves[0]), Stood(curves[1])};
      const std::array<const Prepared *, 2> kept  = {
          stood[0] != nullptr ? &stood[0]->Kept : &curves[0].get(),
          stood[1] != nullptr ? &stood[1]->Kept : &curves[1].get()};
      /**
       *  Whether they lie within reach is told from the curves' own boxes, which hold their arcs.
       */
      const std::optional<double> length = Reached(curves, sane);
      if (!length)
      {
        return;
      }
      const auto met = [](const void *context, const Curves &, const Tolerance &held,
                          Intersections &found) -> void
      {
        const auto &readies = *static_cast<const std::array<const Prepared *, 2> *>(context);
        Intersect(Prepared::Pair{*readies[0], *readies[1]}, held, found);
      };
      Met::Slightly({curves[0].get().Shape(), curves[1].get().Shape()},
                    {.Sane      = sane,
                     .Length    = *length,
                     .Intersect = met,
                     .Deep      = Deeply,
                     .Context   = &kept,
                     .Kept      = {stood[0] != nullptr ? &stood[0]->Chords : nullptr,
                                   stood[1] != nullptr ? &stood[1]->Chords : nullptr}},
                    into);
    }

    /**
     *  Whether two ready curves are of a size and of ways Shrunk leaves as they are, as their kept
     *  sizes tell it: their size not past 2^96 nor under 2^-96, and each ray or endless line one
     *  Prepared keeps a way for, whose largest coordinate lies within them.
     *
     *  @param  curves The two ready curves.
     *  @return        Whether they are; false where it is not known.
     */
    static auto Shrinkless(const Prepared::Pair &curves) -> bool
    {
      const auto paced = [](const Prepared &curve) -> bool
      {
        const Curve &shape = curve.Shape();
        return curve.m_Long > 0.0
            || !(std::holds_alternative<Ray>(shape) || std::holds_alternative<Xline>(shape));
      };
      return paced(curves[0])
          && paced(curves[1])
          && !Outscaled(std::max(curves[0].get().m_Size, curves[1].get().m_Size));
    }

    /**
     *  A ready line as the exact passes read it: what it keeps, its ends at the places 0 and 1.
     */
    struct Bar
    {
      const Prepared *Kept = nullptr; /**< The ready line. */
    };

    /**
     *  A ready ray or endless line as the exact passes read it, as a line whose end lies at an
     *  infinity of places along it, and for an endless line whose start does too: its point at 0,
     *  its way, which is exact, and where it starts.
     */
    struct Track
    {
      const Prepared *Kept  = nullptr; /**< The ready ray or endless line. */
      Point           Start = {};      /**< Its point at 0.                */
      Vector          Way   = {};      /**< Its way.                       */
      double          Begin = 0.0;     /**< Its start's place: 0, or an infinity for an endless
                                            line.                                              */
    };

    /**
     *  Whether a view of a ready straight curve is a line's.
     */
    template <typename StraightT>
    static constexpr bool Barred = std::is_same_v<StraightT, Bar>;

    /**
     *  A ready ray or endless line as the exact passes read it.
     *
     *  @pre             The curve is a ray or an endless line.
     *  @param  straight The ready ray or endless line.
     *  @return          Its track.
     */
    static auto Tracked(const Prepared &straight) -> Track
    {
      if (const auto *ray = std::get_if<Ray>(&straight.Shape()))
      {
        return {.Kept = &straight, .Start = ray->Start, .Way = ray->Way, .Begin = 0.0};
      }
      const auto &xline = std::get<Xline>(straight.Shape());
      return {.Kept  = &straight,
              .Start = xline.Through,
              .Way   = xline.Way,
              .Begin = -std::numeric_limits<double>::infinity()};
    }

    /**
     *  A straight curve's point at 0: a line's start, a ray's, an endless line's point.
     *
     *  @param  bar The ready line.
     *  @return     Its start.
     */
    static auto Origin(Bar bar) -> Point { return bar.Kept->m_Ends[0]; }

    /**
     *  A straight curve's point at 0.
     *
     *  @param  track The ready ray or endless line.
     *  @return       Its point.
     */
    static auto Origin(const Track &track) -> Point { return track.Start; }

    /**
     *  A straight curve's way, rounded: a line's as plain arithmetic gives it.
     *
     *  @param  bar The ready line.
     *  @return     Its way.
     */
    static auto Heading(Bar bar) -> Vector { return Way(*bar.Kept); }

    /**
     *  A straight curve's way, rounded: a ray's or an endless line's as given.
     *
     *  @param  track The ready ray or endless line.
     *  @return       Its way.
     */
    static auto Heading(const Track &track) -> Vector { return track.Way; }

    /**
     *  A ready line's way from its start to its end, exactly: the rounded way and what rounding
     *  left of it, kept.
     *
     *  @param  bar The ready line.
     *  @return     Its way.
     */
    static auto Course(Bar bar) -> Lean::Way
    {
      const Vector way = Heading(bar);
      return {Lean::Double{.High = way.X, .Low = bar.Kept->m_Slip.X},
              Lean::Double{.High = way.Y, .Low = bar.Kept->m_Slip.Y}};
    }

    /**
     *  A ready ray's or endless line's way, exactly: as given, nothing left by rounding.
     *
     *  @param  track The ready ray or endless line.
     *  @return       Its way.
     */
    static auto Course(const Track &track) -> Lean::Way
    {
      return {Lean::Double{.High = track.Way.X, .Low = 0.0},
              Lean::Double{.High = track.Way.Y, .Low = 0.0}};
    }

    /**
     *  Where a line starts, in places along it times a unit: 0.
     *
     *  @return 0.
     */
    static auto First(Bar /*bar*/, double /*unit*/) -> double { return 0.0; }

    /**
     *  Where a ray or an endless line starts, in places along it times a unit above nothing: 0 for
     *  a ray, an infinity for an endless line.
     *
     *  @param  track The ready ray or endless line.
     *  @return       Where.
     */
    static auto First(const Track &track, double /*unit*/) -> double { return track.Begin; }

    /**
     *  Where a line ends, in places along it times a unit: the unit.
     *
     *  @param  unit The unit.
     *  @return      Where.
     */
    static auto Last(Bar /*bar*/, double unit) -> double { return unit; }

    /**
     *  Where a ray or an endless line ends, in places along it times a unit above nothing: an
     *  infinity.
     *
     *  @return An infinity.
     */
    static auto Last(const Track & /*track*/, double /*unit*/) -> double
    {
      return std::numeric_limits<double>::infinity();
    }

    /**
     *  A line's ends.
     *
     *  @param  bar The ready line.
     *  @return     Its start and its end.
     */
    static auto Tips(Bar bar) -> std::array<Point, 2> { return bar.Kept->m_Ends; }

    /**
     *  A ray's or an endless line's ends: its start, or for an endless line, and its end, each a
     *  point whose coordinates change along it an infinity of its way's sign, either way.
     *
     *  @param  track The ready ray or endless line.
     *  @return       Its start and its end.
     */
    static auto Tips(const Track &track) -> std::array<Point, 2>
    {
      const auto beyond = [&track](double sign) -> Point
      {
        const Vector way     = track.Way * sign;
        const double endless = std::numeric_limits<double>::infinity();
        return {.X = Zero(way.X) ? track.Start.X : std::copysign(endless, way.X),
                .Y = Zero(way.Y) ? track.Start.Y : std::copysign(endless, way.Y)};
      };
      return {Zero(track.Begin) ? track.Start : beyond(-1.0), beyond(1.0)};
    }

    /**
     *  Where two ready lines' carriers cross, as Abutting's Crossing gives Intersect's own: where
     *  the exact pass leaves a coordinate.
     *
     *  @param  one   The first ready line.
     *  @param  other The other.
     *  @return       The point; nothing where a coordinate is neither Intersect's quick one nor its
     *                exact way's.
     */
    static auto Abutted(Bar one, Bar other) -> std::optional<Point>
    {
      return Abutting::Crossing(*one.Kept, *other.Kept);
    }

    /**
     *  Where a ray's or an endless line's carrier crosses another straight curve's, where the exact
     *  pass leaves a coordinate: nothing, as Abutting tells two lines alone, and Intersect decides.
     *
     *  @return Nothing.
     */
    template <typename OneT, typename OtherT>
    static auto Abutted(OneT /*one*/, OtherT /*other*/) -> std::optional<Point>
    {
      return std::nullopt;
    }

    /**
     *  What a test of a line's length, or of its shadow on another curve, tells of a line: what the
     *  test tells.
     *
     *  @param  passed Whether it passes: 1, or 0.
     *  @return        The same.
     */
    static auto Lasting(Bar /*bar*/, int passed) -> int { return passed; }

    /**
     *  What a test of a line's length, or of its shadow on another curve, tells of a ray or an
     *  endless line: that it passes, as it runs on without end, and its shadow with it.
     *
     *  @return 1.
     */
    static auto Lasting(const Track & /*track*/, int /*passed*/) -> int { return 1; }

    /**
     *  Whether a point lies on a ready line's carrier exactly, as Axis tells it.
     *
     *  @param  bar   The ready line.
     *  @param  point The point.
     *  @return       Whether it does.
     */
    static auto Carries(Bar bar, Point point) -> bool
    {
      return Axis::On(std::get<Line>(bar.Kept->Shape()), point);
    }

    /**
     *  Whether a crossing of a ray's or an endless line's carrier at zero is made certain here: no,
     *  for Intersect to decide.
     *
     *  @return False.
     */
    static auto Carries(const Track & /*track*/, Point /*point*/) -> bool { return false; }

    /**
     *  Where two ready lines' carriers cross, the double nearest the true point, for two lines
     *  Surely has told cross: their sine above Steep, so their cross product is at least an eighth
     *  of its terms' size.
     *
     *  A coordinate a line holds fixed, as a level or an upright one does, is that line's. Any
     *  other is the first line's start moved along its exact way by t = n / d, n and d the cross
     *  products of the exact gap between the starts and of the first way with the second, each in
     *  lean double-double within 2^-100 of its terms' size, Pn and Pd. Their quotient is within
     *  (Grain Pn + |t| Grain Pd) / |d| of the true one, four times what the cross products leave,
     *  with |d| at least Pd / 8 and so the reciprocal's rounding within a factor 1 + 2^-48, and
     *  within Grain |t| for the quotient's own rounding, eight times what it leaves. Moved along a
     *  step s, that is at most (1 + u) |s| times the place's bound, and the sum's own error at most
     *  a fifth of Grain times its terms' size: under 0.6 of the bound Certain is given. A
     *  coordinate not certain here, or held at zero, is Intersect's own as Abutting's Crossing
     *  gives it for two lines: its quick crossing's, or its exact way's where a coordinate is zero;
     *  for a ray or an endless line it is left to Intersect.
     *
     *  @param  one   The first ready line, ray or endless line, crossing the other as Surely tells
     *                it.
     *  @param  other The other.
     *  @return       The crossing; nothing where a coordinate is certain neither here nor there, as
     *                halfway between two doubles, nor, for a ray or an endless line, here.
     */
    template <typename OneT, typename OtherT>
    static auto Crossed(OneT one, OtherT other) -> std::optional<Point>
    {
      const Point     start  = Origin(one);
      const Point     begin  = Origin(other);
      const Lean::Way way    = Course(one);
      const Lean::Way across = Course(other);
      const auto      held   = [&way, &across](std::size_t axis) -> bool
      {
        return Zero(way.at(axis).High) || Zero(across.at(axis).High);
      };
      const auto fixed = [&way, start, begin](std::size_t axis) -> double
      {
        const bool first = Zero(way.at(axis).High);
        if (axis == 0)
        {
          return first ? start.X : begin.X;
        }
        return first ? start.Y : begin.Y;
      };
      /**
       *  A coordinate held at zero is the zero Intersect's quick crossing gives where it settles
       *  both, and else its exact way's, of the sign of the ways' cross product, as Abutting tells
       *  it.
       */
      if ((held(0) && Zero(fixed(0))) || (held(1) && Zero(fixed(1))))
      {
        return Abutted(one, other);
      }
      if (held(0) && held(1))
      {
        return Point{.X = fixed(0), .Y = fixed(1)};
      }
      const Lean::Way    gap     = Lean::Apart(begin, start);
      const Lean::Sized  turn    = Lean::Cross(way, across);
      const Lean::Sized  reached = Lean::Cross(gap, across);
      const double       inverse = 1.0 / turn.Value.High;
      const Lean::Double along   = Lean::Over(reached.Value, turn.Value, inverse);
      const double       size    = std::abs(along.High);
      const double       bound =
          (Lean::Grain * (reached.Size + (size * turn.Size)) * std::abs(inverse) * Lean::Loose)
          + (Lean::Grain * size);
      std::array<double, 2> point = {};
      for (std::size_t axis = 0; axis < 2; ++axis)
      {
        if (held(axis))
        {
          point.at(axis) = fixed(axis);
          continue;
        }
        const Lean::Double  step  = way.at(axis);
        const Lean::Rounded moved = Lean::Moved(axis == 0 ? start.X : start.Y, step, along);
        if (!Lean::Certain(moved.Value, moved.Rest,
                           (2.0 * std::abs(step.High) * bound) + (Lean::Grain * moved.Size)))
        {
          return Abutted(one, other);
        }
        point.at(axis) = moved.Value;
      }
      return Point{.X = point[0], .Y = point[1]};
    }

    /**
     *  Whether two ready lines cross plainly, or plainly do not, as Straight would tell them, with
     *  no division and no square root, every test worked out and none branched on, so that nothing
     *  waits on a guess: each of Straight's tests on a quotient made a test on a product, passed
     *  only where the quotient Straight rounds passes too.
     *
     *  Each test here is on products of the very doubles Straight rounds its quotients from, clear
     *  of the other side by a factor 1 + 2^-48, 32 times what the four roundings either side can
     *  take:
     *  - the sine above Steep: the same product, as Straight works it;
     *  - a shadow lap / reach under half a tolerance: lap under (length / 2) reach (1 - 2^-48), so
     *    the quotient is under (length / 2)(1 - 30 u)(1 + u); over 32 tolerances likewise;
     *  - a place n / d along a line at least wide from its start: |n| reach over wide |d| (1 +
     *    2^-48), so the rounded place is over wide / reach (1 + 8 u)(1 - u)^2;
     *  - and from its end: |n - d| reach over that and 2^-50 |n| reach more, for the rounding of n
     *    / d before 1 is taken from it, so again over wide / reach (1 + 8 u)(1 - u)^2.
     *
     *  A place so clear of 0 and 1, by Wide tolerances where a tolerance is at least sixteen steps
     *  of a double at the curves' size, as Reached makes it, rounds below 1 exactly where n is
     *  below d. Lines no longer than two Wide tolerances are left to Straight, so that every
     *  product a pair told here is tested on is a normal double, far from underflow.
     *
     *  A ray is met by Intersect's Straightly as a line whose end lies at an infinity of places, an
     *  endless line as one whose start does too, so each test at such an end is n against an
     *  infinity, and passes; and its way, of any length its maker chose, is within 2^-96 and 2^96
     *  as Prepared keeps it, so no product nears underflow whatever its length. A ray's or an
     *  endless line's shadow on the other curve is endless, and Straightly's Aligned finds no
     *  stretch of the other within the tolerance of it: the stretch it weighs runs from the
     *  projection of the ray's start, or from every place, to an end of the other, or on without
     *  end, where Aligned asks for parallel carriers; that end of the other, and the ray's start,
     *  lie past Wide tolerances from the crossing at a sine above Steep, so past eight from the
     *  other's carrier. Only a line's shadow, on whatever it meets, is tested.
     *
     *  @param  one    The first ready line, ray or endless line.
     *  @param  other  The other.
     *  @param  length The tolerance used.
     *  @return        1 where they cross, 0 where they plainly do not, -1 where Straight must
     *                 decide.
     */
    template <typename OneT, typename OtherT>
    static auto Judged(OneT one, OtherT other, double length) -> int
    {
      const Vector way    = Heading(one);
      const Vector across = Heading(other);
      const double reach  = one.Kept->m_Long;
      const double span   = other.Kept->m_Long;
      const double wide   = Wide * length;
      const double cross  = Cross(way, across);
      const int    steep  = Lasting(one, static_cast<int>(reach > 2.0 * wide))
                          * Lasting(other, static_cast<int>(span > 2.0 * wide))
                          * static_cast<int>(std::abs(cross) > Met::Steep * reach * span);
      const double lap    = std::abs(Dot(way, across));
      const double under  = (length / 2.0) * Lean::Tight;
      const double over   = 32.0 * length * Lean::Loose;
      const auto   plain  = [lap, under, over](double size) -> int
      {
        return std::max(static_cast<int>(lap < under * size), static_cast<int>(lap > over * size));
      };
      const double turn   = std::abs(cross);
      const auto   placed = [turn, wide](double place, double size, const auto &track) -> int
      {
        const double edge = wide * turn * Lean::Loose;
        const double far  = std::abs(place) * size;
        return static_cast<int>(std::abs(place - First(track, turn)) * size > edge)
             * static_cast<int>(std::abs(place - Last(track, turn)) * size
                                > edge + (far * 0x1p-50));
      };
      const double sign  = std::copysign(1.0, cross);
      const Vector gap   = Origin(other) - Origin(one);
      const double ahead = sign * Cross(gap, across);
      const double aside = sign * Cross(gap, way);
      const int    sure  = steep
                         * Lasting(other, plain(reach))
                         * Lasting(one, plain(span))
                         * placed(ahead, reach, one)
                         * placed(aside, span, other);
      const int    hit   = static_cast<int>(ahead > First(one, turn))
                         * static_cast<int>(ahead < Last(one, turn))
                         * static_cast<int>(aside > First(other, turn))
                         * static_cast<int>(aside < Last(other, turn));
      return (sure * (hit + 1)) - 1;
    }

    /**
     *  Two lines, rays or endless lines that cross plainly, or plainly do not, told as Straight
     *  tells them, as Judged tells it: the crossing's point the double nearest the true one, from
     *  Crossed. Anything else is left to Straight.
     *
     *  @param  one    The first ready line, ray or endless line.
     *  @param  other  The other.
     *  @param  length The tolerance used.
     *  @return        The crossing, exact, if any; nothing where Straight must decide.
     */
    template <typename OneT, typename OtherT>
    static auto Surely(OneT one, OtherT other, double length) -> std::optional<Sightings>
    {
      const int judged = Judged(one, other, length);
      if (judged < 0)
      {
        return std::nullopt;
      }
      Sightings told;
      if (judged > 0)
      {
        const std::optional<Point> exact = Crossed(one, other);
        if (!exact)
        {
          return std::nullopt;
        }
        told.Add(*exact);
      }
      return told;
    }

    /**
     *  Two lines that cross plainly, or plainly do not: at a sine above Steep, their carriers'
     *  crossing far from every end, a hit where it lies within both. Each line's length is the one
     *  kept, worked out as it would be here.
     *
     *  @param  lines  The two ready lines.
     *  @param  length The tolerance used.
     *  @return        The crossing, roughly, if any; nothing where Intersect must decide.
     */
    static auto Straight(const Prepared::Pair &lines, double length) -> std::optional<Sightings>
    {
      const Vector way   = Way(lines[0]);
      const Vector other = Way(lines[1]);
      const double reach = lines[0].get().m_Long;
      const double span  = lines[1].get().m_Long;
      const double cross = Cross(way, other);
      if (!(std::abs(cross) > Met::Steep * reach * span))
      {
        return std::nullopt;
      }
      /**
       *  Two lines nearly square to each other, the shadow of one on the other about a tolerance
       *  long, are left to Intersect: its thorough way weighs whether they run along each other
       *  there.
       */
      const double lap   = std::abs(Dot(way, other));
      const auto   plain = [length](double shadow) -> bool
      {
        return shadow < length / 2.0 || shadow > 32.0 * length;
      };
      if (!plain(lap / reach) || !plain(lap / span))
      {
        return std::nullopt;
      }
      const Point  start  = lines[0].get().m_Ends[0];
      const Vector gap    = lines[1].get().m_Ends[0] - start;
      const double along  = Cross(gap, other) / cross;
      const int    first  = Seated(along, reach, Wide * length);
      const int    second = Seated(Cross(gap, way) / cross, span, Wide * length);
      if (first == 0 || second == 0)
      {
        return std::nullopt;
      }
      Sightings told;
      if (first > 0 && second > 0)
      {
        told.Add(start + (way * along));
      }
      return told;
    }

    /**
     *  A ready curve's square, kept: a line's way's length squared in lean double-double, within
     *  2^-100 of it, or an arc's radius squared, exactly.
     *
     *  @param  curve The ready line or arc.
     *  @return       Its square.
     */
    static auto Squared(const Prepared &curve) -> Lean::Double
    {
      return {.High = curve.m_Square, .Low = curve.m_Under};
    }

    /**
     *  How a ready arc's crossings are put in order when it is first, as Intersect gives them:
     *  round it the way it sweeps, as Turned tells it.
     *
     *  @param  first The first ready curve.
     *  @param  told  Crossings, in order along the other where the first is an arc.
     *  @return       The same crossings, in order along the first.
     */
    static auto Ordered(const Prepared &first, Sightings told) -> Sightings
    {
      if (!first.m_Round)
      {
        return told;
      }
      const bool  swap  = told.Count == 2 && Turned(first, told.At[0]) > Turned(first, told.At[1]);
      const Point lead  = swap ? told.At[1] : told.At[0];
      const Point trail = swap ? told.At[0] : told.At[1];
      told.At           = {lead, trail};
      return told;
    }

    /**
     *  Whether a point of a ready arc's circle lies well within its sweep, 1, well off it, -1, or
     *  near an end of it, 0, as Swept tells it, every test on the point worked out and none
     *  branched on, so that nothing waits on a guess. A whole circle's seam is near where the point
     *  is within wide of it, its square to the square.
     *
     *  @param  round The ready arc.
     *  @param  point The point.
     *  @param  wide  How far from each end's ray is well away.
     *  @return       Which.
     */
    static auto Sweeps(const Prepared &round, Point point, double wide) -> int
    {
      if (std::abs(round.m_Sweep) >= 2.0 * Pi)
      {
        const Vector seam = point - round.m_Ends[0];
        return static_cast<int>(Dot(seam, seam) > wide * wide);
      }
      const Arc   &arc    = std::get<Arc>(round.Shape());
      const auto   back   = static_cast<std::size_t>(!(round.m_Sweep > 0.0));
      const Vector start  = round.m_Ends.at(back) - arc.Center;
      const Vector end    = round.m_Ends.at(1 - back) - arc.Center;
      const Vector toward = point - arc.Center;
      const double left   = Cross(start, toward);
      const double right  = Cross(toward, end);
      const double edge   = wide * arc.Radius;
      const int    clear =
          static_cast<int>(std::abs(left) > edge) * static_cast<int>(std::abs(right) > edge);
      const int narrow  = static_cast<int>(std::abs(round.m_Sweep) <= Pi);
      const int both    = static_cast<int>(left > 0.0) * static_cast<int>(right > 0.0);
      const int neither = static_cast<int>(left < 0.0) * static_cast<int>(right < 0.0);
      const int inside  = (narrow * both) + ((1 - narrow) * (1 - neither));
      return clear * ((2 * inside) - 1);
    }

    /**
     *  What the two crossings of a chord share, as Chord works them out.
     */
    using Halves = Chord::Halves;

    /**
     *  A crossing worked out, and whether its point is certain, as Chord certifies it.
     */
    using Settle = Chord::Settle;

    /**
     *  What two tests tell of each of a chord's ends, as Seated and Swept tell it: well within, 1,
     *  well off, -1, or near an edge, 0.
     */
    using Tells = std::array<std::array<int, 2>, 2>;

    /**
     *  The hits among a chord's two ends, in the order the ends are given: each end a hit where
     *  both tests tell it well within; nothing where any test tells an end near an edge, or a hit
     *  is not certain. No branch waits on a test: the hits are written in place as counted.
     *
     *  @param  tells What the tests tell of each end.
     *  @param  ends  The ends, exactly where certain.
     *  @return       The hits; nothing where Intersect must decide.
     */
    static auto Kept(const Tells &tells, const std::array<Settle, 2> &ends)
        -> std::optional<Sightings>
    {
      Sightings told;
      int       lost = 0;
      for (std::size_t side = 0; side < 2; ++side)
      {
        const std::array<int, 2> &tell = tells.at(side);
        const Settle             &end  = ends.at(side);
        const int hit           = static_cast<int>(tell[0] > 0) * static_cast<int>(tell[1] > 0);
        lost                   += static_cast<int>(tell[0] == 0)
                                + static_cast<int>(tell[1] == 0)
                                + (hit * static_cast<int>(!end.Sure));
        told.At.at(told.Count)  = end.At;
        told.Count             += static_cast<std::size_t>(hit);
      }
      if (lost > 0)
      {
        return std::nullopt;
      }
      return told;
    }

    /**
     *  Whether a chord's two ends are both told well off, by one test or both, and neither near an
     *  edge: no hit, and nothing more to work out.
     *
     *  @param  tells What the tests tell of each end.
     *  @return       Whether they miss.
     */
    static auto Missed(const Tells &tells) -> bool
    {
      int hits  = 0;
      int edges = 0;
      for (const std::array<int, 2> &tell : tells)
      {
        hits  += static_cast<int>(tell[0] > 0) * static_cast<int>(tell[1] > 0);
        edges += static_cast<int>(tell[0] == 0) + static_cast<int>(tell[1] == 0);
      }
      return hits == 0 && edges == 0;
    }

    /**
     *  Whether any end of a chord is told near an edge, by either test: Intersect must decide, as
     *  Kept tells it, whatever the ends' exact points.
     *
     *  @param  tells What the tests tell of each end.
     *  @return       Whether one is.
     */
    static auto Edged(const Tells &tells) -> bool
    {
      return std::ranges::any_of(tells, [](const std::array<int, 2> &tell) -> bool
                                 { return tell[0] == 0 || tell[1] == 0; });
    }

    /**
     *  A line's carrier and a circle as a quadratic in the place t along the line, P + t d against
     *  |X - O| = r: a = d.d, b = d.(P - O), c = |P - O|^2 - r^2, and the discriminant b^2 - a c,
     *  each in lean double-double with its bound.
     */
    struct Quadratic
    {
      Lean::Double A     = {};  /**< a.                               */
      Lean::Double B     = {};  /**< b.                               */
      Lean::Double Disc  = {};  /**< b^2 - a c.                       */
      double       OffA  = 0.0; /**< How far a may be off.            */
      double       OffB  = 0.0; /**< How far b may be off.            */
      double       OffD  = 0.0; /**< How far the discriminant may be. */
      double       Width = 0.0; /**< Wide tolerances.                 */
    };

    /**
     *  Where a ready line crossing a ready arc's circle at a sine above Steep does so, and which
     *  crossings are hits: the chord form, M = P + d (-b / a) and H = d (sqrt(b^2 - a c) / a),
     *  every crossing M + H or M - H.
     *
     *  -b / a and the root over a are each a quotient by the reciprocal kept with the line, within
     *  (bound of the top + |q| bound of a) / a, four times what their parts leave, and Grain |q|,
     *  eight times what the quotient leaves; the root's own bound is the discriminant's over the
     *  root, and Grain times the root, for what it leaves. A place is told as Lined tells it: well
     *  within the line, well off it, or near an end, by its bound and its length kept, each test
     *  clear by a factor 1 + 2^-48; then, where not near an end, its point within the arc's sweep,
     *  as Swept tells it.
     *
     *  A ray's or an endless line's place is told as a line's whose end, or both ends, lie at an
     *  infinity of places: each test at such an end passes, as Intersect's Roundly, which meets
     *  them as runs over such places, holds no end there. A crossing at zero made certain by Axis
     *  is told of a line alone; a ray's or an endless line's is left to Intersect.
     *
     *  @param  line  The ready line, ray or endless line.
     *  @param  round The ready arc.
     *  @param  shape The quadratic.
     *  @return       The hits, exact, in order along the line; nothing where one is not certain or
     *                a crossing is near an end of either.
     */
    template <typename StraightT>
    static auto Cut(StraightT line, const Prepared &round, const Quadratic &shape)
        -> std::optional<Sightings>
    {
      const Lean::Way    way     = Course(line);
      const Lean::Rooted root    = Lean::Root(shape.Disc);
      const double       inverse = line.Kept->m_Inverse;
      const Lean::Double middle  = Lean::Over(Lean::Minus(shape.B), shape.A, inverse);
      const Lean::Double half    = Lean::Over(root.Value, shape.A, inverse);
      const double       offRoot =
          (shape.OffD * root.Inverse * Lean::Loose) + (Lean::Grain * root.Value.High);
      const double offMiddle =
          ((shape.OffB + (std::abs(middle.High) * shape.OffA)) * inverse * Lean::Loose)
          + (Lean::Grain * std::abs(middle.High));
      const double offHalf = ((offRoot + (half.High * shape.OffA)) * inverse * Lean::Loose)
                           + (Lean::Grain * half.High);
      const double reach   = line.Kept->m_Long * Lean::Tight;
      const double wide    = shape.Width;
      const double off =
          ((offMiddle + offHalf) * Lean::Loose) + (0x1p-50 * (std::abs(middle.High) + half.High));
      const double first = First(line, 1.0);
      const double last  = Last(line, 1.0);
      const auto   lain  = [reach, wide, off, first, last](double place) -> int
      {
        const int start = static_cast<int>((std::abs(place - first) - off) * reach > wide);
        const int end =
            static_cast<int>(((std::abs(place - last) * Lean::Tight) - off) * reach > wide);
        const int inside = static_cast<int>(place > first) * static_cast<int>(place < last);
        return start * end * ((2 * inside) - 1);
      };
      const Point  start = Origin(line);
      const Point  base  = {.X = start.X + (way[0].High * middle.High),
                            .Y = start.Y + (way[1].High * middle.High)};
      const Vector rise  = {.X = way[0].High * half.High, .Y = way[1].High * half.High};
      const auto   tell  = [&](double sign) -> std::array<int, 2>
      {
        const double place = (middle.High + (sign * half.High)) + (middle.Low + (sign * half.Low));
        return {lain(place), Sweeps(round, base + (rise * sign), wide)};
      };
      const Tells tells = {tell(-1.0), tell(1.0)};
      if (Missed(tells))
      {
        return Sightings{};
      }
      if (Edged(tells))
      {
        return std::nullopt;
      }
      Halves halves;
      for (std::size_t axis = 0; axis < 2; ++axis)
      {
        const Lean::Double step = way.at(axis);
        const double       from = axis == 0 ? start.X : start.Y;
        halves.Middle.at(axis)  = Lean::Moved(from, step, middle);
        halves.Half.at(axis)    = Lean::Times(step, half);
        halves.Off.at(axis)     = 2.0 * std::abs(step.High) * (offMiddle + offHalf);
      }
      const Arc &arc  = std::get<Arc>(round.Shape());
      const auto zero = [&line, &arc, &shape](Point spot, double halfway) -> bool
      {
        return static_cast<int>(Axis::Fits(shape.A.High * halfway))
                 * static_cast<int>(Carries(line, spot))
                 * static_cast<int>(Axis::On(arc, spot))
             > 0;
      };
      const std::array<Settle, 2> ends = Chord::Zeroed(Chord::Ends(halves), halves, zero);
      return Kept(tells, {ends[1], ends[0]});
    }

    /**
     *  A level or upright line's crossings with a circle, one whose free coordinate's bound reaches
     *  zero made certain where it is zero: the point of zero and the line's coordinate on the
     *  circle exactly, as Axis tells it, which Intersect rounds to +0, the center's coordinate and
     *  the root's share cancelling exactly, as Chord's Zeroed tells.
     *
     *  @param  line The ready line.
     *  @param  arc  The arc whose circle.
     *  @param  ends Its crossings, each asked where its Naught says so.
     *  @return      The crossings.
     */
    static auto Naught(const Prepared &line, const Arc &arc, std::array<Settle, 2> ends)
        -> std::array<Settle, 2>
    {
      for (Settle &end : ends)
      {
        for (std::size_t axis = 0; axis < 2; ++axis)
        {
          std::array<double, 2> spot   = {end.At.X, end.At.Y};
          const double          centre = axis == 0 ? arc.Center.X : arc.Center.Y;
          spot.at(axis)                = 0.0;
          const Point naught           = {.X = spot[0], .Y = spot[1]};
          if (end.Naught.at(axis)
              && static_cast<int>(Axis::Fits(line.m_Square * centre))
                         * static_cast<int>(Axis::On(arc, naught))
                     > 0)
          {
            end.At   = naught;
            end.Sure = true;
          }
        }
      }
      return ends;
    }

    /**
     *  A ready level or upright line and a ready arc, told as Chorded tells them, the coordinate
     *  the line holds fixed exact: its crossings lie at O +- sqrt(r^2 - g^2) along the other axis,
     *  g = c - O the fixed coordinate c's exact gap from the center's, one exact root.
     *
     *  r^2 - g^2 is worked out in lean double-double, within 2^-102 (r^2 + g^2), and bounded by
     *  Grain (r^2 + g^2); its root within the radicand's bound over the root and Grain times the
     *  root; each crossing one exact sum, rounded twice and twice again, within 10.1 u^2 (|O| +
     *  root), half what Grain gives it. Lined goes on where r^2 - g^2 is above r^2 / 64; such a
     *  line outside the circle is never near enough its box to be asked. A crossing lies well
     *  within the line where it is more than wide from each end along the axis, by its bound and
     *  the half step of a double its rounding may have taken.
     *
     *  A ray or an endless line is told as a line whose end, or both ends, lie at an infinity along
     *  the axis: each test at such an end passes. A crossing at zero made certain by Axis is told
     *  of a line alone; a ray's or an endless line's is left to Intersect. One whose carrier passes
     *  more than eight tolerances outside the circle, g^2 - r^2 above 16 L r + 64 L^2 by its bound
     *  and a factor 1 + 2^-48, meets the arc nowhere, as Chorded tells a slanted one: Intersect's
     *  Roundly finds the carrier past the circle by more than four tolerances, and the ray's start
     *  and the arc's ends farther than that from the other curve, so nothing. A line so far out is
     *  not asked: its box lies more than four tolerances from the circle's, and Reached tells it
     *  far.
     *
     *  @param  line   The ready line, ray or endless line, level or upright.
     *  @param  round  The ready arc, of a radius, sweeping at most once round.
     *  @param  length The tolerance used.
     *  @return        The hits, exact, in order along the line; nothing where Lined must decide.
     */
    template <typename StraightT>
    static auto Levelled(StraightT line, const Prepared &round, double length)
        -> std::optional<Sightings>
    {
      const Arc   &arc  = std::get<Arc>(round.Shape());
      const double wide = Wide * length;
      if (!(line.Kept->m_Square > 0.0) || !(arc.Radius > Wide * wide))
      {
        return std::nullopt;
      }
      const std::array<Point, 2>  tips   = Tips(line);
      const std::array<double, 2> from   = {tips[0].X, tips[0].Y};
      const std::array<double, 2> until  = {tips[1].X, tips[1].Y};
      const std::array<double, 2> centre = {arc.Center.X, arc.Center.Y};
      const std::size_t           fixed  = Zero(until[1] - from[1]) ? 1 : 0;
      const std::size_t           free   = 1 - fixed;
      /**
       *  A line along an axis at -0 crosses where Intersect's quick crossing adds to its -0 a zero
       *  of either sign: left to Intersect. At +0 every way of Intersect's gives +0, as here.
       */
      if (Zero(from.at(fixed)) && std::signbit(from.at(fixed)))
      {
        return std::nullopt;
      }
      const Twofold::Two apart   = Twofold::TwoSum(from.at(fixed), -centre.at(fixed));
      const Lean::Double gap     = {.High = apart.High, .Low = apart.Low};
      const Lean::Double squared = Squared(round);
      const Lean::Double rest    = Lean::Plus(squared, Lean::Minus(Lean::Times(gap, gap)));
      const double       offRest = Lean::Grain * (squared.High + (gap.High * gap.High));
      if (!((rest.High * Lean::Tight) - offRest
            > squared.High * (Met::Steep * Met::Steep) * Lean::Loose))
      {
        if constexpr (!Barred<StraightT>)
        {
          const double outside =
              ((16.0 * length * arc.Radius) + (64.0 * length * length)) * Lean::Loose;
          if ((-rest.High * Lean::Tight) - offRest > outside)
          {
            return Sightings{};
          }
        }
        return std::nullopt;
      }
      const Lean::Rooted root = Lean::Root(rest);
      const double       offRoot =
          (offRest * root.Inverse * Lean::Loose) + (Lean::Grain * root.Value.High);
      const double low  = std::min(from.at(free), until.at(free));
      const double high = std::max(from.at(free), until.at(free));
      const auto   end  = [&](double sign) -> Settle
      {
        const Twofold::Two sum   = Twofold::TwoSum(centre.at(free), sign * root.Value.High);
        const double       part  = sum.Low + (sign * root.Value.Low);
        const double       value = sum.High + part;
        const double       bound =
            (2.0 * offRoot) + (Lean::Grain * (std::abs(centre.at(free)) + root.Value.High));
        const bool            sure = Lean::Certain(value, (sum.High - value) + part, bound);
        std::array<double, 2> spot = {};
        spot.at(fixed)             = from.at(fixed);
        spot.at(free)              = value;
        Settle settle              = {.At = {.X = spot[0], .Y = spot[1]}, .Sure = sure};
        /**
         *  A crossing whose bound reaches zero, the center's coordinate and the root's share of
         *  opposite signs, is asked of Axis by Naught.
         */
        settle.Naught.at(free) =
            static_cast<int>(!sure)
                * static_cast<int>(!(std::abs(value) > bound))
                * static_cast<int>(centre.at(free) * (sign * root.Value.High) < 0.0)
            > 0;
        return settle;
      };
      const auto tell = [&](const Settle &settle) -> std::array<int, 2>
      {
        const double value = free == 0 ? settle.At.X : settle.At.Y;
        const double off   = (2.0 * (offRoot + (Lean::Grain * (std::abs(value) + root.Value.High))))
                           + (0x1p-52 * std::abs(value));
        const int    start = static_cast<int>((std::abs(value - low) * Lean::Tight) - off > wide);
        const int    finish = static_cast<int>((std::abs(high - value) * Lean::Tight) - off > wide);
        const int    inside = static_cast<int>(value > low) * static_cast<int>(value < high);
        return {start * finish * ((2 * inside) - 1), Sweeps(round, settle.At, wide)};
      };
      const bool            rising = until.at(free) > from.at(free);
      std::array<Settle, 2> ends   = {end(rising ? -1.0 : 1.0), end(rising ? 1.0 : -1.0)};
      if constexpr (Barred<StraightT>)
      {
        if (!ends[0].Sure || !ends[1].Sure)
        {
          ends = Naught(*line.Kept, arc, ends);
        }
      }
      return Kept({tell(ends[0]), tell(ends[1])}, ends);
    }

    /**
     *  A ready line and a ready arc that cross plainly, or plainly do not, told as Lined tells them
     *  from the true geometry rather than a rough one: from the quadratic of the line's carrier and
     *  the circle, worked out in lean double-double, with r^2 - off^2 = (b^2 - a c) / a for off the
     *  distance from the center to the carrier.
     *
     *  Lined goes on where the carrier cuts the circle at a sine above Steep, r^2 - off^2 above r^2
     *  / 64, here b^2 - a c above a r^2 / 64. The line misses where the carrier passes more than
     *  eight tolerances outside, off^2 - r^2 above 16 L r + 64 L^2 for the tolerance L, here a c -
     *  b^2 above a times that. All else is left to Lined. Each test here is clear of its bound,
     *  each side's few roundings taken by a factor 1 + 2^-48.
     *
     *  The discriminant's bound is the products' share of their parts' bounds, and Grain times
     *  their size for what the products and their difference leave, four times that.
     *
     *  A ray or an endless line is told as a line whose end, or both ends, lie at an infinity of
     *  places, as Cut tells it.
     *
     *  @param  line   The ready line, ray or endless line.
     *  @param  round  The ready arc, of a radius, sweeping at most once round.
     *  @param  length The tolerance used.
     *  @return        The hits, exact, in order along the line; nothing where Lined must decide.
     */
    template <typename StraightT>
    static auto Chorded(StraightT line, const Prepared &round, double length)
        -> std::optional<Sightings>
    {
      const Arc         &arc  = std::get<Arc>(round.Shape());
      const double       wide = Wide * length;
      const Lean::Double a    = Squared(*line.Kept);
      if (!(a.High > 0.0) || !(arc.Radius > Wide * wide))
      {
        return std::nullopt;
      }
      const Lean::Way    way     = Course(line);
      const Lean::Way    gap     = Lean::Apart(Origin(line), arc.Center);
      const Lean::Sized  b       = Lean::Dot(way, gap);
      const Lean::Sized  away    = Lean::Dot(gap, gap);
      const Lean::Double squared = Squared(round);
      const Lean::Double c       = Lean::Plus(away.Value, Lean::Minus(squared));
      const Lean::Double lead    = Lean::Times(b.Value, b.Value);
      const Lean::Double rest    = Lean::Times(a, c);
      const double       offA    = Lean::Grain * a.High;
      const double       offB    = Lean::Grain * b.Size;
      const double       offC    = Lean::Grain * (away.Size + squared.High);
      const Quadratic    shape   = {.A     = a,
                                    .B     = b.Value,
                                    .Disc  = Lean::Plus(lead, Lean::Minus(rest)),
                                    .OffA  = offA,
                                    .OffB  = offB,
                                    .OffD  = (((2.0 * std::abs(b.Value.High) * offB)
                                               + (offB * offB)
                                               + (a.High * offC)
                                               + (std::abs(c.High) * offA)
                                               + (offA * offC))
                                              * Lean::Loose)
                                           + (Lean::Grain * (lead.High + std::abs(rest.High))),
                                    .Width = wide};
      const double       disc    = shape.Disc.High;
      if ((disc * Lean::Tight) - shape.OffD
          > a.High * squared.High * (Met::Steep * Met::Steep) * Lean::Loose)
      {
        return Cut(line, round, shape);
      }
      const double outside =
          a.High * ((16.0 * length * arc.Radius) + (64.0 * length * length)) * Lean::Loose;
      return (-disc * Lean::Tight) - shape.OffD > outside ? std::optional<Sightings>{Sightings{}}
                                                          : std::nullopt;
    }

    /**
     *  Two ready arcs that cross plainly, or plainly do not, told as Rounds tells them from the
     *  true geometry: from their centers' gap g exact, |g|^2, the chord's k = |g|^2 + r1^2 - r2^2
     *  and big = 4 |g|^2 r1^2 - k^2, each in lean double-double with its bound; every crossing M +
     *  H or M - H, M = O1 + g k / (2 |g|^2) and H = g' sqrt(big) / (2 |g|^2), g' a quarter turn of
     *  g.
     *
     *  Rounds goes on where the circles cross more than four tolerances from touching, outside or
     *  within, |g| between |r1 - r2| + 4 L and r1 + r2 - 4 L, at a sine above Steep, which is big
     *  above r1^2 r2^2 / 16; it tells they miss where |g| is past r1 + r2 + 8 L or short of |r1 -
     *  r2| - 8 L; it leaves all else to Intersect. Each test here is on |g|^2, within its bound,
     *  against squares worked out plainly, within 7 u, clear by a factor 1 + 2^-48. Each crossing
     *  is told then as Rounds tells it: where within both sweeps, or near an end of either, as
     *  Swept tells it.
     *
     *  @param  rounds The two ready arcs, each of a radius, sweeping at most once round.
     *  @param  length The tolerance used.
     *  @return        The hits, exact, in order round the first; nothing where Rounds must decide.
     */
    static auto Circled(const Prepared::Pair &rounds, double length) -> std::optional<Sightings>
    {
      const Prepared &one    = rounds[0];
      const Prepared &other  = rounds[1];
      const Arc      &first  = std::get<Arc>(one.Shape());
      const Arc      &second = std::get<Arc>(other.Shape());
      const double    wide   = Wide * length;
      if (!(std::min(first.Radius, second.Radius) > Wide * wide))
      {
        return std::nullopt;
      }
      const Lean::Way    gap   = Lean::Apart(second.Center, first.Center);
      const Lean::Sized  apart = Lean::Dot(gap, gap);
      const Lean::Double inner = Squared(one);
      const Lean::Double outer = Squared(other);
      const Lean::Double chord = Lean::Plus(Lean::Plus(apart.Value, inner), Lean::Minus(outer));
      const Lean::Double four  = Lean::Times(apart.Value, inner);
      const Lean::Double lead  = {.High = 4.0 * four.High, .Low = 4.0 * four.Low};
      const Lean::Double tail  = Lean::Times(chord, chord);
      const Lean::Double big   = Lean::Plus(lead, Lean::Minus(tail));
      const double       offG  = Lean::Grain * apart.Size;
      const double       offK  = offG + (Lean::Grain * (apart.Size + inner.High + outer.High));
      const double       offB =
          (((4.0 * inner.High * offG) + (2.0 * std::abs(chord.High) * offK) + (offK * offK))
           * Lean::Loose)
          + (Lean::Grain * (lead.High + tail.High));
      const double sum    = first.Radius + second.Radius;
      const double differ = std::abs(first.Radius - second.Radius);
      const double least  = (apart.Value.High * Lean::Tight) - offG;
      const double most   = (apart.Value.High * Lean::Loose) + offG;
      const auto   square = [](double side) -> double
      {
        return side * side;
      };
      const bool crossing =
          most < square(sum - (4.0 * length)) * Lean::Tight
          && least > square(differ + (4.0 * length)) * Lean::Loose
          && (big.High * Lean::Tight) - offB
                 > inner.High * outer.High * (4.0 * Met::Steep * Met::Steep) * Lean::Loose;
      if (!crossing)
      {
        const bool apartWide = least > square(sum + (8.0 * length)) * Lean::Loose;
        const bool within =
            differ > 16.0 * length && most < square(differ - (8.0 * length)) * Lean::Tight;
        return apartWide || within ? std::optional<Sightings>{Sightings{}} : std::nullopt;
      }
      const Lean::Double twice   = {.High = 2.0 * apart.Value.High, .Low = 2.0 * apart.Value.Low};
      const double       inverse = 1.0 / twice.High;
      const Lean::Rooted root    = Lean::Root(big);
      const Lean::Double along   = Lean::Over(chord, twice, inverse);
      const Lean::Double across  = Lean::Over(root.Value, twice, inverse);
      const double offRoot  = (offB * root.Inverse * Lean::Loose) + (Lean::Grain * root.Value.High);
      const double offAlong = ((offK + (std::abs(along.High) * 2.0 * offG)) * inverse * Lean::Loose)
                            + (Lean::Grain * std::abs(along.High));
      const double offAcross = ((offRoot + (across.High * 2.0 * offG)) * inverse * Lean::Loose)
                             + (Lean::Grain * across.High);
      const Point  base      = {.X = first.Center.X + (gap[0].High * along.High),
                                .Y = first.Center.Y + (gap[1].High * along.High)};
      const Vector rise      = {.X = -(gap[1].High * across.High), .Y = gap[0].High * across.High};
      const std::array<Point, 2> rough = {base + rise, base - rise};
      const auto                 tell  = [&](std::size_t side) -> std::array<int, 2>
      {
        return {Sweeps(one, rough.at(side), wide), Sweeps(other, rough.at(side), wide)};
      };
      const Tells tells = {tell(0), tell(1)};
      if (Missed(tells))
      {
        return Sightings{};
      }
      if (Edged(tells))
      {
        return std::nullopt;
      }
      const Halves halves = {
          .Middle = {Lean::Moved(first.Center.X, gap[0], along),
                     Lean::Moved(first.Center.Y, gap[1], along)},
          .Half   = {Lean::Minus(Lean::Times(gap[1], across)), Lean::Times(gap[0], across)},
          .Off    = {
              2.0 * ((std::abs(gap[0].High) * offAlong) + (std::abs(gap[1].High) * offAcross)),
              2.0 * ((std::abs(gap[1].High) * offAlong) + (std::abs(gap[0].High) * offAcross))}};
      const auto zero = [&first, &second, &twice](Point spot, double middle) -> bool
      {
        return static_cast<int>(Axis::Fits(twice.High * middle))
                 * static_cast<int>(Axis::On(first, spot))
                 * static_cast<int>(Axis::On(second, spot))
             > 0;
      };
      const std::array<Settle, 2> ends = Chord::Zeroed(Chord::Ends(halves), halves, zero);
      const bool                  swap = Turned(one, rough[0]) > Turned(one, rough[1]);
      return swap ? Kept({tells[1], tells[0]}, {ends[1], ends[0]}) : Kept(tells, ends);
    }

    /**
     *  A line and a ready arc that cross plainly: the line cutting the circle at a sine above
     *  Steep, each crossing far from the line's ends and the arc's, a hit where it lies within
     *  both. A line plainly outside the circle Chorded tells before; one near it is left to
     *  Intersect.
     *
     *  @param  line   The line.
     *  @param  round  The ready arc, of a radius, sweeping at most once round.
     *  @param  length The tolerance used.
     *  @return        The crossings, roughly, in order along the line; nothing where Intersect must
     *                 decide.
     */
    static auto Lined(const Line &line, const Prepared &round, double length)
        -> std::optional<Sightings>
    {
      const Arc   &arc    = std::get<Arc>(round.Shape());
      const double clear  = 4.0 * length;
      const double wide   = Wide * length;
      const Vector way    = line.End - line.Start;
      const double square = Dot(way, way);
      if (!(square > 0.0) || !(arc.Radius > Wide * wide))
      {
        return std::nullopt;
      }
      const double middle = Dot(arc.Center - line.Start, way) / square;
      const double off    = Distance(line.Start + (way * middle), arc.Center);
      const double ratio  = off / arc.Radius;
      if (!(arc.Radius - off > clear) || !((1.0 - (ratio * ratio)) > Met::Steep * Met::Steep))
      {
        return std::nullopt;
      }
      const double reach = std::sqrt(square);
      const double half  = std::sqrt((arc.Radius * arc.Radius) - (off * off)) / reach;
      Sightings    told;
      for (const double along : {middle - half, middle + half})
      {
        const Point point = line.Start + (way * along);
        const int   lain  = Seated(along, reach, wide);
        const int   swept = lain == 0 ? 0 : Swept(round, point, wide);
        if (swept == 0)
        {
          return std::nullopt;
        }
        if (lain > 0 && swept > 0)
        {
          told.Further = along > middle;
          told.Add(point);
        }
      }
      return told;
    }

    /**
     *  Two ready arcs that cross plainly: their circles crossing at a sine above Steep, each
     *  crossing far from the ends of both, a hit where it lies within both sweeps. Circles plainly
     *  apart, or one plainly within the other, Circled tells before; any others are left to
     *  Intersect.
     *
     *  @param  rounds The two ready arcs, each of a radius, sweeping at most once round.
     *  @param  length The tolerance used.
     *  @return        The crossings, roughly; nothing where Intersect must decide.
     */
    static auto Rounds(const Prepared::Pair &rounds, double length) -> std::optional<Sightings>
    {
      const Arc   &first  = std::get<Arc>(rounds[0].get().Shape());
      const Arc   &second = std::get<Arc>(rounds[1].get().Shape());
      const double clear  = 4.0 * length;
      const double wide   = Wide * length;
      const double apart  = Distance(first.Center, second.Center);
      const double sum    = first.Radius + second.Radius;
      const double gap    = std::abs(first.Radius - second.Radius);
      const double near   = std::min(std::abs(apart - sum), std::abs(apart - gap));
      if (!(std::min(first.Radius, second.Radius) > Wide * wide) || !(near > clear))
      {
        return std::nullopt;
      }
      /**
       *  The crossings are worked from the smaller circle: from the larger, a chord much shorter
       *  than its radius is the difference of two squares all but equal, off by far more than the
       *  tolerance. Circles apart, or one within the other, have no crossing: their cosine lies
       *  past 1 in size by far more than its roundings, as they keep four tolerances from touching,
       *  and the sine test below leaves them to Intersect.
       */
      const Arc   &small = first.Radius <= second.Radius ? first : second;
      const Arc   &large = first.Radius <= second.Radius ? second : first;
      const double inner = small.Radius;
      const double outer = large.Radius;
      const double cosine =
          ((inner * inner) + (outer * outer) - (apart * apart)) / (2.0 * inner * outer);
      if (!((1.0 - (cosine * cosine)) > Met::Steep * Met::Steep))
      {
        return std::nullopt;
      }
      const Vector toward = (large.Center - small.Center) / apart;
      const double along  = ((apart * apart) + (inner * inner) - (outer * outer)) / (2.0 * apart);
      const double height = std::sqrt(std::max(0.0, (inner * inner) - (along * along)));
      const Point  base   = small.Center + (toward * along);
      const Vector rise   = Perpendicular(toward) * height;
      Sightings    told;
      for (const Point point : {base + rise, base - rise})
      {
        const int one   = Swept(rounds[0], point, wide);
        const int other = one == 0 ? 0 : Swept(rounds[1], point, wide);
        if (other == 0)
        {
          return std::nullopt;
        }
        if (one > 0 && other > 0)
        {
          told.Add(point);
        }
      }
      return told;
    }

    /**
     *  Whether two ready curves are ones the passes here tell: plain lines or arcs of a radius,
     *  swept at most once round, at an angle within half of Steep.
     *
     *  @param  curves The two ready curves.
     *  @param  usable The tolerance used, no finer than their grain.
     *  @return        Whether they are.
     */
    static auto Sound(const Prepared::Pair &curves, const Tolerance &usable) -> bool
    {
      const Prepared &one   = curves[0];
      const Prepared &other = curves[1];
      return one.m_Plain
          && other.m_Plain
          && (std::holds_alternative<Line>(one.Shape()) || one.m_Round)
          && (std::holds_alternative<Line>(other.Shape()) || other.m_Round)
          && usable.Radians <= Met::Steep / 2.0;
    }

    /**
     *  Whether two ready curves are ones the exact passes tell with a ray or an endless line among
     *  them: each a ray or an endless line Prepared keeps a way for, a plain line, or a plain arc
     *  of a radius swept at most once round, one of them a ray or an endless line; at an angle
     *  within half of Steep; and of a size Shrunk leaves as it is, so that Intersect meets them as
     *  they are, as Prepared keeps them.
     *
     *  @param  curves The two ready curves.
     *  @param  usable The tolerance used, no finer than their grain.
     *  @return        Whether they are.
     */
    static auto Unbounded(const Prepared::Pair &curves, const Tolerance &usable) -> bool
    {
      const Prepared &one     = curves[0];
      const Prepared &other   = curves[1];
      const auto      endless = [](const Prepared &curve) -> bool
      {
        const Curve &shape = curve.Shape();
        return (std::holds_alternative<Ray>(shape) || std::holds_alternative<Xline>(shape))
            && curve.m_Long > 0.0;
      };
      const auto plain = [](const Prepared &curve) -> bool
      {
        return curve.m_Plain && (std::holds_alternative<Line>(curve.Shape()) || curve.m_Round);
      };
      const bool first  = endless(one);
      const bool second = endless(other);
      return (first || second)
          && (first || plain(one))
          && (second || plain(other))
          && usable.Radians <= Met::Steep / 2.0
          && !Outscaled(std::max(one.m_Size, other.m_Size));
    }

    /**
     *  A ready line, ray or endless line and a ready arc that cross plainly, or plainly do not,
     *  told exactly: a level or upright one by Levelled, any other by Chorded.
     *
     *  @param  line   The ready line, ray or endless line.
     *  @param  round  The ready arc, of a radius, sweeping at most once round.
     *  @param  length The tolerance used.
     *  @return        The hits, exact, in order along the line; nothing where Lined must decide.
     */
    template <typename StraightT>
    static auto Ringed(StraightT line, const Prepared &round, double length)
        -> std::optional<Sightings>
    {
      const Vector way = Heading(line);
      return Zero(way.X) || Zero(way.Y) ? Levelled(line, round, length)
                                        : Chorded(line, round, length);
    }

    /**
     *  The crossings of two ready curves, a ray or an endless line among them, as Unbounded takes
     *  them, that cross plainly, or plainly do not, told exactly: two straight ones by Surely, each
     *  read as its kind is, a line as kept and a ray or an endless line as a track; a ray or an
     *  endless line and an arc by Levelled or Chorded.
     *
     *  @param  curves The two ready curves.
     *  @param  usable The tolerance used, no finer than their grain.
     *  @return        The crossings, exact, in order along the first; nothing where a rough pass or
     *                 Intersect must decide.
     */
    static auto Endlessly(const Prepared::Pair &curves, const Tolerance &usable)
        -> std::optional<Sightings>
    {
      const Prepared &one    = curves[0];
      const Prepared &other  = curves[1];
      const double    length = usable.Length;
      const bool      line   = std::holds_alternative<Line>(one.Shape());
      const bool      across = std::holds_alternative<Line>(other.Shape());
      const bool      round  = std::holds_alternative<Arc>(one.Shape());
      const bool      bent   = std::holds_alternative<Arc>(other.Shape());
      if (round || bent)
      {
        const std::optional<Sightings> sure =
            Ringed(Tracked(round ? other : one), round ? one : other, length);
        return !round || !sure ? sure : Ordered(one, *sure);
      }
      if (line)
      {
        return Surely(Bar{&one}, Tracked(other), length);
      }
      if (across)
      {
        return Surely(Tracked(one), Bar{&other}, length);
      }
      return Surely(Tracked(one), Tracked(other), length);
    }

    /**
     *  The crossings of two ready lines, rays, endless lines or arcs that cross plainly, or plainly
     *  do not, told exactly from the true geometry: by Surely, Levelled, Chorded or Circled.
     *
     *  @param  curves The two ready curves.
     *  @param  usable The tolerance used, no finer than their grain.
     *  @return        The crossings, exact, in order along the first; nothing for any other pair,
     *                 or where a rough pass or Intersect must decide.
     */
    static auto Exacted(const Prepared::Pair &curves, const Tolerance &usable)
        -> std::optional<Sightings>
    {
      if (!Sound(curves, usable))
      {
        return Unbounded(curves, usable) ? Endlessly(curves, usable) : std::nullopt;
      }
      const Prepared &one    = curves[0];
      const Prepared &other  = curves[1];
      const bool      line   = std::holds_alternative<Line>(one.Shape());
      const bool      across = std::holds_alternative<Line>(other.Shape());
      if (line && across)
      {
        return Surely(Bar{&one}, Bar{&other}, usable.Length);
      }
      if (line || across)
      {
        const Bar                      straight = {line ? &one : &other};
        const Prepared                &round    = line ? other : one;
        const std::optional<Sightings> sure     = Ringed(straight, round, usable.Length);
        return line || !sure ? sure : Ordered(one, *sure);
      }
      return Circled(curves, usable.Length);
    }

    /**
     *  The crossings of two ready lines or arcs that cross plainly, or plainly do not, roughly: by
     *  Straight, Lined or Rounds, where the exact passes leave them.
     *
     *  @param  curves The two ready curves.
     *  @param  usable The tolerance used, no finer than their grain.
     *  @return        The crossings, in order along a line that is first; nothing for any other
     *                 pair, for an angle past what Steep allows for, or where Intersect must
     *                 decide.
     */
    static auto Roughly(const Prepared::Pair &curves, const Tolerance &usable)
        -> std::optional<Sightings>
    {
      if (!Sound(curves, usable))
      {
        return std::nullopt;
      }
      const Prepared &one    = curves[0];
      const Prepared &other  = curves[1];
      const auto     *line   = std::get_if<Line>(&one.Shape());
      const auto     *across = std::get_if<Line>(&other.Shape());
      if (line != nullptr && across != nullptr)
      {
        return Straight(curves, usable.Length);
      }
      if (line != nullptr)
      {
        return Lined(*line, other, usable.Length);
      }
      if (across != nullptr)
      {
        return Lined(*across, one, usable.Length);
      }
      return Rounds(curves, usable.Length);
    }

    /**
     *  Crossings told roughly made exact: each the double nearest the true crossing, as Intersect
     *  gives it, the nearer of the carriers' two where they cross twice; in order along the first
     *  curve. A line's one crossing with an arc is worked out alone, the carrier's other left
     *  unworked.
     *
     *  @param  curves The two ready curves.
     *  @param  rough  Their crossings, roughly, in order along a line that is first.
     *  @return        The crossings; nothing where one is not certain, for Intersect to work out.
     */
    static auto Settled(const Prepared::Pair &curves, const Sightings &rough)
        -> std::optional<Sightings>
    {
      if (rough.Count == 0)
      {
        return rough;
      }
      const Prepared &first  = curves[0];
      const Prepared &second = curves[1];
      const auto     *line   = std::get_if<Line>(&first.Shape());
      const auto     *across = std::get_if<Line>(&second.Shape());
      /**
       *  A line that crosses an arc once: that crossing alone is worked out, of the carrier's two.
       */
      if (rough.Count == 1 && (line == nullptr) != (across == nullptr))
      {
        const Line &straight = line != nullptr ? *line : *across;
        const Arc  &arc      = std::get<Arc>((line != nullptr ? second : first).Shape());
        return Crossings::Quickly(straight, arc, rough.Further)
            .transform(
                [](Point point) -> Sightings
                {
                  Sightings exact;
                  exact.Add(point);
                  return exact;
                });
      }
      const Crossings::Roots     roots = Crossings::Quickly({&first.Shape(), &second.Shape()});
      const std::optional<Point> added = roots.At[0];
      const std::optional<Point> taken = roots.At[1];
      const int lost = static_cast<int>(!added.has_value())
                     + (static_cast<int>(roots.Twice) * static_cast<int>(!taken.has_value()));
      if (lost > 0)
      {
        return std::nullopt;
      }
      const Point one   = added.value_or(Point{});
      const Point other = taken.value_or(one);
      Sightings   exact;
      for (const Point near : std::span{rough.At}.first(rough.Count))
      {
        exact.Add(Distance(one, near) <= Distance(other, near) ? one : other);
      }
      return Ordered(first, exact);
    }

    /**
     *  Crossings told exact, each a mark of where they cross, in their order.
     *
     *  @param told The crossings, exact.
     *  @param into The marking they are added to.
     */
    static void Marked(const Sightings &told, Marking &into)
    {
      for (const Point point : std::span{told.At}.first(told.Count))
      {
        into.Marks.push_back({.At = point, .Kind = Contact::Cross});
      }
    }

    /**
     *  Two ready arcs of one circle, its center and radius the same to the bit, met as Intersect's
     *  direct way meets them, as Coincident tells it, into a meeting each thread keeps, and its
     *  hits kept as marks: with none of the ways that leave such arcs first asked.
     *
     *  @param  curves The two ready curves.
     *  @param  usable The tolerance used, no finer than their grain.
     *  @param  into   The marking, empty, that what they meet in is written into.
     *  @return        Whether they are such arcs, and were met so.
     */
    static auto Concyclic(const Prepared::Pair &curves, const Tolerance &usable, Marking &into)
        -> bool
    {
      const auto *first  = std::get_if<Arc>(&curves[0].get().Shape());
      const auto *second = std::get_if<Arc>(&curves[1].get().Shape());
      if (first == nullptr || second == nullptr)
      {
        return false;
      }
      thread_local Intersections whole;
      whole.Hits.clear();
      whole.Overlaps.clear();
      if (!Met::Coincident({first, second}, usable, whole))
      {
        return false;
      }
      for (const Hit &hit : whole.Hits)
      {
        into.Marks.push_back({.At = hit.At, .Kind = hit.Kind});
      }
      into.Overlaps.assign(whole.Overlaps.begin(), whole.Overlaps.end());
      return true;
    }

    /**
     *  Two ready curves met as Intersect meets them, into a meeting each thread keeps, and its hits
     *  kept as marks.
     *
     *  @param curves The two ready curves.
     *  @param sane   The tolerance, as Sane makes it.
     *  @param into   The marking, empty, that what they meet in is written into.
     */
    static void Thoroughly(const Prepared::Pair &curves, const Tolerance &sane, Marking &into)
    {
      thread_local Intersections whole;
      Intersect(curves, sane, whole);
      for (const Hit &hit : whole.Hits)
      {
        into.Marks.push_back({.At = hit.At, .Kind = hit.Kind});
      }
      into.Overlaps.assign(whole.Overlaps.begin(), whole.Overlaps.end());
      into.Spent = whole.Spent;
    }

    /**
     *  Two ready lines met as Intersect meets them made ready, by the same steps: told far apart
     *  from their boxes; else, where both are plain, the length used within Coarse and Fine and the
     *  angle within half of Steep, as Abutting tells a level line and an upright one, or as Surely
     *  tells them, and so with no marking made; else by Straight, by Slanted, and a line of no
     *  length on the other as Dotted tells it, then Intersect's thorough ways, into the marking, as
     *  Intersect goes on. Abutting answers as Intersect does, as Surely does where it answers, so
     *  it goes first: a few comparisons tell such a pair, crossing, at a corner or a joint, or
     *  beside each other, where Surely works out a few dozen products and leaves all but the first.
     *
     *  @pre          Both curves are lines.
     *  @param  lines The two ready lines.
     *  @param  sane  The tolerance, as Sane makes it.
     *  @param  into  The marking, written into where it is not told.
     *  @return       What they meet in, where it is told.
     */
    static auto Plainly(const Prepared::Pair &lines, const Tolerance &sane, Marking &into) -> Plain
    {
      const std::optional<double> length = Reached(lines, sane);
      if (!length)
      {
        return {.Told = true};
      }
      const bool fit = lines[0].get().m_Plain
                    && lines[1].get().m_Plain
                    && !(*length > Coarse)
                    && !(*length < Fine)
                    && sane.Radians <= Met::Steep / 2.0;

      const std::optional<Abutment> square =
          fit ? Abutting::Squarely(lines, *length) : std::nullopt;
      if (square)
      {
        const Mark &mark = square->Marks[0];
        return {.At = mark.At, .Count = square->Count, .Kind = mark.Kind, .Told = true};
      }
      const Bar                  one    = {&lines[0].get()};
      const Bar                  other  = {&lines[1].get()};
      const int                  judged = fit ? Judged(one, other, *length) : -1;
      const std::optional<Point> exact  = judged > 0 ? Crossed(one, other) : std::nullopt;
      if (judged == 0 || exact)
      {
        return {.At = exact.value_or(Point{}), .Count = exact ? 1U : 0U, .Told = true};
      }
      into.Marks.clear();
      into.Overlaps.clear();
      into.Spent                             = false;
      const std::optional<Sightings> rough   = fit ? Straight(lines, *length) : std::nullopt;
      const std::optional<Sightings> settled = rough ? Settled(lines, *rough) : std::nullopt;
      if (settled)
      {
        Marked(*settled, into);
        return {};
      }
      const std::optional<Abutment> slanted =
          fit ? Abutting::Slanted(lines, {.Length = *length, .Radians = sane.Radians})
              : std::nullopt;
      /**
       *  A stretch the lines share is no plain hit: it goes into the marking, as Intersect's would.
       */
      if (slanted && slanted->Stretches > 0)
      {
        into.Overlaps.assign(1, slanted->Stretch);
        return {};
      }
      if (slanted)
      {
        const Mark &mark = slanted->Marks[0];
        return {.At = mark.At, .Count = slanted->Count, .Kind = mark.Kind, .Told = true};
      }
      const std::optional<Abutment> dotted =
          fit ? Abutting::Dotted(lines, sane) : std::optional<Abutment>{};
      if (dotted)
      {
        const Mark &mark = dotted->Marks[0];
        return {.At = mark.At, .Count = dotted->Count, .Kind = mark.Kind, .Told = true};
      }
      Thoroughly(lines, sane, into);
      return {};
    }
  };

  Prepared::Prepared(const Curve &curve)
      : m_Size(Sized(curve)), m_Shape(&curve), m_Legwork(Legworked(curve))
  {
    if (std::optional<Met::Standing> stood = Met::Stood(curve))
    {
      m_Standby = std::make_shared<const Standby>(std::move(*stood));
    }
    /**
     *  A box whose edges are not all finite, as an arc's of numbers near the largest double may be,
     *  is no box: Reachable weighs no such pair by its boxes.
     */
    const std::optional<Box> box = Confines(curve);
    m_Box                        = box.value_or(Box{});
    m_Reach         = std::max({std::abs(m_Box.Min.X), std::abs(m_Box.Max.X), std::abs(m_Box.Min.Y),
                                std::abs(m_Box.Max.Y)});
    m_Boxed         = box.has_value() && std::isfinite(m_Reach);
    const auto *arc = std::get_if<Arc>(&curve);
    const auto *line  = std::get_if<Line>(&curve);
    m_Plain           = m_Boxed && (arc != nullptr || line != nullptr);
    const auto *ray   = std::get_if<Ray>(&curve);
    const auto *xline = std::get_if<Xline>(&curve);
    if (ray != nullptr || xline != nullptr)
    {
      /**
       *  A ray or an endless line the exact passes take: its point and its way finite, the way's
       *  largest coordinate within 2^-96 and 2^96, so that Shrunk leaves its pace as it is. Its way
       *  is exact, nothing left by rounding; its length and square are kept as a line's are.
       */
      const Point  start = ray != nullptr ? ray->Start : xline->Through;
      const Vector way   = ray != nullptr ? ray->Way : xline->Way;
      const double reach = std::max(std::abs(way.X), std::abs(way.Y));
      if (!std::isfinite(start.X + start.Y + reach) || !(reach > 0.0) || Outscaled(reach))
      {
        return;
      }
      const Lean::Way    exact  = {Lean::Double{.High = way.X, .Low = 0.0},
                                   Lean::Double{.High = way.Y, .Low = 0.0}};
      const Lean::Double square = Lean::Dot(exact, exact).Value;
      m_Ends                    = {start, start};
      m_Long                    = Length(way);
      m_Square                  = square.High;
      m_Under                   = square.Low;
      m_Inverse                 = 1.0 / square.High;
      return;
    }
    if (line != nullptr)
    {
      const Lean::Way    way    = Lean::Apart(line->End, line->Start);
      const Lean::Double square = Lean::Dot(way, way).Value;
      m_Slip                    = {.X = way[0].Low, .Y = way[1].Low};
      m_Long                    = Length(Vector{.X = way[0].High, .Y = way[1].High});
      m_Square                  = square.High;
      m_Under                   = square.Low;
      m_Inverse                 = square.High > 0.0 ? 1.0 / square.High : 0.0;
      m_Ends                    = {line->Start, line->End};
    }
    if (arc == nullptr || !m_Plain)
    {
      return;
    }
    /**
     *  An arc held as points, as Intersect's direct ways hold it: its ends by the library's sine
     *  and cosine, within a few steps of a double of the exact ones, asked only which side of them
     *  points far off lie.
     */
    const double sweep = arc->Sweep.InRadians();
    const auto   cut   = [arc](Vector way) -> Point
    {
      return arc->Center + (way * arc->Radius);
    };
    const Twofold::Two square = Twofold::TwoProduct(arc->Radius, arc->Radius);
    m_Round  = arc->Radius > 0.0 && std::abs(sweep) > 0.0 && std::abs(sweep) <= 2.0 * Pi;
    m_Square = square.High;
    m_Under  = square.Low;
    m_Sweep  = sweep;
    m_Ends   = {cut(Trig::Way(arc->Start)), cut(Trig::Way(arc->Start, arc->Sweep))};
  }

  void Intersect(const Prepared::Pair &curves, const Tolerance &tolerance, Intersections &into)
  {
    const Tolerance sane = Sane(tolerance);
    into.Hits.clear();
    into.Overlaps.clear();
    into.Spent = false;
    /**
     *  A polyline with remote arcs is met through what stands for it, kept ready beside it: nothing
     *  kept for the polyline itself is worked out from such an arc's far center.
     */
    if (Preparing::Stood(curves[0]) != nullptr || Preparing::Stood(curves[1]) != nullptr)
    {
      Preparing::Slightly(curves, sane, into);
      return;
    }
    const std::optional<double> length = Preparing::Reached(curves, sane);
    if (length)
    {
      const Met::Vetting vetted;
      Entered({.Pair   = {.Shapes = {&curves[0].get().Shape(), &curves[1].get().Shape()},
                          .Legs   = {Preparing::Worked(curves[0]), Preparing::Worked(curves[1])}},
               .Length = *length,
               .Kept   = &curves},
              sane, into);
    }
  }

  void Intersect(const Prepared::Pair &curves, const Tolerance &tolerance, Marking &into)
  {
    const Tolerance sane = Sane(tolerance);
    into.Marks.clear();
    into.Overlaps.clear();
    into.Spent = false;
    if (Preparing::Stood(curves[0]) != nullptr || Preparing::Stood(curves[1]) != nullptr)
    {
      Preparing::Thoroughly(curves, sane, into);
      return;
    }
    const std::optional<double> length = Preparing::Reached(curves, sane);
    if (!length)
    {
      return;
    }
    /**
     *  The passes, each told as Intersect tells it or not at all, so their order is only how fast:
     *  first the exact ones for crossings at a steep angle, then the bands where curves nearly
     *  touch or cross shallowly, then the rough ones for steep crossings left, then ends on the
     *  other curve, then steep crossings with an end near as Roundly or Circularly tells them, and
     *  lines stopping short; last arcs of one circle, which every other pass leaves at once.
     */
    if (!(*length > Coarse) && !(*length < Fine))
    {
      const Tolerance                usable = {.Length = *length, .Radians = sane.Radians};
      const std::optional<Sightings> sure   = Preparing::Exacted(curves, usable);
      if (sure)
      {
        Preparing::Marked(*sure, into);
        return;
      }
      std::optional<Abutment> told = Glancing::Met(curves, usable);
      if (!told)
      {
        const std::optional<Sightings> rough = Preparing::Roughly(curves, usable);
        const std::optional<Sightings> exact =
            rough ? Preparing::Settled(curves, *rough) : std::nullopt;
        if (exact)
        {
          Preparing::Marked(*exact, into);
          return;
        }
        told = Abutting::Met(curves, usable);
      }
      if (!told)
      {
        told = Abutting::Dotted(curves, sane);
      }
      if (!told)
      {
        told = Glancing::Steeply(curves, usable);
      }
      if (told)
      {
        const std::span<const Mark>    marks  = std::span{told->Marks}.first(told->Count);
        const std::span<const Overlap> shared = std::span{&told->Stretch, told->Stretches};
        into.Marks.assign(marks.begin(), marks.end());
        into.Overlaps.assign(shared.begin(), shared.end());
        return;
      }
      if (Preparing::Concyclic(curves, usable, into))
      {
        return;
      }
    }
    Preparing::Thoroughly(curves, sane, into);
  }

  auto Plainly(const Prepared::Pair &lines, const Tolerance &sane, Marking &into) -> Plain
  {
    return Preparing::Plainly(lines, sane, into);
  }

  auto Shrinkless(const Prepared::Pair &kept) -> bool
  {
    return Preparing::Shrinkless(kept);
  }
} // namespace Geoxy
