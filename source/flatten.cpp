/***************************************************************************************************
 *  @file      source/flatten.cpp
 *  @brief     Curves as chords within a tolerance.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "geoxy/flatten.hpp"

#include "basis.hpp"
#include "bow.hpp"
#include "kernel/piece.hpp"
#include "scratch.hpp"
#include "swift/section.hpp"
#include "trig.hpp"
#include "walk.hpp"

#include "geoxy/accuracy.hpp"
#include "geoxy/angle.hpp"
#include "geoxy/arc.hpp"
#include "geoxy/at.hpp"
#include "geoxy/bezier/cubic.hpp"
#include "geoxy/bezier/quadratic.hpp"
#include "geoxy/curve.hpp"
#include "geoxy/domain.hpp"
#include "geoxy/ellipse.hpp"
#include "geoxy/line.hpp"
#include "geoxy/matrix.hpp"
#include "geoxy/point.hpp"
#include "geoxy/polyline.hpp"
#include "geoxy/ray.hpp"
#include "geoxy/segment.hpp"
#include "geoxy/spline.hpp"
#include "geoxy/vector.hpp"
#include "geoxy/xline.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <numbers>
#include <optional>
#include <span>
#include <type_traits>
#include <utility>
#include <variant>
#include <vector>

namespace Geoxy
{
  namespace
  {
    /**
     *  How deep a piece is halved at most: 2^16 chords, MostChords.
     */
    constexpr int Deepest = 16;

    /**
     *  How many points of an arc or a bulge are worked out each from the one before, by a turn,
     *  after one worked out anew: so few that what the turns' roundings add up to stays within
     *  Wheeling's slack.
     */
    constexpr std::size_t Relay = 32;

    /**
     *  A circular sweep to be flattened, and how far rounding moves its points.
     */
    struct Wheeling
    {
      double Radius = 0.0; /**< The radius.                                                  */
      double Sweep  = 0.0; /**< The sweep, in radians.                                       */
      double Slack  = 0.0; /**< How far a point may lie from the place it stands for.        */
      double Skew   = 0.0; /**< How far a step may sweep past an even one, in radians, more. */
    };

    /**
     *  How many even chords a circular sweep takes to stay within a tolerance, and a quarter turn
     *  each.
     *
     *  A chord of sweep s strays r (1 - cos(s / 2)) = 2 r sin²(s / 4) from its arc, and the chord
     *  between two points each within e of a chord's ends strays no more than e farther: the points
     *  of the one at a share of its length are within e of the other's. So each step, with its skew
     *  k, needs 2 r sin²((s + k) / 4) no more than the tolerance t less the slack e, and s is 4
     *  asin(sqrt(w / 2r)) less k for w, the share 1 - 2^-20 of t less e, which leaves room for the
     *  roundings of the square root and the arcsine; asin as acos(1 - x) is 0 for x under 1e-17.
     *  The sweep is taken 2^-40 longer, for the even step's own rounding.
     *
     *  @param  sweeping  The sweep.
     *  @param  tolerance How far the arc may stray.
     *  @return           The count, at least 1 and at most MostChords: MostChords where the
     *                    tolerance is not past the slack and the skew.
     */
    auto Chords(const Wheeling &sweeping, double tolerance) -> std::size_t
    {
      const double within = (tolerance - sweeping.Slack) * (1.0 - 0x1p-20);
      if (!(within > 0.0))
      {
        return MostChords;
      }
      const double span = 2.0 * std::abs(sweeping.Radius);
      const double step =
          (4.0 * Trig::Asin(std::sqrt(within < span ? within / span : 1.0))) - sweeping.Skew;
      const double swept = std::abs(sweeping.Sweep) * (1.0 + 0x1p-40);
      const auto   most  = static_cast<double>(MostChords);
      if (!(swept < step * most))
      {
        return MostChords;
      }
      const double needed = std::max(std::ceil(std::abs(sweeping.Sweep) / (std::numbers::pi / 2.0)),
                                     std::ceil(swept / step));
      return static_cast<std::size_t>(std::clamp(needed, 1.0, most));
    }

    /**
     *  A unit vector, or a displacement, turned by a unit vector's angle, rounded as written.
     *
     *  @param  way  What is turned.
     *  @param  turn The turn, a unit vector.
     *  @return      It turned.
     */
    auto Rolled(Vector way, Vector turn) -> Vector
    {
      return {.X = (way.X * turn.X) - (way.Y * turn.Y), .Y = (way.Y * turn.X) + (way.X * turn.Y)};
    }

    /**
     *  How far rounding a point's coordinates to the nearest doubles moves it, at most, where they
     *  are no larger than given: half the step between the doubles at each, together.
     *
     *  @param  most The largest each coordinate may be, in magnitude.
     *  @return      The distance; infinite for a bound past the doubles.
     */
    auto Halfsteps(Vector most) -> double
    {
      const auto half = [](double bound) -> double
      {
        return std::ldexp(1.0, std::ilogb(std::max(0x1p-1022, bound * (1.0 + 0x1p-40))) - 53);
      };
      return Trig::Hypot(half(most.X), half(most.Y)) * (1.0 + 0x1p-40);
    }

    /**
     *  How an arc's chords are sized: its points, as At gives them or turned from one that is, lie
     *  within Halfsteps of its coordinates' bounds, |c| + |r|, and |r| (2^-45 + 2^-52 (|start| +
     *  |sweep|)) of the circle's points at the angles they stand for, and its steps sweep no more
     *  than 2^-51 (|start| + 2 |sweep|) past an even one (tests/proof/flatten.py).
     *
     *  @param  arc The arc.
     *  @return     Its sweep and the rounding.
     */
    auto Leeways(const Arc &arc) -> Wheeling
    {
      const double radius = std::abs(arc.Radius);
      const double start  = std::abs(arc.Start.InRadians());
      const double sweep  = std::abs(arc.Sweep.InRadians());
      const Vector most   = {.X = std::abs(arc.Center.X) + radius,
                             .Y = std::abs(arc.Center.Y) + radius};
      return {.Radius = radius,
              .Sweep  = sweep,
              .Slack  = Halfsteps(most) + (radius * (0x1p-45 + (0x1p-52 * (start + sweep)))),
              .Skew   = 0x1p-51 * (start + (2.0 * sweep))};
    }

    /**
     *  The powers of a turn, from the nought-th: each the one before turned by it, as many as a run
     *  of points takes.
     *
     *  @param  turn  The turn, a unit vector.
     *  @param  count How many, at most Relay.
     *  @return       The powers, the first count of them.
     */
    auto Powers(Vector turn, std::size_t count) -> std::array<Vector, Relay>
    {
      std::array<Vector, Relay> powers = {};
      powers[0]                        = {.X = 1.0, .Y = 0.0};
      for (std::size_t index = 1; index < count; ++index)
      {
        powers.at(index) = Rolled(powers.at(index - 1), turn);
      }
      return powers;
    }

    /**
     *  An arc's points at even steps of its sweep: every Relay-th as At gives it, at its share of
     *  the sweep, and those after it that point's direction turned by the step's powers, each
     *  scaled by the radius and moved to the center as At moves it; its end as At gives it. Each
     *  point is one turn from its run's first: none waits on the one before.
     *
     *  @param  arc    The arc.
     *  @param  chords How many steps, at least one.
     *  @return        The points.
     */
    auto Wheeled(const Arc &arc, std::size_t chords) -> std::vector<Point>
    {
      std::vector<Point> points;
      points.reserve(chords + 1);
      const auto                      steps = static_cast<double>(chords);
      const std::array<Vector, Relay> powers =
          Powers(Direction(arc.Sweep / steps), std::min(Relay, chords));
      for (std::size_t first = 0; first < chords; first += Relay)
      {
        const Vector way =
            Direction(arc.Start + (arc.Sweep * (static_cast<double>(first) / steps)));
        points.push_back(arc.Center + (way * arc.Radius));
        const std::size_t count = std::min(Relay, chords - first);
        for (std::size_t index = 1; index < count; ++index)
        {
          points.push_back(arc.Center + (Rolled(way, powers.at(index)) * arc.Radius));
        }
      }
      points.push_back(At(arc, 1.0));
      return points;
    }

    /**
     *  How a bulge's chords are sized: its points lie within Halfsteps of its coordinates' bounds,
     *  |s| + D, and 2^-40 D of its circle's, s its start and D the farthest any lies from it, its
     *  chord's length or, past a half turn, its circle's width; and its steps sweep no more than
     *  2^-50 of its sweep past an even one. Along's roundings and the turns' are a few hundred
     *  steps of a double of D; its sweep, and so its circle, rests on the library's arc tangent and
     *  sine.
     *
     *  @param  bow The bulge.
     *  @return     Its sweep and the rounding.
     */
    auto Leeways(const Bight &bow) -> Wheeling
    {
      const double chord = Length(bow.Chord);
      const double half  = std::abs(bow.Half);
      const double reach =
          half < 1.0 && std::abs(bow.Sweep) > std::numbers::pi ? chord / half : chord;
      const Vector most = {.X = std::abs(bow.Start.X) + reach, .Y = std::abs(bow.Start.Y) + reach};
      return {.Radius = chord / (2.0 * std::abs(Trig::Sin(bow.Sweep / 2.0))),
              .Sweep  = bow.Sweep,
              .Slack  = Halfsteps(most) + (0x1p-40 * reach),
              .Skew   = 0x1p-50 * std::abs(bow.Sweep)};
    }

    /**
     *  A bulge's points between its ends at even steps of its sweep: every Relay-th as Along gives
     *  it, at its share of the sweep, and those after it by the chords of the steps from it, added
     *  to its way from the start, which is then added to the start: each point rounded once at the
     *  size of its coordinates, and none waiting on the one before.
     *
     *  The chord of the sweep from a to b is the whole chord turned by (a + b - S) / 2 and scaled
     *  by sin((b - a) / 2) / sin(S / 2), with no point worked out from the circle's center, which a
     *  small bulge puts far off; each step's chord is the one before turned by the step, so the
     *  first k of a run's, from its first point, are its first chord turned by the sum of the
     *  step's first k powers.
     *
     *  @param points Where the points are added.
     *  @param bow    The bulge.
     *  @param chords How many steps.
     */
    void Strided(std::vector<Point> &points, const Bight &bow, std::size_t chords)
    {
      const auto                steps = static_cast<double>(chords);
      const double              step  = bow.Sweep / steps;
      const double              scale = Trig::Sin(step / 2.0) / bow.Half;
      const std::size_t         most  = std::min(Relay, chords - 1);
      std::array<Vector, Relay> sums  = Powers(Direction(Radians(step)), most);
      Vector                    sum   = {};
      for (std::size_t index = 0; index < most; ++index)
      {
        sum            = sum + sums.at(index);
        sums.at(index) = sum;
      }
      for (std::size_t first = 1; first < chords; first += Relay)
      {
        const double swept = bow.Sweep * static_cast<double>(first) / steps;
        const Vector lead  = Turned(bow.Chord * (Trig::Sin(swept / 2.0) / bow.Half),
                                    Radians((swept - bow.Sweep) / 2.0));
        const Vector chord = Turned(bow.Chord * scale, Radians(swept + ((step - bow.Sweep) / 2.0)));
        points.push_back(bow.Start + lead);
        const std::size_t count = std::min(Relay, chords - first);
        for (std::size_t index = 1; index < count; ++index)
        {
          points.push_back(bow.Start + (lead + Rolled(chord, sums.at(index - 1))));
        }
      }
    }

    /**
     *  A piece's chords appended, halving it until each half is flat, the ends of the halves in
     *  order.
     *
     *  @param points    The points so far, ending at the piece's start.
     *  @param piece     The piece.
     *  @param tolerance How far it may stray.
     */
    void Halve(std::vector<Point> &points, const Kernel::Piece &piece, double tolerance)
    {
      struct Pending
      {
        Kernel::Piece Piece = {}; /**< The piece left to flatten. */
        int           Depth = 0;  /**< How often it was halved.   */
      };

      std::vector<Pending> stack = {{.Piece = piece, .Depth = 0}};
      while (!stack.empty())
      {
        const Pending top = stack.back();
        stack.pop_back();
        if (top.Depth == Deepest || Kernel::Flat(top.Piece, tolerance))
        {
          points.push_back(Kernel::End(top.Piece));
          continue;
        }
        const std::array<Kernel::Piece, 2> halves = Kernel::Split(top.Piece, 0.5);
        stack.push_back({.Piece = halves[1], .Depth = top.Depth + 1});
        stack.push_back({.Piece = halves[0], .Depth = top.Depth + 1});
      }
    }

    /**
     *  The points of a curve at even steps of its parameter from 0 to 1.
     *
     *  @tparam CurveT The curve's kind.
     *  @param  curve  The curve.
     *  @param  chords How many steps.
     *  @return        The points, its ends as At gives them.
     */
    template <typename CurveT>
    auto Even(const CurveT &curve, std::size_t chords) -> std::vector<Point>
    {
      std::vector<Point> points;
      points.reserve(chords + 1);
      for (std::size_t index = 0; index <= chords; ++index)
      {
        points.push_back(At(curve, static_cast<double>(index) / static_cast<double>(chords)));
      }
      return points;
    }

    /**
     *  What sizes an ellipse's chords, its map divided by a power of two near its size.
     *
     *  The speed of its parameter t is |M τ(t)|, τ(t) = (−sin t, cos t), whose square is b² + (a² −
     *  b²) sin²(t − tip): least, b, at a tip of its major axis, most, a, a quarter turn on. The
     *  chord of the parameter from m − s / 2 to m + s / 2 strays from the curve by exactly |det M|
     *  (1 − cos(s / 2)) / |M τ(m)|, the unit circle's sagitta times their ratio, as the distance
     *  from the chord's line is the unit circle's scaled by it, and largest at the middle. Between
     *  two tips' quarters the speed is monotone, so a chord whose middle is between two places is
     *  held to the larger ratio of the two; and there the arc never runs past its chord's ends, so
     *  its distance from the chord's line is its distance from the chord.
     */
    struct Sizing
    {
      Matrix Map    = {};  /**< The map, as given: the points are its own.          */
      double Tip    = 0.0; /**< The parameter of a tip of the major axis.           */
      Vector Toward = {};  /**< The unit circle's point there.                      */
      double Least  = 0.0; /**< b², divided twice.                                  */
      double Rise   = 0.0; /**< a² − b², divided twice.                             */
      double Widest = 0.0; /**< a, divided: the ratio at its largest.               */
      double Area   = 0.0; /**< |det M|, divided twice, and its rounding over.      */
      double Within = 0.0; /**< The tolerance, divided, less what rounding may add. */
    };

    /**
     *  What sizes an ellipse's chords to a tolerance.
     *
     *  @param  map       The ellipse's map.
     *  @param  tolerance How far the curve may stray from its chords.
     *  @return           The sizing; nothing where the tolerance is not past the rounding of the
     *                    points, or not a number, or the map is not finite.
     */
    auto Calibrated(const Matrix &map, double tolerance) -> std::optional<Sizing>
    {
      const double most =
          std::max({std::abs(map.X.X), std::abs(map.X.Y), std::abs(map.Y.X), std::abs(map.Y.Y)});
      int power = 0;
      std::frexp(most, &power);
      const double scale  = std::ldexp(1.0, power);
      const Vector x      = map.X / scale;
      const Vector y      = map.Y / scale;
      const double mean   = (Dot(x, x) + Dot(y, y)) / 2.0;
      const double apart  = (Dot(x, x) - Dot(y, y)) / 2.0;
      const double across = Dot(x, y);
      const double half   = Trig::Hypot(apart, across);
      const double square = mean + half;
      const double area   = std::abs(Cross(x, y));
      const double slack =
          (0x1p-42 * 4.0 * most) + (0x1p-50 * (std::abs(map.Origin.X) + std::abs(map.Origin.Y)));
      const double within = (tolerance - slack) * (1.0 - 0x1p-16) / scale;
      if (!(within > 0.0) || !std::isfinite(square))
      {
        return std::nullopt;
      }
      const double tip = Trig::Atan2(across, apart) / 2.0;
      return Sizing{.Map    = map,
                    .Tip    = tip,
                    .Toward = Trig::Way(tip),
                    .Least  = square > 0.0 ? area * area / square : 0.0,
                    .Rise   = 2.0 * half,
                    .Widest = std::sqrt(square),
                    .Area   = (area * (1.0 + 0x1p-49)) + (0x1p-49 * square),
                    .Within = within};
    }

    /**
     *  How far a chord strays at a place of an ellipse's unit circle for each of the unit circle's:
     *  |det M| / |M τ|, at most a.
     *
     *  The square of the speed is told by the sine of the turn from the tip, so it keeps its digits
     *  by a flat ellipse's tip, and the area is padded by the rounding that leaves.
     *
     *  @param  sizing The sizing.
     *  @param  unit   The place on the unit circle.
     *  @return        The ratio.
     */
    auto Ratio(const Sizing &sizing, Vector unit) -> double
    {
      const double sine  = Cross(sizing.Toward, unit);
      const double speed = std::sqrt(sizing.Least + (sizing.Rise * sine * sine));
      if (!(sizing.Area < sizing.Widest * speed))
      {
        return sizing.Widest;
      }
      return sizing.Area / speed;
    }

    /**
     *  The longest step of an ellipse's parameter whose chord strays no more than the tolerance
     *  where a chord strays a ratio of the unit circle's, a quarter turn at most.
     *
     *  The unit circle's chord of a step s strays 1 − cos(s / 2) = 2 sin²(s / 4), so the step is 4
     *  asin(sqrt(k)) for k = tolerance / (2 ratio), a quarter turn where k is sin²(pi / 8); 4
     *  sqrt(k) (1 + k / 6) is short of it by less than a part in 600, and never over.
     *
     *  @param  sizing The sizing.
     *  @param  ratio  The ratio.
     *  @return        The step, in radians.
     */
    auto Stepped(const Sizing &sizing, double ratio) -> double
    {
      constexpr double Quartering = 0.14644660940672624;
      if (!(2.0 * ratio * Quartering > sizing.Within))
      {
        return std::numbers::pi / 2.0;
      }
      const double share = sizing.Within / (2.0 * ratio);
      return 4.0 * std::sqrt(share) * (1.0 + (share / 6.0));
    }

    /**
     *  How far a run of even chords may reach from a place on an ellipse, toward its minor axis:
     *  sixteen steps, or less where the speed doubles before them, where steps that size would
     *  waste chords; two steps at least, as each run costs its end and its turn.
     *
     *  Where the speed doubles is found from the sine of the turn from the tip by the tangent of
     *  arcsine there, a little short of it.
     *
     *  @param  sizing The sizing.
     *  @param  unit   The place on the unit circle.
     *  @param  step   The step the run takes.
     *  @return        The reach, in radians.
     */
    auto Doubling(const Sizing &sizing, Vector unit, double step) -> double
    {
      const double sine   = std::abs(Cross(sizing.Toward, unit));
      const double cosine = std::abs(Dot(sizing.Toward, unit));
      const double grown  = std::min(
          sizing.Rise, (4.0 * (sizing.Least + (sizing.Rise * sine * sine))) - sizing.Least);
      const double root = std::sqrt(sizing.Rise);
      const double rise = std::max(0.0, std::sqrt(grown) - (sine * root));
      const double run  = cosine * root;
      const double most = 16.0 * step;
      if (!(rise < most * run))
      {
        return most;
      }
      return std::max(2.0 * step, rise / run);
    }

    /**
     *  A place of an ellipse's unit circle, and its parameter.
     */
    struct Anchor
    {
      double Angle = 0.0; /**< The parameter, in radians.     */
      Vector Unit  = {};  /**< The unit circle's point there. */
    };

    /**
     *  Where a quarter of an ellipse between the directions of its axes starts: its unit circle's
     *  point, the tip's turned by quarter turns, exactly.
     *
     *  @param  sizing The sizing.
     *  @param  index  Which quarter, a whole number, the one after the tip's 0.
     *  @return        Its start.
     */
    auto Corner(const Sizing &sizing, double index) -> Anchor
    {
      const Vector                toward = sizing.Toward;
      const std::array<Vector, 4> turned = {toward, Vector{.X = -toward.Y, .Y = toward.X},
                                            Vector{.X = -toward.X, .Y = -toward.Y},
                                            Vector{.X = toward.Y, .Y = -toward.X}};
      const double                turn   = index - (4.0 * std::floor(index / 4.0));
      return {.Angle = sizing.Tip + (index * (std::numbers::pi / 2.0)),
              .Unit  = turned.at(static_cast<std::size_t>(turn))};
    }

    /**
     *  A point of an ellipse at a place of its unit circle.
     *
     *  @param  map  The ellipse's map.
     *  @param  unit The place.
     *  @return      The point.
     */
    auto Ellipsed(const Matrix &map, Vector unit) -> Point
    {
      return map.Origin + ((map.X * unit.X) + (map.Y * unit.Y));
    }

    /**
     *  The points of a piece of an ellipse within one quarter, from the end nearer a tip of its
     *  major axis to the other, both ends among them: runs of even chords, each sized by the ratio
     *  at its start, the larger, its points turned one from the next and its end found anew.
     *
     *  @param  sizing The sizing.
     *  @param  span   The end nearer the tip, then the other.
     *  @param  points Where the points are added.
     *  @return        Whether they are no more than MostChords and one in all.
     */
    auto Graded(const Sizing &sizing, std::array<Anchor, 2> span, std::vector<Point> &points)
        -> bool
    {
      Anchor       from  = span[0];
      const Anchor until = span[1];
      const double way   = until.Angle > from.Angle ? 1.0 : -1.0;
      double       ratio = Ratio(sizing, from.Unit);
      points.push_back(Ellipsed(sizing.Map, from.Unit));
      while (way * (until.Angle - from.Angle) > 0.0)
      {
        const double step   = Stepped(sizing, ratio);
        const double ahead  = from.Angle + (way * Doubling(sizing, from.Unit, step));
        const Anchor next   = way * (until.Angle - ahead) > 0.0
                                ? Anchor{.Angle = ahead, .Unit = Trig::Way(ahead)}
                                : until;
        const double turn   = next.Angle - from.Angle;
        const double chords = std::max(1.0, std::ceil((std::abs(turn) / step) - 0x1p-20));
        if (!(static_cast<double>(points.size()) + chords <= static_cast<double>(MostChords + 1)))
        {
          return false;
        }
        const double cosine = Trig::Cos(turn / chords);
        const double sine   = Trig::Sin(turn / chords);
        Vector       unit   = from.Unit;
        const auto   count  = static_cast<std::size_t>(chords);
        for (std::size_t index = 1; index < count; ++index)
        {
          unit = {.X = (unit.X * cosine) - (unit.Y * sine),
                  .Y = (unit.Y * cosine) + (unit.X * sine)};
          points.push_back(Ellipsed(sizing.Map, unit));
        }
        points.push_back(Ellipsed(sizing.Map, next.Unit));
        from  = next;
        ratio = Ratio(sizing, next.Unit);
      }
      return true;
    }

    /**
     *  The points of a piece of an ellipse within one quarter, after those before it and in the
     *  order of its sweep: walked from the end nearer a tip of its major axis, and turned round
     *  where that is the piece's end. The last point before it, its start, is found anew.
     *
     *  @param  sizing The sizing.
     *  @param  piece  Its start and its end.
     *  @param  points Where the points are added.
     *  @return        Whether they are no more than MostChords and one in all.
     */
    auto Pieced(const Sizing &sizing, std::array<Anchor, 2> piece, std::vector<Point> &points)
        -> bool
    {
      if (!points.empty())
      {
        points.pop_back();
      }
      const auto mark = static_cast<std::ptrdiff_t>(points.size());
      if (!(Ratio(sizing, piece[0].Unit) < Ratio(sizing, piece[1].Unit)))
      {
        return Graded(sizing, piece, points);
      }
      if (!Graded(sizing, {piece[1], piece[0]}, points))
      {
        return false;
      }
      std::reverse(points.begin() + mark, points.end());
      return true;
    }

    /**
     *  An ellipse's points within a tolerance, quarter by quarter between the directions of its
     *  axes. A piece of a quarter of a sweep under 2^-20 so short that the curve along it is within
     *  2^-17 of the tolerance, its sweep times a, is left to the next: an end a hair from a tip
     *  took a chord of no length. Its end, At's, stands for the next's, and the chord from it
     *  strays no more than that farther, which the 2^-16 the sizing leaves out holds.
     *
     *  @param  sizing  The sizing.
     *  @param  ellipse The ellipse.
     *  @return         The points, its ends as At gives them; nothing where they would be more than
     *                  MostChords and one, or its sweep would take more quarters, or it starts
     *                  nowhere.
     */
    auto Quartered(const Sizing &sizing, const Ellipse &ellipse)
        -> std::optional<std::vector<Point>>
    {
      constexpr double Quarter = std::numbers::pi / 2.0;
      const double     start   = ellipse.Start.InRadians();
      const double     sweep   = ellipse.Sweep.InRadians();
      const double     end     = start + sweep;
      if (!std::isfinite(start) || !(std::abs(sweep) <= Quarter * static_cast<double>(MostChords)))
      {
        return std::nullopt;
      }
      const double       first   = std::floor((start - sizing.Tip) / Quarter);
      const double       last    = std::floor((end - sizing.Tip) / Quarter);
      const double       way     = sweep > 0.0 ? 1.0 : -1.0;
      const auto         count   = static_cast<std::size_t>(std::abs(last - first));
      const Anchor       opening = {.Angle = start, .Unit = Trig::Way(start)};
      const Anchor       closing = {.Angle = end, .Unit = Trig::Way(end)};
      std::vector<Point> points;
      points.reserve(64);
      for (std::size_t index = 0; index <= count; ++index)
      {
        const double quarter = first + (way * static_cast<double>(index));
        const Anchor entry   = Corner(sizing, way > 0.0 ? quarter : quarter + 1.0);
        const Anchor exit    = Corner(sizing, way > 0.0 ? quarter + 1.0 : quarter);
        const Anchor from    = way * (start - entry.Angle) > 0.0 ? opening : entry;
        const Anchor until   = way * (exit.Angle - end) > 0.0 ? closing : exit;
        const double spans   = way * (until.Angle - from.Angle);
        if ((spans > 0x1p-20 || spans * sizing.Widest > sizing.Within * 0x1p-17)
            && !Pieced(sizing, {from, until}, points))
        {
          return std::nullopt;
        }
      }
      if (points.empty())
      {
        return std::vector<Point>{At(ellipse, 0.0), At(ellipse, 1.0)};
      }
      points.front() = At(ellipse, 0.0);
      points.back()  = At(ellipse, 1.0);
      return points;
    }

    /**
     *  A place of a Bézier curve and its derivative there.
     */
    struct Station
    {
      double Along = 0.0; /**< The parameter.                   */
      Point  At    = {};  /**< The point, as At gives it.       */
      Vector Way   = {};  /**< The derivative by the parameter. */
    };

    /**
     *  A quadratic Bézier curve's place, by de Casteljau's algorithm as At takes it.
     *
     *  @param  quadratic The curve.
     *  @param  along     The parameter.
     *  @return           The place.
     */
    auto Stationed(const Bezier::Quadratic &quadratic, double along) -> Station
    {
      const Point start = Lerp(quadratic.Start, quadratic.Control, along);
      const Point end   = Lerp(quadratic.Control, quadratic.End, along);
      return {.Along = along, .At = Lerp(start, end, along), .Way = (end - start) * 2.0};
    }

    /**
     *  A cubic Bézier curve's place, by de Casteljau's algorithm as At takes it.
     *
     *  @param  cubic The curve.
     *  @param  along The parameter.
     *  @return       The place.
     */
    auto Stationed(const Bezier::Cubic &cubic, double along) -> Station
    {
      const Point first  = Lerp(cubic.Start, cubic.Out, along);
      const Point second = Lerp(cubic.Out, cubic.In, along);
      const Point third  = Lerp(cubic.In, cubic.End, along);
      const Point start  = Lerp(first, second, along);
      const Point end    = Lerp(second, third, along);
      return {.Along = along, .At = Lerp(start, end, along), .Way = (end - start) * 3.0};
    }

    /**
     *  A curve's derivative as a polynomial: its coefficients of the square, the parameter and one.
     */
    using Velocities = std::array<Vector, 3>;

    /**
     *  A quadratic Bézier curve's derivative.
     *
     *  @param  quadratic The curve.
     *  @return           Its coefficients.
     */
    auto Derived(const Bezier::Quadratic &quadratic) -> Velocities
    {
      const Vector out    = (quadratic.Control - quadratic.Start) * 2.0;
      const Vector inward = (quadratic.End - quadratic.Control) * 2.0;
      return {Vector{}, inward - out, out};
    }

    /**
     *  A cubic Bézier curve's derivative.
     *
     *  @param  cubic The curve.
     *  @return       Its coefficients.
     */
    auto Derived(const Bezier::Cubic &cubic) -> Velocities
    {
      const Vector out    = (cubic.Out - cubic.Start) * 3.0;
      const Vector middle = (cubic.In - cubic.Out) * 3.0;
      const Vector inward = (cubic.End - cubic.In) * 3.0;
      return {out - (middle * 2.0) + inward, (middle - out) * 2.0, out};
    }

    /**
     *  Where a Bézier curve's chords are cut: at 0, at 1, and between them at the real parts of the
     *  roots of its derivative taken as a complex polynomial, where it turns hardest, through any
     *  cusp or turn back, so that no chord's stretch of the curve doubles back over it.
     *
     *  @param  hodograph Its derivative.
     *  @return           The cuts in order, and how many past the first; the rest are 1.
     */
    auto Cut(const Velocities &hodograph) -> std::pair<std::array<double, 4>, std::size_t>
    {
      const Vector          square   = hodograph[0];
      const Vector          linear   = hodograph[1];
      const Vector          constant = hodograph[2];
      const double          leading  = Dot(square, square);
      std::array<double, 2> roots    = {-1.0, -1.0};
      if (leading > 0.0)
      {
        const Vector apart = {.X = ((linear.X * linear.X) - (linear.Y * linear.Y))
                                 - (4.0 * ((square.X * constant.X) - (square.Y * constant.Y))),
                              .Y = (2.0 * linear.X * linear.Y)
                                 - (4.0 * ((square.X * constant.Y) + (square.Y * constant.X)))};
        const double size  = Length(apart);
        Vector       root  = {.X = std::sqrt((size + apart.X) / 2.0),
                              .Y = std::copysign(std::sqrt((size - apart.X) / 2.0), apart.Y)};
        if (Dot(linear, root) < 0.0)
        {
          root = root * -1.0;
        }
        const Vector half = (linear + root) * -0.5;
        roots[0]          = Dot(half, square) / leading;
        const double last = Dot(half, half);
        if (last > 0.0)
        {
          roots[1] = Dot(constant, half) / last;
        }
      }
      else if (Dot(linear, linear) > 0.0)
      {
        roots[0] = -Dot(constant, linear) / Dot(linear, linear);
      }
      std::array<double, 4> cuts  = {0.0, 1.0, 1.0, 1.0};
      std::size_t           count = 1;
      for (const double root : roots)
      {
        if (root > 0.0 && root < 1.0)
        {
          cuts.at(count++) = root;
        }
      }
      /**
       *  The cuts past count are 1, no less than any root kept: sorting all four leaves the first
       *  count as sorting those alone would, and GCC sees the array's size.
       */
      std::ranges::sort(cuts);
      return {cuts, count};
    }

    /**
     *  The least and the most of a cubic polynomial over 0 to 1, from its Bernstein coefficients:
     *  at its ends, or where its derivative, a quadratic, is nothing, by the quotients of the
     *  stable formula, their divisors kept off nothing so that a root that is not there falls at an
     *  end.
     *
     *  @param  bernstein The coefficients.
     *  @return           The least, then the most.
     */
    auto Bounds(std::array<double, 4> bernstein) -> std::array<double, 2>
    {
      constexpr double Tiny   = 0x1p-1000;
      const double     first  = bernstein[1] - bernstein[0];
      const double     second = bernstein[2] - bernstein[1];
      const double     third  = bernstein[3] - bernstein[2];
      const double     square = first - (2.0 * second) + third;
      const double     linear = 2.0 * (second - first);
      const double     apart  = std::max(0.0, (linear * linear) - (4.0 * square * first));
      const double     half   = -(linear + std::copysign(std::sqrt(apart), linear)) / 2.0;
      const auto       valued = [&bernstein](double along) -> double
      {
        const double rest = 1.0 - along;
        return (bernstein[0] * rest * rest * rest)
             + (3.0 * bernstein[1] * rest * rest * along)
             + (3.0 * bernstein[2] * rest * along * along)
             + (bernstein[3] * along * along * along);
      };
      const double one = valued(
          std::clamp(half / std::copysign(std::max(std::abs(square), Tiny), square), 0.0, 1.0));
      const double other =
          valued(std::clamp(first / std::copysign(std::max(std::abs(half), Tiny), half), 0.0, 1.0));
      return {std::min({bernstein[0], bernstein[3], one, other}),
              std::max({bernstein[0], bernstein[3], one, other})};
    }

    /**
     *  How far a Bézier curve's stretch between two places strays from the chord between them, at
     *  most, exactly but for rounding: the stretch, a cubic or less, is the Bézier curve of the two
     *  places and a third of its parameter's span along each tangent, so its distance from the
     *  chord's line and its place along the chord are cubics of known Bernstein coefficients. Where
     *  the place's inner coefficients are within the chord the stretch is too, and its distance
     *  from the chord is from the line; a point of it past an end is as far as it is from the line
     *  and past the end together.
     *
     *  @param  from  The place the chord starts.
     *  @param  until The place it ends.
     *  @return       The distance; infinite for a chord of no length.
     */
    auto Strays(const Station &from, const Station &until) -> double
    {
      const Vector chord  = until.At - from.At;
      const double square = Dot(chord, chord);
      if (!(square > 0.0))
      {
        return std::numeric_limits<double>::infinity();
      }
      const double                third  = (until.Along - from.Along) / 3.0;
      const double                length = std::sqrt(square);
      const std::array<double, 2> aside  = Bounds({0.0, third * Cross(from.Way, chord) / length,
                                                   -third * Cross(until.Way, chord) / length, 0.0});
      const double                away   = std::max(-aside[0], aside[1]);
      const double                leave  = third * Dot(from.Way, chord) / square;
      const double                arrive = 1.0 - (third * Dot(until.Way, chord) / square);
      if (std::max({-leave, leave - 1.0, -arrive, arrive - 1.0}) <= 0.0)
      {
        return away;
      }
      const std::array<double, 2> along = Bounds({0.0, leave, arrive, 1.0});
      return Trig::Hypot(away, std::max({0.0, -along[0], along[1] - 1.0}) * length);
    }

    /**
     *  Where a Bézier curve's chords are tried to, along one stretch between cuts.
     */
    struct Stride
    {
      double End    = 0.0; /**< Where the stretch ends.                                */
      double Span   = 0.0; /**< The span of the parameter tried first.                 */
      double Within = 0.0; /**< How far a chord may stray, less what rounding may add. */
    };

    /**
     *  The longest chord tried from a place of a Bézier curve whose stretch strays no more than the
     *  tolerance: tries grown or shrunk by the square root of how far the last strayed, as a short
     *  chord strays as its span squared; a chord of no length is too short to tell, and grown.
     *
     *  @tparam CurveT The curve's kind.
     *  @param  curve  The curve.
     *  @param  from   The place.
     *  @param  stride The stretch, its span grown or shrunk by the tries.
     *  @return        Where the chord ends; nothing if forty tries found none.
     */
    template <typename CurveT>
    auto Longest(const CurveT &curve, const Station &from, Stride &stride) -> std::optional<Station>
    {
      std::optional<Station> best;
      for (int tries = 0; tries < 40; ++tries)
      {
        const double  reach = from.Along + stride.Span;
        const Station until =
            Stationed(curve, stride.End - reach <= stride.Span * 0x1p-20 ? stride.End : reach);
        const double away = Strays(from, until);
        if (away <= stride.Within)
        {
          best = until;
          if (!(until.Along < stride.End) || away > 0.7 * stride.Within)
          {
            break;
          }
          stride.Span *=
              std::min(2.0, 0.95 * std::sqrt(stride.Within / std::max(away, stride.Within / 16.0)));
          continue;
        }
        if (best)
        {
          break;
        }
        const bool blank = std::isinf(away) && until.Along < stride.End;
        stride.Span *= blank ? 2.0 : std::clamp(0.95 * std::sqrt(stride.Within / away), 0.1, 0.9);
      }
      return best;
    }

    /**
     *  How a Bézier curve's chords are tried: how far they may stray, and where each stretch's
     *  first is tried to.
     */
    struct Tempo
    {
      double Within = 0.0;   /**< How far a chord may stray, less what rounding may add. */
      double Span   = 0.0;   /**< The span of the last chord no cut cut short.           */
      bool   Carry  = false; /**< Whether a stretch's first chord is tried to that span,
                                  not to the stretch's end: the curve one piece of many,
                                  as fast as the one before.                             */
    };

    /**
     *  A Bézier curve's chords, greedy: from each place the longest chord tried, the chords ending
     *  at the cuts.
     *
     *  @tparam CurveT The curve's kind.
     *  @param  curve  The curve.
     *  @param  pace   How the chords are tried, the span of the last no cut cut short kept.
     *  @param  points Where the points go, the curve's start among them.
     *  @return        Whether every chord was found and they are no more than MostChords.
     */
    template <typename CurveT>
    auto Greedy(const CurveT &curve, Tempo &pace, std::vector<Point> &points) -> bool
    {
      const auto [cuts, count] = Cut(Derived(curve));
      Station from             = Stationed(curve, 0.0);
      points.push_back(from.At);
      for (std::size_t index = 0; index < count; ++index)
      {
        Stride stride = {.End = cuts.at(index + 1), .Span = 0.0, .Within = pace.Within};
        stride.Span   = stride.End - from.Along;
        if (pace.Carry && pace.Span > 0.0)
        {
          stride.Span = std::min(stride.Span, pace.Span);
        }
        while (from.Along < stride.End)
        {
          const std::optional<Station> best = Longest(curve, from, stride);
          if (!best || points.size() > MostChords)
          {
            return false;
          }
          stride.Span = best->Along - from.Along;
          if (best->Along < stride.End)
          {
            pace.Span = stride.Span;
          }
          from = *best;
          points.push_back(from.At);
        }
      }
      return true;
    }

    /**
     *  A Bézier curve as chords: greedy where the tolerance is past rounding and the chords are
     *  found; else halved until each half is flat.
     *
     *  @tparam CurveT    The curve's kind.
     *  @param  curve     The curve.
     *  @param  points    Its control points.
     *  @param  tolerance How far it may stray.
     *  @return           Its points.
     */
    template <typename CurveT>
    auto Chorded(const CurveT &curve, std::span<const Point> points, double tolerance)
        -> std::vector<Point>
    {
      double most = 0.0;
      for (const Point point : points)
      {
        most = std::max({most, std::abs(point.X), std::abs(point.Y)});
      }
      Tempo pace = {
          .Within = (tolerance - (0x1p-44 * most)) * (1.0 - 0x1p-20), .Span = 0.0, .Carry = false};
      std::vector<Point> chords;
      if (pace.Within > 0.0 && Greedy(curve, pace, chords))
      {
        return chords;
      }
      chords = {points.front()};
      Halve(chords, Kernel::Polynomial(points), tolerance);
      return chords;
    }

    /**
     *  A piece of a chain of chords appended to it: greedy where the tolerance is past rounding and
     *  the chords are found, else halved until each half is flat; its start, the chain's last point
     *  all but to the bit, left out.
     *
     *  @tparam CurveT The piece's kind.
     *  @param  curve  The piece, as a Bézier curve.
     *  @param  pace   How far the chain may stray, the span its last chord no cut cut short took
     *                 kept, each stretch's first chord tried to it.
     *  @param  points The chain, ending where the piece starts.
     */
    template <typename CurveT>
    void Tacked(const CurveT &curve, Tempo &pace, std::vector<Point> &points)
    {
      const Kernel::Piece piece = Kernel::Polynomial(curve);
      double              most  = 0.0;
      for (const Homogeneous &pole : Kernel::Poles(piece))
      {
        most = std::max({most, std::abs(pole.X), std::abs(pole.Y)});
      }
      Tempo             held = {.Within = (pace.Within - (0x1p-44 * most)) * (1.0 - 0x1p-20),
                                .Span   = pace.Span,
                                .Carry  = true};
      const std::size_t mark = points.size();
      if (held.Within > 0.0 && Greedy(curve, held, points))
      {
        points.erase(points.begin() + static_cast<std::ptrdiff_t>(mark));
        pace.Span = held.Span;
        return;
      }
      points.resize(mark);
      Halve(points, piece, pace.Within);
    }

    /**
     *  A polynomial spline of degree two or three as chords, span by span, each span's Bézier piece
     *  greedy as a Bézier curve's: its control points within a few steps of a double of the true
     *  ones, which the tolerance less what rounding may add leaves room for many times over. Its
     *  ends are At's, a few steps of a double from the pieces'.
     *
     *  @param  spline    The spline.
     *  @param  tolerance How far it may stray, above nothing.
     *  @return           Its points; nothing for a spline the kernels do not cut into such pieces.
     */
    auto Chordwise(const Spline &spline, double tolerance) -> std::optional<std::vector<Point>>
    {
      const std::optional<List<Swift::Section>> sections = Swift::Sectioned(spline);
      if (!sections)
      {
        return std::nullopt;
      }
      std::vector<Point> points;
      points.reserve(32 * sections->size());
      points.push_back(sections->front().Poles[0]);
      Tempo pace = {.Within = tolerance, .Span = 0.0, .Carry = true};
      for (const Swift::Section &section : *sections)
      {
        const std::array<Point, 4> &poles = section.Poles;
        if (section.Degree == 2)
        {
          Tacked(Bezier::Quadratic{.Start = poles[0], .Control = poles[1], .End = poles[2]}, pace,
                 points);
        }
        else
        {
          Tacked(Bezier::Cubic{.Start = poles[0], .Out = poles[1], .In = poles[2], .End = poles[3]},
                 pace, points);
        }
      }
      const Range domain = Domain(spline);
      points.front()     = At(spline, domain.Start);
      points.back()      = At(spline, domain.End);
      return points;
    }

    /**
     *  The power of two that brings a Bézier curve's control points near 1, where they lie past
     *  2^400 or below 2^-400: there the square of a chord's deviation, of the tolerance or of a
     *  coordinate overflows or underflows, and a curve of 1e-160 strayed past its tolerance a
     *  thousand times.
     *
     *  @param  points The control points.
     *  @return        The power, 0 where they need none, are all nothing, or are not finite.
     */
    auto LeveledPower(std::span<const Point> points) -> int
    {
      double most = 0.0;
      for (const Point point : points)
      {
        most = std::max({most, std::abs(point.X), std::abs(point.Y)});
      }
      const bool moderate = most >= 0x1p-400 && most <= 0x1p400;
      if (moderate || !(most > 0.0) || !std::isfinite(most))
      {
        return 0;
      }
      return -(std::ilogb(most) + 1);
    }

    /**
     *  A point times a power of two, exactly where it stays a normal double.
     *
     *  @param  point The point.
     *  @param  power The power of two.
     *  @return       The point scaled.
     */
    auto Unshrunk(Point point, int power) -> Point
    {
      return {.X = std::ldexp(point.X, power), .Y = std::ldexp(point.Y, power)};
    }

    /**
     *  Chords worked out on a curve scaled by a power of two, scaled back: as exact as the chords
     *  of the curve itself, its every sum, product and quotient scaled alike, but with none of them
     *  past the doubles.
     *
     *  @param  chords The chords of the scaled curve.
     *  @param  power  The power the curve was scaled by.
     *  @return        The chords of the curve.
     */
    auto Unscaled(std::vector<Point> chords, int power) -> std::vector<Point>
    {
      for (Point &point : chords)
      {
        point = Unshrunk(point, -power);
      }
      return chords;
    }

    /**
     *  A rational piece's homogeneous blossom: de Casteljau's rounds, each at its own parameter.
     *
     *  At (a, ..., a, b, ..., b), k of them b, it is the k-th control point of the piece's stretch
     *  from a to b; at one parameter repeated, the piece's point there.
     *
     *  @param  piece The piece, of degree three at most.
     *  @param  args  The parameters, the first Degree of them.
     *  @return       The blossom.
     */
    auto Homoblossom(const Kernel::Piece &piece, std::array<double, 3> args) -> Homogeneous
    {
      std::array<Homogeneous, 4> cells = {};
      for (std::size_t index = 0; index <= piece.Degree; ++index)
      {
        cells.at(index) = piece.Poles.at(index);
      }
      for (std::size_t round = 1; round <= piece.Degree; ++round)
      {
        const double along = args.at(round - 1);
        const double rest  = 1.0 - along;
        for (std::size_t index = 0; index + round <= piece.Degree; ++index)
        {
          const Homogeneous lhs = cells.at(index);
          const Homogeneous rhs = cells.at(index + 1);
          cells.at(index)       = {.X = (lhs.X * rest) + (rhs.X * along),
                                   .Y = (lhs.Y * rest) + (rhs.Y * along),
                                   .W = (lhs.W * rest) + (rhs.W * along)};
        }
      }
      return cells[0];
    }

    /**
     *  How far a rational piece's stretch strays from the chord between its ends, at most.
     *
     *  The stretch's control points are blossoms, over their weights, which are above nothing: the
     *  stretch lies within their hull, so no farther from the chord than the farthest of them.
     *  Where each lies within the chord's ends along it, the stretch does too, and its distance
     *  from the chord's line is the polynomial sum of B_k w_k d_k over the weight's, sum B_k w_k,
     *  d_k the control points' distances from the line: the first at its most, as the cubic
     *  Bernstein bounds tell it, over the least w_k, which is no more than the weight's sum.
     *
     *  @param  piece The piece, of degree two or three, its weights above nothing.
     *  @param  from  Where the stretch starts.
     *  @param  until Where it ends.
     *  @return       The distance; infinite for a chord of no length.
     */
    auto Strayed(const Kernel::Piece &piece, double from, double until) -> double
    {
      const std::size_t     degree  = piece.Degree;
      std::array<Point, 4>  points  = {};
      std::array<double, 4> weights = {};
      for (std::size_t index = 0; index <= degree; ++index)
      {
        std::array<double, 3> args = {};
        for (std::size_t round = 0; round < degree; ++round)
        {
          args.at(round) = round + index < degree ? from : until;
        }
        const Homogeneous point = Homoblossom(piece, args);
        points.at(index)        = Projected(point);
        weights.at(index)       = point.W;
      }
      const Point  start  = points[0];
      const Vector chord  = points.at(degree) - start;
      const double square = Dot(chord, chord);
      if (!(square > 0.0))
      {
        return std::numeric_limits<double>::infinity();
      }
      const double          length = std::sqrt(square);
      double                hull   = 0.0;
      double                least  = weights[0];
      bool                  inside = true;
      std::array<double, 4> aside  = {};
      for (std::size_t index = 0; index <= degree; ++index)
      {
        const Vector away  = points.at(index) - start;
        const double along = Dot(away, chord) / square;
        const double off   = Cross(chord, away) / length;
        const double past  = std::max({0.0, -along, along - 1.0}) * length;
        hull               = std::max(hull, Trig::Hypot(off, past));
        inside             = inside && past <= 0.0;
        least              = std::min(least, weights.at(index));
        aside.at(index)    = weights.at(index) * off;
      }
      if (!inside)
      {
        return hull;
      }
      const std::array<double, 4> cubic =
          degree == 2 ? std::array{aside[0], (aside[0] + (2.0 * aside[1])) / 3.0,
                                   ((2.0 * aside[1]) + aside[2]) / 3.0, aside[2]}
                      : aside;
      const std::array<double, 2> bounds = Bounds(cubic);
      return std::min(hull, std::max(-bounds[0], bounds[1]) / least);
    }

    /**
     *  A rational spline of degree two or three as chords, span by span, greedy: from each place
     *  the chord tried to a span grown or shrunk by the square root of how far the last strayed, as
     *  a short chord strays as its span squared, kept where its stretch strays no farther than the
     *  tolerance less what rounding may add, the span carried into the next stretch.
     *
     *  Each span's piece is Kernel::Section's, its points and its stretches' control points within
     *  some (4 d + 2) (w_max / w_min) steps of a double of its poles' largest coordinate of the
     *  true ones: 64 (d + 1) such steps are kept from the tolerance. Its ends are At's.
     *
     *  @param  spline    The spline, rational.
     *  @param  tolerance How far it may stray, above nothing.
     *  @return           Its points; nothing for one not valid, of a degree other than two or
     *                    three, a tolerance rounding leaves no room within, or more tries than four
     *                    times the most chords. A chord of no length is shrunk, as one too long.
     */
    auto Rationed(const Spline &spline, double tolerance) -> std::optional<std::vector<Point>>
    {
      const auto degree = static_cast<std::size_t>(spline.Degree);
      if (degree < 2 || degree > 3 || !Valid(spline))
      {
        return std::nullopt;
      }
      const Range        domain = Domain(spline);
      std::vector<Point> points = {At(spline, domain.Start)};
      double             span   = 1.0;
      std::size_t        tries  = 0;
      for (std::size_t knot = degree; knot < spline.Poles.size(); ++knot)
      {
        if (!(spline.Knots[knot] < spline.Knots[knot + 1]))
        {
          continue;
        }
        double most  = 0.0;
        double light = std::numeric_limits<double>::infinity();
        double heavy = 0.0;
        for (std::size_t index = knot - degree; index <= knot; ++index)
        {
          const Point pole = spline.Poles[index];
          most             = std::max({most, std::abs(pole.X), std::abs(pole.Y)});
          light            = std::min(light, spline.Weights[index]);
          heavy            = std::max(heavy, spline.Weights[index]);
        }
        const double slack =
            64.0 * static_cast<double>(degree + 1) * (heavy / light) * 0x1p-53 * most;
        const double within = (tolerance - slack) * (1.0 - 0x1p-20);
        if (!(within > 0.0))
        {
          return std::nullopt;
        }
        const Kernel::Piece piece = Kernel::Section(spline, knot);
        double              from  = 0.0;
        while (from < 1.0)
        {
          if (++tries > 4 * MostChords)
          {
            return std::nullopt;
          }
          const double until = 1.0 - from <= span * (1.0 + 0x1p-20) ? 1.0 : from + span;
          const double away  = Strayed(piece, from, until);
          const double grown = 0.95 * std::sqrt(within / std::max(away, within / 16.0));
          if (away <= within)
          {
            points.push_back(Projected(Homoblossom(piece, {until, until, until})));
            span = (until - from) * std::min(2.0, grown);
            from = until;
            continue;
          }
          span = (until - from) * std::clamp(grown, 0.1, 0.9);
        }
      }
      points.back() = At(spline, domain.End);
      return points;
    }
  } // namespace

  auto Flatten(const Line &line, Accuracy /*accuracy*/) -> std::vector<Point>
  {
    return {line.Start, line.End};
  }

  auto Flatten(const Arc &arc, Accuracy accuracy) -> std::vector<Point>
  {
    return Wheeled(arc, Chords(Leeways(arc), accuracy.Error));
  }

  auto Flatten(const Ellipse &ellipse, Accuracy accuracy) -> std::vector<Point>
  {
    if (const std::optional<Sizing> sizing = Calibrated(ellipse.Map, accuracy.Error))
    {
      if (std::optional<std::vector<Point>> points = Quartered(*sizing, ellipse))
      {
        return *std::move(points);
      }
    }
    return Even(ellipse, MostChords);
  }

  auto Flatten(const Polyline &polyline, Accuracy accuracy) -> std::vector<Point>
  {
    if (polyline.Vertices.empty())
    {
      return {};
    }

    /**
     *  Each bulge sized first, so that the points are reserved once.
     */
    struct Sized
    {
      Bight       Bow    = {}; /**< The bulge.                    */
      std::size_t Chords = 0;  /**< Its chords; 1 for a straight. */
    };

    const std::span<const Vertex> vertices = polyline.Vertices;
    const std::size_t             count    = Count(polyline);
    List<Sized>                   sized(count);
    std::size_t                   total = 1;
    for (std::size_t index = 0; index < count; ++index)
    {
      const Step step = {.Start    = vertices[index].Position,
                         .End      = vertices[(index + 1) % vertices.size()].Position,
                         .Bulge    = vertices[index].Bulge,
                         .Fraction = 0.0};
      if (!Straightness(step))
      {
        const Bight bow = Bowed(step.Start, step.End, step.Bulge);
        sized[index]    = {.Bow = bow, .Chords = Chords(Leeways(bow), accuracy.Error)};
      }
      else
      {
        sized[index] = {.Bow = {}, .Chords = 1};
      }
      total += sized[index].Chords;
    }
    std::vector<Point> points;
    points.reserve(total);
    points.push_back(vertices.front().Position);
    for (std::size_t index = 0; index < count; ++index)
    {
      if (sized[index].Chords > 1)
      {
        Strided(points, sized[index].Bow, sized[index].Chords);
      }
      points.push_back(vertices[(index + 1) % vertices.size()].Position);
    }
    return points;
  }

  auto Flatten(const Spline &spline, Accuracy accuracy) -> std::vector<Point>
  {
    if (accuracy.Error > 0.0)
    {
      if (std::optional<std::vector<Point>> spanned = Rational(spline)
                                                        ? Rationed(spline, accuracy.Error)
                                                        : Chordwise(spline, accuracy.Error))
      {
        return *std::move(spanned);
      }
    }
    std::vector<Point> points;
    Kernel::Sections(spline,
                     [&points, accuracy](const Kernel::Piece &piece, Range /*span*/) -> void
                     {
                       if (points.empty())
                       {
                         points.push_back(Kernel::Start(piece));
                       }
                       Halve(points, piece, accuracy.Error);
                     });
    return points;
  }

  auto Flatten(const Bezier::Quadratic &quadratic, Accuracy accuracy) -> std::vector<Point>
  {
    const std::array points = {quadratic.Start, quadratic.Control, quadratic.End};
    if (const int power = LeveledPower(points); power != 0)
    {
      const Bezier::Quadratic scaled = {.Start   = Unshrunk(quadratic.Start, power),
                                        .Control = Unshrunk(quadratic.Control, power),
                                        .End     = Unshrunk(quadratic.End, power)};
      return Unscaled(Flatten(scaled, {.Error = std::ldexp(accuracy.Error, power)}), power);
    }
    return Chorded(quadratic, points, accuracy.Error);
  }

  auto Flatten(const Bezier::Cubic &cubic, Accuracy accuracy) -> std::vector<Point>
  {
    const std::array points = {cubic.Start, cubic.Out, cubic.In, cubic.End};
    if (const int power = LeveledPower(points); power != 0)
    {
      const Bezier::Cubic scaled = {.Start = Unshrunk(cubic.Start, power),
                                    .Out   = Unshrunk(cubic.Out, power),
                                    .In    = Unshrunk(cubic.In, power),
                                    .End   = Unshrunk(cubic.End, power)};
      return Unscaled(Flatten(scaled, {.Error = std::ldexp(accuracy.Error, power)}), power);
    }
    return Chorded(cubic, points, accuracy.Error);
  }

  auto Flatten(const Curve &curve, Accuracy accuracy) -> std::vector<Point>
  {
    return std::visit(
        [accuracy](const auto &shape) -> std::vector<Point>
        {
          using Kind = std::decay_t<decltype(shape)>;
          if constexpr (std::is_same_v<Kind, Ray> || std::is_same_v<Kind, Xline>)
          {
            return {};
          }
          else
          {
            return Flatten(shape, accuracy);
          }
        },
        curve);
  }
} // namespace Geoxy
