/***************************************************************************************************
 *  @file      source/swift/oval.cpp
 *  @brief     A line, a ray or an endless line met with an ellipse by its chord: Lagrange's
 *             identity makes the discriminant of the fourth degree, worked in double-double, and
 *             each crossing rounded once and proved the double nearest the true point.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "elliptic.hpp"
#include "real.hpp"
#include "sure.hpp"
#include "swift/lane.hpp"
#include "trig.hpp"

#include "geoxy/angle.hpp"
#include "geoxy/at.hpp"
#include "geoxy/ellipse.hpp"
#include "geoxy/matrix.hpp"
#include "geoxy/point.hpp"
#include "geoxy/tolerance.hpp"
#include "geoxy/vector.hpp"

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <numbers>
#include <optional>
#include <utility>

namespace Geoxy::Swift
{
  namespace
  {
    using Sure::Number;

    /**
     *  An ellipse's map as the kernel weighs it: its axes and center, the determinant D = X x Y in
     *  doubles with a bound on its rounding, and the least it stretches any length by, sigma_min >=
     *  |D| / |M|_F, as sigma_min sigma_max = |D| and sigma_max <= |M|_F.
     */
    struct Map
    {
      Vector X      = {};  /**< Where the unit x axis goes.            */
      Vector Y      = {};  /**< Where the unit y axis goes.            */
      Point  Origin = {};  /**< The center.                            */
      double Det    = 0.0; /**< X x Y, rounded.                        */
      double Slack  = 0.0; /**< How far Det may be from the exact one. */
      double Least  = 0.0; /**< sigma_min of the map, at the least.    */
    };

    /**
     *  An ellipse's map, as the kernel weighs it. D = X x Y rounded is within a step of each
     *  product and of their difference, under 2^-51 of the products' sizes; a map Inverse takes is
     *  far from flat, its determinant some 2^-40 of its axes' product at the least, so the least
     *  stretch is well above zero.
     *
     *  @param  ellipse The ellipse.
     *  @return         Its map.
     */
    auto Mapped(const Ellipse &ellipse) -> Map
    {
      const Vector axis  = ellipse.Map.X;
      const Vector other = ellipse.Map.Y;
      const double det   = Cross(axis, other);
      const double frob  = std::sqrt(Dot(axis, axis) + Dot(other, other)) * (1.0 + 0x1p-50);
      const double slack = 0x1p-51 * (std::abs(axis.X * other.Y) + std::abs(axis.Y * other.X));
      return {.X      = axis,
              .Y      = other,
              .Origin = ellipse.Map.Origin,
              .Det    = det,
              .Slack  = slack,
              .Least  = (std::abs(det) - slack) / frob * Sure::Shrink};
    }

    /**
     *  How far a point is from an ellipse's carrier, at the least: its distance from the unit
     *  circle in the circle's frame, | |adj(M) (p - o)| - |D| | / |D|, times the least stretch,
     *  each part less what its rounding may take off.
     *
     *  @param  map   The ellipse's map.
     *  @param  point The point.
     *  @return       The distance, at the least.
     */
    auto Offness(const Map &map, Point point) -> double
    {
      const Vector gap    = point - map.Origin;
      const double across = (gap.X * map.Y.Y) - (gap.Y * map.Y.X);
      const double upward = (map.X.X * gap.Y) - (map.X.Y * gap.X);
      const double slack  = 0x1p-50
                          * (std::abs(gap.X * map.Y.Y)
                             + std::abs(gap.Y * map.Y.X)
                             + std::abs(map.X.X * gap.Y)
                             + std::abs(map.X.Y * gap.X));
      const double size   = std::abs(map.Det);
      const double off    = std::abs(Trig::Hypot(across, upward) - size) * (1.0 - 0x1p-50);
      return map.Least * (off - (2.0 * slack) - map.Slack) / ((size + map.Slack) * (1.0 + 0x1p-50));
    }

    /**
     *  How a carrier passes an ellipse, as plain doubles tell it.
     */
    enum class Passing : std::uint8_t
    {
      Misses, /**< It misses the whole ellipse by more than the reach. */
      Cuts,   /**< It cuts into it, certainly.                         */
      Near,   /**< Neither is certain.                                 */
    };

    /**
     *  How a lane's carrier passes an ellipse, in doubles, as most pairs miss by far. With e =
     *  adj(M) d and w = g x d, g = P - o, the carrier lies |w| / |e| from the center in the
     *  circle's frame. d is rounded once at the most, g once, and each product and sum once more:
     *  e's coordinates and w are each within 2^-50 of their terms' sizes, and A = |e|^2 within
     *  2^-50 of itself and twice e's slack. Missing the whole ellipse, the carrier lies |w| /
     *  sqrt(A) - 1 off it in the circle's frame, sigma_min times that at the least.
     *
     *  @param  lane  The lane.
     *  @param  map   The ellipse's map.
     *  @param  reach How far it must miss by.
     *  @return       How it passes; near where A is all but nothing, as for a tiny lane.
     */
    auto Passed(const Lane &lane, const Map &map, double reach) -> Passing
    {
      const Vector way    = lane.Rough;
      const Vector gap    = lane.Start - map.Origin;
      const double aheadX = (way.X * map.Y.Y) - (way.Y * map.Y.X);
      const double aheadY = (map.X.X * way.Y) - (map.X.Y * way.X);
      const double slackX = 0x1p-50 * (std::abs(way.X * map.Y.Y) + std::abs(way.Y * map.Y.X));
      const double slackY = 0x1p-50 * (std::abs(map.X.X * way.Y) + std::abs(map.X.Y * way.X));
      const double square = (aheadX * aheadX) + (aheadY * aheadY);
      const double spread = (0x1p-50 * square)
                          + (2.0 * ((std::abs(aheadX) * slackX) + (std::abs(aheadY) * slackY)))
                          + (slackX * slackX)
                          + (slackY * slackY);
      const double cross  = Cross(gap, way);
      const double wobble = 0x1p-50 * (std::abs(gap.X * way.Y) + std::abs(gap.Y * way.X));
      const double low    = (square - spread) * Sure::Shrink;
      const double high   = (square + spread) * Sure::Grow;
      const double least  = (std::abs(cross) - wobble) * Sure::Shrink;
      const double most   = (std::abs(cross) + wobble) * Sure::Grow;
      if (!(low > Narrowest * Narrowest))
      {
        return Passing::Near;
      }
      if (least > 0.0 && least * least > high * Sure::Grow)
      {
        const double root   = std::sqrt(high) * Sure::Grow;
        const double missed = map.Least * ((least - root) / root) * Sure::Shrink;
        return missed > reach ? Passing::Misses : Passing::Near;
      }
      return most * most < low * Sure::Shrink ? Passing::Cuts : Passing::Near;
    }

    /**
     *  What a carrier's two crossings with an ellipse share, in double-double: c = adj(M) g, e =
     *  adj(M) d, A = e.e, B = c.e and w = g x d. Lagrange's identity, (c.e)^2 - |c|^2 |e|^2 = -(c x
     *  e)^2 with c x e = D w, makes B^2 - A C = D^2 (A - w^2): the crossings are t = (-B -+ |D|
     *  sqrt(A - w^2)) / A, the middle -B / A and the half |D| sqrt(A - w^2) / A. The shallower cap
     *  the carrier cuts off is 1 - |w| / sqrt(A) = (A - w^2) / (sqrt(A) (sqrt(A) + |w|)) high in
     *  the circle's frame, sigma_min times that at the least.
     */
    struct Chord
    {
      std::array<Number, 2> Offset = {};  /**< c, the lane's start in the circle's frame
                                               times D.                                  */
      std::array<Number, 2> Ahead  = {};  /**< e, its way so.                        */
      Number                Middle = {};  /**< -B / A.                               */
      Number                Half   = {};  /**< |D| sqrt(A - w^2) / A.                */
      double                Cap    = 0.0; /**< The shallower cap's height, at least. */
      double                Facing = 1.0; /**< D's sign.                             */
    };

    /**
     *  What a lane's two crossings with an ellipse share.
     *
     *  @param  lane The lane, its carrier cutting into the ellipse.
     *  @param  map  The ellipse's map.
     *  @return      The chord.
     */
    auto Chorded(const Lane &lane, const Map &map) -> Chord
    {
      const Number dx      = lane.Way[0];
      const Number dy      = lane.Way[1];
      const Number gapX    = Sure::Apart(lane.Start.X, map.Origin.X);
      const Number gapY    = Sure::Apart(lane.Start.Y, map.Origin.Y);
      const Number axisX   = Sure::Exact(map.X.X);
      const Number axisY   = Sure::Exact(map.X.Y);
      const Number otherX  = Sure::Exact(map.Y.X);
      const Number otherY  = Sure::Exact(map.Y.Y);
      const Number aheadX  = (dx * otherY) - (dy * otherX);
      const Number aheadY  = (dy * axisX) - (dx * axisY);
      const Number offsetX = (gapX * otherY) - (gapY * otherX);
      const Number offsetY = (gapY * axisX) - (gapX * axisY);
      const Number square  = (aheadX * aheadX) + (aheadY * aheadY);
      const Number along   = (offsetX * aheadX) + (offsetY * aheadY);
      const Number cross   = (gapX * dy) - (gapY * dx);
      const Number inside  = square - (cross * cross);
      const double facing  = map.Det > 0.0 ? 1.0 : -1.0;
      const Number size =
          (Sure::Product(map.X.X, map.Y.Y) - Sure::Product(map.X.Y, map.Y.X)) * Sure::Exact(facing);
      const double top = std::sqrt(Sure::Most(square)) * Sure::Grow;
      return {
          .Offset = {offsetX, offsetY},
          .Ahead  = {aheadX, aheadY},
          .Middle = Sure::Quotient(-along, square),
          .Half   = Sure::Quotient(size * Sure::Root(inside), square),
          .Cap = map.Least * Sure::Least(inside) / (top * (top + Sure::Most(cross))) * Sure::Shrink,
          .Facing = facing};
    }

    /**
     *  A lane and an ellipse being met, and what is worked out once about them.
     */
    struct OvalScene
    {
      const Lane      *Straight  = nullptr; /**< The lane.                     */
      const Ellipse   *Shape     = nullptr; /**< The ellipse.                  */
      Map              Frame     = {};      /**< Its map.                      */
      Chord            Cut       = {};      /**< Their chord.                  */
      Geoxy::Tolerance Tolerance = {};      /**< The tolerance.                */
      bool             Whole     = false;   /**< Whether the ellipse is whole. */
    };

    /**
     *  Whether a crossing's angle on the unit circle lies within an arc's sweep: by the signs of
     *  its cross products with the sweep's ends, each decided where it is clear of what rounding
     *  could make of it.
     *
     *  @param  ellipse The ellipse, swept no more than once round.
     *  @param  seen    The crossing's point on the unit circle, times D's size.
     *  @param  slack   How far that may be from the exact one.
     *  @return         Inside or out; near where a sign is not clear.
     */
    auto Arcwise(const Ellipse &ellipse, Vector seen, double slack) -> Side
    {
      const double sweep   = ellipse.Sweep.InRadians();
      const bool   forward = sweep > 0.0;
      const Vector first   = Direction(ellipse.Start);
      const Vector last    = Direction(ellipse.Start + ellipse.Sweep);
      const Vector start   = forward ? first : last;
      const Vector end     = forward ? last : first;
      const double clear   = (0x1p-48 * Length(seen)) + (2.0 * slack);
      const double enter   = Cross(start, seen);
      const double leave   = Cross(seen, end);
      if (!(std::abs(enter) > clear) || !(std::abs(leave) > clear))
      {
        return Side::Near;
      }
      const bool inside = std::abs(sweep) <= std::numbers::pi ? enter > 0.0 && leave > 0.0
                                                              : enter > 0.0 || leave > 0.0;
      return inside ? Side::In : Side::Out;
    }

    /**
     *  The sine of the angle between a line and an ellipse's tangent where it crosses.
     *
     *  @param  map  The ellipse's map.
     *  @param  seen The crossing's point on the unit circle, times D.
     *  @param  way  The line's way.
     *  @return      The sine.
     */
    auto Sine(const Map &map, Vector seen, Vector way) -> double
    {
      const Vector tangent = (map.Y * seen.X) - (map.X * seen.Y);
      return std::abs(Cross(way, tangent)) / (Length(way) * Length(tangent));
    }

    /**
     *  One of a lane's crossings with an ellipse: where it lies on the lane, then on the arc, by
     *  its point on the unit circle c + t e; kept where it lies on both, steep, the cap deeper than
     *  the reach, and its point settled.
     *
     *  @param  scene The lane and the ellipse.
     *  @param  sign  Which root: -1 or 1.
     *  @param  spots Where it goes.
     *  @return       Whether it was decided.
     */
    auto Rooted(const OvalScene &scene, double sign, Spots &spots) -> bool
    {
      const Lane  &lane   = *scene.Straight;
      const Chord &chord  = scene.Cut;
      const double reach  = Reach(scene.Tolerance);
      const Number along  = sign < 0.0 ? chord.Middle - chord.Half : chord.Middle + chord.Half;
      const double margin = (reach / (lane.Speed * Sure::Shrink))
                          + (along.Bound * Sure::Grow)
                          + (0x1p-50 * std::abs(along.Value.High));
      const Side   place  = Within(lane.Along, along.Value.High, margin);
      if (place != Side::In)
      {
        return place == Side::Out;
      }
      const Number across = chord.Offset[0] + (along * chord.Ahead[0]);
      const Number upward = chord.Offset[1] + (along * chord.Ahead[1]);
      const Vector seen   = {.X = chord.Facing * across.Value.High,
                             .Y = chord.Facing * upward.Value.High};
      const double off =
          (across.Bound + std::abs(across.Value.Low) + upward.Bound + std::abs(upward.Value.Low))
          * Sure::Grow;
      const Side arc = scene.Whole ? Side::In : Arcwise(*scene.Shape, seen, off);
      if (arc != Side::In)
      {
        return arc == Side::Out;
      }
      const auto point =
          Sure::Rounded(std::array{Sure::Exact(lane.Start.X) + (lane.Way[0] * along),
                                   Sure::Exact(lane.Start.Y) + (lane.Way[1] * along)});
      const bool steep = Steep(Sine(scene.Frame, seen, lane.Rough), scene.Tolerance);
      if (!(chord.Cap > reach) || !steep || !point)
      {
        return false;
      }
      spots.push_back({.At    = {.X = (*point)[0], .Y = (*point)[1]},
                       .Along = {along.Value.High, 0.0},
                       .Own   = false});
      return true;
    }

    /**
     *  Whether every end of either curve lies far from the other, a whole ellipse's seam as an end:
     *  so that no end is a hit, and no crossing lies near one.
     *
     *  @param  scene The lane and the ellipse.
     *  @return       Whether they do.
     */
    auto Ended(const OvalScene &scene) -> bool
    {
      const Lane  &lane  = *scene.Straight;
      const double reach = Reach(scene.Tolerance);
      for (std::size_t index = 0; index < lane.Count; ++index)
      {
        if (!(Offness(scene.Frame, lane.Ends.at(index)) > reach))
        {
          return false;
        }
      }
      const bool start = Away(lane, At(*scene.Shape, 0.0)) > reach;
      const bool end   = scene.Whole || Away(lane, At(*scene.Shape, 1.0)) > reach;
      return start && end;
    }
  } // namespace

  auto Oval(const Lane &lane, const Ellipse &ellipse, const Tolerance &tolerance)
      -> std::optional<Spots>
  {
    const Vector axis  = ellipse.Map.X;
    const Vector other = ellipse.Map.Y;
    const Point  home  = ellipse.Map.Origin;
    const double sweep = ellipse.Sweep.InRadians();
    if (!Held({axis.X, axis.Y, other.X, other.Y, home.X, home.Y})
        || !(std::abs(Cross(axis, other)) > Narrowest)
        || Zero(sweep)
        || Coiled(ellipse))
    {
      return std::nullopt;
    }
    const Map     map     = Mapped(ellipse);
    const Passing passing = Passed(lane, map, Reach(tolerance));
    if (passing != Passing::Cuts)
    {
      return passing == Passing::Misses ? std::optional<Spots>{std::in_place} : std::nullopt;
    }
    const OvalScene scene = {.Straight  = &lane,
                             .Shape     = &ellipse,
                             .Frame     = map,
                             .Cut       = Chorded(lane, map),
                             .Tolerance = tolerance,
                             .Whole     = !(std::abs(sweep) < Turn.InRadians())};
    Spots           spots;
    for (const double sign : {-1.0, 1.0})
    {
      if (!Rooted(scene, sign, spots))
      {
        return std::nullopt;
      }
    }
    if (!Ended(scene))
    {
      return std::nullopt;
    }
    return spots;
  }
} // namespace Geoxy::Swift
