/***************************************************************************************************
 *  @file      source/abut.hpp
 *  @brief     Two ready curves where an end of one lies on the other, a corner, a joint, a door's
 *             swing on its wall, told from what they keep as Intersect tells them, to the bit.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#ifndef GEOXY_ABUT_HPP
#define GEOXY_ABUT_HPP

#include "met/linear.hpp"
#include "real.hpp"
#include "thorough.hpp"

#include "geoxy/intersect.hpp"
#include "geoxy/line.hpp"
#include "geoxy/point.hpp"
#include "geoxy/prepared.hpp"
#include "geoxy/tolerance.hpp"
#include "geoxy/vector.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <optional>

namespace Geoxy
{
  /**
   *  What two ready curves meet in where an end of one lies on the other, told: two marks at the
   *  most, in order along the first.
   */
  struct Abutment
  {
    std::array<Mark, 2> Marks     = {}; /**< The marks.                           */
    std::size_t         Count     = 0;  /**< How many there are.                  */
    Overlap             Stretch   = {}; /**< A stretch they share, where they do. */
    std::size_t         Stretches = 0;  /**< How many there are: none, or one.    */
  };

  /**
   *  A line as Jointed reads it, worked out once for each joint it is in: its way, the largest
   *  magnitude of its way's coordinates and of its ends', and whether its coordinates sum to a
   *  finite number; and for Seamless, sixteen steps of a double at its ends' largest magnitude, its
   *  start moved along its way, and whether half its end is exact.
   */
  struct Limb
  {
    const Line *Shape   = nullptr; /**< The line.                                       */
    Vector      Way     = {};      /**< Its end less its start, rounded.                */
    double      Largest = 0.0;     /**< The largest magnitude of its way's coordinates. */
    double      Reach   = 0.0;     /**< The largest magnitude of its ends' coordinates. */
    double      Grain   = 0.0;     /**< Sixteen steps of a double at its reach, or an
                                        infinity where its coordinates' sum is not
                                        finite.                                        */
    Point Beyond = {};             /**< Its start moved along its way, rounded.      */
    bool  Finite = false;          /**< Whether its four coordinates' sum is finite. */
    bool  Halved = false;          /**< Whether each coordinate of its end, halved and
                                        doubled, is itself.                            */
  };

  /**
   *  Two ready curves met where an end of one lies on the other, as Intersect meets them, from what
   *  they keep: the pairs a drawing's corners, joints and door swings make, which no pass that
   *  tells crossings far from every end can tell, and which Intersect's direct ways tell after a
   *  long way round.
   *
   *  Each test here is on exact differences, or on numbers clear of what Intersect compares them
   *  with by more than its roundings can close; a pair any test is not clear for is left to
   *  Intersect.
   */
  struct Abutting
  {
    /**
     *  A level line and an upright one, either first: nothing, one crossing, or one hit at an end,
     *  all at their carriers' crossing (x, y), x the upright's and y the level's, which is a
     *  double.
     *
     *  Each end's gap from the other's carrier is its own coordinate less that carrier's, exact in
     *  its sign and in whether it is zero, and rounded once. A line whose ends both lie past four
     *  tolerances of the other's carrier, to one side, is beside it, and they meet in nothing. Else
     *  every end must lie on the other's carrier exactly, or past four tolerances from it: then
     *  each carrier's crossing lies on the other line, and the ends on it make it a hit at an end,
     *  the rest leaving it a crossing. The proof that Intersect's direct way answers so, its every
     *  rounding taken, is with the definition.
     *
     *  @param  lines  The two ready lines.
     *  @param  length The length they are met within, as Reached gives it.
     *  @return        What they meet in; nothing where Intersect must decide: an end near but not
     *                 on the other's carrier, a line no longer than four tolerances, a tolerance
     *                 finer than 2^-47 of the lines' numbers, or a crossing on an axis, whose zero
     *                 Intersect gives a sign to.
     */
    static auto Squarely(const Prepared::Pair &lines, double length) -> std::optional<Abutment>;

    /**
     *  A level or upright line and an arc of a radius, either first, where an end of either lies at
     *  a crossing of the line's carrier with the arc's circle: a door's swing ending on its wall, a
     *  wall ending on an arc.
     *
     *  The carrier crosses the circle at c +- sqrt(r^2 - g^2) along the line, c the center's
     *  coordinate there and g the gap from the center to the carrier, one exact root, each crossing
     *  certified the double nearest it as Levelled certifies it, one on an axis by Axis. A line at
     *  +0 crosses at +0, as every way of Intersect's gives it. Every end of either curve must lie
     *  within half the tolerance of a crossing, or farther than four tolerances from the other
     *  curve, and each crossing on the line or off it, on the arc or off it, each by more than
     *  Intersect's roundings can close; where not, Intersect decides. A crossing on both is a hit,
     *  at an end where an end near a crossing lies within the tolerance of it, and placed on both
     *  curves as Intersect places it.
     *
     *  @param  curves The ready line and arc.
     *  @param  length The length they are met within, as Reached gives it.
     *  @return        What they meet in; nothing where Intersect must decide: a slanted line, an
     *                 arc whole or no wider than 1024 tolerances, a carrier that cuts the circle at
     *                 a sine of 1 / sqrt(32) or less, a tolerance finer than 2^-42 of the curves'
     *                 numbers, or a line on an axis at -0.
     */
    static auto Rimmed(const Prepared::Pair &curves, double length) -> std::optional<Abutment>;

    /**
     *  Two ready lines at any angle where an end of one lies near the other, or where they cross
     *  near an end, as Intersect's direct way meets them: a corner, a joint, two sides of a polygon
     *  end to end, a crossing near an end that no plain pass tells.
     *
     *  Every test Intersect's Straightly makes on its rough numbers is made here on the same
     *  numbers, worked out by the same steps, so that it comes out the same, or on squares that
     *  tell it where they lie clear of its threshold by more than the roundings of either side can
     *  close. Where an end lies near the other line, Endly's hit at the carriers' crossing is
     *  Intersect's own quick crossing, placed on both lines as Intersect places it, or where that
     *  does not settle its exact way's, a coordinate on an axis told by Axis; where Endly gives the
     *  pair up, as at a narrow angle, Butt's at two ends that meet is halfway between them. A
     *  crossing farther than four tolerances from every end is found as Endly's is.
     *
     *  @param  lines  The two ready lines.
     *  @param  usable The tolerance they are met within: the length as Reached gives it, within
     *                 Coarse and Fine, and the angle as Sane makes it.
     *  @return        What they meet in; nothing where Intersect must decide: a stretch they share,
     *                 an angle within the tolerance's, a crossing whose coordinate neither settles
     *                 nor lies on an axis, or an end Intersect's ways leave to its thorough one.
     */
    static auto Slanted(const Prepared::Pair &lines, const Tolerance &usable)
        -> std::optional<Abutment>;

    /**
     *  Whether two lines joined end to end, the first's end the very point of the second's start or
     *  its start of the second's end, as a loop's sides are, meet in that point alone as Intersect
     *  meets them made ready: one hit there, which a cut of either line at its own end would make
     *  nothing of.
     *
     *  Intersect's passes before Abutting's tell such lines nothing, and Met tells them by Squarely
     *  or by Slanted, which ask here only what the joint does not settle exactly: the steps that
     *  decide are Slanted's own, on the same numbers, each asked first of a bound it is never below
     *  where that settles it.
     *
     *  @param  lines The two lines, of finite coordinates.
     *  @param  sane  The tolerance, as Sane makes it.
     *  @return       Whether they do; false where they are not so joined, or Intersect may answer
     *                anything else: a stretch along the other near the joint, an end near the other
     *                line, or a hit not at the joint's very point.
     */
    static auto Jointed(const std::array<const Line *, 2> &lines, const Tolerance &sane) -> bool;

    /**
     *  Whether two lines joined end to end meet in that point alone as Intersect meets them made
     *  ready, as the overload of two lines tells it, from what each keeps as a limb: the numbers
     *  every joint of a line asks of it, worked out once.
     *
     *  @param  limbs The two lines as limbs, as Limbed makes them.
     *  @param  sane  The tolerance, as Sane makes it.
     *  @return       Whether they do, as the overload of two lines tells it.
     */
    static auto Jointed(const std::array<Limb, 2> &limbs, const Tolerance &sane) -> bool;

    /**
     *  A line as Jointed reads it.
     *
     *  @param  line The line.
     *  @return      Its limb, pointing at it.
     */
    static auto Limbed(const Line &line) -> Limb;

    /**
     *  Where two ready lines' carriers cross, as Intersect's Exactly finds the point: its quick
     *  crossing where that settles both coordinates, by the same steps; else its exact way's, a
     *  coordinate the line's own where a line holds it, the quick one's where that settles it, and
     *  zero where the crossing lies on the axis, of the ways' cross product's sign.
     *
     *  @param  lines The two ready lines, crossing at a sine past Steep.
     *  @return       The point; nothing where a coordinate is neither of those.
     */
    static auto Crossing(const Prepared::Pair &lines) -> std::optional<Point>;

    /**
     *  Where two ready lines' carriers cross, as the overload of a pair tells it, the lines given
     *  one by one, so that a caller holding no pair of them makes none.
     *
     *  @param  one   The first ready line.
     *  @param  other The other, crossing it at a sine past Steep.
     *  @return       The point, as the overload of a pair gives it.
     */
    static auto Crossing(const Prepared &one, const Prepared &other) -> std::optional<Point>;

    /**
     *  Two ready curves where an end of one lies on the other, as Intersect meets them, each told
     *  by its pass: a level line and an upright one by Squarely, else two lines by Slanted, a level
     *  or upright line and an arc by Rimmed.
     *
     *  @param  curves The two ready curves.
     *  @param  usable The tolerance they are met within: the length as Reached gives it, within
     *                 Coarse and Fine, and the angle as Sane makes it.
     *  @return        What they meet in; nothing for any other pair, or where Intersect must
     *                 decide.
     */
    static auto Met(const Prepared::Pair &curves, const Tolerance &usable)
        -> std::optional<Abutment>;

    /**
     *  Two ready lines, one of no length, a point, and the other many tolerances long, as Intersect
     *  meets them: its thorough way's, as no direct way tells them. The point is a hit at an end
     *  where it lies within the tolerance the thorough way holds the lines to, Usable's, of its
     *  nearest point on the other line, halfway between the two; nothing where it lies farther.
     *  Where it lies within a few steps of a double of that tolerance, Intersect decides.
     *
     *  @param  lines The two ready lines.
     *  @param  sane  The tolerance, as Sane makes it.
     *  @return       What they meet in; nothing for any other pair, or where Intersect must decide.
     */
    static auto Dotted(const Prepared::Pair &lines, const Tolerance &sane)
        -> std::optional<Abutment>;
  };

  inline auto Abutting::Limbed(const Line &line) -> Limb
  {
    const Vector way    = line.End - line.Start;
    const double reach  = std::max({std::abs(line.Start.X), std::abs(line.Start.Y),
                                    std::abs(line.End.X), std::abs(line.End.Y)});
    const bool   finite = std::isfinite(line.Start.X + line.Start.Y + line.End.X + line.End.Y);
    return {.Shape   = &line,
            .Way     = way,
            .Largest = Largest(way),
            .Reach   = reach,
            .Grain   = finite ? 16.0 * Above(reach) : std::numeric_limits<double>::infinity(),
            .Beyond  = line.Start + way,
            .Finite  = finite,
            .Halved  = Exactly((line.End.X * 0.5) + (line.End.X * 0.5), line.End.X)
                    && Exactly((line.End.Y * 0.5) + (line.End.Y * 0.5), line.End.Y)};
  }

  /**
   *  Whether two lines joined end to end, the first's end the very point of the second's start, are
   *  told by Jointed the way a polygon's sides most often are, every test that way passing, and so
   *  meet in their joint alone: worked out on the numbers each limb keeps, inline, so that a run of
   *  a loop's joints goes through at a few nanoseconds each. Where it says no, Jointed may yet say
   *  yes.
   *
   *  It asks what Jointed and Propped ask, by the same steps on the same numbers, where they answer
   *  yes:
   *  - the joint as Selfsame tells it, the other two ends apart;
   *  - Joinable's length, the largest of the two limbs' grains being sixteen steps at the largest
   *    of the eight coordinates, as a step of a double grows with its magnitude; or an infinity
   *    where Joinable finds a line not plain, which no length within Coarse is;
   *  - the ways not below half the length at their largest, which Asunder takes first of the joint
   *    and the near end, the first line's way being their difference, and of the joint and the far
   *    end, whose difference the second's way is in magnitude, the joint being its start but for
   *    the sign of a zero, which no difference with it keeps past its magnitude;
   *  - the ways' dot product not below zero, so that Aligned is not asked;
   *  - Endly's sine not past Steep, by the bound Jointed takes first or by the lengths it takes
   *    then, so that Propped is asked;
   *  - Asunder of the near and far ends;
   *  - each other end's gap from the box of the other line's run, which Remote asks first, the
   *    run's last point being its limb's beyond;
   *  - and the halfway point of the two lines' ends at the joint as Propped works it out: those
   *    ends are the joint itself but for the sign of a zero, each way being finite within Coarse,
   *    so their halves are one and their sum the joint where half of it is exact, as the first's
   *    limb tells.
   *
   *  Squarer, asked first by Jointed, answers yes or leaves the rest to answer, and the rest
   *  answers yes here.
   *
   *  @param  one   The first line as a limb, its end at the second's start.
   *  @param  other The second.
   *  @param  sane  The tolerance, as Sane makes it.
   *  @return       Whether every test passes; then Jointed tells them.
   */
  inline auto Seamless(const Limb &one, const Limb &other, const Tolerance &sane) -> bool
  {
    const Line  &first  = *one.Shape;
    const Line  &second = *other.Shape;
    const double length = std::max({sane.Length, one.Grain, other.Grain});
    const double half   = length / 2;
    const double clear  = 4.0 * length;
    const double least  = one.Largest * other.Largest;
    const double cross  = Cross(one.Way, other.Way);
    const auto   gap    = [](double value, double from, double until) -> double
    {
      const double low  = std::min(from, until);
      const double high = std::max(from, until);
      if (value < low)
      {
        return low - value;
      }
      return value > high ? value - high : 0.0;
    };
    const bool joined =
        Exactly(first.End.X, second.Start.X)
        && Exactly(first.End.Y, second.Start.Y)
        && !(Exactly(first.Start.X, second.End.X) && Exactly(first.Start.Y, second.End.Y));
    const bool fit   = !(length > Coarse) && !(length < Fine) && sane.Radians <= Met::Steep / 2.0;
    const bool apart = one.Largest > half
                    && other.Largest > half
                    && Dot(one.Way, other.Way) >= 0.0
                    && Largest(second.End - first.Start) > half
                    && std::max(gap(first.Start.X, second.Start.X, other.Beyond.X),
                                gap(first.Start.Y, second.Start.Y, other.Beyond.Y))
                           > clear
                    && std::max(gap(second.End.X, first.Start.X, one.Beyond.X),
                                gap(second.End.Y, first.Start.Y, one.Beyond.Y))
                           > clear
                    && one.Halved;
    if (!joined || !fit || !apart)
    {
      return false;
    }
    /**
     *  Each way's largest coordinate is past half the length, no finer than Fine, so their product
     *  is past 2^-1000, as is the product of their lengths, which no way within Coarse takes past
     *  every double: Jointed's tests of both pass.
     */
    if (std::abs(cross) <= Met::Steep * least)
    {
      return true;
    }
    return !(std::abs(cross / (Length(one.Way) * Length(other.Way))) > Met::Steep);
  }
} // namespace Geoxy

#endif /* GEOXY_ABUT_HPP */
