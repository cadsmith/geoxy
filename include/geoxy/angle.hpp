/***************************************************************************************************
 *  @file      geoxy/angle.hpp
 *  @brief     Angles, held in the unit they were made in, so that the unit is never a guess and an
 *             angle given in degrees is worked exactly.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#ifndef GEOXY_ANGLE_HPP
#define GEOXY_ANGLE_HPP

#include "geoxy/api.hpp"
#include "geoxy/vector.hpp"

#include <compare>
#include <functional>
#include <numbers>

namespace Geoxy
{
  /**
   *  The double nearest pi.
   */
  inline constexpr double Pi = std::numbers::pi;

  /**
   *  An angle, held in the unit it was made in: degrees by Degrees, radians by Radians.
   *
   *  DXF keeps an ARC's angles in degrees and an ELLIPSE's parameters in radians, and each comes
   *  back exactly as given. In degrees a turn, its quarters and eighths and every whole degree are
   *  doubles, so an angle held in degrees wraps and sweeps exactly: the reductions are proved in
   *  tests/proof/angle.v. In radians none of them is a double, so nothing reduces exactly, but its
   *  directions are reduced to 2^-201 and come out as correctly rounded as in degrees. It may run
   *  past a turn, as a sweep does.
   */
  class Angle
  {
  public:
    /**
     *  No angle, in radians.
     */
    constexpr Angle() = default;

    /**
     *  The angle in radians.
     *
     *  @return The radians: as held, or the degrees converted to the nearest double, so 90, 180 and
     *          360 degrees give pi / 2, pi and 2 pi.
     */
    [[nodiscard]] constexpr auto InRadians() const -> double
    {
      return m_InDegrees ? Converted(m_Value, true) : m_Value;
    }

    /**
     *  The angle in degrees.
     *
     *  @return The degrees: as held, or the radians converted to the nearest double, so pi / 2, pi
     *          and 2 pi give 90, 180 and 360 degrees.
     */
    [[nodiscard]] constexpr auto InDegrees() const -> double
    {
      return m_InDegrees ? m_Value : Converted(m_Value, false);
    }

    /**
     *  Whether the angle is held in degrees, where it is worked exactly.
     *
     *  @return Whether it was made by Degrees, or from angles that were.
     */
    [[nodiscard]] constexpr auto HeldInDegrees() const -> bool { return m_InDegrees; }

    /**
     *  Whether two angles are equal: held in one unit, exactly; in two, in degrees.
     *
     *  Mixed, the radians are rounded into degrees first (InDegrees, to the nearest double) and
     *  compared there: Radians(Pi) == Degrees(180), though pi's double is not 180 degrees.
     *
     *  @param  lhs One angle.
     *  @param  rhs The other.
     *  @return     Whether they are equal.
     */
    [[nodiscard]] friend constexpr auto operator==(const Angle &lhs, const Angle &rhs) -> bool
    {
      if (lhs.m_InDegrees == rhs.m_InDegrees)
      {
        return std::equal_to<double>{}(lhs.m_Value, rhs.m_Value);
      }
      return std::equal_to<double>{}(lhs.InDegrees(), rhs.InDegrees());
    }

    /**
     *  The order of two angles: held in one unit, exactly; in two, in degrees.
     *
     *  Mixed, the radians are rounded into degrees first, as for ==, so two angles a rounding apart
     *  can compare equal.
     *
     *  @param  lhs One angle.
     *  @param  rhs The other.
     *  @return     The order, unordered when either is NaN.
     */
    [[nodiscard]] friend constexpr auto operator<=>(const Angle &lhs, const Angle &rhs)
        -> std::partial_ordering
    {
      if (lhs.m_InDegrees == rhs.m_InDegrees)
      {
        return lhs.m_Value <=> rhs.m_Value;
      }
      return lhs.InDegrees() <=> rhs.InDegrees();
    }

    /**
     *  The angle of so many radians.
     *
     *  @param  radians The radians, as given.
     *  @return         The angle.
     */
    friend constexpr auto Radians(double radians) -> Angle;

    /**
     *  The angle of so many degrees.
     *
     *  @param  degrees The degrees, as given.
     *  @return         The angle.
     */
    friend constexpr auto Degrees(double degrees) -> Angle;

  private:
    /**
     *  The angle of a value in a unit, made only through Radians and Degrees.
     *
     *  @param value     The value.
     *  @param inDegrees Whether it is degrees.
     */
    constexpr Angle(double value, bool inDegrees) : m_Value(value), m_InDegrees(inDegrees) {}

    /**
     *  Degrees to radians or back: the product with a constant held in two doubles, and where that
     *  leaves the rounding in doubt, in three, its rounding then decided by an exact sign.
     *
     *  @param  value     The value.
     *  @param  toRadians Whether it is degrees to radians.
     *  @return           The value in the other unit, the nearest double, always:
     *                    tests/proof/convert.py.
     */
    [[nodiscard]] GEOXY_API static auto Converted(double value, bool toRadians) -> double;

    double m_Value     = 0.0;   /**< The value, as given.          */
    bool   m_InDegrees = false; /**< Whether the value is degrees. */
  };

  /**
   *  The angle of so many radians.
   *
   *  @param  radians The radians, as given, past a turn too.
   *  @return         The angle, held in radians.
   */
  [[nodiscard]] constexpr auto Radians(double radians) -> Angle
  {
    return {radians, false};
  }

  /**
   *  The angle of so many degrees.
   *
   *  @param  degrees The degrees, as given, past a turn too.
   *  @return         The angle, held in degrees exactly.
   */
  [[nodiscard]] constexpr auto Degrees(double degrees) -> Angle
  {
    return {degrees, true};
  }

  /**
   *  A whole turn, 360 degrees: 2 pi in radians.
   */
  inline constexpr Angle Turn = Degrees(360.0);

  /**
   *  The sum of two angles: in their unit when they share one, else in degrees.
   *
   *  Mixed, the radians are rounded into degrees first (InDegrees, to the nearest double), so the
   *  sum rounds twice; Wrapped and Within on it are then exact on those degrees.
   *
   *  @param  lhs One angle.
   *  @param  rhs The other.
   *  @return     Their sum, rounded once.
   */
  [[nodiscard]] constexpr auto operator+(Angle lhs, Angle rhs) -> Angle
  {
    if (lhs.HeldInDegrees() || rhs.HeldInDegrees())
    {
      return Degrees(lhs.InDegrees() + rhs.InDegrees());
    }
    return Radians(lhs.InRadians() + rhs.InRadians());
  }

  /**
   *  The difference of two angles: in their unit when they share one, else in degrees.
   *
   *  Mixed, the radians are rounded into degrees first, as for +.
   *
   *  @param  lhs The angle subtracted from.
   *  @param  rhs The angle subtracted.
   *  @return     Their difference, rounded once.
   */
  [[nodiscard]] constexpr auto operator-(Angle lhs, Angle rhs) -> Angle
  {
    if (lhs.HeldInDegrees() || rhs.HeldInDegrees())
    {
      return Degrees(lhs.InDegrees() - rhs.InDegrees());
    }
    return Radians(lhs.InRadians() - rhs.InRadians());
  }

  /**
   *  The angle the other way.
   *
   *  @param  angle The angle.
   *  @return       Its negative, exactly, in its unit.
   */
  [[nodiscard]] constexpr auto operator-(Angle angle) -> Angle
  {
    return angle.HeldInDegrees() ? Degrees(-angle.InDegrees()) : Radians(-angle.InRadians());
  }

  /**
   *  An angle scaled.
   *
   *  @param  lhs The angle.
   *  @param  rhs The scale.
   *  @return     The angle times the scale, rounded once, in its unit.
   */
  [[nodiscard]] constexpr auto operator*(Angle lhs, double rhs) -> Angle
  {
    return lhs.HeldInDegrees() ? Degrees(lhs.InDegrees() * rhs) : Radians(lhs.InRadians() * rhs);
  }

  /**
   *  An angle divided.
   *
   *  @param  lhs The angle.
   *  @param  rhs The divisor.
   *  @return     The angle over the divisor, rounded once, in its unit.
   */
  [[nodiscard]] constexpr auto operator/(Angle lhs, double rhs) -> Angle
  {
    return lhs.HeldInDegrees() ? Degrees(lhs.InDegrees() / rhs) : Radians(lhs.InRadians() / rhs);
  }

  /**
   *  The same direction within (-180, 180] degrees.
   *
   *  Not [0, 360): -1e-20 degrees has no double there but 0 and 360 - 2^-44, while (-180, 180]
   *  holds every wrapped double.
   *
   *  @param  angle The angle.
   *  @return       In degrees: the angle less whole turns, exactly when it is held in degrees; +0
   *                for a whole number of turns.
   */
  [[nodiscard]] GEOXY_API auto Wrapped(Angle angle) -> Angle;

  /**
   *  Whether a sweep from a start passes a direction, its ends included.
   *
   *  The sweep runs counterclockwise when positive; a turn or more passes every direction. The end
   *  is the exact start + sweep, not its rounding, so the answer never flips near an end. Held in
   *  degrees, all three are decided exactly; an angle held in radians is taken in degrees to 2^-104
   *  first.
   *
   *  @param  direction The direction, any turn of it.
   *  @param  start     Where the sweep starts, any turn of it.
   *  @param  sweep     The sweep.
   *  @return           Whether it passes the direction; false when any is not finite.
   */
  [[nodiscard]] GEOXY_API auto Within(Angle direction, Angle start, Angle sweep) -> bool;

  /**
   *  The unit vector at an angle from the x axis, counterclockwise.
   *
   *  Held in degrees, the angle is reduced exactly; in radians, by Payne and Hanek's reduction to
   *  2^-201 of a quarter turn, past any number of turns. The double nearest a multiple of a quarter
   *  turn, below 16 radians, is that quarter turn: 2 pi as a double is a whole turn.
   *
   *  Each component is rounded by Ziv's strategy: a fast pair within 2^-68 of the exact value
   *  rounds wherever it is far enough from a halfway point between doubles, as it is about 9,999
   *  times in 10,000, and the rest are worked to as many bits as it takes, never on a halfway
   *  point, by Lindemann's theorem in radians and Niven's in degrees.
   *
   *  @param  angle The angle, past a turn too.
   *  @return       Each component the double nearest the exact value, always, so exact wherever
   *                that is a double (at multiples of 30 degrees and at quarter turns); a zero
   *                component +0.
   */
  [[nodiscard]] GEOXY_API auto Direction(Angle angle) -> Vector;

  /**
   *  The unit vector part of the way along a sweep.
   *
   *  @param  start The start.
   *  @param  sweep The sweep.
   *  @param  along The part of the sweep, 0 at the start and 1 at its end.
   *  @return       The direction of start + sweep * along: held in degrees, that sum and product
   *                not rounded, so the end of a sweep is its exact end; else as Direction of it.
   */
  [[nodiscard]] GEOXY_API auto Direction(Angle start, Angle sweep, double along) -> Vector;

  /**
   *  The angle of a displacement from the x axis.
   *
   *  A mirrored -0 does not turn pi into -pi. Rounded by Ziv's strategy, as Direction is.
   *
   *  @param  vector The displacement.
   *  @return        In radians, in [-pi, pi], the double nearest the exact angle, always; so on the
   *                 axes the double nearest a multiple of pi / 2, and 45 degrees exactly on a
   *                 diagonal; 0 for no displacement. -pi only where the exact angle is just above
   *                 it, y < 0 and x < 0 with y below an ulp of x.
   */
  [[nodiscard]] GEOXY_API auto Heading(Vector vector) -> Angle;

  /**
   *  The signed angle from one displacement to another.
   *
   *  Worked from the cross and dot products, exactly, so it keeps its precision near 0 and near 180
   *  degrees, and is 0 or 180 exactly for parallel displacements.
   *
   *  @pre          The nonzero components of each displacement are within 2^960 of each other.
   *  @param  start The displacement turned from.
   *  @param  end   The displacement turned to.
   *  @return       In radians, in [-pi, pi], positive counterclockwise, as Heading is accurate.
   */
  [[nodiscard]] GEOXY_API auto Between(Vector start, Vector end) -> Angle;

  /**
   *  The order of two displacements by their exact headings in (-180, 180], from signs alone.
   *
   *  @pre        The nonzero components of each displacement are within 2^960 of each other.
   *  @param  lhs One displacement.
   *  @param  rhs The other.
   *  @return     Their order; no displacement is at 0 with the positive x axis.
   */
  [[nodiscard]] GEOXY_API auto Order(Vector lhs, Vector rhs) -> std::strong_ordering;
} // namespace Geoxy

#endif /* GEOXY_ANGLE_HPP */
