/***************************************************************************************************
 *  @file      source/secular.cpp
 *  @brief     The place of an ellipse nearest a point, by Eberly's secular equation in the
 *             ellipse's principal frame, polished on its own map.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "secular.hpp"

#include "real.hpp"
#include "trig.hpp"

#include "geoxy/angle.hpp"
#include "geoxy/at.hpp"
#include "geoxy/ellipse.hpp"
#include "geoxy/matrix.hpp"
#include "geoxy/place.hpp"
#include "geoxy/point.hpp"
#include "geoxy/vector.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <numbers>
#include <optional>
#include <span>

namespace Geoxy
{
  namespace
  {
    /**
     *  An ellipse's map divided by a power of two near its size, and its principal frame: its unit
     *  circle turned by phi goes to its longest semi-axis, a, and turned a quarter more to its
     *  shortest, b.
     *
     *  Found in doubles, the two semi-axes are square to a rounding of the ellipse's size, so the
     *  frame is the ellipse's to that rounding, and what is found in it is polished on the map.
     */
    struct PrincipalFrame
    {
      Vector X      = {};  /**< The map's first axis, divided.          */
      Vector Y      = {};  /**< Its second axis, divided.               */
      Vector Major  = {};  /**< The longest semi-axis, divided.         */
      Vector Minor  = {};  /**< The shortest semi-axis, divided.        */
      double Cosine = 1.0; /**< The cosine of phi.                      */
      double Sine   = 0.0; /**< Its sine.                               */
      double Spread = 0.0; /**< a² − b², divided twice.                 */
      double Size   = 0.0; /**< The squares of the axes, summed.        */
      double Scale  = 1.0; /**< The power of two the map is divided by. */
    };

    /**
     *  An ellipse's principal frame.
     *
     *  The square of the length its map gives the unit circle's point at an angle t is (xx + yy) /
     *  2 + ((xx − yy) / 2) cos 2t + xy sin 2t, largest where 2t is the angle of (xx − yy, 2 xy),
     *  whose length is a² − b²; phi's cosine and sine come from that angle's by the half-angle
     *  formulas, the larger of the two by a root, the other by a quotient.
     *
     *  @param  map The ellipse's map.
     *  @return     The frame; nothing for a circle, or a map that is not finite.
     */
    auto Oriented(const Matrix &map) -> std::optional<PrincipalFrame>
    {
      const double most =
          std::max({std::abs(map.X.X), std::abs(map.X.Y), std::abs(map.Y.X), std::abs(map.Y.Y)});
      int power = 0;
      std::frexp(most, &power);
      const double scale   = std::ldexp(1.0, power);
      const Vector x       = map.X / scale;
      const Vector y       = map.Y / scale;
      const double squareX = Dot(x, x);
      const double squareY = Dot(y, y);
      const double apart   = squareX - squareY;
      const double twice   = 2.0 * Dot(x, y);
      const double spread  = std::sqrt((apart * apart) + (twice * twice));
      if (!(spread > 0.0))
      {
        return std::nullopt;
      }
      double cosine = 0.0;
      double sine   = 0.0;
      if (apart >= 0.0)
      {
        cosine = std::sqrt((spread + apart) / (2.0 * spread));
        sine   = twice / (2.0 * spread * cosine);
      }
      else
      {
        sine   = std::sqrt((spread - apart) / (2.0 * spread));
        cosine = twice / (2.0 * spread * sine);
      }
      return PrincipalFrame{.X      = x,
                            .Y      = y,
                            .Major  = (x * cosine) + (y * sine),
                            .Minor  = (y * cosine) - (x * sine),
                            .Cosine = cosine,
                            .Sine   = sine,
                            .Spread = spread,
                            .Size   = squareX + squareY,
                            .Scale  = scale};
    }

    /**
     *  A point's terms in the secular equation of one side of an ellipse's major axis, the point at
     *  (x, y) in the principal frame, both taken not below nothing.
     *
     *  On the point's own side its foot is (ax / (s + a² − b²), by / s) on the unit circle, for the
     *  root s past nothing of G(s) = (ax / (s + a² − b²))² + (by / s)² − 1. Across the major axis
     *  it is (ax / (a² − b² − s), −by / s), the root of G with the first term's s turned negative,
     *  from nothing to a² − b².
     */
    struct Secular
    {
      double Along  = 0.0; /**< ax, a times the point's place along the major axis. */
      double Across = 0.0; /**< by, b times its place across it.                    */
      double Spread = 0.0; /**< a² − b².                                            */
      double Side   = 1.0; /**< 1 on the point's side of the major axis, −1 across. */
    };

    /**
     *  Where Newton's method starts on the point's own side: below the root, where G is not below
     *  nothing, as near it as bounds tell.
     *
     *  The root is past by, past ax − (a² − b²), and past the length of (ax, by) less a² − b²,
     *  where G would be nothing were its first denominator its second. Near the evolute's cusp on
     *  the major axis, ax near a² − b², those leave it far below when by is small, and Newton's
     *  method would creep up a half at a step: G there is at least (by / s)² − 2 (s + gap) / (a² −
     *  b²), for gap = a² − b² − ax, which is not below nothing where s² (s + gap) is at most (by)²
     *  (a² − b²) / 2, and so at the s given, by a cube root, a square root or a quotient.
     *
     *  @param  secular The point's terms, on its own side.
     *  @return         The start.
     */
    auto Start(const Secular &secular) -> double
    {
      const double along  = secular.Along;
      const double across = secular.Across;
      const double spread = secular.Spread;
      const double least  = std::max(across, Trig::Hypot(along, across) - spread);
      if (!(4.0 * along > 3.0 * spread) || !(4.0 * across < spread))
      {
        return least;
      }
      const double gap   = spread - along;
      const double rise  = std::abs(gap);
      const double third = Trig::Cbrt(across);
      const double cube  = third * third * Trig::Cbrt(spread / (gap > 0.0 ? 4.0 : 8.0));
      double       cusp  = cube + std::max(-gap, 0.0);
      if (rise > cube)
      {
        cusp = gap > 0.0 ? across * std::sqrt(spread / (4.0 * gap))
                         : rise + ((across / gap) * (across / gap) * spread / 8.0);
      }
      return std::max(least, cusp);
    }

    /**
     *  The root of the secular equation of one side, by Newton's method from below it, where G is
     *  not below nothing.
     *
     *  G falls and bends up there, so no step passes the root but by rounding; a step of less than
     *  a 2^-26th of where it is leaves the place good to the last bits, a step of Newton's method
     *  on the map then making up any rounding.
     *
     *  @param  secular The point's terms.
     *  @param  start   Where to start.
     *  @return         The root; across the major axis, below nothing where there is none, G rising
     *                  again before it comes down to nothing, or the step passing a² − b².
     */
    auto Newtoned(const Secular &secular, double start) -> double
    {
      const double pole =
          secular.Side > 0.0 ? std::numeric_limits<double>::infinity() : secular.Spread;
      double root = start;
      for (;;)
      {
        const double near   = 1.0 / (secular.Spread + (secular.Side * root));
        const double far    = 1.0 / root;
        const double along  = secular.Along * near;
        const double across = secular.Across * far;
        const double excess = (along * along) + (across * across) - 1.0;
        if (!(excess > 0.0))
        {
          return root;
        }
        const double fall = (secular.Side * along * along * near) + (across * across * far);
        if (!(fall > 0.0))
        {
          return -1.0;
        }
        const double step = excess / (2.0 * fall);
        const double next = root + step;
        if (!(next < pole))
        {
          return -1.0;
        }
        if (!(step > 0x1p-26 * root))
        {
          return next;
        }
        root = next;
      }
    }

    /**
     *  A place on a frame's unit circle, as the map's own unit circle has it.
     *
     *  @param  frame The frame.
     *  @param  unit  The place on the frame's unit circle.
     *  @return       The place on the map's.
     */
    auto Unturned(const PrincipalFrame &frame, Vector unit) -> Vector
    {
      return {.X = (frame.Cosine * unit.X) - (frame.Sine * unit.Y),
              .Y = (frame.Sine * unit.X) + (frame.Cosine * unit.Y)};
    }

    /**
     *  A place on the map's unit circle made of length one, then polished by a step of Newton's
     *  method on (C − p) · C', C the map's own point there and p the point.
     *
     *  @param  frame The frame.
     *  @param  away  The point from the ellipse's center, divided as the map is.
     *  @param  unit  The place, of a length near one.
     *  @return       The place polished; only of length one where the distance does not bend up, at
     *                the evolute and past it, where a step could run off.
     */
    auto Polished(const PrincipalFrame &frame, Vector away, Vector unit) -> Vector
    {
      const Vector turned = unit / std::sqrt(Dot(unit, unit));
      const Vector here   = (frame.X * turned.X) + (frame.Y * turned.Y);
      const Vector off    = here - away;
      const Vector way    = (frame.Y * turned.X) - (frame.X * turned.Y);
      const double bend   = Dot(way, way) - Dot(off, here);
      if (!(bend > 0.0))
      {
        return turned;
      }
      const double turn = -Dot(off, way) / bend;
      return Vector{.X = turned.X - (turn * turned.Y), .Y = turned.Y + (turn * turned.X)}
           / std::sqrt(1.0 + (turn * turn));
    }

    /**
     *  How far a point is from the place of a map's unit circle.
     *
     *  @param  frame The frame.
     *  @param  away  The point from the ellipse's center, divided as the map is.
     *  @param  unit  The place, of length one.
     *  @return       The distance, divided as the map is.
     */
    auto Remoteness(const PrincipalFrame &frame, Vector away, Vector unit) -> double
    {
      return Length((frame.X * unit.X) + (frame.Y * unit.Y) - away);
    }

    /**
     *  The parameter at which an ellipse swept round no more than once passes an angle of its unit
     *  circle.
     *
     *  @param  ellipse The ellipse.
     *  @param  angle   The angle, in radians.
     *  @return         The parameter, from 0 to 1; nothing where the sweep does not pass it.
     */
    auto Unwound(const Ellipse &ellipse, double angle) -> std::optional<double>
    {
      constexpr double Whole  = 2.0 * std::numbers::pi;
      const double     sweep  = ellipse.Sweep.InRadians();
      const double     turned = angle - ellipse.Start.InRadians();
      const double     turns  = turned / Whole;
      const double     along =
          (turned - (Whole * (sweep > 0.0 ? std::floor(turns) : std::ceil(turns)))) / sweep;
      if (!(along <= 1.0))
      {
        return std::nullopt;
      }
      return std::max(along, 0.0);
    }

    /**
     *  A place on the map's unit circle made of length one, then polished by a step of Newton's
     *  method on (C − p) · C', where the distance is least or most alike: a foot of any kind, of a
     *  root that is not double, so that the step's slope is not nothing.
     *
     *  @param  frame The frame.
     *  @param  away  The point from the ellipse's center, divided as the map is.
     *  @param  unit  The place, of a length near one.
     *  @return       The place polished.
     */
    auto Trued(const PrincipalFrame &frame, Vector away, Vector unit) -> Vector
    {
      const Vector turned = unit / std::sqrt(Dot(unit, unit));
      const Vector here   = (frame.X * turned.X) + (frame.Y * turned.Y);
      const Vector off    = here - away;
      const Vector way    = (frame.Y * turned.X) - (frame.X * turned.Y);
      const double turn   = -Dot(off, way) / (Dot(way, way) - Dot(off, here));
      return Vector{.X = turned.X - (turn * turned.Y), .Y = turned.Y + (turn * turned.X)}
           / std::sqrt(1.0 + (turn * turn));
    }

    /**
     *  A place weighed for the nearest: its parameter and its distance from the point.
     */
    struct Candidate
    {
      double Along = 0.0; /**< The parameter.                       */
      double Away  = 0.0; /**< The distance, divided as the map is. */
    };
  } // namespace

  auto Closest(const Ellipse &ellipse, Point point) -> std::optional<Place>
  {
    const std::optional<PrincipalFrame> framed = Oriented(ellipse.Map);
    if (!framed || Zero(ellipse.Sweep.InRadians()))
    {
      return std::nullopt;
    }
    const PrincipalFrame &frame  = *framed;
    const Vector          away   = (point - ellipse.Map.Origin) / frame.Scale;
    const double          along  = Dot(away, frame.Major);
    const double          across = Dot(away, frame.Minor);
    const Secular         near   = {
        .Along = std::abs(along), .Across = std::abs(across), .Spread = frame.Spread, .Side = 1.0};
    const auto unturned = [&frame, along, across](Vector foot, double side) -> Vector
    {
      return Unturned(
          frame, {.X = std::copysign(foot.X, along), .Y = side * std::copysign(foot.Y, across)});
    };

    /**
     *  The nearest foot, on the frame's unit circle in the point's quarter: on the major axis
     *  inside the evolute, where G has no root, the two feet either side of it, at a cosine of ax /
     *  (a² − b²), as near.
     */
    double root = 0.0;
    Vector foot = {};
    if (near.Across > 0.0 || near.Along > near.Spread)
    {
      root = Newtoned(near, Start(near));
      foot = {.X = near.Along / (near.Spread + root), .Y = near.Across / root};
    }
    else
    {
      const double cosine = near.Along / near.Spread;
      foot                = {.X = cosine, .Y = std::sqrt(1.0 - (cosine * cosine))};
    }
    const Vector                first = Polished(frame, away, unturned(foot, 1.0));
    const std::optional<double> swept = Unwound(ellipse, Trig::Atan2(first.Y, first.X));
    const double                tie   = 0x1p-44
                                      * (std::abs(away.X) + std::abs(away.Y))
                                      * frame.Size
                                      * (1.0 + (frame.Size / frame.Spread));
    if (swept && near.Across > tie)
    {
      return Place{.At = At(ellipse, *swept), .Along = *swept};
    }

    /**
     *  Left out by the sweep, or so near the major axis that the frame's rounding may have put it
     *  on the wrong side: the second least across the axis, if the point is inside the evolute, and
     *  the ends where the first is left out, weighed with it.
     */
    std::array<Candidate, 4> candidates = {};
    std::size_t              count      = 0;
    if (swept)
    {
      candidates.at(count++) = {.Along = *swept, .Away = Remoteness(frame, away, first)};
    }
    std::optional<Vector> other;
    if (!(near.Across > 0.0))
    {
      if (foot.Y > 0.0)
      {
        other = foot;
      }
    }
    else if (root < near.Spread)
    {
      const double beyond = Newtoned(
          {.Along = near.Along, .Across = near.Across, .Spread = near.Spread, .Side = -1.0}, root);
      if (beyond > 0.0)
      {
        other = Vector{.X = near.Along / (near.Spread - beyond), .Y = near.Across / beyond};
      }
    }
    if (other)
    {
      const Vector second = Polished(frame, away, unturned(*other, -1.0));
      if (const std::optional<double> beside = Unwound(ellipse, Trig::Atan2(second.Y, second.X)))
      {
        candidates.at(count++) = {.Along = *beside, .Away = Remoteness(frame, away, second)};
      }
    }
    if (!swept)
    {
      const double start     = ellipse.Start.InRadians();
      const double end       = start + ellipse.Sweep.InRadians();
      candidates.at(count++) = {.Along = 0.0, .Away = Remoteness(frame, away, Trig::Way(start))};
      candidates.at(count++) = {.Along = 1.0, .Away = Remoteness(frame, away, Trig::Way(end))};
    }
    const std::span<Candidate> weighed = std::span{candidates}.first(count);
    std::ranges::sort(weighed, {}, &Candidate::Along);
    const double best = std::ranges::min(weighed, {}, &Candidate::Away).Along;
    return Place{.At = At(ellipse, best), .Along = best};
  }

  auto Plumbs(const Ellipse &ellipse, Point point) -> std::optional<Footholds>
  {
    const std::optional<PrincipalFrame> framed = Oriented(ellipse.Map);
    if (!framed)
    {
      return std::nullopt;
    }
    const PrincipalFrame &frame  = *framed;
    const Vector          away   = (point - ellipse.Map.Origin) / frame.Scale;
    const double          along  = Dot(away, frame.Major);
    const double          across = Dot(away, frame.Minor);
    const double          spread = frame.Spread;
    const double          reach  = (std::abs(away.X) + std::abs(away.Y)) * frame.Size;
    /**
     *  The point's place against the evolute, (cbrt((ax)²) + cbrt((by)²))³ against (a² − b²)².
     */
    const double third  = Trig::Cbrt(along * along) + Trig::Cbrt(across * across);
    const double cusped = third * third * third;
    const double square = spread * spread;
    const bool   thin   = !(frame.Size - spread > 0x1p-40 * frame.Size);
    if (thin
        || !(std::abs(along) > 0x1p-40 * reach)
        || !(std::abs(across) > 0x1p-40 * reach)
        || !(std::abs(cusped - square) > 0x1p-30 * square))
    {
      return std::nullopt;
    }
    const Secular near = {
        .Along = std::abs(along), .Across = std::abs(across), .Spread = spread, .Side = 1.0};
    const Secular far = {
        .Along = std::abs(across), .Across = std::abs(along), .Spread = spread, .Side = 1.0};
    std::array<Vector, 4> feet  = {};
    std::size_t           count = 0;
    const double          root  = Newtoned(near, Start(near));
    const double          back  = Newtoned(far, Start(far));
    feet.at(count++)            = {.X = near.Along / (spread + root), .Y = near.Across / root};
    feet.at(count++)            = {.X = -near.Along / back, .Y = -near.Across / (spread + back)};
    if (cusped < square)
    {
      /**
       *  Within the evolute both roots between the poles, each from the side of its own pole: at
       *  by, below the first, where G is ax² / (a² − b² − by)² and falls, as (cbrt((ax)²) +
       *  cbrt((by)²))³ below (a² − b²)² makes ax² by below (a² − b² − by)³; and at ax, measured
       *  from the other pole, alike. G bends up, so neither run passes its root, and the margin
       *  from the evolute keeps the two roots apart and simple.
       */
      const double first =
          Newtoned({.Along = near.Along, .Across = near.Across, .Spread = spread, .Side = -1.0},
                   near.Across);
      const double second = Newtoned(
          {.Along = near.Across, .Across = near.Along, .Spread = spread, .Side = -1.0}, near.Along);
      feet.at(count++) = {.X = near.Along / (spread - first), .Y = -near.Across / first};
      feet.at(count++) = {.X = near.Along / second, .Y = -near.Across / (spread - second)};
    }
    Footholds holds = {.Angles = {}, .Count = count};
    for (std::size_t index = 0; index < count; ++index)
    {
      const Vector foot      = feet.at(index);
      const Vector placed    = {.X = foot.X * std::copysign(1.0, along),
                                .Y = foot.Y * std::copysign(1.0, across)};
      const Vector unit      = Trued(frame, away, Unturned(frame, placed));
      holds.Angles.at(index) = Trig::Atan2(unit.Y, unit.X);
    }
    return holds;
  }
} // namespace Geoxy
