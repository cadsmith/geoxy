/***************************************************************************************************
 *  @file      source/angle.cpp
 *  @brief     Angles: degrees and radians converted, and the two exact tests on degrees, wrapping
 *             and sweeping.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "geoxy/angle.hpp"

#include "degrees.hpp"
#include "radians.hpp"
#include "real.hpp"
#include "turning.hpp"
#include "twofold.hpp"

#include <array>
#include <bit>
#include <cmath>
#include <cstdint>
#include <optional>

namespace Geoxy
{
  namespace
  {
    constexpr double Assurance = 0x1p-40; /**< Beyond it the rounded sweep test cannot be wrong. */

    /**
     *  An angle in degrees as a pair: exactly when held in degrees; when in radians, less whole
     *  turns and within 2^-104.
     *
     *  @param  angle The angle.
     *  @return       Its degrees, as a pair.
     */
    auto AsDegrees(Angle angle) -> Twofold::Two
    {
      if (angle.HeldInDegrees())
      {
        return {.High = angle.InDegrees(), .Low = 0.0};
      }
      return Turning::Reduced(angle.InRadians());
    }

    /**
     *  A sweep in degrees as a pair: exactly when held in degrees, within 2^-104 when in radians.
     *
     *  @param  sweep The sweep.
     *  @return       Its degrees, as a pair, not reduced.
     */
    auto Degreed(Angle sweep) -> Twofold::Two
    {
      if (sweep.HeldInDegrees())
      {
        return {.High = sweep.InDegrees(), .Low = 0.0};
      }
      const double radians = sweep.InRadians();
      /**
       *  The double nearest a multiple of a quarter turn is that quarter turn, below 16 radians, as
       *  in Direction: a sweep of 2 pi as a double is a whole turn.
       */
      if (std::abs(radians) < 16.0)
      {
        const double       quarters = std::rint(radians / Table::HalfPi.High);
        const Twofold::Two product  = Twofold::TwoProduct(quarters, Table::HalfPi.High);
        const double       apart =
            ((radians - product.High) - product.Low) - (quarters * Table::HalfPi.Low);
        int exponent = 0;
        static_cast<void>(std::frexp(radians, &exponent));
        if (!Zero(quarters) && std::abs(apart) < std::ldexp(1.0, exponent - 54))
        {
          return {.High = 90.0 * quarters, .Low = 0.0};
        }
      }
      const double high = radians * Table::ToDegrees.High;
      return {.High = high,
              .Low  = std::fma(radians, Table::ToDegrees.High, -high)
                    + (radians * Table::ToDegrees.Low)};
    }
  } // namespace

  auto Angle::Converted(double value, bool toRadians) -> double
  {
    /**
     *  The double nearest value C, C = pi / 180 or 180 / pi, for every double value.
     *
     *  tests/proof/convert.py proves no double's product with C comes nearer a halfway point than
     *  2^-61 ulp. Here C is held in three doubles to 2^-155, value C as an exact sum of four
     *  doubles and a fifth rounded, and the halfway point is passed or not by an exact sign: the
     *  answer is off by below 2^-100 ulp from the exact product, so it is always on the right side.
     */
    constexpr double   Underflows = 0x1p-960;
    constexpr double   Lift       = 0x1p120;
    constexpr double   Drop       = 0x1p-120;
    const Twofold::Two factor     = toRadians ? Table::ToRadians : Table::ToDegrees;
    const double       third      = toRadians ? Table::ToRadiansThird : Table::ToDegreesThird;
    const double       size       = std::abs(value);
    const double       sign       = std::copysign(1.0, value);
    /**
     *  First the product as a pair, of C's two leading doubles: within 2^-104 of the exact one, the
     *  two doubles off C by 2^-108 and the fused add's rounding 2^-105. Ziv's test, given 2^-80 to
     *  keep its own roundings inside, settles all but one input in 2^27 or so; those, and any
     *  outside the range where nothing underflows or overflows, go the exact way below.
     */
    if (size > 0x1p-900 && size < 0x1p900)
    {
      const Twofold::Two product = Twofold::TwoProduct(size, factor.High);
      const double       low     = std::fma(size, factor.Low, product.Low);
      if (const std::optional<double> settled = Turning::Settled(product.High, low, 0x1p-80))
          [[likely]]
      {
        return sign * *settled;
      }
    }
    /**
     *  A value so small its products' roundoffs would underflow is lifted by 2^120 and dropped
     *  after, exactly; its result is then rounded to the subnormal grid where it falls there.
     */
    const bool   lifted = size < Underflows;
    const double taken  = lifted ? size * Lift : size;

    const Twofold::Two first = Twofold::TwoProduct(taken, factor.High);
    if (!std::isfinite(first.High) || Zero(first.High))
    {
      return sign * first.High;
    }
    const Twofold::Two second = Twofold::TwoProduct(taken, factor.Low);
    const Twofold::Two middle = Twofold::TwoSum(first.Low, second.High);
    const double       small  = (middle.Low + second.Low) + (taken * third);
    const Twofold::Two lead   = Twofold::TwoSum(first.High, middle.High);

    /**
     *  The grid the result is rounded to: the lead's own ulp, or the subnormal grid, 2^-1074
     *  dropped, where the lifted result is below 2^-902.
     */
    constexpr double   Grid = 0x1p-954;
    const bool         fine = lifted && lead.High < 0x1p-902;
    const double       base = fine ? std::rint(lead.High / Grid) * Grid : lead.High;
    const Twofold::Two rest = Twofold::TwoSum(lead.High - base, lead.Low);
    const double       unit =
        fine ? Grid
             : std::bit_cast<double>(std::bit_cast<std::uint64_t>(base) & 0x7ff0000000000000U)
                   * 0x1p-52;
    /**
     *  Below a power of two the grid is twice as fine; the rest is below half of it either way.
     */
    /**
     *  Arithmetic, not branches: the step to a neighbour is taken only within 2^-52 ulp of a
     *  halfway point, which no ordinary input reaches, and tests reach by tests/proof's hardest.
     */
    const double way   = std::copysign(1.0, rest.High + (rest.Low + small));
    const bool   power = (std::bit_cast<std::uint64_t>(base) & 0xfffffffffffffU) == 0U;
    const auto   finer = static_cast<double>(static_cast<unsigned>(way < 0.0)
                                             & static_cast<unsigned>(!fine)
                                             & static_cast<unsigned>(power));
    const double step  = std::copysign(unit * (1.0 - (finer / 2.0)), way);
    /**
     *  Past the halfway point to the neighbour, by the exact sign of rest - step / 2.
     */
    const Twofold::Two apart   = Twofold::TwoSum(rest.High, -step / 2.0);
    const double       beyond  = apart.High + ((apart.Low + rest.Low) + small);
    const auto         passes  = static_cast<double>(beyond * way > 0.0);
    double             nearest = base + (passes * step);
    if (lifted)
    {
      nearest *= Drop;
    }
    return sign * nearest;
  }

  auto Wrapped(Angle angle) -> Angle
  {
    if (angle.HeldInDegrees() || !std::isfinite(angle.InRadians()))
    {
      return Degrees(Turning::Wrap(angle.InDegrees()) + 0.0);
    }
    const Twofold::Two degrees = Turning::Reduced(angle.InRadians());
    return Degrees(Turning::Wrap(degrees.High + degrees.Low) + 0.0);
  }

  auto Within(Angle direction, Angle start, Angle sweep) -> bool
  {
    if (!direction.HeldInDegrees() && !start.HeldInDegrees() && !sweep.HeldInDegrees())
    {
      /**
       *  In radians the ends are the only places the answer can turn on exactly, 2 pi being no
       *  double: there it is decided exactly, and elsewhere in degrees to 2^-104.
       */
      const double toward = direction.InRadians();
      const double from   = start.InRadians();
      const double run    = sweep.InRadians();
      if (std::isfinite(toward)
          && std::isfinite(from)
          && std::isfinite(run)
          && (Turning::SumSign(std::array{toward, -from}) == 0
              || Turning::SumSign(std::array{toward, -from, -run}) == 0))
      {
        return true;
      }
    }
    const auto held = [](Angle angle) -> double
    {
      return angle.HeldInDegrees() ? angle.InDegrees() : angle.InRadians();
    };
    if (!std::isfinite(held(direction))
        || !std::isfinite(held(start))
        || !std::isfinite(held(sweep)))
    {
      return false;
    }
    Twofold::Two toward = AsDegrees(direction);
    Twofold::Two from   = AsDegrees(start);
    Twofold::Two run    = Degreed(sweep);
    /**
     *  A turn or more passes every direction. The big part reaches 360 only where the sweep does:
     *  in degrees it is the sweep; in radians the double taken for 2 pi is 360 exactly, and no
     *  other comes within an ulp of 360 from below.
     */
    if (std::abs(run.High) >= 360.0)
    {
      return true;
    }
    /**
     *  A clockwise sweep is the counterclockwise one mirrored; negating is exact.
     */
    const double way = run.High < 0.0 ? -1.0 : 1.0;
    toward           = {.High = Turning::Wrap(toward.High * way), .Low = toward.Low * way};
    from             = {.High = Turning::Wrap(from.High * way), .Low = from.Low * way};
    run              = {.High = run.High * way, .Low = run.Low * way};

    /**
     *  Rounded first, each step off by at most an ulp of 360: away by 2^-40 from where the turn
     *  added or the comparison could change, the rounded answer is the exact one.
     */
    const double apart   = (toward.High - from.High) + (toward.Low - from.Low);
    const double rounded = apart < 0.0 ? apart + 360.0 : apart;
    if (std::abs(rounded - run.High) > Assurance && std::abs(apart) > Assurance)
    {
      return rounded <= run.High;
    }

    /**
     *  Exactly: the whole turns that bring the difference into [0, 360), found by exact signs, so
     *  the test is (toward - from) mod 360 <= run. The difference is below a turn already: both are
     *  within (-180, 180] but for a part below half an ulp of 180, and those come to less than the
     *  ulp between 180 and the double below it.
     */
    double turns = 0.0;
    while (
        Turning::SumSign(std::array{toward.High, -from.High, toward.Low, -from.Low, turns * 360.0})
        < 0)
    {
      turns += 1.0;
    }
    return Turning::SumSign(std::array{toward.High, -from.High, toward.Low, -from.Low,
                                       turns * 360.0, -run.High, -run.Low})
        <= 0;
  }
} // namespace Geoxy
