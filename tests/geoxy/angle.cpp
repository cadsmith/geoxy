/***************************************************************************************************
 *  @file      tests/geoxy/angle.cpp
 *  @brief     Tests of angles: exact answers against whole numbers, rounded ones against MPFR's,
 *             worked out once by build/lab/angle and written here.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "geoxy/angle.hpp"

#include "draws.hpp"
#include "precise.hpp"
#include "radians.hpp"
#include "real.hpp"

#include "geoxy/vector.hpp"

#include <doctest/doctest.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <compare>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <initializer_list>
#include <limits>
#include <optional>
#include <ranges>
#include <span>

namespace Geoxy
{
  static_assert(std::equal_to<double>{}(Radians(1.5).InRadians(), 1.5));
  static_assert(std::equal_to<double>{}(Degrees(30.0).InDegrees(), 30.0));
  static_assert(std::equal_to<double>{}(Angle{}.InRadians(), 0.0));
  static_assert(Degrees(1) + Degrees(2) == Degrees(3));
  static_assert(Radians(3) - Radians(2) == Radians(1));
  static_assert(-Degrees(1)             == Degrees(-1));
  static_assert(Degrees(1) * 2.0        == Degrees(2));
  static_assert(Radians(2) / 2.0        == Radians(1));
  static_assert(Degrees(1) < Degrees(2));
  static_assert(Turn.HeldInDegrees());
  static_assert(!Radians(1).HeldInDegrees());

  namespace
  {
    /**
     *  A turn in units of 2^-44 degree.
     */
    constexpr std::int64_t WholeTurn = std::int64_t{360} << 44U;

    /**
     *  Degrees as a whole number of 2^-44 degree, when they are one and below 2^11.
     *
     *  @param  degrees The degrees.
     *  @return         The whole number, or none.
     */
    auto Fixed(double degrees) -> std::optional<std::int64_t>
    {
      const double scaled = std::ldexp(degrees, 44);
      if (!(std::abs(degrees) < 2048.0) || !std::equal_to<double>{}(scaled, std::trunc(scaled)))
      {
        return std::nullopt;
      }
      return static_cast<std::int64_t>(scaled);
    }

    /**
     *  A whole number reduced into [0, a turn).
     *
     *  @param  value The whole number.
     *  @return       The value less whole turns.
     */
    auto Reduced(std::int64_t value) -> std::int64_t
    {
      const std::int64_t rest = value % WholeTurn;
      return rest < 0 ? rest + WholeTurn : rest;
    }

    /**
     *  Within, worked in whole numbers.
     *
     *  @param  direction The direction.
     *  @param  start     The start.
     *  @param  sweep     The sweep.
     *  @return           Whether the sweep passes the direction.
     */
    auto Truly(std::int64_t direction, std::int64_t start, std::int64_t sweep) -> bool
    {
      if (sweep >= WholeTurn || sweep <= -WholeTurn)
      {
        return true;
      }
      if (sweep < 0)
      {
        return Reduced(start - direction) <= -sweep;
      }
      return Reduced(direction - start) <= sweep;
    }

    /**
     *  The double so many doubles away from another.
     *
     *  @param  value The double.
     *  @param  steps How many doubles away, either way.
     *  @return       The double there.
     */
    auto Nudged(double value, int steps) -> double
    {
      for (int step = 0; step < std::abs(steps); ++step)
      {
        value = std::nextafter(value, steps < 0 ? -Infinite : Infinite);
      }
      return value;
    }

    /**
     *  Whether two doubles are the same double, the sign of a zero included.
     *
     *  @param  lhs One.
     *  @param  rhs The other.
     *  @return     Whether they are.
     */
    auto Same(double lhs, double rhs) -> bool
    {
      return std::equal_to<double>{}(lhs, rhs) && std::signbit(lhs) == std::signbit(rhs);
    }

    /**
     *  The sign of a whole number.
     *
     *  @param  value The number.
     *  @return       1, 0 or -1.
     */
    auto Signed(int value) -> int
    {
      return static_cast<int>(value > 0) - static_cast<int>(value < 0);
    }

    /**
     *  Which stretch of headings a displacement is in, in order: (-180, 0), 0, (0, 180), 180.
     *
     *  @param  x Along the x axis.
     *  @param  y Along the y axis.
     *  @return   0 to 3; no displacement is 1.
     */
    auto Half(int x, int y) -> int
    {
      if (y != 0)
      {
        return y < 0 ? 0 : 2;
      }
      return x < 0 ? 3 : 1;
    }

    /**
     *  cos and sin of an angle, each the double nearest the exact value, from MPFR.
     */
    struct CosSin
    {
      double Angle  = 0.0; /**< The angle.  */
      double Cosine = 0.0; /**< Its cosine. */
      double Sine   = 0.0; /**< Its sine.   */
    };

    /**
     *  In degrees.
     */
    constexpr std::array InDegrees{
        CosSin{.Angle = 0x1p+0, .Cosine = 0x1.ffec097f5af8ap-1, .Sine = 0x1.1df0b2b89dd1ep-6},
        CosSin{.Angle  = 0x1.0624dd2f1a9fcp-10,
               .Cosine = 0x1.fffffffeb111dp-1,
               .Sine   = 0x1.24d15211536dap-16},
        CosSin{.Angle  = 0x1.67fef9db22d0ep+6,
               .Cosine = 0x1.24d15211596efp-16,
               .Sine   = 0x1.fffffffeb111dp-1},
        CosSin{.Angle = 0x1.b7cdfd9d7bdbbp-34, .Cosine = 0x1p+0, .Sine = 0x1.eb443a0930a1fp-40},
        CosSin{.Angle  = 0x1.edd2f1a9fbe77p+6,
               .Cosine = -0x1.1a43870d371f9p-1,
               .Sine   = 0x1.ab2a9b25a9ab9p-1},
        CosSin{.Angle = -0x1.0f8p+8, .Cosine = 0x1.ace214390ca91p-6, .Sine = 0x1.ffd315bbf4275p-1},
        CosSin{
            .Angle = 0x1.e84804p+19, .Cosine = 0x1.6807e1489cb33p-3, .Sine = -0x1.f806c25684ea8p-1},
        CosSin{.Angle = 0x1.68p+5, .Cosine = 0x1.6a09e667f3bcdp-1, .Sine = 0x1.6a09e667f3bcdp-1},
    };

    /**
     *  In radians, the last at the double nearest a multiple of a quarter turn of all, which is
     *  past 16 and so not taken for one.
     */
    constexpr std::array InRadians{
        CosSin{.Angle =    0x1p+0, .Cosine =  0x1.14a280fb5068cp-1, .Sine =  0x1.aed548f090ceep-1},
        CosSin{.Angle =    0x1p-1, .Cosine =   0x1.c1528065b7d5p-1, .Sine =   0x1.eaee8744b05fp-2},
        CosSin{.Angle = -0x1.8p+1, .Cosine = -0x1.fae04be85e5d2p-1, .Sine = -0x1.210386db6d55bp-3},
        CosSin{
            .Angle = 0x1.e848p+19, .Cosine = 0x1.df9df9906d32cp-1, .Sine = -0x1.6664b2568d867p-2},
        CosSin{.Angle  = 0x1.0f0cf064dd592p+73,
               .Cosine = 0x1.0be2cef01c8f4p-1,
               .Sine   = -0x1.b453ab76bf397p-1},
        CosSin{.Angle = 0x1.6ac5b262ca1ffp+849, .Cosine = -0x1.14ae72e6ba22fp-61, .Sine = 0x1p+0},
    };

    /**
     *  A heading in radians, the double nearest the exact one, from MPFR.
     */
    struct Headed
    {
      double X       = 0.0; /**< Along the x axis. */
      double Y       = 0.0; /**< Along the y axis. */
      double Radians = 0.0; /**< The heading.      */
    };

    constexpr std::array Headings{
        Headed{.X = 0x1.8p+1, .Y =                0x1p+2, .Radians =  0x1.dac670561bb4fp-1},
        Headed{.X =   0x1p+0, .Y = 0x1.19799812dea11p-40, .Radians = 0x1.19799812dea11p-40},
        // NOLINTNEXTLINE(modernize-use-std-numbers): just short of pi, not pi
        Headed{.X =   -0x1p+0, .Y = 0x1.19799812dea11p-40, .Radians =  0x1.921fb5444244cp+1},
        Headed{.X = -0x1.4p+2, .Y =             -0x1.8p+3, .Radians = -0x1.f730bd281f69bp+0},
        Headed{.X       = 0x1.56e1fc2f8f359p-997,
               .Y       = 0x1.01297d23ab683p-995,
               .Radians = 0x1.3fc176b7a856p+0},
        Headed{.X = 0x1p+0, .Y = 0x1p+0, .Radians = 0x1.921fb54442d18p-1},
    };

    /**
     *  A double and the double nearest its product with pi / 180 or 180 / pi, from exact rationals:
     *  the first the nearest a halfway point any double comes (tests/proof/convert.py), the last
     *  rounded to the subnormal grid.
     */
    struct Converted
    {
      double Given   = 0.0; /**< The double.                 */
      double Nearest = 0.0; /**< The nearest to its product. */
    };

    constexpr std::array ToRadians{
        Converted{.Given =    0x1.96bdf4aa9cd3bp+0, .Nearest =    0x1.c655cf14d66cbp-6},
        Converted{.Given =    0x1.f2fcf8a18a34ep+0, .Nearest =    0x1.16b00466ffef2p-5},
        Converted{.Given =    0x1.a860000000000p+0, .Nearest =    0x1.da080b755b362p-6},
        Converted{.Given =    0x1.ce00000000000p+0, .Nearest =    0x1.02078bc788bdfp-5},
        Converted{.Given = 0x0.00000000065afp-1022, .Nearest = 0x0.00000000001c6p-1022},
        Converted{.Given = 0x0.000000000001ep-1022, .Nearest = 0x0.0000000000001p-1022},
        Converted{.Given = 0x0.0000000000001p-1022, .Nearest =                0x0.0p+0},
        Converted{.Given = 0x1.0000000000000p-1022, .Nearest = 0x0.0477d1a894a75p-1022},
        Converted{.Given = 0x1.8000000000000p-1000, .Nearest = 0x1.acee9f37bebd6p-1006},
    };

    constexpr std::array ToDegrees{
        Converted{.Given =    0x1.db0fb301078bep+0, .Nearest =    0x1.a94c09279849fp+6},
        Converted{.Given =    0x1.b19beb28a54d9p+0, .Nearest =    0x1.842fe02da15d7p+6},
        Converted{.Given =    0x1.88282350430f4p+0, .Nearest =    0x1.5f13b733aa70ep+6},
        Converted{.Given =    0x1.2eb176c1f58fcp+0, .Nearest =    0x1.0efc2cf3bf3dcp+6},
        Converted{.Given = 0x0.00000000065afp-1022, .Nearest = 0x0.000000016c20ap-1022},
        Converted{.Given = 0x0.000000000001ep-1022, .Nearest = 0x0.00000000006b7p-1022},
        Converted{.Given = 0x0.0000000000001p-1022, .Nearest = 0x0.0000000000039p-1022},
        Converted{.Given = 0x1.0000000000000p-1022, .Nearest = 0x1.ca5dc1a63c1f8p-1017},
        Converted{.Given = 0x1.8000000000000p-1000, .Nearest =  0x1.57c6513cad17ap-994},
    };

    /**
     *  Degrees whose cosine or sine is within 2^-17 of a halfway point between doubles, and their
     *  negatives: the fast path cannot round them, from MPFR. The first five give one angle's
     *  cosine and sine, turned.
     */
    constexpr std::array HardDegrees{
        CosSin{.Angle  = 0x1.fab126e978d5p+6,
               .Cosine = -0x1.31ca64547ad52p-1,
               .Sine   = 0x1.9aa73d164efb7p-1},
        CosSin{.Angle  = -0x1.fab126e978d5p+6,
               .Cosine = -0x1.31ca64547ad52p-1,
               .Sine   = -0x1.9aa73d164efb7p-1},
        CosSin{.Angle  = 0x1.cc229567f734ep+6,
               .Cosine = -0x1.b14eddb32c778p-2,
               .Sine   = 0x1.cfe6e6e39cda6p-1},
        CosSin{.Angle  = -0x1.cc229567f734ep+6,
               .Cosine = -0x1.b14eddb32c778p-2,
               .Sine   = -0x1.cfe6e6e39cda6p-1},
        CosSin{.Angle  = 0x1.32ac49ba5e354p+8,
               .Cosine = 0x1.31ca64547ad52p-1,
               .Sine   = -0x1.9aa73d164efb7p-1},
        CosSin{.Angle  = -0x1.32ac49ba5e354p+8,
               .Cosine = 0x1.31ca64547ad52p-1,
               .Sine   = 0x1.9aa73d164efb7p-1},
        CosSin{.Angle  = -0x1.b0ffc425f238cp+7,
               .Cosine = -0x1.9b93bf6f23a84p-1,
               .Sine   = 0x1.308bc5c82622p-1},
        CosSin{.Angle  = 0x1.b0ffc425f238cp+7,
               .Cosine = -0x1.9b93bf6f23a84p-1,
               .Sine   = -0x1.308bc5c82622p-1},
        CosSin{.Angle  = 0x1.d2a76c8b43958p+7,
               .Cosine = -0x1.31ca64547ad52p-1,
               .Sine   = -0x1.9aa73d164efb7p-1},
        CosSin{.Angle  = -0x1.d2a76c8b43958p+7,
               .Cosine = -0x1.31ca64547ad52p-1,
               .Sine   = 0x1.9aa73d164efb7p-1},
        CosSin{.Angle  = 0x1.b1589374bc6a8p+7,
               .Cosine = -0x1.9aa73d164efb7p-1,
               .Sine   = -0x1.31ca64547ad52p-1},
        CosSin{.Angle  = -0x1.b1589374bc6a8p+7,
               .Cosine = -0x1.9aa73d164efb7p-1,
               .Sine   = 0x1.31ca64547ad52p-1},
        CosSin{.Angle  = 0x1.c8c49ba5e353fp+2,
               .Cosine = 0x1.fc0871544396ap-1,
               .Sine   = 0x1.fce576a6b33dap-4},
        CosSin{.Angle  = -0x1.c8c49ba5e353fp+2,
               .Cosine = 0x1.fc0871544396ap-1,
               .Sine   = -0x1.fce576a6b33dap-4},
        CosSin{.Angle  = -0x1.4353b645a1cacp+8,
               .Cosine = 0x1.9aa73d164efb7p-1,
               .Sine   = 0x1.31ca64547ad52p-1},
        CosSin{.Angle  = 0x1.4353b645a1cacp+8,
               .Cosine = 0x1.9aa73d164efb7p-1,
               .Sine   = -0x1.31ca64547ad52p-1},
        CosSin{.Angle  = 0x1.a11a9fbe76c8bp+7,
               .Cosine = -0x1.c1bb81f8bb3cap-1,
               .Sine   = -0x1.e96d5673e1807p-2},
        CosSin{.Angle  = -0x1.a11a9fbe76c8bp+7,
               .Cosine = -0x1.c1bb81f8bb3cap-1,
               .Sine   = 0x1.e96d5673e1807p-2},
        CosSin{.Angle  = 0x1.57127057d3aa9p+8,
               .Cosine = 0x1.e9d0da313fa13p-1,
               .Sine   = -0x1.2a2848f56d2dbp-2},
        CosSin{.Angle  = -0x1.57127057d3aa9p+8,
               .Cosine = 0x1.e9d0da313fa13p-1,
               .Sine   = 0x1.2a2848f56d2dbp-2},
    };

    /**
     *  Radians likewise, past 2^50 three of them.
     */
    constexpr std::array HardRadians{
        CosSin{.Angle  = -0x1.1ea7c1abd20eep+1,
               .Cosine = -0x1.3d6c1016965dcp-1,
               .Sine   = -0x1.91baf1f66b7f1p-1},
        CosSin{.Angle  = 0x1.1ea7c1abd20eep+1,
               .Cosine = -0x1.3d6c1016965dcp-1,
               .Sine   = 0x1.91baf1f66b7f1p-1},
        CosSin{.Angle  = 0x1.e2ff856178049p+2,
               .Cosine = 0x1.3595f02010587p-2,
               .Sine   = 0x1.e80a4659e5485p-1},
        CosSin{.Angle  = -0x1.e2ff856178049p+2,
               .Cosine = 0x1.3595f02010587p-2,
               .Sine   = -0x1.e80a4659e5485p-1},
        CosSin{.Angle  = -0x1.c2e667b0e3e7ap+2,
               .Cosine = 0x1.725d618a9ad4dp-1,
               .Sine   = -0x1.61843b59f3b04p-1},
        CosSin{.Angle  = 0x1.c2e667b0e3e7ap+2,
               .Cosine = 0x1.725d618a9ad4dp-1,
               .Sine   = 0x1.61843b59f3b04p-1},
        CosSin{.Angle  = -0x1.1cabf5cf9b41bp+3,
               .Cosine = -0x1.ba12287f58751p-1,
               .Sine   = -0x1.024c02225e2bp-1},
        CosSin{.Angle  = 0x1.1cabf5cf9b41bp+3,
               .Cosine = -0x1.ba12287f58751p-1,
               .Sine   = 0x1.024c02225e2bp-1},
        CosSin{.Angle  = 0x1.15368f8ca94b6p+3,
               .Cosine = -0x1.7274e37847a56p-1,
               .Sine   = 0x1.616b98f986ceep-1},
        CosSin{.Angle  = -0x1.15368f8ca94b6p+3,
               .Cosine = -0x1.7274e37847a56p-1,
               .Sine   = -0x1.616b98f986ceep-1},
        CosSin{.Angle  = -0x1.574f4c22d402ep+2,
               .Cosine = 0x1.3699633ac851dp-1,
               .Sine   = 0x1.9707483296038p-1},
        CosSin{.Angle  = 0x1.574f4c22d402ep+2,
               .Cosine = 0x1.3699633ac851dp-1,
               .Sine   = -0x1.9707483296038p-1},
        CosSin{.Angle  = 0x1.2ed5dcaaaa447p+1,
               .Cosine = -0x1.6d893025b0beap-1,
               .Sine   = 0x1.6681e0ffcb768p-1},
        CosSin{.Angle  = -0x1.2ed5dcaaaa447p+1,
               .Cosine = -0x1.6d893025b0beap-1,
               .Sine   = -0x1.6681e0ffcb768p-1},
        CosSin{.Angle  = -0x1.3d383ef88d69ap+3,
               .Cosine = -0x1.c4277e8f71c32p-1,
               .Sine   = 0x1.e06b2a6d5e482p-2},
        CosSin{.Angle  = 0x1.3d383ef88d69ap+3,
               .Cosine = -0x1.c4277e8f71c32p-1,
               .Sine   = -0x1.e06b2a6d5e482p-2},
        CosSin{.Angle  = 0x1.9e8e3b48d83f6p+2,
               .Cosine = 0x1.f65eff6dd65bap-1,
               .Sine   = 0x1.8b5186bfd176p-3},
        CosSin{.Angle  = -0x1.9e8e3b48d83f6p+2,
               .Cosine = 0x1.f65eff6dd65bap-1,
               .Sine   = -0x1.8b5186bfd176p-3},
        CosSin{.Angle  = 0x1.28e6b30562b36p+989,
               .Cosine = -0x1.f305335cf3d5ap-1,
               .Sine   = -0x1.ca37c56b045dp-3},
        CosSin{.Angle  = -0x1.28e6b30562b36p+989,
               .Cosine = -0x1.f305335cf3d5ap-1,
               .Sine   = 0x1.ca37c56b045dp-3},
        CosSin{.Angle  = 0x1.4423fb76c3ce6p+215,
               .Cosine = -0x1.54a357978654p-4,
               .Sine   = -0x1.fe39f463d13b1p-1},
        CosSin{.Angle  = -0x1.4423fb76c3ce6p+215,
               .Cosine = -0x1.54a357978654p-4,
               .Sine   = 0x1.fe39f463d13b1p-1},
        CosSin{.Angle  = 0x1.823592c86b3dp+664,
               .Cosine = 0x1.f8a080ef4266dp-1,
               .Sine   = -0x1.5a4f200713486p-3},
        CosSin{.Angle  = -0x1.823592c86b3dp+664,
               .Cosine = 0x1.f8a080ef4266dp-1,
               .Sine   = 0x1.5a4f200713486p-3},
        /**
         *  From the suites, its sine within 2^-16.4 of a halfway point.
         */
        CosSin{.Angle  = 0x1.8db958179047ep+1,
               .Cosine = -0x1.ffb292948d298p-1,
               .Sine   = 0x1.198919387535p-5},
    };

    /**
     *  Headings likewise, two in each octant, then four from the suites: 1 - 2^-50 across, the
     *  heading within 2^-20.8 of a halfway point.
     */
    constexpr std::array HardHeadings{
        Headed{
            .X = 0x1.781ed6cd4e78p-1, .Y = 0x1.27ee5ba2782b8p-2, .Radians = 0x1.7fcc6b438fb11p-2},
        Headed{
            .X = 0x1.781ed6cd4e78p-1, .Y = -0x1.27ee5ba2782b8p-2, .Radians = -0x1.7fcc6b438fb11p-2},
        Headed{
            .X = 0x1.3d3a2d26daa5p-3, .Y = -0x1.04e2c14d954fp-1, .Radians = -0x1.469320a2421f2p+0},
        Headed{.X = 0x1.3d3a2d26daa5p-3, .Y = 0x1.04e2c14d954fp-1, .Radians = 0x1.469320a2421f2p+0},
        Headed{
            .X = -0x1.b9e1cc82409bp-2, .Y = -0x1.f95db3d052eap-5, .Radians = -0x1.7ff2c5fd22a78p+1},
        Headed{
            .X = -0x1.b9e1cc82409bp-2, .Y = 0x1.f95db3d052eap-5, .Radians = 0x1.7ff2c5fd22a78p+1},
        Headed{.X       = -0x1.5864ee10ae628p-3,
               .Y       = -0x1.598f13b24db83p-1,
               .Radians = -0x1.d0a2a60d9e803p+0},
        Headed{
            .X = -0x1.5864ee10ae628p-3, .Y = 0x1.598f13b24db83p-1, .Radians = 0x1.d0a2a60d9e803p+0},
        Headed{
            .X = -0x1.b24e0d7936099p-1, .Y = -0x1.d91cc671b8c2p-5, .Radians = -0x1.896c2555b806p+1},
        Headed{
            .X = -0x1.b24e0d7936099p-1, .Y = 0x1.d91cc671b8c2p-5, .Radians = 0x1.896c2555b806p+1},
        Headed{.X       = -0x1.0f1a5f8dde5f2p-1,
               .Y       = -0x1.fb0ad9b0fccc2p-1,
               .Radians = -0x1.07e8fdbdb6c94p+1},
        Headed{
            .X = -0x1.0f1a5f8dde5f2p-1, .Y = 0x1.fb0ad9b0fccc2p-1, .Radians = 0x1.07e8fdbdb6c94p+1},
        Headed{.X       = 0x1.ffffffffffff8p-1,
               .Y       = 0x1.7ffffe33094c9p-25,
               .Radians = 0x1.7ffffe33094cbp-25},
        Headed{.X       = 0x1.ffffffffffff8p-1,
               .Y       = -0x1.7ffffe33094c9p-25,
               .Radians = -0x1.7ffffe33094cbp-25},
        Headed{.X       = 0x1.ffffffffffff8p-1,
               .Y       = 0x1.7ffffe469898cp-25,
               .Radians = 0x1.7ffffe469898ep-25},
        Headed{.X       = 0x1.ffffffffffff8p-1,
               .Y       = -0x1.7ffffe469898cp-25,
               .Radians = -0x1.7ffffe469898ep-25},
    };

    /**
     *  Subnormal headings, which the fast path leaves to the slow one: the first within 2^-17 of a
     *  halfway point, the rest anywhere; then three just above the subnormal range, whose pairs sit
     *  on a halfway point, which a test whose error bound underflowed let through wrong.
     */
    constexpr std::array SubnormalHeadings{
        Headed{
            .X = 0x1.e89b001c34dd7p+2, .Y = 0x1.77b5562bac8p-1031, .Radians = 0x1.89b284c3fbp-1034},
        Headed{.X       = 0x1.e89b001c34dd7p+2,
               .Y       = -0x1.77b5562bac8p-1031,
               .Radians = -0x1.89b284c3fbp-1034},
        Headed{.X = 0x1.10ac0b5d2e5e3p+1, .Y =   -0x1p-1074, .Radians =      -0x0p+0},
        Headed{.X = 0x1.2d56082460113p+1, .Y = 0x1.e6p-1064, .Radians = 0x1.9dp-1065},
        Headed{
            .X = 0x1.bb2b7a2f0b906p+1, .Y = -0x1.8786716ep-1037, .Radians = -0x1.c455903cap-1039},
        Headed{.X =  0x1.9a3cddb5a713p+0, .Y =     0x1.5ba1p-1057, .Radians =    0x1.b1dcp-1058},
        Headed{.X = 0x1.e7a3d636b8a83p+0, .Y =       -0x1.1p-1070, .Radians =      -0x1.2p-1071},
        Headed{.X = 0x1.3e3089e4cc001p+0, .Y = 0x1.c89dd6fcp-1044, .Radians = 0x1.6f5f247p-1044},
        Headed{.X       = 0x1.55c9c9b5be36bp+2,
               .Y       = -0x1.d13af74192d1cp-1019,
               .Radians = -0x1.5c7556c2657e9p-1021},
        Headed{.X       = 0x1.2c21e531c67p+1,
               .Y       = 0x1.1c17bb5a1715p-1019,
               .Radians = 0x1.e4a33d3e524adp-1021},
        Headed{.X       = 0x1.30a4224d10d9fp+0,
               .Y       = 0x1.39feb6a14a3cp-1021,
               .Radians = 0x1.07dc3e7812dd5p-1021},
    };

    /**
     *  The heading of one displacement from another, likewise, and with the two swapped.
     */
    struct Apart
    {
      Vector First   = {};  /**< The first displacement.  */
      Vector Second  = {};  /**< The second displacement. */
      double Radians = 0.0; /**< The heading between.     */
    };

    constexpr std::array HardBetweens{
        Apart{.First   = {.X = -0x1.8b10729e9d4a6p-1, .Y = 0x1.7e9a2dca46f1ep-1},
              .Second  = {.X = 0x1.76b084901ca3p-1, .Y = -0x1.c0321ca6f9fc8p-4},
              .Radians = -0x1.42a45b55164cep+1},
        Apart{.First   = {.X = 0x1.76b084901ca3p-1, .Y = -0x1.c0321ca6f9fc8p-4},
              .Second  = {.X = -0x1.8b10729e9d4a6p-1, .Y = 0x1.7e9a2dca46f1ep-1},
              .Radians = 0x1.42a45b55164cep+1},
        Apart{.First   = {.X = 0x1.bf4975a37bfcp-4, .Y = -0x1.52557f2b5a1b9p-1},
              .Second  = {.X = -0x1.d2330c375c066p-2, .Y = -0x1.77a656f1b53d2p-2},
              .Radians = -0x1.0e6afcf54bee2p+0},
        Apart{.First   = {.X = -0x1.d2330c375c066p-2, .Y = -0x1.77a656f1b53d2p-2},
              .Second  = {.X = 0x1.bf4975a37bfcp-4, .Y = -0x1.52557f2b5a1b9p-1},
              .Radians = 0x1.0e6afcf54bee2p+0},
        Apart{.First   = {.X = -0x1.5c3e27f7c6d6p-2, .Y = -0x1.04a4c88c0cbc4p-1},
              .Second  = {.X = -0x1.c4fa2e88e51bdp-1, .Y = 0x1.303d5d2eb25f4p-1},
              .Radians = -0x1.92c298fdf0b1bp+0},
        Apart{.First   = {.X = -0x1.c4fa2e88e51bdp-1, .Y = 0x1.303d5d2eb25f4p-1},
              .Second  = {.X = -0x1.5c3e27f7c6d6p-2, .Y = -0x1.04a4c88c0cbc4p-1},
              .Radians = 0x1.92c298fdf0b1bp+0},
        Apart{.First   = {.X = -0x1.a5615e3155606p-1, .Y = 0x1.d8ca5b08a863p-4},
              .Second  = {.X = 0x1.cbcfed06560fp-1, .Y = 0x1.8af7e94042138p-2},
              .Radians = -0x1.4c5ce13d26daep+1},
        Apart{.First   = {.X = 0x1.cbcfed06560fp-1, .Y = 0x1.8af7e94042138p-2},
              .Second  = {.X = -0x1.a5615e3155606p-1, .Y = 0x1.d8ca5b08a863p-4},
              .Radians = 0x1.4c5ce13d26daep+1},
        Apart{.First   = {.X = -0x1.aa0a5a9e5ea2cp-1, .Y = -0x1.2cb24f80627ffp-1},
              .Second  = {.X = 0x1.c680cadbbb9e4p-1, .Y = 0x1.3f916dfbcf422p-1},
              .Radians = 0x1.91e50eef3cb58p+1},
        Apart{.First   = {.X = 0x1.c680cadbbb9e4p-1, .Y = 0x1.3f916dfbcf422p-1},
              .Second  = {.X = -0x1.aa0a5a9e5ea2cp-1, .Y = -0x1.2cb24f80627ffp-1},
              .Radians = -0x1.91e50eef3cb58p+1},
        Apart{.First   = {.X = -0x1.7e0f087f06b3p-5, .Y = 0x1.bde1e08c324cp-5},
              .Second  = {.X = 0x1.614fc3b0c53ecp-2, .Y = -0x1.153b8b141dbf2p-1},
              .Radians = 0x1.800f22f42ced9p+1},
        Apart{.First   = {.X = 0x1.614fc3b0c53ecp-2, .Y = -0x1.153b8b141dbf2p-1},
              .Second  = {.X = -0x1.7e0f087f06b3p-5, .Y = 0x1.bde1e08c324cp-5},
              .Radians = -0x1.800f22f42ced9p+1},
    };

    /**
     *  A number at a precision: the sum of 2^place for every place in the ranges, truncated.
     *
     *  @param  bits     The precision.
     *  @param  ranges   Each [low, high) of places.
     *  @param  negative Whether the number is below zero.
     *  @return          The number.
     */
    auto Ones(int bits, std::initializer_list<std::array<int, 2>> ranges, bool negative = false)
        -> Precise::Fixed
    {
      Precise::Fixed value = {.Negative = negative};
      for (const auto &[low, high] : ranges)
      {
        for (int place = std::max(low + bits, 0); place < high + bits; ++place)
        {
          const auto bit = static_cast<std::size_t>(place);
          value.Words.resize(std::max(value.Words.size(), (bit / 32) + 1));
          value.Words[bit / 32] |= 1U << (bit % 32);
        }
      }
      return value;
    }
  } // namespace

  TEST_CASE("angle: the conversions are the nearest double, at the hardest cases and subnormal")
  {
    for (const Converted &converted : ToRadians)
    {
      CHECK(Degrees(converted.Given).InRadians()  == converted.Nearest);
      CHECK(Degrees(-converted.Given).InRadians() == -converted.Nearest);
    }
    for (const Converted &converted : ToDegrees)
    {
      CHECK(Radians(converted.Given).InDegrees()  == converted.Nearest);
      CHECK(Radians(-converted.Given).InDegrees() == -converted.Nearest);
    }
    CHECK(Degrees(Infinite).InRadians()  == Infinite);
    CHECK(Radians(-Infinite).InDegrees() == -Infinite);
    CHECK(std::isnan(Degrees(std::numeric_limits<double>::quiet_NaN()).InRadians()));
    CHECK(Same(Radians(-0.0).InDegrees(), -0.0));
    CHECK(Radians(1e308).InDegrees() == Infinite);
    /**
     *  Just below a power of two the grid is twice as fine.
     */
    CHECK(Degrees(Radians(1).InDegrees()).InRadians() == 1.0);
  }

  TEST_CASE("angle: each unit is held as given, and converted correctly rounded")
  {
    int kept = 0;
    for (int whole = -720; whole <= 720; ++whole)
    {
      kept += Same(Degrees(whole).InDegrees(), whole) ? 1 : 0;
    }
    CHECK(kept == 1441);
    CHECK(Degrees(90).InRadians()       == Pi / 2.0);
    CHECK(Degrees(180).InRadians()      == Pi);
    CHECK(Degrees(360).InRadians()      == 2.0 * Pi);
    CHECK(Degrees(-90).InRadians()      == -Pi / 2.0);
    CHECK(Radians(Pi).InDegrees()       == 180.0);
    CHECK(Radians(Pi / 2.0).InDegrees() == 90.0);
    CHECK(Radians(2.0 * Pi).InDegrees() == 360.0);
    CHECK(Radians(1).InDegrees()        == 0x1.ca5dc1a63c1f8p+5);
    CHECK(Degrees(1).InRadians()        == 0x1.1df46a2529d39p-6);
    CHECK(Radians(0.25).InRadians()     == 0.25);
    CHECK(Turn             == Radians(2.0 * Pi));
    CHECK(Turn.InRadians() == 2.0 * Pi);
  }

  TEST_CASE("angle: the arithmetic keeps a shared unit, and mixes into degrees")
  {
    CHECK((Degrees(30) + Degrees(60)).InDegrees() == 90.0);
    CHECK((Radians(1) + Radians(2)).InRadians()   == 3.0);
    CHECK_FALSE((Radians(1) + Radians(2)).HeldInDegrees());
    CHECK((Degrees(90) + Radians(Pi / 2.0)) == Degrees(180));
    CHECK((Degrees(90) + Radians(Pi / 2.0)).HeldInDegrees());
    CHECK((-Radians(1)).InRadians() == -1.0);
    CHECK(Radians(1) < Degrees(58));
    CHECK(Degrees(57) < Radians(1));
  }

  TEST_CASE("angle: angles are ordered by their size, in one unit or across two")
  {
    const Angle less = Degrees(10);
    const Angle more = Degrees(20);
    CHECK(less < more);
    CHECK(Radians(1) < Radians(2));
    CHECK(Degrees(10) < Radians(1));
    CHECK_FALSE(more < less);
  }

  TEST_CASE("angle: wrapped within (-180, 180] exactly in degrees, nearest from radians")
  {
    CHECK(Wrapped(Degrees(-180)).InDegrees()   == 180.0);
    CHECK(Wrapped(Degrees(540)).InDegrees()    == 180.0);
    CHECK(Wrapped(Degrees(-540)).InDegrees()   == 180.0);
    CHECK(Wrapped(Degrees(270)).InDegrees()    == -90.0);
    CHECK(Wrapped(Degrees(-1e-20)).InDegrees() == -1e-20);
    CHECK(Same(Wrapped(Degrees(-360)).InDegrees(), 0.0));
    CHECK(Wrapped(Radians(Pi)).InDegrees()        == 180.0);
    CHECK(Wrapped(Radians(-Pi / 2.0)).InDegrees() == -90.0);
    CHECK(Wrapped(Radians(4.0 * Pi)).InDegrees()  == 0.0);
    CHECK(Wrapped(Radians(7.0)).InDegrees()
          == doctest::Approx((7.0 * 180.0 / Pi) - 360.0).epsilon(1e-15));
    const double rest = std::fmod(1e300, 360.0);
    CHECK(Wrapped(Degrees(1e300)).InDegrees() == (rest > 180.0 ? rest - 360.0 : rest));

    Draws draws{1};
    int   wrong = 0;
    for (int index = 0; index < 200'000; ++index)
    {
      const double degrees = std::ldexp(std::round(draws.Near(0, 0x1p50)), -40);
      const double wrapped = Wrapped(Degrees(degrees)).InDegrees();
      const auto   given   = Fixed(degrees);
      const auto   held    = Fixed(wrapped);
      const bool   good =
          given && held && Reduced(*given - *held) == 0 && wrapped > -180.0 && wrapped <= 180.0;
      wrong += good ? 0 : 1;
    }
    CHECK(wrong == 0);
  }

  TEST_CASE("angle: within a sweep, exactly, the exact end included")
  {
    CHECK(Within(Degrees(30), Degrees(10), Degrees(20)));
    CHECK_FALSE(Within(Degrees(Nudged(30, 1)), Degrees(10), Degrees(20)));
    CHECK(Within(Degrees(350), Degrees(10), Degrees(-20)));
    CHECK(Within(Degrees(370), Degrees(10), Degrees(-0.0)));
    CHECK(Within(Degrees(123), Degrees(4), Degrees(-360)));
    CHECK_FALSE(Within(Degrees(std::numeric_limits<double>::quiet_NaN()), Degrees(0), Degrees(10)));
    CHECK(Within(Radians(Pi / 2.0), Degrees(0), Degrees(90)));
    CHECK(Within(Degrees(0.1 + 89.9), Degrees(0.1), Degrees(89.9)));
    CHECK(Within(Radians(0.25 + 1.25), Radians(0.25), Radians(1.25)));
    CHECK_FALSE(Within(Radians(Nudged(1.5, 1)), Radians(0.25), Radians(1.25)));
    CHECK_FALSE(Within(Radians(0.1 + 1.4), Radians(0.1), Radians(1.4))); // 1.5, past the exact end
    CHECK(Within(Radians(0.3), Radians(0.3), Radians(-0.0)));
    CHECK(Within(Radians(1e22), Radians(1e22), Radians(1e-300)));
    CHECK(Within(Radians(-1e-17), Radians(0), Radians(2.0 * Pi)));
    CHECK(Within(Radians(5), Radians(1), Radians(-Pi / 2.0 * 3.0)));

    Draws draws{2};
    int   wrong   = 0;
    int   checked = 0;
    for (int index = 0; index < 300'000; ++index)
    {
      const double start = index % 3 == 0 ? std::round(draws.Near(0, 720)) : draws.Near(0, 720);
      const double sweep = index % 5 == 0 ? std::round(draws.Near(0, 400)) : draws.Near(0, 400);
      const int    ulps  = static_cast<int>(draws.Next() % 9U) - 4;
      const double turns = std::round(draws.Near(0, 2.5));
      const double direction =
          index % 2 == 0 ? Nudged(start + sweep, ulps) + (360.0 * turns) : draws.Near(0, 720);
      const auto toward = Fixed(direction);
      const auto from   = Fixed(start);
      const auto by     = Fixed(sweep);
      if (!toward || !from || !by)
      {
        continue;
      }
      ++checked;
      const bool answer  = Within(Degrees(direction), Degrees(start), Degrees(sweep));
      wrong             += answer == Truly(*toward, *from, *by) ? 0 : 1;
    }
    CHECK(checked > 100'000);
    CHECK(wrong == 0);
  }

  TEST_CASE("angle: a direction in degrees is correctly rounded")
  {
    for (const CosSin &turned : InDegrees)
    {
      const Vector way = Direction(Degrees(turned.Angle));
      CHECK(way.X == turned.Cosine);
      CHECK(way.Y == turned.Sine);
    }
    for (int turn = -4; turn <= 4; ++turn)
    {
      const double base = 360.0 * turn;
      CHECK(Direction(Degrees(base))           == Vector{.X = 1, .Y = 0});
      CHECK(Direction(Degrees(base + 90.0))    == Vector{.X = 0, .Y = 1});
      CHECK(Direction(Degrees(base + 180.0))   == Vector{.X = -1, .Y = 0});
      CHECK(Direction(Degrees(base + 270.0))   == Vector{.X = 0, .Y = -1});
      CHECK(Direction(Degrees(base + 30.0)).Y  == 0.5);
      CHECK(Direction(Degrees(base + 60.0)).X  == 0.5);
      CHECK(Direction(Degrees(base + 150.0)).Y == 0.5);
      CHECK(Direction(Degrees(base + 240.0)).X == -0.5);
    }
    CHECK(Same(Direction(Degrees(180)).Y, 0.0));
    CHECK(std::isnan(Direction(Degrees(Infinite)).X));
  }

  TEST_CASE("angle: a direction in degrees keeps its symmetries bit for bit")
  {
    Draws draws{3};
    int   broken = 0;
    for (int index = 0; index < 200'000; ++index)
    {
      const double degrees     = draws.Near(0, 720);
      const Vector way         = Direction(Degrees(degrees));
      const Vector back        = Direction(Degrees(-degrees));
      broken                  += Same(back.X, way.X) && Same(back.Y, -way.Y + 0.0) ? 0 : 1;
      const double complement  = 90.0 - degrees;
      if (std::equal_to<double>{}(90.0 - complement, degrees))
      {
        broken += Same(Direction(Degrees(complement)).Y, way.X) ? 0 : 1;
      }
      broken += std::abs(std::hypot(way.X, way.Y) - 1.0) <= 0x1p-51 ? 0 : 1;
    }
    CHECK(broken == 0);
  }

  TEST_CASE("angle: along a sweep in degrees, the end is the exact end")
  {
    CHECK(Direction(Degrees(10.25), Degrees(79.75), 1.0) == Vector{.X = 0, .Y = 1});
    CHECK(Direction(Degrees(10), Degrees(20), 1.0).Y     == 0.5);
    CHECK(Direction(Degrees(1e20), Degrees(-1e20), 1.0)  == Vector{.X = 1, .Y = 0});
    CHECK(Direction(Degrees(0.1), Degrees(89.9), 1.0).X < 0.0);
    CHECK(Direction(Radians(1), Radians(1), 0.5)  == Direction(Radians(1.5)));
    CHECK(Direction(Degrees(10), Radians(1), 0.5) == Direction(Degrees(10) + Radians(0.5)));
  }

  TEST_CASE("angle: a direction in radians is correctly rounded, past any number of turns")
  {
    for (const CosSin &turned : InRadians)
    {
      const Vector way = Direction(Radians(turned.Angle));
      CHECK(way.X == turned.Cosine);
      CHECK(way.Y == turned.Sine);
    }
    CHECK(Direction(Radians(0)) == Vector{.X = 1, .Y = 0});
    CHECK(Same(Direction(Radians(-0.0)).Y, 0.0));
    CHECK(Direction(Radians(1e-300)) == Vector{.X = 1, .Y = 1e-300});
    CHECK(std::isnan(Direction(Radians(Infinite)).Y));

    Draws draws{4};
    int   broken = 0;
    for (int index = 0; index < 100'000; ++index)
    {
      const double radians  = draws.Near(0, 1e3);
      const Vector way      = Direction(Radians(radians));
      const Vector back     = Direction(Radians(-radians));
      broken               += Same(back.X, way.X) && Same(back.Y, -way.Y + 0.0) ? 0 : 1;
      broken               += std::abs(way.X - std::cos(radians)) <= 0x1p-52 ? 0 : 1;
      broken               += std::abs(way.Y - std::sin(radians)) <= 0x1p-52 ? 0 : 1;
    }
    CHECK(broken == 0);
  }

  TEST_CASE("angle: in radians the double nearest a quarter turn is that quarter turn")
  {
    for (int power = 0; power <= 3; ++power)
    {
      for (const double sign : {1.0, -1.0})
      {
        const double radians = std::ldexp(sign * Pi / 2.0, power);
        const Vector axis    = Direction(Degrees(std::ldexp(sign * 90.0, power)));
        CHECK(Direction(Radians(radians)) == axis);
        CHECK(Direction(Radians(Nudged(radians, 1))) != axis);
      }
    }
    CHECK(Direction(Radians(1.5 * Pi)) == Vector{.X = 0, .Y = -1});
    CHECK(Direction(Radians(2.0 * Pi)) == Vector{.X = 1, .Y = 0});
    CHECK(Direction(Radians(Pi))       == Vector{.X = -1, .Y = 0});
    CHECK(Direction(Radians(-Pi / 2))  == Vector{.X = 0, .Y = -1});
  }

  TEST_CASE("angle: a heading is correctly rounded in radians, exact on axes and diagonals")
  {
    for (const Headed &headed : Headings)
    {
      CHECK(Heading({.X = headed.X, .Y = headed.Y}).InRadians() == headed.Radians);
    }
    for (const double size : {1e-300, 0.1, 1.0, 3.0, 1e300})
    {
      CHECK(Heading({.X = size, .Y = 0}).InRadians()      == 0.0);
      CHECK(Heading({.X = -size, .Y = -0.0}).InRadians()  == Pi);
      CHECK(Heading({.X = 0, .Y = size}).InRadians()      == Pi / 2.0);
      CHECK(Heading({.X = 0, .Y = -size}).InRadians()     == -Pi / 2.0);
      CHECK(Heading({.X = size, .Y = size}).InDegrees()   == 45.0);
      CHECK(Heading({.X = -size, .Y = size}).InDegrees()  == 135.0);
      CHECK(Heading({.X = -size, .Y = -size}).InDegrees() == -135.0);
    }
    CHECK(Heading({}).InRadians() == 0.0);
    CHECK_FALSE(Heading({.X = 1, .Y = 2}).HeldInDegrees());
  }

  TEST_CASE("angle: between two displacements, from their exact cross and dot products")
  {
    CHECK(Between({.X = 1, .Y = 0}, {.X = 0, .Y = 1}).InRadians()   == Pi / 2.0);
    CHECK(Between({.X = 1, .Y = 0}, {.X = 0, .Y = -1}).InRadians()  == -Pi / 2.0);
    CHECK(Between({.X = 1, .Y = 0}, {.X = -1, .Y = 0}).InRadians()  == Pi);
    CHECK(Between({.X = 1, .Y = 3}, {.X = 2, .Y = 6}).InRadians()   == 0.0);
    CHECK(Between({.X = 1, .Y = 3}, {.X = -2, .Y = -6}).InRadians() == Pi);
    CHECK(Between({.X = 1, .Y = 2}, {.X = 1.0000000000000002, .Y = 2}).InRadians()
          == -0x1.9999999999999p-54);
    CHECK(Between({.X = 1, .Y = 0}, {.X = 1, .Y = 1e-12}).InRadians()
          == doctest::Approx(1e-12).epsilon(1e-15));
  }

  TEST_CASE("angle: the order of headings, exactly, as Heading orders them")
  {
    int wrong = 0;
    for (int ax = -9; ax <= 9; ++ax)
    {
      for (int ay = -9; ay <= 9; ++ay)
      {
        for (int bx = -9; bx <= 9; ++bx)
        {
          for (int by = -9; by <= 9; ++by)
          {
            const Vector lhs   = {.X = static_cast<double>(ax), .Y = static_cast<double>(ay)};
            const Vector rhs   = {.X = static_cast<double>(bx), .Y = static_cast<double>(by)};
            const int    left  = Half(ax, ay);
            const int    right = Half(bx, by);
            const bool   same  = left == right && (left == 0 || left == 2);
            const int    truth = same ? Signed((ay * bx) - (ax * by)) : Signed(left - right);
            const std::strong_ordering order = Order(lhs, rhs);
            const int got        = static_cast<int>(order > 0) - static_cast<int>(order < 0);
            wrong               += got == truth ? 0 : 1;
            const double first   = Heading(lhs).InRadians();
            const double second  = Heading(rhs).InRadians();
            wrong += (truth < 0 && first > second) || (truth > 0 && first < second) ? 1 : 0;
          }
        }
      }
    }
    CHECK(wrong == 0);
  }

  TEST_CASE("angle: the arithmetic mixes units at run time too")
  {
    CHECK((Radians(Pi) - Degrees(90)) == Degrees(90));
    CHECK((Radians(Pi) - Degrees(90)).HeldInDegrees());
    CHECK((Radians(3.0) / 2.0).InRadians() == 1.5);
    CHECK_FALSE((Radians(3.0) / 2.0).HeldInDegrees());
  }

  TEST_CASE("angle: headings past the ordinary: infinite, NaN, tiny and lifted")
  {
    constexpr double NotANumber = std::numeric_limits<double>::quiet_NaN();
    CHECK(Heading({.X = Infinite, .Y = 1}).InRadians()          == 0.0);
    CHECK(Heading({.X = 1, .Y = Infinite}).InRadians()          == Pi / 2.0);
    CHECK(Heading({.X = -Infinite, .Y = -Infinite}).InDegrees() == -135.0);
    CHECK(std::isnan(Heading({.X = NotANumber, .Y = 1}).InRadians()));
    CHECK(std::isnan(Heading({.X = 1, .Y = NotANumber}).InRadians()));
    /**
     *  atan t is t, nearest, for t below 2^-400; a subnormal one is lifted first.
     */
    CHECK(Heading({.X = 1, .Y = 1e-200}).InRadians()  == 1e-200);
    CHECK(Heading({.X = 1, .Y = -1e-200}).InRadians() == -1e-200);
    CHECK(Heading({.X = 1, .Y = 1e-310}).InRadians()  == 1e-310);
  }

  TEST_CASE("angle: between and order past the ordinary: scaled, infinite and none")
  {
    constexpr double NotANumber = std::numeric_limits<double>::quiet_NaN();
    CHECK(std::isnan(Between({.X = Infinite, .Y = 0}, {.X = 1, .Y = 0}).InRadians()));
    CHECK(std::isnan(Between({.X = 1, .Y = 0}, {.X = NotANumber, .Y = 0}).InRadians()));
    CHECK(Between({}, {.X = 1, .Y = 0}).InRadians() == 0.0);
    CHECK(Between({.X = 1, .Y = 0}, {}).InRadians() == 0.0);
    /**
     *  Components past 2^480 or below 2^-480 are scaled by powers of two first, exactly.
     */
    CHECK(Between({.X = 1e-200, .Y = 1}, {.X = 1, .Y = 1e-200}).InRadians()
          == Between({.X = 0, .Y = 1}, {.X = 1, .Y = 0}).InRadians());
    CHECK(Between({.X = 1, .Y = 1e200}, {.X = 1, .Y = 0}).InRadians()        == -Pi / 2.0);
    CHECK(Between({.X = 1e305, .Y = 0}, {.X = 0, .Y = 1e305}).InRadians()    == Pi / 2.0);
    CHECK(Between({.X = 1e-305, .Y = 0}, {.X = 0, .Y = -2e-305}).InRadians() == -Pi / 2.0);
    /**
     *  Products below 2^-500 are lifted before the division.
     */
    CHECK(Between({.X = 1e-140, .Y = 0}, {.X = 0, .Y = 1e-140}).InRadians() == Pi / 2.0);
    /**
     *  A tiny angle between, below the axis.
     */
    CHECK(Between({.X = 1, .Y = 0}, {.X = 1, .Y = -1e-200}).InRadians() == -1e-200);

    CHECK((Order({.X = 1e300, .Y = 1}, {.X = 1e300, .Y = 2}) < 0));
    CHECK((Order({.X = 1, .Y = 1}, {.X = 1, .Y = 1e300}) < 0));
    CHECK((Order({.X = 1, .Y = 1e-300}, {.X = 1, .Y = 1}) < 0));
  }

  TEST_CASE("angle: directions past the ordinary: huge, near a quarter turn, not finite")
  {
    constexpr double NotANumber = std::numeric_limits<double>::quiet_NaN();
    CHECK(std::isnan(Direction(Degrees(Infinite), Degrees(1), 0.5).X));
    /**
     *  Past 2^45 degrees the reduction takes fmod first, exactly.
     */
    CHECK(Direction(Degrees(1e20)) == Direction(Degrees(std::fmod(1e20, 360.0))));
    CHECK(Direction(Degrees(1e300), Degrees(0), 0.5)
          == Direction(Degrees(std::fmod(1e300, 360.0))));
    /**
     *  Within 2^-30 of a quarter turn past 16 radians, Payne and Hanek's reduction decides.
     */
    const Vector near = Direction(Radians((Pi / 2.0) * 11.0));
    CHECK(std::abs(near.X) < 1e-14);
    CHECK(near.Y == -1.0);
    CHECK(Wrapped(Radians((Pi / 2.0) * 11.0)).InDegrees() == doctest::Approx(-90.0).epsilon(1e-12));
    /**
     *  Radians reduced into each quarter turn, and past 16.
     */
    CHECK(Wrapped(Radians(3.2)).InDegrees()
          == doctest::Approx((3.2 * 180.0 / Pi) - 360.0).epsilon(1e-15));
    CHECK(Wrapped(Radians(4.7)).InDegrees()
          == doctest::Approx((4.7 * 180.0 / Pi) - 360.0).epsilon(1e-15));
    CHECK(Wrapped(Radians(100.0)).InDegrees()
          == doctest::Approx((100.0 * 180.0 / Pi) - (16.0 * 360.0)).epsilon(1e-13));
    CHECK(Wrapped(Radians(0x1.6ac5b262ca1ffp+849)).InDegrees()
          == doctest::Approx(90.0).epsilon(1e-15));
    CHECK(std::isnan(Wrapped(Radians(Infinite)).InDegrees()));
    CHECK(std::isnan(Wrapped(Radians(NotANumber)).InDegrees()));
  }

  TEST_CASE("angle: within past the ordinary: radian sweeps, not finite")
  {
    constexpr double NotANumber = std::numeric_limits<double>::quiet_NaN();
    CHECK(Within(Degrees(10), Degrees(0), Radians(20.0)));
    CHECK(Within(Degrees(10), Degrees(0), Radians(0.5)));
    CHECK_FALSE(Within(Degrees(40), Degrees(0), Radians(0.5)));
    CHECK(Within(Radians(1), Radians(0), Radians(1e308)));
    CHECK_FALSE(Within(Radians(NotANumber), Radians(0), Radians(1)));
    CHECK_FALSE(Within(Radians(0.5), Radians(NotANumber), Radians(1)));
    CHECK_FALSE(Within(Radians(0.5), Radians(0), Radians(NotANumber)));
    CHECK_FALSE(Within(Degrees(0), Degrees(NotANumber), Degrees(1)));
    CHECK_FALSE(Within(Degrees(0), Degrees(1), Degrees(NotANumber)));
  }

  TEST_CASE("angle: the hard cases are correctly rounded, the slow way")
  {
    for (const CosSin &turned : HardDegrees)
    {
      const std::uint64_t before = Precise::Slowed();
      const Vector        way    = Direction(Degrees(turned.Angle));
      CHECK(way.X == turned.Cosine);
      CHECK(way.Y == turned.Sine);
      CHECK(Precise::Slowed() > before);
    }
    for (const CosSin &turned : HardRadians)
    {
      const std::uint64_t before = Precise::Slowed();
      const Vector        way    = Direction(Radians(turned.Angle));
      CHECK(way.X == turned.Cosine);
      CHECK(way.Y == turned.Sine);
      CHECK(Precise::Slowed() > before);
    }
    for (const auto &rows :
         {std::span<const Headed>{HardHeadings}, std::span<const Headed>{SubnormalHeadings}})
    {
      for (const Headed &headed : rows)
      {
        const std::uint64_t before = Precise::Slowed();
        CHECK(Heading({.X = headed.X, .Y = headed.Y}).InRadians() == headed.Radians);
        /**
         *  But one that rounds to zero: its pair is zeros, the quotient below 2^-1075.
         */
        CHECK((Precise::Slowed() > before || std::equal_to<double>{}(headed.Radians, 0.0)));
      }
    }
    for (const Apart &apart : HardBetweens)
    {
      const std::uint64_t before = Precise::Slowed();
      CHECK(Between(apart.First, apart.Second).InRadians() == apart.Radians);
      CHECK(Precise::Slowed() > before);
    }
  }

  TEST_CASE("angle: the hard cases along a sweep, held as a pair of degrees")
  {
    /**
     *  A hair past a hard angle rounds as it does, and the slow way too.
     */
    const std::uint64_t before = Precise::Slowed();
    CHECK(Direction(Degrees(0x1.fab126e978d5p+6), Degrees(0x1p-80), 1.0)
          == Vector{.X = -0x1.31ca64547ad52p-1, .Y = 0x1.9aa73d164efb7p-1});
    CHECK(Direction(Degrees(0x1.b1589374bc6a8p+7), Degrees(-0x1p-80), 1.0)
          == Vector{.X = -0x1.9aa73d164efb7p-1, .Y = -0x1.31ca64547ad52p-1});
    CHECK(Precise::Slowed() >= before + 2);
    /**
     *  Tiny ones, whose sine is their radians, from the pair.
     */
    CHECK(Direction(Degrees(0x1p-520), Degrees(0x1p-600), 1.0)
          == Vector{.X = 1, .Y = 0x1.1df46a2529d39p-526});
    CHECK(Direction(Degrees(-0x1p-520), Degrees(0x1.8p-590), 1.0)
          == Vector{.X = 1, .Y = -0x1.1df46a2529d39p-526});
  }

  TEST_CASE("angle: exact and near-exact inputs never take the slow path")
  {
    const std::uint64_t before = Precise::Slowed();
    for (int fifteen = -48; fifteen <= 48; ++fifteen)
    {
      /**
       *  Multiples of 15 degrees, exactly and in radians, within a few ulps of 4 pi.
       */
      const double degrees = 15.0 * fifteen;
      const Vector exact   = Direction(Degrees(degrees));
      const Vector near    = Direction(Radians(fifteen * (Pi / 12.0)));
      CHECK(std::abs(near.X - exact.X) <= 0x1p-48);
      CHECK(std::abs(near.Y - exact.Y) <= 0x1p-48);
      CHECK(Direction(Degrees(degrees), Degrees(15), 0.5) == Direction(Degrees(degrees + 7.5)));
      /**
       *  Their headings, from the rounded directions.
       */
      CHECK(Heading(exact).InRadians() == doctest::Approx(Wrapped(Degrees(degrees)).InRadians()));
      CHECK(Between(exact, Direction(Degrees(degrees + 30.0))).InRadians()
            == doctest::Approx(Pi / 6.0));
    }
    for (const double size : {0x1p-1000, 1e-300, 0.1, 1.0, 3.0, 1e300})
    {
      /**
       *  Axes and diagonals, and a 3-4-5 triangle's.
       */
      for (const double x : {-size, 0.0, size})
      {
        for (const double y : {-size, 0.0, size})
        {
          static_cast<void>(Heading({.X = x, .Y = y}));
          static_cast<void>(Between({.X = x, .Y = y}, {.X = size, .Y = 0}));
        }
      }
      static_cast<void>(Heading({.X = 3.0 * size, .Y = 4.0 * size}));
      static_cast<void>(Between({.X = 3, .Y = 4}, {.X = -4, .Y = 3}));
    }
    CHECK(Precise::Slowed() == before);
  }

  TEST_CASE("angle: the slow path's parts: pi past its table, a power of two, the subnormal grid")
  {
    /**
     *  pi from its table, and by Machin's formula past it, which agrees to the table's last bit.
     */
    CHECK(std::ranges::equal(Precise::Pi(4096).Words, Table::PiWords));
    CHECK(
        std::ranges::equal(Precise::Pi(4096 + 1024).Words | std::views::drop(32), Table::PiWords));

    /**
     *  The arctangent on a diagonal the slow way, its weight 1/2 exactly, on either side.
     */
    const std::array<double, 1> one      = {1.0};
    const auto                  diagonal = [&one](int bits) -> Precise::Fixed
    {
      return Precise::Arc(one, one, bits);
    };
    CHECK(Precise::Nearest({.High = Pi / 4.0, .Low = 0x1p-60}, 0x1p-40, diagonal)  == Pi / 4.0);
    CHECK(Precise::Nearest({.High = Pi / 4.0, .Low = -0x1p-60}, 0x1p-40, diagonal) == Pi / 4.0);

    /**
     *  Just below 1 the doubles are twice as close: 1 - 2^-54 is the halfway point.
     */
    const auto below = [](int bits) -> Precise::Fixed
    {
      return Ones(bits, {{-80, -54}, {-53, 0}});
    };
    const auto above = [](int bits) -> Precise::Fixed
    {
      return Ones(bits, {{-80, -79}, {-54, 0}});
    };
    CHECK(Precise::Nearest({.High = 1.0, .Low = -0x1p-54}, 0x1p-40, below)
          == std::nextafter(1.0, 0.0));
    CHECK(Precise::Nearest({.High = 1.0, .Low = -0x1p-54}, 0x1p-40, above) == 1.0);

    /**
     *  On the subnormal grid, either way of 1.5 2^-1074, either sign.
     */
    for (const bool negative : {false, true})
    {
      const double sign   = negative ? -1.0 : 1.0;
      const auto   upward = [negative](int bits) -> Precise::Fixed
      {
        return Ones(bits, {{-1100, -1099}, {-1075, -1073}}, negative);
      };
      const auto downward = [negative](int bits) -> Precise::Fixed
      {
        return Ones(bits, {{-1076, -1075}, {-1074, -1073}}, negative);
      };
      CHECK(Precise::Nearest({.High = sign * 0x1p-1074}, 0x1p-40, upward)   == sign * 0x1p-1073);
      CHECK(Precise::Nearest({.High = sign * 0x1p-1074}, 0x1p-40, downward) == sign * 0x1p-1074);
    }
  }
} // namespace Geoxy
