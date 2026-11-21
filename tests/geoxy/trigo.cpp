/***************************************************************************************************
 *  @file      tests/geoxy/trigo.cpp
 *  @brief     Tests for the library's own sines, cosines, tangents, arc functions and lengths.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "trigo.hpp"

#include "draws.hpp"
#include "golden.hpp"

#include "geoxy/angle.hpp"
#include "geoxy/vector.hpp"

#include <doctest/doctest.h>

#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <numbers>
#include <sstream>
#include <string>

namespace Geoxy
{
  namespace
  {
    constexpr double TrigoNan = std::numeric_limits<double>::quiet_NaN();
    constexpr double TrigoInf = std::numeric_limits<double>::infinity();

    /**
     *  Whether two doubles have the same bits, the sign of zero included; any two not a numbers
     *  pass.
     */
    auto Same(double got, double want) -> bool
    {
      return (std::isnan(got) && std::isnan(want))
          || std::bit_cast<std::uint64_t>(got) == std::bit_cast<std::uint64_t>(want);
    }

    /**
     *  How far a double is from an exact value given as a pair, in steps of a double of the value.
     */
    auto TrigoError(double got, double truth, double rest) -> double
    {
      const int power = std::max(std::ilogb(truth), -1022);
      return std::abs(((got - truth) - rest) / std::ldexp(1.0, power - 52));
    }

    /**
     *  A function of the golden vectors at its arguments, and the most it may be off.
     */
    struct TrigoCase
    {
      double Got   = 0.0; /**< Its value.           */
      double Bound = 0.0; /**< Its bound, in steps. */
    };

    auto TrigoEvaluated(const std::string &name, double first, double second) -> TrigoCase
    {
      if (name == "sin")
      {
        return {.Got = Trigo::Sin(first), .Bound = 0.56};
      }
      if (name == "cos")
      {
        return {.Got = Trigo::Cos(first), .Bound = 0.56};
      }
      if (name == "tan")
      {
        return {.Got = Trigo::Tan(first), .Bound = 0.59};
      }
      if (name == "sind")
      {
        return {.Got = Trigo::WayDegrees(first).Y, .Bound = 0.56};
      }
      if (name == "cosd")
      {
        return {.Got = Trigo::WayDegrees(first).X, .Bound = 0.56};
      }
      if (name == "atan")
      {
        return {.Got = Trigo::Atan(first), .Bound = 0.51};
      }
      if (name == "atan2")
      {
        return {.Got = Trigo::Atan2(first, second), .Bound = 0.51};
      }
      if (name == "asin")
      {
        return {.Got = Trigo::Asin(first), .Bound = 0.53};
      }
      return {.Got = Trigo::Acos(first), .Bound = 0.52};
    }

    /**
     *  FNV-1a over the bits of each value.
     */
    class TrigoChecksum
    {
    public:
      void Add(double value)
      {
        const auto bits = std::bit_cast<std::uint64_t>(value);
        for (unsigned shift = 0; shift < 64; shift += 8)
        {
          m_Hash ^= (bits >> shift) & 0xFFU;
          m_Hash *= 0x100000001B3U;
        }
      }

      [[nodiscard]] auto Value() const -> std::uint64_t { return m_Hash; }

    private:
      std::uint64_t m_Hash = 0xCBF29CE484222325U;
    };
  } // namespace

  TEST_CASE("trigo: each function is within its bound of mpmath's values")
  {
    std::istringstream lines{std::string{Tests::Golden}};
    std::string        line;
    int                count = 0;
    while (std::getline(lines, line))
    {
      if (line.empty())
      {
        continue;
      }
      std::istringstream fields{line};
      std::string        name;
      std::string        first;
      std::string        second;
      std::string        third;
      std::string        fourth;
      fields >> name >> first >> second >> third;
      const bool      two   = name == "atan2" && static_cast<bool>(fields >> fourth);
      const double    x     = std::stod(first);
      const double    y     = two ? std::stod(second) : 0.0;
      const double    truth = std::stod(two ? third : second);
      const double    rest  = std::stod(two ? fourth : third);
      const TrigoCase found = TrigoEvaluated(name, x, y);
      INFO(line);
      CHECK(TrigoError(found.Got, truth, rest) <= found.Bound);
      ++count;
    }
    CHECK(count == 794);
  }

  TEST_CASE("trigo: an angle held in degrees is reduced exactly")
  {
    CHECK(Same(Trigo::WayDegrees(30.0).Y, 0.5));
    CHECK(Same(Trigo::WayDegrees(150.0).Y, 0.5));
    CHECK(Same(Trigo::WayDegrees(-30.0).Y, -0.5));
    CHECK(Same(Trigo::WayDegrees(60.0).X, 0.5));
    CHECK(Same(Trigo::WayDegrees(120.0).X, -0.5));
    CHECK(Same(Trigo::WayDegrees(90.0).X, 0.0));
    CHECK(Same(Trigo::WayDegrees(90.0).Y, 1.0));
    CHECK(Same(Trigo::WayDegrees(-90.0).X, 0.0));
    CHECK(Same(Trigo::WayDegrees(-90.0).Y, -1.0));
    CHECK(Same(Trigo::WayDegrees(180.0).X, -1.0));
    CHECK(Same(Trigo::WayDegrees(180.0).Y, 0.0));
    CHECK(Same(Trigo::WayDegrees(-180.0).Y, -0.0));
    CHECK(Same(Trigo::WayDegrees(270.0).X, 0.0));
    CHECK(Same(Trigo::WayDegrees(270.0).Y, -1.0));
    CHECK(Same(Trigo::WayDegrees(360.0).Y, 0.0));
    CHECK(Same(Trigo::WayDegrees(-360.0).Y, -0.0));
    CHECK(Same(Trigo::WayDegrees(540.0).X, -1.0));
    CHECK(Same(Trigo::WayDegrees(0.0).Y, 0.0));
    CHECK(Same(Trigo::WayDegrees(-0.0).Y, -0.0));
    CHECK(Same(Trigo::WayDegrees(-0.0).X, 1.0));

    /**
     *  Past 2^40, whole turns are taken off in whole numbers: below 2^52, from 2^52 to 2^55, and
     *  beyond, through the cycle of 2's powers modulo 45.
     */
    for (const double degrees : {0x1p40 + 0.5, 0x1.8p51 + 90.25, 0x1p53 + 2.0, 0x1.8p54 + 8.0,
                                 1e300, -1e300, 0x1p60, -0x1.fffffffffffffp1023, 0x1.23456789p70})
    {
      INFO(degrees);
      const double turned = std::fmod(degrees, 360.0);
      CHECK(Same(Trigo::WayDegrees(degrees).X, Trigo::WayDegrees(turned).X));
      CHECK(Same(Trigo::WayDegrees(degrees).Y, Trigo::WayDegrees(turned).Y));
    }
    CHECK(Same(Trigo::WayDegrees(1e300).Y, 0.0));
    CHECK(Same(Trigo::WayDegrees(1e300).X, 1.0));

    /**
     *  Below 2^-900 the radians are worked on their own, and rounded by hand below the normals: the
     *  values are MPFR's, for mpmath's float() rounds twice there.
     */
    struct Tiny
    {
      double Degrees = 0.0;
      double Radians = 0.0;
    };

    for (const Tiny tiny : {
             Tiny{.Degrees =                 0x1p-950, .Radians =   0x1.1df46a2529d39p-956},
             Tiny{.Degrees =  -0x1.b333333333333p-950, .Radians =  -0x1.e61f813f2d814p-956},
             Tiny{.Degrees =  -0x1.fffffffffffffp-901, .Radians =  -0x1.1df46a2529d39p-906},
             Tiny{.Degrees =                0x1p-1020, .Radians =  0x0.11df46a2529d4p-1022},
             Tiny{.Degrees =             0x0.01p-1022, .Radians =  0x0.000477d1a894ap-1022},
             Tiny{.Degrees =    -0x0.000000000cp-1022, .Radians = -0x0.000000000035ap-1022},
             Tiny{.Degrees =  0x0.000000012d687p-1022, .Radians =  0x0.000000000542bp-1022},
             Tiny{.Degrees =  0x0.00004c0000000p-1022, .Radians =  0x0.00000153923e1p-1022},
             Tiny{.Degrees =  0x0.0000000000007p-1022, .Radians =                      0.0},
             Tiny{.Degrees = -0x1.868af36535aebp-1019, .Radians = -0x1.b43d8f55e3e68p-1025},
             Tiny{.Degrees =  0x1.02bdef050f976p-1019, .Radians =  0x1.21047b68a7e08p-1025},
             Tiny{.Degrees = -0x1.bae2b9cbad0fcp-1021, .Radians =  -0x1.eeb5409301dap-1027},
         })
    {
      INFO(tiny.Degrees);
      CHECK(Same(Trigo::WayDegrees(tiny.Degrees).Y, tiny.Radians));
      CHECK(Same(Trigo::WayDegrees(tiny.Degrees).X, 1.0));
    }
    CHECK(std::isnan(Trigo::WayDegrees(TrigoInf).X));
    CHECK(std::isnan(Trigo::WayDegrees(-TrigoInf).Y));
    CHECK(std::isnan(Trigo::WayDegrees(TrigoNan).Y));
  }

  TEST_CASE("trigo: sines of zeros, tiny angles, infinities and not a number")
  {
    CHECK(Same(Trigo::Sin(0.0), 0.0));
    CHECK(Same(Trigo::Sin(-0.0), -0.0));
    CHECK(Same(Trigo::Cos(-0.0), 1.0));
    CHECK(Same(Trigo::Tan(-0.0), -0.0));
    CHECK(Same(Trigo::Sin(0x1p-1074), 0x1p-1074));
    CHECK(Same(Trigo::Tan(-0x1p-600), -0x1p-600));
    CHECK(Same(Trigo::Cos(0x1.fffffffffffffp-28), 1.0));
    CHECK(Same(Trigo::Way(-0x1p-30).X, 1.0));
    CHECK(Same(Trigo::Way(-0x1p-30).Y, -0x1p-30));
    for (const double angle : {TrigoInf, -TrigoInf, TrigoNan})
    {
      CHECK(std::isnan(Trigo::Sin(angle)));
      CHECK(std::isnan(Trigo::Cos(angle)));
      CHECK(std::isnan(Trigo::Tan(angle)));
      CHECK(std::isnan(Trigo::Way(angle).X));
      CHECK(std::isnan(Trigo::Way(angle).Y));
    }
  }

  TEST_CASE("trigo: an angle of any size is reduced, from 2^30 by Payne and Hanek")
  {
    Draws draws{0x7472'6967U};
    for (int round = 0; round < Rounds(4000); ++round)
    {
      const double angle = std::ldexp(draws.Near(1.5, 0.5), 30 + (round % 994));
      const Vector way   = Trigo::Way(angle);
      INFO(angle);
      CHECK(Same(way.X, Trigo::Cos(angle)));
      CHECK(Same(way.Y, Trigo::Sin(angle)));
      CHECK(Same(Trigo::Sin(-angle), -way.Y));
      CHECK(Same(Trigo::Cos(-angle), way.X));
      CHECK(Same(Trigo::Tan(-angle), -Trigo::Tan(angle)));
      CHECK(std::abs(((way.X * way.X) + (way.Y * way.Y)) - 1.0) < 0x1p-50);
    }
  }

  TEST_CASE("trigo: the way of an angle is its cosine and sine together")
  {
    Draws draws{0x7761'7973U};
    for (int round = 0; round < Rounds(20000); ++round)
    {
      const double angle = draws.Near(0.0, round % 2 == 0 ? 7.0 : 1e6);
      const Vector way   = Trigo::Way(angle);
      INFO(angle);
      CHECK(Same(way.X, Trigo::Cos(angle)));
      CHECK(Same(way.Y, Trigo::Sin(angle)));
    }
  }

  TEST_CASE("trigo: the arc tangent of 2^27 and past it is near a quarter turn")
  {
    /**
     *  trigo's own atan returned 2^27 itself at exactly 2^27: its quick path stopped there and its
     *  slow path took only what lay past it.
     */
    constexpr double Quarter = std::numbers::pi / 2.0;
    CHECK(Same(Trigo::Atan(0x1p27), 0x1.921fb52442d18p+0));
    CHECK(Same(Trigo::Atan(-0x1p27), -0x1.921fb52442d18p+0));
    CHECK(Same(Trigo::Atan(0x1p70), Quarter));
    CHECK(Same(Trigo::Atan(-TrigoInf), -Quarter));
    CHECK(Same(Trigo::Atan(0x1p40), 0x1.921fb54441d18p+0));
    CHECK(Same(Trigo::Atan(-0x1p-30), -0x1p-30));
    CHECK(Same(Trigo::Atan(-0.0), -0.0));
    CHECK(std::isnan(Trigo::Atan(TrigoNan)));
  }

  TEST_CASE("trigo: the angle of a displacement is IEEE 754's at zeros and infinities")
  {
    CHECK(Same(Trigo::Atan2(0.0, 1.0), 0.0));
    CHECK(Same(Trigo::Atan2(-0.0, 1.0), -0.0));
    CHECK(Same(Trigo::Atan2(0.0, -0.0), Pi));
    CHECK(Same(Trigo::Atan2(-0.0, -1.0), -Pi));
    CHECK(Same(Trigo::Atan2(2.0, 0.0), Pi / 2.0));
    CHECK(Same(Trigo::Atan2(-2.0, -0.0), -Pi / 2.0));
    CHECK(Same(Trigo::Atan2(TrigoInf, TrigoInf), Pi / 4.0));
    CHECK(Same(Trigo::Atan2(TrigoInf, -TrigoInf), 0x1.2d97c7f3321d2p+1));
    CHECK(Same(Trigo::Atan2(-TrigoInf, 5.0), -Pi / 2.0));
    CHECK(Same(Trigo::Atan2(5.0, TrigoInf), 0.0));
    CHECK(Same(Trigo::Atan2(-5.0, -TrigoInf), -Pi));
    CHECK(std::isnan(Trigo::Atan2(TrigoNan, 1.0)));
    CHECK(std::isnan(Trigo::Atan2(1.0, TrigoNan)));

    /**
     *  Sides far apart in size: the angle is the ratio, off a multiple of pi / 4.
     */
    CHECK(Same(Trigo::Atan2(0x1p-70, 3.0), 0x1p-70 / 3.0));
    CHECK(Same(Trigo::Atan2(-0x1p-70, -3.0), -Pi));
    CHECK(Same(Trigo::Atan2(3.0, -0x1p-70), Pi / 2.0));
    CHECK(Same(Trigo::Atan2(-3.0, 0x1p-70), -Pi / 2.0));
    CHECK(Same(Trigo::Atan2(0x1p-1000, -0x1p-900), Pi));

    /**
     *  Sides beyond [2^-960, 2^1000): scaled by a power of two first, which changes no angle.
     */
    Draws draws{0x6174'616EU};
    for (int round = 0; round < Rounds(2000); ++round)
    {
      const double rise  = draws.Near(0.0, 1.0);
      const double run   = draws.Near(0.0, 1.0);
      const int    power = round % 2 == 0 ? 1010 : -1010;
      INFO(rise);
      INFO(run);
      CHECK(Same(Trigo::Atan2(std::ldexp(rise, power), std::ldexp(run, power)),
                 Trigo::Atan2(rise, run)));
    }
  }

  TEST_CASE("trigo: arc sines and cosines at their ends")
  {
    CHECK(Same(Trigo::Asin(1.0), Pi / 2.0));
    CHECK(Same(Trigo::Asin(-1.0), -Pi / 2.0));
    CHECK(Same(Trigo::Acos(-1.0), Pi));
    CHECK(Same(Trigo::Acos(1.0), 0.0));
    CHECK(Same(Trigo::Asin(-0.0), -0.0));
    CHECK(Same(Trigo::Asin(0x1p-40), 0x1p-40));
    CHECK(Same(Trigo::Acos(0.0), Pi / 2.0));
    CHECK(Same(Trigo::Acos(-0x1p-40), 0x1.921fb54443d18p+0));
    CHECK(Same(Trigo::Asin(0.5), 0x1.0c152382d7366p-1));
    CHECK(Same(Trigo::Acos(0.5), 0x1.0c152382d7366p+0));
    for (const double outside : {1.0000000000000002, -2.0, TrigoInf, TrigoNan})
    {
      CHECK(std::isnan(Trigo::Asin(outside)));
      CHECK(std::isnan(Trigo::Acos(outside)));
    }
  }

  TEST_CASE("trigo: the length of a displacement, scaled where its squares would not hold")
  {
    CHECK(Same(Trigo::Hypot(3.0, -4.0), 5.0));
    CHECK(Same(Trigo::Hypot(-0.0, 0.0), 0.0));
    CHECK(Same(Trigo::Hypot(0x1p-1074, 0.0), 0x1p-1074));
    CHECK(Same(Trigo::Hypot(3e300, 4e300), 5e300));
    CHECK(Same(Trigo::Hypot(0x1.fffffffffffffp1023, 1.0), 0x1.fffffffffffffp1023));
    CHECK(std::isinf(Trigo::Hypot(0x1.fffffffffffffp1023, 0x1.fffffffffffffp1023)));
    CHECK(std::isinf(Trigo::Hypot(TrigoNan, -TrigoInf)));
    CHECK(std::isinf(Trigo::Hypot(TrigoInf, TrigoNan)));
    CHECK(std::isnan(Trigo::Hypot(TrigoNan, 1.0)));
    CHECK(std::isnan(Trigo::Hypot(1.0, TrigoNan)));
    CHECK(std::isnan(Trigo::Hypot(0.0, TrigoNan)));

    /**
     *  Pythagorean triples at every scale: exact lengths, so the nearest double.
     */
    Draws draws{0x6879'706FU};
    for (int round = 0; round < Rounds(4000); ++round)
    {
      const double larger  = std::floor(draws.Near(5000.0, 4990.0));
      const double smaller = std::floor(draws.Near(larger / 2.0, (larger / 2.0) - 1.0));
      const int    power   = (round % 2000) - 1000;
      const double run     = std::ldexp((larger * larger) - (smaller * smaller), power);
      const double rise    = std::ldexp(2.0 * larger * smaller, power);
      INFO(run);
      INFO(rise);
      CHECK(Same(Trigo::Hypot(run, rise),
                 std::ldexp((larger * larger) + (smaller * smaller), power)));
    }
  }

  TEST_CASE("trigo: logarithms, powers, inverse hyperbolic sines and cube roots")
  {
    struct Worked
    {
      int    Function = 0;   /**< 0 log2, 1 exp2, 2 pow, 3 asinh, 4 cbrt. */
      double First    = 0.0; /**< The argument, or the base.              */
      double Second   = 0.0; /**< The power.                              */
      double Truth    = 0.0; /**< mpmath's value, nearest.                */
      double Rest     = 0.0; /**< What that leaves out, nearest.          */
    };

    /**
     *  From mpmath at 400 bits.
     */
    for (const Worked worked : {
             Worked{.Function = 0,
                    .First    = 0x1.8p+0,
                    .Second   = 0.0,
                    .Truth    = 0x1.2b803473f7ad1p-1,
                    .Rest     = -0x1.817fd3b7d7e5dp-58},
             Worked{.Function = 0,
                    .First    = std::numbers::sqrt2,
                    .Second   = 0.0,
                    .Truth    = 0x1.0000000000001p-1,
                    .Rest     = -0x1.c991a39819edfp-57},
             Worked{.Function = 0,
                    .First    = 0x1.fffffffffffffp-1,
                    .Second   = 0.0,
                    .Truth    = -0x1.71547652b82fep-53,
                    .Rest     = -0x1.cf33ba5220647p-107},
             Worked{.Function = 0,
                    .First    = 0x1.0000000000001p+0,
                    .Second   = 0.0,
                    .Truth    = 0x1.71547652b82fdp-52,
                    .Rest     = 0x1.7b365759f7d4fp-106},
             Worked{.Function = 0,
                    .First    = 0x1.3c5p+700,
                    .Second   = 0.0,
                    .Truth    = 0x1.5e271103dc3afp+9,
                    .Rest     = 0x1.5983d192306a7p-45},
             Worked{.Function = 0,
                    .First    = 0x1.e1p-1030,
                    .Second   = 0.0,
                    .Truth    = -0x1.0145c44fcc968p+10,
                    .Rest     = 0x1.f86493917b407p-44},
             Worked{.Function = 0,
                    .First    = 0x1.92p-3,
                    .Second   = 0.0,
                    .Truth    = -0x1.2caa569330c4fp+1,
                    .Rest     = 0x1.29b7bfe661cfdp-55},
             Worked{.Function = 1,
                    .First    = 0x1p-1,
                    .Second   = 0.0,
                    .Truth    = std::numbers::sqrt2,
                    .Rest     = -0x1.bdd3413b26456p-54},
             Worked{.Function = 1,
                    .First    = -0x1p-1,
                    .Second   = 0.0,
                    .Truth    = 0x1.6a09e667f3bcdp-1,
                    .Rest     = -0x1.bdd3413b26456p-55},
             Worked{.Function = 1,
                    .First    = 0x1.fd7p+9,
                    .Second   = 0.0,
                    .Truth    = 0x1.d5818dcfba487p+1018,
                    .Rest     = 0x1.2ed02d75b3707p+963},
             Worked{.Function = 1,
                    .First    = 0x1.999999999999ap-4,
                    .Second   = 0.0,
                    .Truth    = 0x1.125fbee250664p+0,
                    .Rest     = 0x1.575e72323eb32p-55},
             Worked{.Function = 2,
                    .First    = 0x1.4p+5,
                    .Second   = 0x1.ccccccccccccdp-1,
                    .Truth    = 0x1.ba8fd5778503fp+4,
                    .Rest     = 0x1.9c0f3d8d8b7d3p-51},
             Worked{.Function = 2,
                    .First    = 0x1.2aa1127df7c41p+0,
                    .Second   = 0x1.ccccccccccccdp-1,
                    .Truth    = 0x1.261095fcd13c8p+0,
                    .Rest     = -0x1.5482e1b147e4fp-54},
             Worked{.Function = 2,
                    .First    = 0x1.8p+2,
                    .Second   = 0x1.8p-1,
                    .Truth    = 0x1.eab5536a3a356p+1,
                    .Rest     = 0x1.bb6edf8af5ae5p-55},
             Worked{.Function = 2,
                    .First    = 0x1.fp+5,
                    .Second   = 0x1p-1,
                    .Truth    = 0x1.f7efbeb8d4f12p+2,
                    .Rest     = -0x1.e84e8a6ff3256p-54},
             Worked{.Function = 3,
                    .First    = 0x1.e853448214359p-4,
                    .Second   = 0.0,
                    .Truth    = 0x1.e72d019df1f12p-4,
                    .Rest     = 0x1.ffe5f8f43ca64p-60},
             Worked{.Function = 3,
                    .First    = -0x1.fff15f3865de9p-7,
                    .Second   = 0.0,
                    .Truth    = -0x1.ffec0a7e72eb3p-7,
                    .Rest     = -0x1.8d468d320979ep-61},
             Worked{.Function = 3,
                    .First    = 0x1.8p+1,
                    .Second   = 0.0,
                    .Truth    = 0x1.d185b507edc0ep+0,
                    .Rest     = -0x1.04d61d27f1036p-59},
             Worked{.Function = 3,
                    .First    = -0x1.5p+20,
                    .Second   = 0.0,
                    .Truth    = -0x1.da7f2d405c985p+3,
                    .Rest     = 0x1.2faa2e176c81fp-51},
             Worked{.Function = 3,
                    .First    = 0x1.4p+40,
                    .Second   = 0.0,
                    .Truth    = 0x1.ca465c63fcf6fp+4,
                    .Rest     = -0x1.00917b5f572e7p-50},
             Worked{.Function = 3,
                    .First    = 0x1.0000000000001p-26,
                    .Second   = 0.0,
                    .Truth    = 0x1.0000000000001p-26,
                    .Rest     = -0x1.5555555555559p-81},
             Worked{.Function = 4,
                    .First    = -0x1.4p+3,
                    .Second   = 0.0,
                    .Truth    = -0x1.13c484138704fp+1,
                    .Rest     = 0x1.abaee946d90b0p-54},
             Worked{.Function = 4,
                    .First    = 0x1.5555555555555p-2,
                    .Second   = 0.0,
                    .Truth    = 0x1.63003fbb4c375p-1,
                    .Rest     = 0x1.28b0486c7699fp-55},
             Worked{.Function = 4,
                    .First    = 0x1.23456789abcdep+900,
                    .Second   = 0.0,
                    .Truth    = 0x1.0b4145573271cp+300,
                    .Rest     = 0x1.3dce6b2bc3fb5p+245},
             Worked{.Function = 4,
                    .First    = -0x1.8p-1050,
                    .Second   = 0.0,
                    .Truth    = -0x1.250bfe1b082f5p-350,
                    .Rest     = 0x1.91cacc5e8d601p-408},
         })
    {
      const std::array<double, 5> got = {Trigo::Log2(worked.First), Trigo::Exp2(worked.First),
                                         Trigo::Pow(worked.First, worked.Second),
                                         Trigo::Asinh(worked.First), Trigo::Cbrt(worked.First)};
      INFO(worked.First);
      CHECK(TrigoError(got.at(static_cast<std::size_t>(worked.Function)), worked.Truth, worked.Rest)
            <= 0.8);
    }

    /**
     *  Exact at powers of two, whole powers and cubes.
     */
    for (int power = -1074; power <= 1023; power += 7)
    {
      CHECK(Same(Trigo::Log2(std::ldexp(1.0, power)), static_cast<double>(power)));
      CHECK(Same(Trigo::Exp2(static_cast<double>(power)), std::ldexp(1.0, power)));
    }
    CHECK(Same(Trigo::Cbrt(27.0), 3.0));
    CHECK(Same(Trigo::Cbrt(-0x1p-1074), -0x1p-358));
    CHECK(Same(Trigo::Exp2(-1076.0), 0.0));
    CHECK(Same(Trigo::Pow(0x1.8p-1000, 1.5), 0.0));
    CHECK(Same(Trigo::Pow(0.0, 0.5), 0.0));
    CHECK(Same(Trigo::Pow(0x1p+800, 2.5), TrigoInf));
    CHECK(Same(Trigo::Pow(TrigoInf, 0.9), TrigoInf));

    /**
     *  Zeros, infinities and not a number.
     */
    CHECK(Same(Trigo::Log2(0.0), -TrigoInf));
    CHECK(Same(Trigo::Log2(TrigoInf), TrigoInf));
    CHECK(std::isnan(Trigo::Log2(-1.0)));
    CHECK(std::isnan(Trigo::Log2(TrigoNan)));
    CHECK(Same(Trigo::Exp2(1e9), TrigoInf));
    CHECK(Same(Trigo::Exp2(-TrigoInf), 0.0));
    CHECK(std::isnan(Trigo::Exp2(TrigoNan)));
    CHECK(std::isnan(Trigo::Pow(-2.0, 0.5)));
    CHECK(Same(Trigo::Asinh(-0.0), -0.0));
    CHECK(Same(Trigo::Asinh(-TrigoInf), -TrigoInf));
    CHECK(std::isnan(Trigo::Asinh(TrigoNan)));
    CHECK(Same(Trigo::Asinh(0x1p-30), 0x1p-30));
    CHECK(Same(Trigo::Cbrt(-0.0), -0.0));
    CHECK(Same(Trigo::Cbrt(TrigoInf), TrigoInf));
    CHECK(std::isnan(Trigo::Cbrt(TrigoNan)));

    /**
     *  Odd where they should be, and inverse to one another within a few steps.
     */
    Draws draws{0x6C6F'6732U};
    for (int round = 0; round < Rounds(4000); ++round)
    {
      const double value = std::ldexp(draws.Near(1.5, 0.5), (round % 400) - 200);
      INFO(value);
      CHECK(Same(Trigo::Asinh(-value), -Trigo::Asinh(value)));
      CHECK(Same(Trigo::Cbrt(-value), -Trigo::Cbrt(value)));
      const double root = Trigo::Cbrt(value);
      CHECK(std::abs((root * root * root) - value) <= 0x1p-50 * value);
      CHECK(std::abs(Trigo::Exp2(Trigo::Log2(value)) - value) <= 0x1p-44 * value);
    }
  }

  /**
   *  The functions are worked from operations IEEE 754 rounds one way, so every machine gives the
   *  same bits; a failure on one platform means it computes them otherwise, not that the code
   *  changed.
   */
  TEST_CASE("trigo: every function gives the same bits on every platform")
  {
    Draws         draws{0x7369'6E65U};
    TrigoChecksum sum;
    for (int index = 0; index < 100'000; ++index)
    {
      const double angle = draws.Near(0.0, index % 4 == 0 ? 1e4 : 4.0);
      const double ratio = draws.Near(0.0, index % 3 == 0 ? 50.0 : 1.0);
      const Vector way   = Trigo::Way(angle);
      const Vector turn  = Trigo::WayDegrees(angle * 60.0);
      sum.Add(Trigo::Sin(angle));
      sum.Add(Trigo::Tan(angle));
      sum.Add(way.X);
      sum.Add(turn.X);
      sum.Add(turn.Y);
      sum.Add(Trigo::Atan(ratio));
      sum.Add(Trigo::Atan2(ratio, angle));
      sum.Add(Trigo::Asin(ratio / 50.0));
      sum.Add(Trigo::Acos(ratio / 50.0));
      sum.Add(Trigo::Hypot(angle, ratio));
      sum.Add(Trigo::Log2(std::abs(angle)));
      sum.Add(Trigo::Exp2(ratio));
      sum.Add(Trigo::Pow(std::abs(angle), 0.9));
      sum.Add(Trigo::Asinh(ratio));
      sum.Add(Trigo::Cbrt(angle));
    }
    CHECK(sum.Value() == 8280383835261003768U);
  }
} // namespace Geoxy
