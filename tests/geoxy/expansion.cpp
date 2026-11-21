/***************************************************************************************************
 *  @file      tests/geoxy/expansion.cpp
 *  @brief     Tests for numbers held exactly as sums of doubles, and the double nearest a quotient.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "expansion.hpp"

#include "real.hpp"

#include "geoxy/angle.hpp"
#include "geoxy/curve.hpp"
#include "geoxy/ellipse.hpp"
#include "geoxy/intersect.hpp"
#include "geoxy/line.hpp"
#include "geoxy/point.hpp"

#include <doctest/doctest.h>

#include <algorithm>
#include <bit>
#include <cmath>
#include <cstdint>
#include <limits>
#include <numbers>
#include <utility>

namespace Geoxy::Expansion
{
  TEST_CASE("expansion: sums and products of doubles held exactly")
  {
    const double big   = 0x1p53;
    const Number sum   = Sum(Of(big), Of(1.0));
    const Number again = Difference(sum, Of(big));
    CHECK(Approximate(again)         == 1.0);
    CHECK(Sign(Difference(sum, sum)) == 0);
    CHECK(Sign(Difference(Of(1.0), Of(2.0))) == -1);

    const double third   = 1.0 / 3.0;
    const Number squared = Product(Of(third), Of(third));
    CHECK(squared.Parts.Count() == 2);
    CHECK(Approximate(Difference(squared, Scaled(Of(third), third))) == 0.0);
    CHECK(Sign(Scaled(Of(1.0), 0.0)) == 0);
    CHECK(Approximate(Scaled(Sum(Of(1.0), Of(1024.0)), 3.0)) == 3075.0);

    /**
     *  More limbs than are held in place: the rest in the thread's pool, the number the same.
     */
    Number many;
    double total = 0.0;
    for (int index = 0; index < 34; ++index)
    {
      many  = Sum(many, Of(std::ldexp(1.0, (60 * index) - 1000)));
      total = total + std::ldexp(1.0, (60 * index) - 1000);
    }
    CHECK(many.Parts.Count() > Parts::Held);
    CHECK(Approximate(many)             == total);
    CHECK(Sign(Difference(many, many))  == 0);
    CHECK(Sign(Product(many, Of(-1.0))) == -1);
    CHECK(Product(many, many).Parts.Count() > Parts::Held);
    CHECK(Sign(Product(Number{}, Of(2.0)))        == 0);
    CHECK(Approximate(Difference(Number{}, many)) == -total);
  }

  TEST_CASE("expansion: limbs copied and moved, held in place or in the pool")
  {
    Number few = Of(-0x1.8p-1070);
    Number lots;
    for (int index = 0; index < 20; ++index)
    {
      lots = Sum(lots, Of(std::ldexp(1.0, 64 * index)));
    }
    REQUIRE(lots.Parts.Count() > Parts::Held);
    const double  value = Approximate(lots);
    const Number &same  = lots;
    Number        copy  = same;
    CHECK(Approximate(copy) == value);
    copy = few;
    CHECK(Approximate(copy) == -0x1.8p-1070);
    copy = same;
    CHECK(Approximate(copy) == value);
    copy = copy.Parts.Count() > 0 ? copy : few;
    CHECK(Approximate(copy) == value);
    Number moved = std::move(copy);
    CHECK(Approximate(moved) == value);
    Number held = Number{few};
    CHECK(Approximate(held) == -0x1.8p-1070);
    held = std::move(moved);
    CHECK(Approximate(held) == value);
    held = Number{few};
    CHECK(Approximate(held) == -0x1.8p-1070);
    Number &alias = held;
    held          = alias;
    held          = std::move(alias);
    CHECK(Approximate(held) == -0x1.8p-1070);
    /**
     *  A number of many limbs that comes out of few: back in place.
     */
    const Number shrunk = Difference(lots, Difference(lots, Of(5.0)));
    CHECK(shrunk.Parts.Count() == 1);
    CHECK(Approximate(shrunk)  == 5.0);
  }

  TEST_CASE("expansion: the double nearest a number, ties to even, subnormals and infinities")
  {
    const double most = std::numeric_limits<double>::max();
    const double tiny = std::numeric_limits<double>::denorm_min();
    CHECK(Approximate(Of(tiny))  == tiny);
    CHECK(Approximate(Of(-most)) == -most);
    CHECK(Approximate(Product(Of(tiny), Of(0.5))) == 0.0);
    CHECK(Approximate(Product(Of(tiny), Of(0.75))) == tiny);
    CHECK(Approximate(Product(Of(tiny), Of(1.5))) == 2 * tiny);
    CHECK(Approximate(Product(Of(-tiny), Of(0x1p-60))) == 0.0);
    CHECK(Approximate(Product(Of(0x1p-1022), Of(0.5))) == 0x1p-1023);
    CHECK(Approximate(Product(Of(most), Of(2.0))) == std::numeric_limits<double>::infinity());
    CHECK(Approximate(Product(Of(-most), Of(2.0))) == -std::numeric_limits<double>::infinity());
    CHECK(Approximate(Sum(Of(most), Of(0x1p970))) == std::numeric_limits<double>::infinity());
    CHECK(Approximate(Sum(Of(most), Of(0x1.fp969))) == most);
    CHECK(Approximate(Sum(Of(1.0), Of(0x1p-53))) == 1.0);
    CHECK(Approximate(Sum(Of(1.0 + 0x1p-52), Of(0x1p-53))) == 1.0 + 0x1p-51);
    CHECK(Approximate(Sum(Sum(Of(1.0), Of(0x1p-53)), Of(0x1p-400))) == 1.0 + 0x1p-52);
    CHECK(Approximate(Sum(Of(-1.0), Of(-0x1p-53))) == -1.0);
    CHECK(Approximate(Of(std::ldexp(1.0, 63)))     == 0x1p63);
    CHECK(Approximate(Of(-std::ldexp(1.0, 63)))    == -0x1p63);
    CHECK(Sign(Of(-std::ldexp(1.0, 63)))           == -1);
    CHECK(Approximate(Product(Of(0x1p-600), Of(0x1p-600))) == 0.0);
    CHECK(Sign(Product(Of(0x1p-600), Of(-0x1p-600)))       == -1);
  }

  TEST_CASE("expansion: the double nearest a quotient, a square root in it or not")
  {
    /**
     *  The square root and the quotient are rounded correctly by the hardware: the same doubles.
     */
    CHECK(Rounded({.Top = {}, .Bottom = Of(1.0), .Root = Of(2.0), .Below = Of(1.0)}, 1)
          == std::sqrt(2.0));
    CHECK(Rounded({.Top = {}, .Bottom = Of(1.0), .Root = Of(2.0), .Below = Of(1.0)}, -1)
          == -std::sqrt(2.0));
    CHECK(Rounded({.Top = Of(1.0), .Bottom = {}, .Root = {}, .Below = Of(3.0)}, 1)   == 1.0 / 3.0);
    CHECK(Rounded({.Top = Of(-1.0), .Bottom = {}, .Root = {}, .Below = Of(-3.0)}, 1) == 1.0 / 3.0);
    CHECK(Rounded({.Top = Of(7.0), .Bottom = {}, .Root = {}, .Below = Of(1.0)}, 1)   == 7.0);

    /**
     *  2^53 + 1 - 2^53 cancels to 1, which a sum of rounded parts would lose.
     */
    const Number lost = Difference(Sum(Of(0x1p53), Of(1.0)), Of(0x1p53));
    CHECK(Rounded({.Top = lost, .Bottom = {}, .Root = {}, .Below = Of(1.0)}, 1) == 1.0);

    /**
     *  sqrt(2) - 1.4142135623730951 is tiny: the guess is far off, and is walked back.
     */
    const double root = std::numbers::sqrt2;
    const double tiny =
        Rounded({.Top = Of(-root), .Bottom = Of(1.0), .Root = Of(2.0), .Below = Of(1.0)}, 1);
    CHECK(tiny < 0.0);
    CHECK(std::abs(tiny) < 2e-16);

    /**
     *  A tie: 3/2 lies halfway between 1 and 2 only for a format of one bit, never for doubles; but
     *  (1 + 2^-53) is halfway between 1 and its next double, and rounds to the even, 1.
     */
    const Number half = Sum(Of(1.0), Of(0x1p-53));
    CHECK(Rounded({.Top = half, .Bottom = {}, .Root = {}, .Below = Of(1.0)}, 1) == 1.0);
    const Number oddHalf = Sum(Of(1.0 + 0x1p-52), Of(0x1p-53));
    CHECK(Rounded({.Top = oddHalf, .Bottom = {}, .Root = {}, .Below = Of(1.0)}, 1)
          == 1.0 + 0x1p-51);

    /**
     *  Just short of that halfway point, by less than double-double holds: the estimate is the even
     *  double above, yet the double nearest is the odd one below.
     */
    const Number shy = Sum(oddHalf, Of(-0x1p-130));
    CHECK(Rounded({.Top = shy, .Bottom = {}, .Root = {}, .Below = Of(1.0)}, 1) == 1.0 + 0x1p-52);
    const Number past = Sum(oddHalf, Of(0x1p-130));
    CHECK(Rounded({.Top = past, .Bottom = {}, .Root = {}, .Below = Of(1.0)}, 1) == 1.0 + 0x1p-51);

    /**
     *  Within a step of a halfway point, 1 + 2^-53 + 2^-120 sqrt(2): beyond double-double's reach,
     *  settled by exact comparisons, on the halfway point itself and either side of it.
     */
    const Number   near   = Sum(Of(1.0), Of(0x1p-53));
    const Quotient beyond = {
        .Top = near, .Bottom = Of(0x1p-120), .Root = Of(2.0), .Below = Of(1.0)};
    CHECK(Rounded(beyond, 1)  == 1.0 + 0x1p-52);
    CHECK(Rounded(beyond, -1) == 1.0);

    /**
     *  sqrt(4) - 2 is 0 exactly: the estimate is the value, found so in one comparison.
     */
    CHECK(Rounded({.Top = Of(-2.0), .Bottom = Of(1.0), .Root = Of(4.0), .Below = Of(1.0)}, 1)
          == 0.0);

    /**
     *  max - max, or past the largest double: past exact sums, the estimate.
     */
    const double most = std::numeric_limits<double>::max();
    CHECK(Rounded({.Top = Of(most), .Bottom = Of(-most), .Root = Of(1.0), .Below = Of(1.0)}, 1)
          == 0.0);
    CHECK(
        std::isinf(Rounded({.Top = Of(most), .Bottom = {}, .Root = {}, .Below = Of(0x1p-10)}, 1)));
  }

  namespace
  {
    /**
     *  Splitmix64: numbers that look random, the same on every run and every machine.
     */
    struct Mixer
    {
      std::uint64_t State = 0x9E3779B97F4A7C15; /**< Where the sequence is. */

      /**
       *  The next number.
       *
       *  @return The number.
       */
      auto Next() -> std::uint64_t
      {
        State             += 0x9E3779B97F4A7C15;
        std::uint64_t mix  = State;
        mix                = (mix ^ (mix >> 30U)) * 0xBF58476D1CE4E5B9;
        mix                = (mix ^ (mix >> 27U)) * 0x94D049BB133111EB;
        return mix ^ (mix >> 31U);
      }

      /**
       *  A double of any size between 2^-60 and 2^61.
       *
       *  @return The double.
       */
      auto Double() -> double
      {
        const double unit = 1.0 + (static_cast<double>(Next() >> 11U) * 0x1p-53);
        return std::ldexp(unit, static_cast<int>(Next() % 121) - 60);
      }
    };

    /**
     *  A double d, its step s, and d + s/2 + a hair: past the halfway point the next double, short
     *  of it d, on it the even of the two. The same value spelled as a sum, with a square root that
     *  cancels a part of it, and over a divisor, either sign: each must round alike.
     *
     *  @param mixer Where the numbers come from.
     *  @param flip  1 or -1: the value's sign.
     */
    void Halved(Mixer &mixer, double flip)
    {
      const double d     = mixer.Double();
      const double step  = std::nextafter(d, std::numeric_limits<double>::infinity()) - d;
      const int    sign  = static_cast<int>(mixer.Next() % 3) - 1;
      const double hair  = sign * std::ldexp(step, -2 - static_cast<int>(mixer.Next() % 159));
      const double above = d + step;
      double       want  = (std::bit_cast<std::uint64_t>(d) & 1U) == 0 ? d : above;
      if (sign != 0)
      {
        want = sign > 0 ? above : d;
      }
      const Number value = Sum(Sum(Of(d), Of(step / 2)), Of(hair));
      const double r     = mixer.Double();
      const auto   below = static_cast<double>((2 * (mixer.Next() % 7)) + 1);
      const Number over  = Of(flip * below);
      const Number root  = Product(Of(r), Of(r));
      CHECK(
          Rounded(
              {.Top = Scaled(value, flip * below), .Bottom = {}, .Root = {}, .Below = Of(below)}, 1)
          == flip * want);
      CHECK(Rounded({.Top    = Scaled(Sum(value, Of(-r)), flip * below),
                     .Bottom = over,
                     .Root   = root,
                     .Below  = Of(below)},
                    1)
            == flip * want);
      CHECK(Rounded({.Top    = Scaled(Sum(value, Of(r)), flip * below),
                     .Bottom = over,
                     .Root   = root,
                     .Below  = Of(below)},
                    -1)
            == flip * want);
    }
  } // namespace

  TEST_CASE("expansion: on a halfway point or a hair either side, every way a quotient is spelled")
  {
    Mixer mixer;
    for (int round = 0; round < 2000; ++round)
    {
      Halved(mixer, 1.0);
      Halved(mixer, -1.0);
    }
  }

  TEST_CASE("expansion: a quotient whose parts cancel, past where its conjugate's square holds")
  {
    /**
     *  2x - x, x near 2^913 or 2^-950: the conjugate's numerator overflows or underflows, so the
     *  plain double-double estimate decides; well inside the range, the conjugate does.
     */
    const Number big = Of(0x1.3p913);
    CHECK(Rounded({.Top = Scaled(big, 2.0), .Bottom = big, .Root = Of(1.0), .Below = Of(1.0)}, -1)
          == 0x1.3p913);
    const Number tiny = Of(0x1.3p-950);
    CHECK(Rounded({.Top = Scaled(tiny, 2.0), .Bottom = tiny, .Root = Of(1.0), .Below = Of(1.0)}, -1)
          == 0x1.3p-950);
    /**
     *  1 from 2^400 + 1 less 2^400, and 2^-950 from 1 + 2^-950 less 1, the second too small for the
     *  conjugate's estimate: the plain estimate's bound too wide to settle, the exact search
     *  finding the value a double at once.
     */
    CHECK(Rounded({.Top    = Sum(Of(0x1p400), Of(1.0)),
                   .Bottom = Of(0x1p400),
                   .Root   = Of(1.0),
                   .Below  = Of(1.0)},
                  -1)
          == 1.0);
    CHECK(Rounded({.Top    = Sum(Of(1.0), Of(0x1p-950)),
                   .Bottom = Of(1.0),
                   .Root   = Of(1.0),
                   .Below  = Of(1.0)},
                  -1)
          == 0x1p-950);
    const Number small = Of(3.0);
    CHECK(
        Rounded({.Top = Scaled(small, 2.0), .Bottom = small, .Root = Of(1.0), .Below = Of(1.0)}, -1)
        == 3.0);
  }

  TEST_CASE(
      "expansion: the step to a double's neighbours, read from its bits, as nextafter gives it")
  {
    const auto stepped = [](double value) -> double
    {
      constexpr double Endless = std::numeric_limits<double>::infinity();
      return std::min(std::nextafter(value, Endless) - value,
                      value - std::nextafter(value, -Endless));
    };
    Mixer mixer;
    for (int round = 0; round < 100000; ++round)
    {
      const auto value = std::bit_cast<double>(mixer.Next());
      if (std::isfinite(value))
      {
        CHECK(Spacing(value) == stepped(value));
      }
    }
    for (int power = -1074; power <= 1023; ++power)
    {
      const double two = std::ldexp(1.0, power);
      CHECK(Spacing(two)  == stepped(two));
      CHECK(Spacing(-two) == stepped(-two));
    }
    for (const double edge : {0.0, -0.0, 0x1p-1022, std::numeric_limits<double>::max(),
                              std::numeric_limits<double>::denorm_min(), 0x1.fffffffffffffp-1023})
    {
      CHECK(Spacing(edge) == stepped(edge));
    }
    CHECK(std::isnan(Spacing(std::numeric_limits<double>::infinity())));
    CHECK(std::isnan(Spacing(std::numeric_limits<double>::quiet_NaN())));
  }

  TEST_CASE("expansion: the step up from a magnitude, read from its bits, as nextafter gives it")
  {
    const auto stepped = [](double magnitude) -> double
    {
      return std::nextafter(magnitude, std::numeric_limits<double>::infinity()) - magnitude;
    };
    Mixer mixer;
    for (int round = 0; round < 100000; ++round)
    {
      const double magnitude = std::abs(std::bit_cast<double>(mixer.Next()));
      if (std::isfinite(magnitude))
      {
        CHECK(Above(magnitude) == stepped(magnitude));
      }
    }
    for (int power = -1074; power <= 1023; ++power)
    {
      const double two = std::ldexp(1.0, power);
      CHECK(Above(two) == stepped(two));
    }
    for (const double edge : {0.0, 0x1p-1022, std::numeric_limits<double>::max(),
                              std::numeric_limits<double>::denorm_min(), 0x1.fffffffffffffp-1023})
    {
      CHECK(Above(edge) == stepped(edge));
    }
    CHECK(std::isnan(Above(std::numeric_limits<double>::infinity())));
    CHECK(std::isnan(Above(std::numeric_limits<double>::quiet_NaN())));
  }

  /**
   *  A term below 2^-900 leaves the estimate no bound, and every double is searched; the halving
   *  took the difference of the two ends' places, which passed the largest whole number for a value
   *  of 1.5 or more: (3 + 2^-950) / 1 came out as 9.000000000000002, and -(3 + 2^-950) as
   *  3.0000000000000004.
   */
  TEST_CASE("expansion: a tiny term beside a value of 1.5 or more is searched over every double")
  {
    const double most = std::numeric_limits<double>::max();
    for (const double top : {3.0, 1.5, 0.75, -3.0, 0x1p500, most, -most})
    {
      CAPTURE(top);
      const Quotient parts = {
          .Top = Of(top), .Bottom = Of(0x1p-950), .Root = Of(1.0), .Below = Of(1.0)};
      CHECK(Rounded(parts, 1)  == top);
      CHECK(Rounded(parts, -1) == top);
    }
    /**
     *  Past the largest double's halfway point either way, and a term no double holds: the search
     *  finds the infinities, and the double nearest.
     */
    CHECK(Rounded({.Top = Of(-most), .Bottom = {}, .Root = {}, .Below = Of(0x1p-10)}, 1)
          == -std::numeric_limits<double>::infinity());
    const Number lost = Product(Of(0x1p-600), Of(0x1p-600));
    CHECK(Rounded({.Top = Sum(Of(1.0), lost), .Bottom = lost, .Root = Of(1.0), .Below = Of(2.0)}, 1)
          == 0.5);
    CHECK(Rounded({.Top = Of(1.0), .Bottom = Of(1.0), .Root = lost, .Below = Of(4.0)}, -1) == 0.25);
    CHECK(Rounded({.Top = lost, .Bottom = {}, .Root = {}, .Below = lost}, 1) == 1.0);
    CHECK(Rounded({.Top = Of(1.0), .Bottom = {}, .Root = {}, .Below = lost}, 1)
          == std::numeric_limits<double>::infinity());
    /**
     *  2^-400 - 2^-401 sqrt(4), too small for the conjugate's estimate and left to the plain way,
     *  whose bound is the terms' size: zero, a double, found so in one comparison.
     */
    CHECK(Rounded({.Top = Of(0x1p-400), .Bottom = Of(0x1p-401), .Root = Of(4.0), .Below = Of(1.0)},
                  -1)
          == 0.0);
  }

  /**
   *  The top and the root's share cancel, each near 2^-520, their squares' difference 2.4e-319
   *  among the subnormals though the value, near 2^-161, is not: read with a few bits, it gave the
   *  estimate a bound it did not keep, and a double 2^-17 of the value off. Such a numerator leaves
   *  the estimate to the plain way.
   */
  TEST_CASE("expansion: a cancelling quotient whose numerator falls among the subnormals")
  {
    const double   top   = std::numbers::sqrt3 * 0x1p-520 * (1 + 0x1p-21);
    const Quotient parts = {
        .Top = Of(top), .Bottom = Of(0x1p-510), .Root = Of(3 * 0x1p-20), .Below = Of(0x1p-380)};
    CHECK(Rounded(parts, -1) == 0x1.bb67ae85189b5p-161);
    CHECK(Rounded(parts, 1)  == 0x1.bb67b5732384cp-139);
    /**
     *  Its divisor among the subnormals, its numerator not: the plain way too.
     */
    const double   near  = std::numbers::sqrt3 * 0x1p-100 * (1 + 0x1p-30);
    const Quotient lower = {
        .Top = Of(near), .Bottom = Of(0x1p-100), .Root = Of(3.0), .Below = Of(0x1p-850)};
    CHECK(Rounded(lower, -1) == 0x1.bb67ae3136a2fp+720);
    CHECK(Rounded(lower, 1)  == 0x1.bb67ae88fb9a0p+751);
    /**
     *  A top past the largest double over a divisor that brings it back: the largest double.
     */
    const double most = std::numeric_limits<double>::max();
    CHECK(Rounded({.Top = Scaled(Of(most), 4.0), .Bottom = {}, .Root = {}, .Below = Of(4.0)}, 1)
          == most);
  }

  /**
   *  A line 1e-160 long from a point of an ellipse: where the crossing's quotient cancels, its
   *  divisor, the gap times what all is divided by, came to nothing below the doubles' range and
   *  was divided by, which the sanitize build stops on. The estimate is left to the plain way, as a
   *  numerator lost so is, and the hit is the same: the line's start, a quarter turn short of the
   *  ellipse's end.
   */
  TEST_CASE("expansion: a line 1e-160 long met with an ellipse divides by nothing nowhere")
  {
    const Line    line    = {.Start = {.X = 0, .Y = 0}, .End = {.X = 1e-160, .Y = 1e-160}};
    const Ellipse ellipse = {
        .Map   = {.X = {.X = 2, .Y = 0}, .Y = {.X = 0, .Y = 1}, .Origin = {.X = 0, .Y = 1}},
        .Start = Degrees(0),
        .Sweep = Degrees(360)};
    for (const bool swapped : {false, true})
    {
      CAPTURE(swapped);
      const Intersections met =
          swapped ? Intersect(Curve{ellipse}, Curve{line}, {.Length = 1e-9, .Radians = 1e-9})
                  : Intersect(Curve{line}, Curve{ellipse}, {.Length = 1e-9, .Radians = 1e-9});
      REQUIRE(met.Hits.size() == 1);
      CHECK(met.Hits[0].At    == Point{});
      CHECK((swapped ? met.Hits[0].Second : met.Hits[0].First).Along == 0);
      CHECK((swapped ? met.Hits[0].First : met.Hits[0].Second).Along == 0.75);
    }
  }
} // namespace Geoxy::Expansion
