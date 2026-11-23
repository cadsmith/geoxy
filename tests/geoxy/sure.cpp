/***************************************************************************************************
 *  @file      tests/geoxy/sure.cpp
 *  @brief     Tests for numbers in double-double with a bound: every bound holds against the exact
 *             value, and the certificate rounds to the nearest double or says it cannot.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "sure.hpp"

#include "draws.hpp"
#include "expansion.hpp"

#include <doctest/doctest.h>

#include <array>
#include <cmath>
#include <limits>
#include <utility>

namespace Geoxy::Sure
{
  namespace
  {
    /**
     *  A number held exactly: its two doubles summed.
     *
     *  @param  number The number.
     *  @return        Its value, exactly.
     */
    auto Held(const Number &number) -> Expansion::Number
    {
      return Expansion::Sum(Expansion::Of(number.Value.High), Expansion::Of(number.Value.Low));
    }

    /**
     *  Whether an exact value lies within a number's bound of it: value - (number - bound) and
     *  (number + bound) - value not below zero.
     *
     *  @param  number The number.
     *  @param  value  The exact value.
     *  @return        Whether it does.
     */
    auto Within(const Number &number, const Expansion::Number &value) -> bool
    {
      const Expansion::Number low =
          Expansion::Difference(Held(number), Expansion::Of(number.Bound));
      const Expansion::Number high = Expansion::Sum(Held(number), Expansion::Of(number.Bound));
      return Expansion::Sign(Expansion::Difference(value, low))  >= 0
          && Expansion::Sign(Expansion::Difference(high, value)) >= 0;
    }

    /**
     *  A number with a bound, from a double and a part of it as its error, the error made exact.
     *
     *  @param  draws Where the numbers come from.
     *  @return       The number, and the exact value it stands for.
     */
    auto Drawn(Draws &draws) -> std::pair<Number, Expansion::Number>
    {
      const double highFraction = draws.Near(0.0, 1.0);
      const double high  = std::ldexp(highFraction, static_cast<int>(draws.Next() % 64) - 32);
      const double error = high * std::ldexp(draws.Near(0.0, 1.0), -60);
      return {Number{.Value = {.High = high, .Low = 0.0}, .Bound = std::abs(error)},
              Expansion::Sum(Expansion::Of(high), Expansion::Of(error))};
    }
  } // namespace

  TEST_CASE("sure: exact differences and products carry no bound, but below the normal doubles")
  {
    const Number apart = Apart(1.0, 0x1p-60);
    CHECK(apart.Value.High == 1.0);
    CHECK(apart.Value.Low  == -0x1p-60);
    CHECK(apart.Bound      == 0.0);
    const Number product = Product(1.0 + 0x1p-52, 1.0 + 0x1p-52);
    CHECK(product.Value.High == 1.0 + 0x1p-51);
    CHECK(product.Value.Low  == 0x1p-104);
    CHECK(product.Bound      == 0.0);
    CHECK(Product(0x1p-600, 0x1p-600).Bound == Tiny);
    CHECK(Exact(3.0).Bound         == 0.0);
    CHECK((-Exact(3.0)).Value.High == -3.0);
  }

  TEST_CASE("sure: sums, products, quotients and roots within their bounds of the exact values")
  {
    Draws     draws{0x5EED};
    const int rounds = Rounds(20000);
    int       missed = 0;
    for (int round = 0; round < rounds; ++round)
    {
      const auto [one, first]     = Drawn(draws);
      const auto [other, second]  = Drawn(draws);
      missed                     += Within(one + other, Expansion::Sum(first, second)) ? 0 : 1;
      missed += Within(one - other, Expansion::Difference(first, second)) ? 0 : 1;
      missed += Within(one * other, Expansion::Product(first, second)) ? 0 : 1;
      /**
       *  A quotient q of x by y holds the exact one within b where q - b and q + b put q y - x on
       *  either side of zero, y above zero; a root r of x where (r - b)^2 and (r + b)^2 put x
       *  between them.
       */
      const Number            bottom   = {.Value = {.High = std::abs(other.Value.High), .Low = 0.0},
                                          .Bound = other.Bound};
      const Expansion::Number below    = Expansion::Of(bottom.Value.High);
      const Number            quotient = Quotient(one, bottom);
      const auto              times    = [&below, &first](double side, const Number &value) -> int
      {
        const Expansion::Number edge =
            Expansion::Sum(Held(value), Expansion::Of(side * value.Bound));
        return Expansion::Sign(Expansion::Difference(Expansion::Product(edge, below), first));
      };
      missed += times(-1.0, quotient) <= 0 && times(1.0, quotient) >= 0 ? 0 : 1;
      const Number            square = one * one;
      const Number            root   = Root(square);
      const Expansion::Number exact  = Expansion::Product(first, first);
      const auto              power  = [&exact](double side, const Number &value) -> int
      {
        const Expansion::Number edge =
            Expansion::Sum(Held(value), Expansion::Of(side * value.Bound));
        return Expansion::Sign(Expansion::Difference(Expansion::Product(edge, edge), exact));
      };
      missed += power(-1.0, root) <= 0 && power(1.0, root) >= 0 ? 0 : 1;
    }
    CHECK(missed == 0);
  }

  TEST_CASE("sure: a quotient or a root whose bound reaches zero settles nothing")
  {
    /**
     *  Its bound is infinite or not a number, so that no decision made from it holds and the
     *  certificate settles none of its doubles.
     */
    const Number zero = {.Value = {.High = 1e-20, .Low = 0.0}, .Bound = 1e-19};
    CHECK(!Rounded(Quotient(Exact(1.0), zero)));
    CHECK(!(Least(Quotient(Exact(1.0), zero)) > 0.0));
    CHECK(!Rounded(Root(zero)));
    CHECK(!Rounded(Root(Exact(-1.0))));
    CHECK(Least(zero) < 0.0);
    CHECK(Most(zero) > 1e-19);
    CHECK(Rounded(Root(Exact(4.0))) == 2.0);
    CHECK(Rounded(Quotient(Exact(1.0), Exact(4.0)))   == 0.25);
    CHECK(Rounded(std::array{Exact(1.0), Exact(2.0)}) == std::array{1.0, 2.0});
    CHECK(!Rounded(std::array{Exact(1.0), zero}));
  }

  TEST_CASE("sure: the certificate rounds to the nearest double, never across a halfway point")
  {
    /**
     *  Exact values round to themselves; values whose bound leaves their rounding clear, to the
     *  double nearest; a value exactly halfway between two doubles, or one whose bound reaches a
     *  halfway point, is not settled, ties to even or not.
     */
    CHECK(Rounded(Exact(1.5)) == 1.5);
    const Number near = {.Value = {.High = 1.0, .Low = 0x1p-54}, .Bound = 0x1p-80};
    CHECK(Rounded(near) == 1.0);
    const Number halfway = {.Value = {.High = 1.0, .Low = 0x1p-53}, .Bound = 0.0};
    CHECK(!Rounded(halfway));
    const Number wide = {.Value = {.High = 1.0, .Low = 0x1p-54}, .Bound = 0x1p-54};
    CHECK(!Rounded(wide));
    const Number below = {.Value = {.High = 1.0, .Low = -0x1p-55}, .Bound = 0x1p-80};
    CHECK(Rounded(below) == 1.0);
    const Number over = {.Value = {.High = 1.0 + 0x1p-52, .Low = -0x1p-53}, .Bound = 0.0};
    CHECK(!Rounded(over));
    const double infinity = std::numeric_limits<double>::infinity();
    CHECK(!Rounded(Exact(infinity)));
    CHECK(!Rounded({.Value = {.High = 1.0, .Low = 0.0}, .Bound = infinity}));
    CHECK(!Rounded(
        {.Value = {.High = 1.0, .Low = 0.0}, .Bound = std::numeric_limits<double>::quiet_NaN()}));
    CHECK(!Rounded({.Value = {.High = 0.0, .Low = 0.0}, .Bound = Tiny}));
  }
} // namespace Geoxy::Sure
