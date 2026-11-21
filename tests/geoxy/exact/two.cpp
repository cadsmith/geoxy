/***************************************************************************************************
 *  @file      tests/geoxy/exact/two.cpp
 *  @brief     Tests of exact sums and products of two doubles, against Shewchuk's own arithmetic.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "exact/two.hpp"

#include "../draws.hpp"
#include "original.hpp"
#include "same.hpp"

#include <doctest/doctest.h>

#include <array>
#include <cmath>
#include <cstddef>
#include <vector>

namespace Geoxy::Exact
{
  namespace
  {
    constexpr int Count = 20'000;

    /**
     *  A double with its own sign and exponent, from 2^-60 to 2^60.
     */
    auto Any(Draws &draws) -> double
    {
      const double mantissa = draws.Near(1.5, 0.5);
      const auto   exponent = static_cast<int>(draws.Next() % 121U) - 60;
      const double value    = std::ldexp(mantissa, exponent);
      return (draws.Next() % 2U == 0) ? value : -value;
    }

    /**
     *  An expansion, as the original reads one: a copy one longer than it is, since it reads one
     *  past the end of each input.
     */
    auto Padded(std::vector<double> expansion) -> std::vector<double>
    {
      expansion.push_back(0.0);
      return expansion;
    }

    /**
     *  The exact sum of two expansions, by the original.
     */
    auto OriginalSum(const std::vector<double> &first, const std::vector<double> &second)
        -> std::vector<double>
    {
      std::vector<double> one = Padded(first);
      std::vector<double> two = Padded(second);
      std::vector<double> sum(first.size() + second.size());
      const int           length =
          fast_expansion_sum_zeroelim(static_cast<int>(first.size()), one.data(),
                                      static_cast<int>(second.size()), two.data(), sum.data());
      sum.resize(static_cast<std::size_t>(length));
      return sum;
    }

    /**
     *  The exact product of an expansion and a double, by the original.
     */
    auto OriginalScale(const std::vector<double> &expansion, double factor) -> std::vector<double>
    {
      std::vector<double> one = Padded(expansion);
      std::vector<double> product(2 * expansion.size());
      const int length = scale_expansion_zeroelim(static_cast<int>(expansion.size()), one.data(),
                                                  factor, product.data());
      product.resize(static_cast<std::size_t>(length));
      return product;
    }

    /**
     *  Whether two expansions hold the same value: their difference, by the original, is one zero.
     */
    auto Equal(const std::vector<double> &ours, std::vector<double> original) -> bool
    {
      for (double &component : original)
      {
        component = -component;
      }
      const std::vector<double> difference = OriginalSum(ours, original);
      return difference.size() == 1 && IsZero(difference.front());
    }

    auto Of(Two value) -> std::vector<double>
    {
      return {value.Small, value.Big};
    }

    template <std::size_t SizeT>
    auto Of(const std::array<double, SizeT> &value) -> std::vector<double>
    {
      return {value.begin(), value.end()};
    }

    struct TwoFixture
    {
      TwoFixture() { exactinit(); }
    };
  } // namespace

  TEST_CASE("two: zero and sign")
  {
    CHECK(IsZero(0.0));
    CHECK(IsZero(-0.0));
    CHECK(!IsZero(0x1.0p-1074));
    CHECK(!IsZero(-1.0));
    CHECK(Sign(2.0)  == 1);
    CHECK(Sign(-2.0) == -1);
    CHECK(Sign(0.0)  == 0);
    CHECK(Sign(-0.0) == 0);
  }

  TEST_CASE_FIXTURE(TwoFixture,
                    "two: sums and differences are exact, their parts rounded as summed")
  {
    Draws draws{0x7375'6D73U};

    int exact = 0;
    for (int index = 0; index < Count; ++index)
    {
      const double a = Any(draws);
      const double b = Any(draws);

      const Two sum        = TwoSum(a, b);
      const Two difference = TwoDiff(a, b);

      exact += static_cast<int>(Bits(sum.Big) == Bits(a + b)
                                && Bits(difference.Big) == Bits(a - b)
                                && Equal(Of(sum), OriginalSum({a}, {b}))
                                && Equal(Of(difference), OriginalSum({a}, {-b})));
    }
    CHECK(exact == Count);
  }

  TEST_CASE_FIXTURE(TwoFixture, "two: products and squares are exact")
  {
    Draws draws{0x7072'6F64U};

    int exact = 0;
    for (int index = 0; index < Count; ++index)
    {
      const double a = Any(draws);
      const double b = Any(draws);

      const Two product  = TwoProduct(a, b);
      const Two square   = Square(a);
      exact             += static_cast<int>(Bits(product.Big) == Bits(a * b)
                                            && Equal(Of(product), OriginalScale({a}, b))
                                            && Bits(square.Big) == Bits(TwoProduct(a, a).Big)
                                            && Same(square.Small, TwoProduct(a, a).Small));
    }
    CHECK(exact == Count);
  }

  TEST_CASE("two: the roundoff of a difference that is exact is zero")
  {
    CHECK(IsZero(TwoDiff(1.5, 0.25).Small));
    CHECK(IsZero(TwoDiffTail(3.0, 1.0, 2.0)));
    CHECK(Bits(TwoDiff(1.0, 0x1.0p-60).Small) == Bits(-0x1.0p-60));
  }
} // namespace Geoxy::Exact
