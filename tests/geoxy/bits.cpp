/***************************************************************************************************
 *  @file      tests/geoxy/bits.cpp
 *  @brief     The checksum that holds every platform to the same bits.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "draws.hpp"

#include <doctest/doctest.h>

#include <bit>
#include <cmath>
#include <cstdint>

namespace Geoxy
{
  namespace
  {
    /**
     *  FNV-1a over the bits of each value.
     */
    class Checksum
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

    constexpr int BitSamples = 1'000'000;
  } // namespace

  /**
   *  The checksum is one number on every platform only while no compiler fuses a multiply and an
   *  add, and the basic operations round as IEEE 754 says. A failure here means a build flag
   *  changed, not the code.
   */
  TEST_CASE("bits: the basic operations give the same bits on every platform")
  {
    Draws draws{0x6765'6F6DU};

    Checksum sum;
    for (int index = 0; index < BitSamples; ++index)
    {
      const double ax = draws.Near(5e6, 1e3);
      const double ay = draws.Near(-3e6, 1e3);
      const double bx = draws.Near(0, 1);
      const double by = draws.Near(0, 1e-3);
      sum.Add((ax * by) - (ay * bx));
      sum.Add((ax * bx) + (ay * by));
      sum.Add(ax / by);
      sum.Add(std::sqrt((bx * bx) + (by * by)));
      sum.Add((ax - ay) * (bx + by));
    }
    CHECK(sum.Value() == 2682263569248122725U);
  }
} // namespace Geoxy
