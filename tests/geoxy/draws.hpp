/***************************************************************************************************
 *  @file      tests/geoxy/draws.hpp
 *  @brief     Seeded draws, the same on every platform.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#ifndef GEOXY_TESTS_DRAWS_HPP
#define GEOXY_TESTS_DRAWS_HPP

#include <algorithm>
#include <cstdint>

namespace Geoxy
{
  /**
   *  How many rounds a randomized test runs: all of them in a build that is optimized, as the
   *  release, sanitize and coverage presets are; one in sixteen where it is not, as the debug
   *  preset is, where a round costs seventeen times as much and the others run them all anyway.
   *
   *  @param  full The rounds an optimized build runs.
   *  @return      The rounds this build runs, at least one.
   */
  constexpr auto Rounds(int full) -> int
  {
#if defined(__OPTIMIZE__) || (defined(_MSC_VER) && defined(NDEBUG))
    return full;
#else
    return std::max(1, full / 16);
#endif
  }

  /**
   *  SplitMix64: the same draws on every platform, unlike the distributions of <random>.
   */
  class Draws
  {
  public:
    explicit Draws(std::uint64_t seed) : m_State(seed) {}

    auto Next() -> std::uint64_t
    {
      m_State             += 0x9E3779B97F4A7C15U;
      std::uint64_t mixed  = m_State;
      mixed                = (mixed ^ (mixed >> 30U)) * 0xBF58476D1CE4E5B9U;
      mixed                = (mixed ^ (mixed >> 27U)) * 0x94D049BB133111EBU;
      return mixed ^ (mixed >> 31U);
    }

    /**
     *  A double in [offset - spread, offset + spread), from 53 random bits, exactly.
     */
    auto Near(double offset, double spread) -> double
    {
      constexpr double Unit = 0x1.0p-53;

      const auto bits = static_cast<double>(Next() >> 11U);
      return offset + (((2 * bits * Unit) - 1) * spread);
    }

  private:
    std::uint64_t m_State;
  };
} // namespace Geoxy

#endif /* GEOXY_TESTS_DRAWS_HPP */
