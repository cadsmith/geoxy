/***************************************************************************************************
 *  @file      tests/fuzz/oracle.hpp
 *  @brief     What an input of the fuzzer breaks of the exact predicates, against Shewchuk's own.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#ifndef GEOXY_TESTS_FUZZ_ORACLE_HPP
#define GEOXY_TESTS_FUZZ_ORACLE_HPP

#include <array>
#include <cstdint>
#include <span>
#include <string>
#include <string_view>

namespace Geoxy::Fuzz
{
  /**
   *  What a fuzz target drives, each a target geoxy_fuzz_<name> and a folder of inputs.
   */
  enum class Target : std::uint8_t
  {
    Batch,  /**< Exact::OrientSigns, point by point against Exact::Orient. */
    Circle, /**< Exact::InCircle.                                          */
    Orient, /**< Exact::Orient.                                            */
    Two,    /**< The exact sums and products of two doubles.               */
  };

  /**
   *  Every target, each with its name.
   */
  struct Named
  {
    Target           Which; /**< The target.                         */
    std::string_view Name;  /**< Its name, of its binary and folder. */
  };

  /**
   *  The targets, in the order of their names.
   */
  inline constexpr std::array<Named, 4> Targets{
      Named{.Which = Target::Batch , .Name = "batch" },
      Named{.Which = Target::Circle, .Name = "circle"},
      Named{.Which = Target::Orient, .Name = "orient"},
      Named{.Which = Target::Two   , .Name = "two"   },
  };

  /**
   *  What an input breaks: a sign that is not the original's, or a law of the determinant that does
   *  not hold.
   *
   *  Every input runs. The signs of finite places are checked against the original's, scaled into
   *  its range by a power of two where they lie outside it, which keeps both signs; where no
   *  scaling brings them in, against the whole number sums. A NaN or infinite coordinate must give
   *  0. The sums and products of two doubles keep their precondition: both between 2^-400 and
   *  2^400.
   *
   *  @param  target What the input drives.
   *  @param  bytes  The input.
   *  @return        What broke, with the input decoded; empty when nothing did.
   */
  auto Broken(Target target, std::span<const std::uint8_t> bytes) -> std::string;
} // namespace Geoxy::Fuzz

#endif /* GEOXY_TESTS_FUZZ_ORACLE_HPP */
