/***************************************************************************************************
 *  @file      tests/fuzz/two.cpp
 *  @brief     The libFuzzer target geoxy_fuzz_two: the exact sums and products of two doubles
 *             against Shewchuk's routines.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "oracle.hpp"

#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <print>
#include <span>
#include <string>

/**
 *  Runs one input, and stops the fuzzer on what it breaks, telling it.
 *
 *  @param  data The input.
 *  @param  size Its length.
 *  @return      0, the input kept or not as libFuzzer decides.
 */
extern "C" auto LLVMFuzzerTestOneInput(const std::uint8_t *data, std::size_t size) -> int
{
  const std::string broken = Geoxy::Fuzz::Broken(Geoxy::Fuzz::Target::Two, {data, size});
  if (!broken.empty())
  {
    std::println(stderr, "{}", broken);
    std::abort();
  }
  return 0;
}
