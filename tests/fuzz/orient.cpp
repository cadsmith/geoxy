/***************************************************************************************************
 *  @file      tests/fuzz/orient.cpp
 *  @brief     The libFuzzer target geoxy_fuzz_orient: Exact::Orient against Shewchuk's orient2d,
 *             and the laws of its determinant.
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
  const std::string broken = Geoxy::Fuzz::Broken(Geoxy::Fuzz::Target::Orient, {data, size});
  if (!broken.empty())
  {
    std::println(stderr, "{}", broken);
    std::abort();
  }
  return 0;
}
