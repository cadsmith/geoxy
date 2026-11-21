/***************************************************************************************************
 *  @file      source/exact/strict.hpp
 *  @brief     The floating point arithmetic the exact predicates are proved for, enforced.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#ifndef GEOXY_EXACT_STRICT_HPP
#define GEOXY_EXACT_STRICT_HPP

#include <cfloat>
#include <limits>

/**
 *  The proofs of tests/proof take every operation of the predicates to be IEEE 754 binary64, each
 *  product and sum rounded to nearest on its own, in the order written. A compiler may break that
 *  three ways, each changing a sign the proofs promise: fusing a product into a sum, one rounding
 *  where two are written, as Clang does on arm64 by default; reassociating a sum; or holding a
 *  value in more precision than a double, as the x87 unit does. Fast math does all of them and
 *  assumes no infinity or NaN arrives, which the predicates handle.
 *
 *  So the exact sources hold to it themselves, whatever the flags of the build: they include this
 *  header, through two.hpp and bounds.hpp, which turns fusing and reassociation off for the rest of
 *  the file, and stops the build where the arithmetic is not IEEE 754 binary64 without excess
 *  precision, or where fast math is on. Clang's -ffp-contract=fast alone overrides the pragma;
 *  geoxy's CMake passes -ffp-contract=off after a user's flags, and tests/geoxy/bits.cpp fails on
 *  any contraction. tests/proof/machine.py reads the machine code the compiler emits, on arm64 and
 *  x86-64, and finds no fused product the source did not ask for.
 */

#if defined(__FAST_MATH__) || (defined(__FINITE_MATH_ONLY__) && __FINITE_MATH_ONLY__ > 0)
#  error "geoxy's exact predicates need IEEE 754 arithmetic: build without fast or finite math"
#endif

#ifdef __clang__
#  pragma clang fp contract(off)
#  pragma clang fp reassociate(off)
#elifdef __GNUC__
#  pragma GCC optimize("fp-contract=off")
#elifdef _MSC_VER
#  pragma fp_contract(off)
#  pragma float_control(precise, on)
#endif

static_assert(std::numeric_limits<double>::is_iec559,
              "geoxy's exact predicates need IEEE 754 binary64 doubles");
static_assert(std::numeric_limits<double>::digits == 53,
              "geoxy's exact predicates need IEEE 754 binary64 doubles");
static_assert(FLT_EVAL_METHOD == 0,
              "geoxy's exact predicates need each double operation rounded to a double");

#endif /* GEOXY_EXACT_STRICT_HPP */
