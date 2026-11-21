/***************************************************************************************************
 *  @file      tests/geoxy/exact/original.hpp
 *  @brief     Shewchuk's own predicates.c, built in C alongside the tests: the oracle every answer
 *             of source/exact must equal, to the bit but for the sign of a zero (exact/two.hpp).
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#ifndef GEOXY_TESTS_EXACT_ORIGINAL_HPP
#define GEOXY_TESTS_EXACT_ORIGINAL_HPP

/**
 *  The names and the types are the C library's own, so they keep its spelling, its global variables
 *  and its leading return types.
 */
// NOLINTBEGIN(readability-identifier-naming)
// NOLINTBEGIN(cppcoreguidelines-avoid-non-const-global-variables)
// NOLINTBEGIN(modernize-use-trailing-return-type)
// NOLINTBEGIN(readability-identifier-length)
extern "C"
{
  void          exactinit();
  double        orient2d(double *pa, double *pb, double *pc);
  double        incircle(double *pa, double *pb, double *pc, double *pd);
  int           fast_expansion_sum_zeroelim(int elen, double *e, int flen, double *f, double *h);
  int           scale_expansion_zeroelim(int elen, double *e, double b, double *h);
  extern double epsilon;
  extern double splitter;
  extern double resulterrbound;
  extern double ccwerrboundA;
  extern double iccerrboundA;
  extern double iccerrboundB;
  extern double iccerrboundC;
}
// NOLINTEND(readability-identifier-length)
// NOLINTEND(modernize-use-trailing-return-type)
// NOLINTEND(cppcoreguidelines-avoid-non-const-global-variables)
// NOLINTEND(readability-identifier-naming)

#endif /* GEOXY_TESTS_EXACT_ORIGINAL_HPP */
