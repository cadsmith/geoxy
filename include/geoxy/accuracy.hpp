/***************************************************************************************************
 *  @file      geoxy/accuracy.hpp
 *  @brief     How close an approximation must come, for what is worked out only approximately.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#ifndef GEOXY_ACCURACY_HPP
#define GEOXY_ACCURACY_HPP

namespace Geoxy
{
  /**
   *  How close an approximation must come: the most it may be off, in drawing units.
   *
   *  What is worked out exactly takes none. What is approximated, chords along a curve or the
   *  length of a spline, takes one, and does more work the finer it is, within a budget that bounds
   *  the work whatever is asked. It is not a tolerance: it says how well to answer, not when two
   *  things count as one.
   */
  struct Accuracy
  {
    double Error = 0.0; /**< The most the answer may be off; at or below 0, as near as can be. */
  };
} // namespace Geoxy

#endif /* GEOXY_ACCURACY_HPP */
