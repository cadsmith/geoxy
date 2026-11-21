/***************************************************************************************************
 *  @file      geoxy/error.hpp
 *  @brief     Why an operation that can fail gave no answer.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#ifndef GEOXY_ERROR_HPP
#define GEOXY_ERROR_HPP

#include <cstdint>

namespace Geoxy
{
  /**
   *  Why an operation that can fail gave no answer: a short list, shared by every one of them, so
   *  that a caller tells the user why in one place.
   */
  enum class Error : std::uint8_t
  {
    Degenerate,    /**< What was given cannot have an answer: a loop that does not close. */
    NotFinite,     /**< A number given is not finite.                                     */
    TooLarge,      /**< A size asked for is past what the shapes allow.                   */
    NoConvergence, /**< An approximation did not reach its accuracy.                      */
    Spent,         /**< A work budget ran out, so a part of the answer may be missing.    */
  };
} // namespace Geoxy

#endif /* GEOXY_ERROR_HPP */
