/***************************************************************************************************
 *  @file      source/quadrature.hpp
 *  @brief     How hard an integral may be worked at.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#ifndef GEOXY_QUADRATURE_HPP
#define GEOXY_QUADRATURE_HPP

#include "geoxy/accuracy.hpp"

#include <cstddef>

namespace Geoxy
{
  /**
   *  How hard an integral may be worked at.
   */
  struct Effort
  {
    Accuracy    Within   = {};     /**< How far off the integral may be.      */
    std::size_t Halvings = 65'536; /**< How many halvings may be spent on it. */
  };
} // namespace Geoxy

#endif /* GEOXY_QUADRATURE_HPP */
