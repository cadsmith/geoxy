/***************************************************************************************************
 *  @file      geoxy/version.hpp
 *  @brief     The release of geoxy a library was built from.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#ifndef GEOXY_VERSION_HPP
#define GEOXY_VERSION_HPP

#include "geoxy/api.hpp"

#include <string_view>

namespace Geoxy
{
  /**
   *  The release of geoxy this library was built from.
   *
   *  @return The version, such as "0.1.0".
   */
  [[nodiscard]] GEOXY_API auto Version() -> std::string_view;
} // namespace Geoxy

#endif /* GEOXY_VERSION_HPP */
