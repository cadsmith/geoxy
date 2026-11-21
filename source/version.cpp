/***************************************************************************************************
 *  @file      source/version.cpp
 *  @brief     The release of geoxy, from the project's version.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "geoxy/version.hpp"

#include <string_view>

namespace Geoxy
{
  auto Version() -> std::string_view
  {
    return GEOXY_VERSION;
  }
} // namespace Geoxy
