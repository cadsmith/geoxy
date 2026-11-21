/***************************************************************************************************
 *  @file      geoxy/api.hpp
 *  @brief     GEOXY_API, the mark on what geoxy exports.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#ifndef GEOXY_API_HPP
#define GEOXY_API_HPP

/**
 *  @def GEOXY_API
 *  The mark on every function geoxy exports.
 *
 *  CMake defines GEOXY_STATIC when it builds geoxy as a static library, which exports nothing, so
 *  there is nothing to mark. A shared library marks with `__declspec(dllexport)` while CMake builds
 *  it as `geoxy_EXPORTS` on Windows, and `__declspec(dllimport)` for a program using it; elsewhere
 *  it hides every symbol but those marked with default visibility. A compiler with no symbol
 *  visibility marks nothing.
 */

#ifdef GEOXY_STATIC
#  define GEOXY_API
#else
#  if defined(_WIN32) || defined(__CYGWIN__)
#    ifdef geoxy_EXPORTS
#      define GEOXY_API __declspec(dllexport)
#    else
#      define GEOXY_API __declspec(dllimport)
#    endif
#  else
#    if defined(__GNUC__) && __GNUC__ >= 4
#      define GEOXY_API __attribute__((visibility("default")))
#    else
#      define GEOXY_API
#    endif
#  endif
#endif

#endif /* GEOXY_API_HPP */
