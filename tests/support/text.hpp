/***************************************************************************************************
 *  @file      tests/support/text.hpp
 *  @brief     Curves as one line of text, exactly, to write a case down and replay it.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#ifndef GEOXY_TEXT_HPP
#define GEOXY_TEXT_HPP

#include "geoxy/curve.hpp"

#include <optional>
#include <string>
#include <string_view>

namespace Geoxy
{
  /**
   *  A curve as one line of text, its kind then its numbers, that reads back as the same curve to
   *  the bit.
   *
   *  Each number is the shortest decimal that reads back as its double; an angle is in radians, or
   *  in degrees with a d after it, as it is held:
   *
   *  - `line x y x y`, its start then its end;
   *  - `ray x y x y` and `xline x y x y`, its point then its way;
   *  - `arc x y radius start sweep`;
   *  - `ellipse x y x y x y start sweep`, its map's x axis, y axis and origin;
   *  - `polyline open|closed count x y bulge...`;
   *  - `spline degree count knot... count x y... count weight...`;
   *  - `quadratic x y x y x y` and `cubic x y x y x y x y`, its points in order.
   *
   *  @param  curve The curve.
   *  @return       The text.
   */
  [[nodiscard]] auto ToText(const Curve &curve) -> std::string;

  /**
   *  The curve a line of text holds, as ToText writes one.
   *
   *  @param  text The text, its words apart by spaces.
   *  @return      The curve, or none when the text is not one: an unknown kind, a word that is not
   *               a number, a count past the words, or words left over.
   */
  [[nodiscard]] auto FromText(std::string_view text) -> std::optional<Curve>;
} // namespace Geoxy

#endif /* GEOXY_TEXT_HPP */
