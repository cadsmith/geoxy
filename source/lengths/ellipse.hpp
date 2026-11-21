/***************************************************************************************************
 *  @file      source/lengths/ellipse.hpp
 *  @brief     The length along an ellipse, by Carlson's elliptic integrals: each whole turn of it,
 *             and what is left.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#ifndef GEOXY_LENGTHS_ELLIPSE_HPP
#define GEOXY_LENGTHS_ELLIPSE_HPP

#include "geoxy/domain.hpp"
#include "geoxy/ellipse.hpp"

#include <cstddef>

namespace Geoxy::Lengths
{
  /**
   *  The length of a stretch of an ellipse, a times E(u) at its end less E(u) at its start, by
   *  Carlson's integrals in double-double, rounded once.
   *
   *  Its ends are its angles held as pairs of degrees: the start's whole turns taken off, the
   *  sweep's kept, as each is held. Each sum and product of pairs, each root and quotient, is
   *  within some 2^-102 of itself, and the cosine and sine of an end within 2^-68; the half turns
   *  of E(m) are counted, and what is left between the ends is one integral by the addition
   *  theorem, so the length is within some 2^-66 of itself and is the double nearest it but within
   *  that of a halfway point. A map of no area runs along a segment and back, m = 1, where E(v) =
   *  sin v.
   *
   *  @param  ellipse The ellipse.
   *  @param  stretch The stretch of its parameter.
   *  @return         Its length, positive whichever way it runs; nothing for a map of nothing or
   *                  one not finite; not a number for an angle that is not finite.
   */
  auto Arclength(const Ellipse &ellipse, Range stretch) -> double;
} // namespace Geoxy::Lengths

#endif /* GEOXY_LENGTHS_ELLIPSE_HPP */
