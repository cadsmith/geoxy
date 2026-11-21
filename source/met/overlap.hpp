/***************************************************************************************************
 *  @file      source/met/overlap.hpp
 *  @brief     The stretches two curves share: snapped to the ends of parts, settled as kept or as a
 *             band about a touch, and joined where one runs on into the next.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#ifndef GEOXY_MET_OVERLAP_HPP
#define GEOXY_MET_OVERLAP_HPP

#include "met/judge.hpp"
#include "met/parts.hpp"
#include "scratch.hpp"

#include "geoxy/intersect.hpp"

#include <array>
#include <cstddef>
#include <span>

namespace Geoxy::Met
{
  /**
   *  Shared stretches that run on into each other joined, in order along the first curve.
   *
   *  @param  overlaps The stretches.
   *  @return          The stretches joined.
   */
  auto Joined(List<Overlap> overlaps) -> List<Overlap>;

  /**
   *  Where each curve's parts end: as parameters of the first, those of the second's that lie
   *  within the tolerance of it taken to its nearest point; and infinity either way, so that every
   *  place has a boundary beyond it.
   *
   *  @param  scene The curves being met.
   *  @param  parts Each curve's parts.
   *  @return       The parameters, in order.
   */
  auto Boundaries(const Scene &scene, std::array<const List<Part> *, 2> parts) -> List<double>;

  /**
   *  A shared stretch with each end taken out to the nearest boundary beyond it, where the curves
   *  stay within the tolerance all the way: stretches found piece by piece stop short of where a
   *  piece ends by as much as the part that was too flat to halve again.
   *
   *  @param  scene   The curves being met.
   *  @param  bounds  Where the curves' parts end, as parameters of the first, in order.
   *  @param  overlap The stretch.
   *  @return         The stretch, its ends taken out where they reach.
   */
  auto Snapped(const Scene &scene, const List<double> &bounds, Overlap overlap) -> Overlap;

  /**
   *  Whether a shared stretch is one: each of its ends at a boundary, where a part of either curve
   *  ends, a curve's end, a vertex, a knot or a closed curve's seam, or out at infinity; or where
   *  either curve all but stalls, as where one turns back on itself. A stretch that ends where the
   *  curves part is the band where they stay within the tolerance about a touch or a shallow
   *  crossing, and is a hit; whatever the curves are.
   *
   *  @param  scene   The curves being met.
   *  @param  bounds  Where the curves' parts end, as parameters of the first, in order.
   *  @param  overlap The stretch.
   *  @return         Whether it is shared.
   */
  auto Settled(const Scene &scene, const List<double> &bounds, const Overlap &overlap) -> bool;

  /**
   *  The hits within a band where the curves are tangent made touches: within its stretch of the
   *  first curve, or within the tolerance of either end, as a hit polished a step of a double past
   *  it is.
   *
   *  @param scene The curves being met.
   *  @param bands The bands.
   *  @param hits  The hits.
   */
  void Touched(const Scene &scene, std::span<const Overlap> bands, std::span<Hit> hits);
} // namespace Geoxy::Met

#endif /* GEOXY_MET_OVERLAP_HPP */
