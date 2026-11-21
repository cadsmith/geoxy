/***************************************************************************************************
 *  @file      source/traced.hpp
 *  @brief     How tracing orders the halves of an arrangement's edges about a vertex, for what
 *             reads its loops with no face traced.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#ifndef GEOXY_TRACED_HPP
#define GEOXY_TRACED_HPP

#include "arranged.hpp"
#include "scratch.hpp"

#include "geoxy/vector.hpp"

#include <cstddef>
#include <limits>
#include <span>

namespace Geoxy
{
  /**
   *  The least angle within which two tangents are taken as one, some steps of a double.
   */
  inline constexpr double LeastAngle = 64.0 * std::numeric_limits<double>::epsilon();

  /**
   *  A half leaving a vertex: the way it leaves, and how fast it turns left.
   */
  struct Spoke
  {
    Vector      Way  = {};  /**< Its tangent where it leaves.            */
    double      Curl = 0.0; /**< Its curvature there, positive leftward. */
    std::size_t Half = 0;   /**< The half.                               */
  };

  /**
   *  How two ways compare counterclockwise, as Order compares them: by their halves of the plane,
   *  told here, and in one half by Order, which no call is made to for most pairs.
   *
   *  @param  lhs One way.
   *  @param  rhs The other.
   *  @return     Below 0 where the first comes first, above where the second does, else 0.
   */
  auto Preceding(Vector lhs, Vector rhs) -> int;

  /**
   *  The halves leaving a vertex of three or more in counterclockwise order: by their exact
   *  tangents and their curvature, as Fanned orders them, which orders them wholly, so that the
   *  order they are given in is no matter.
   *
   *  @param arranged The arrangement.
   *  @param around   The vertex's halves leaving it, put in order.
   *  @param spokes   Room for their spokes.
   */
  void Ordered(const Arranged &arranged, std::span<std::size_t> around, List<Spoke> &spokes);
} // namespace Geoxy

#endif /* GEOXY_TRACED_HPP */
