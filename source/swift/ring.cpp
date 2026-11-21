/***************************************************************************************************
 *  @file      source/swift/ring.cpp
 *  @brief     An arc or an ellipse met with a curve: its circle's equation along each piece of the
 *             other.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "swift/ring.hpp"

#include "scratch.hpp"
#include "swift/circuit.hpp"
#include "swift/gather.hpp"
#include "swift/lane.hpp"
#include "swift/rim.hpp"

#include "geoxy/arc.hpp"
#include "geoxy/curve.hpp"
#include "geoxy/intersect.hpp"
#include "geoxy/place.hpp"
#include "geoxy/tolerance.hpp"

#include <algorithm>
#include <array>
#include <cstddef>
#include <optional>
#include <variant>

namespace Geoxy::Swift
{
  auto Ringed(const Curves &curves, const Tolerance &tolerance, Intersections &into)
      -> std::optional<bool>
  {
    const std::array<bool, 2> conics  = {Conic(curves[0]), Conic(curves[1])};
    const std::array<bool, 2> beziers = {Polynomial(curves[0]), Polynomial(curves[1])};
    const bool                arcs    = std::holds_alternative<Arc>(curves[0].get())
                                     && std::holds_alternative<Arc>(curves[1].get());
    const bool paired = (conics[0] && (conics[1] || beziers[1])) || (conics[1] && beziers[0]);
    if (!paired || arcs)
    {
      return std::nullopt;
    }
    /**
     *  The ring is a conic, the second where both are, unless only the first is an arc, whose
     *  equation is a circle's.
     */
    const std::size_t which =
        !conics[1] || (conics[0] && std::holds_alternative<Arc>(curves[0].get())) ? 0 : 1;
    const std::optional<Ring> carrier = Rimmed(curves.at(which));
    if (!carrier
        || !Rung(*carrier, {.Pair = &curves, .Side = 1 - which, .Tolerance = tolerance}, into))
    {
      into.Hits.clear();
      return false;
    }
    std::ranges::sort(into.Hits, {}, [](const Hit &hit) -> double { return hit.First.Along; });
    return true;
  }

  auto Humped(const Hump &hump, const Curve &other, const Tolerance &tolerance)
      -> std::optional<Spots>
  {
    const std::optional<Ring>        carrier = Encircled(hump);
    const std::optional<List<Found>> found =
        carrier ? Gathered(*carrier, other, tolerance) : std::nullopt;
    if (!found)
    {
      return std::nullopt;
    }
    Spots spots;
    for (const Found &crossing : *found)
    {
      spots.push_back({.At    = crossing.At,
                       .Along = {static_cast<double>(hump.Index) + 0.5, crossing.Along},
                       .Own   = crossing.Own});
    }
    return spots;
  }
} // namespace Geoxy::Swift
