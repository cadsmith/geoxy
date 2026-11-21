/***************************************************************************************************
 *  @file      source/region.cpp
 *  @brief     Regions: where a point lies against one, its area, and the booleans of two, read from
 *             the faces of their loops arranged together.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "geoxy/region.hpp"

#include "arranged.hpp"
#include "exact/two.hpp"
#include "scratch.hpp"

#include "geoxy/accuracy.hpp"
#include "geoxy/at.hpp"
#include "geoxy/curve.hpp"
#include "geoxy/domain.hpp"
#include "geoxy/error.hpp"
#include "geoxy/loop.hpp"
#include "geoxy/point.hpp"
#include "geoxy/tolerance.hpp"

#include <array>
#include <cstddef>
#include <cstdint>
#include <expected>
#include <optional>
#include <utility>
#include <vector>

namespace Geoxy
{
  namespace
  {
    /**
     *  Two regions, the first first.
     */
    using Pair = std::array<const Region *, 2>;

    /**
     *  Whether each edge of a loop ends within a tolerance of where the next starts: what Closed
     *  tells of a loop whose edges are known to have ends and finite numbers, as those of loops an
     *  arrangement was made of without a fault have.
     *
     *  @param  loop      The loop.
     *  @param  tolerance The tolerance.
     *  @return           Whether they do, none for a loop of no edge.
     */
    auto Sealed(const Loop &loop, const Tolerance &tolerance) -> bool
    {
      const std::size_t count = loop.Edges.size();
      for (std::size_t index = 0; index < count; ++index)
      {
        const Curve &edge = loop.Edges[index];
        const Curve &next = loop.Edges[(index + 1) % count];
        if (!Same(At(edge, Domain(edge).End), At(next, Domain(next).Start), tolerance))
        {
          return false;
        }
      }
      return count > 0;
    }

    /**
     *  A boolean of two regions: their loops arranged together, each weighing one along its run in
     *  its own operand, and the faces kept that the rule keeps by where each operand winds about
     *  them; whether the loops close asked only once they are known to be lines and arcs.
     *
     *  @param  regions   The two regions.
     *  @param  keeps     The rule.
     *  @param  tolerance The tolerance.
     *  @return           The region, or why there is none.
     */
    auto Combined(const Pair &regions, Keeps keeps, const Tolerance &tolerance)
        -> std::expected<Region, Error>
    {
      Inputs inputs = {.Runs    = {},
                       .Tallies = {},
                       .Memory  = &Scratchpad::Shared(),
                       .Chained = true,
                       .Spare   = true};
      for (std::size_t operand = 0; operand < regions.size(); ++operand)
      {
        const Tally tally = {operand == 0 ? 1 : 0, operand == 0 ? 0 : 1};
        for (const Loop &loop : regions.at(operand)->Loops)
        {
          inputs.Runs.emplace_back(loop.Edges);
          inputs.Tallies.push_back(tally);
        }
      }
      Arranged arranged = Arrange(inputs, tolerance);
      if (arranged.Fault)
      {
        return std::unexpected(*arranged.Fault);
      }
      for (const Region *region : regions)
      {
        for (const Loop &loop : region->Loops)
        {
          if (!Sealed(loop, tolerance))
          {
            return std::unexpected(Error::Degenerate);
          }
        }
      }
      /**
       *  Read with no face traced where the windings ferried agree; where they disagree the faces
       *  are no plane's, and the split arrangement is read so where it has one and its windings
       *  agree, as Mended would choose it. Else the faces traced and mended.
       *
       *  A cluster of vertices merged by distance through which two lines pass that meet nowhere in
       *  it leaves faces no plane's, whose windings are not a function of the face: the arrangement
       *  with every such cluster split into its points, made from the same cuts, is kept where it
       *  is a plane's. Where it is not, or none is twisted, the first with the pieces that cross
       *  between two vertices made one, as Melded makes them, where that one is.
       */
      Region      region = {};
      const Ferry ferry  = Ferried(arranged, keeps, region);
      if (ferry == Ferry::Done
          || (ferry == Ferry::Twisted
              && arranged.Split
              && Ferried(*arranged.Split, keeps, region) == Ferry::Done))
      {
        return region;
      }
      Trace(arranged);
      const Arranged           chosen   = Mended(std::move(arranged), false);
      const std::vector<Tally> windings = Windings(chosen);
      std::vector<bool>        kept(windings.size());
      for (std::size_t index = 0; index < windings.size(); ++index)
      {
        kept[index] = keeps(windings[index][0] != 0, windings[index][1] != 0);
      }
      return Stitched(chosen, kept);
    }
  } // namespace

  auto Locate(const Region &region, Point point, const Tolerance &tolerance) -> Location
  {
    std::int64_t total = 0;
    for (const Loop &loop : region.Loops)
    {
      const std::optional<int> winding = Winding(loop, point, tolerance);
      if (!winding)
      {
        return Location::On;
      }
      total += *winding;
    }
    return total != 0 ? Location::Inside : Location::Outside;
  }

  auto Area(const Region &region, Accuracy accuracy) -> double
  {
    double big   = 0.0;
    double small = 0.0;
    for (const Loop &loop : region.Loops)
    {
      const Exact::Two sum  = Exact::TwoSum(big, Area(loop, accuracy));
      big                   = sum.Big;
      small                += sum.Small;
    }
    return big + small;
  }

  auto Filled(const Region &region, const Tolerance &tolerance) -> std::expected<Region, Error>
  {
    const Region none = {};
    return Combined({&region, &none}, [](bool first, bool) -> bool { return first; }, tolerance);
  }

  auto Union(const Region &first, const Region &second, const Tolerance &tolerance)
      -> std::expected<Region, Error>
  {
    return Combined(
        {&first, &second}, [](bool inFirst, bool inSecond) -> bool { return inFirst || inSecond; },
        tolerance);
  }

  auto Difference(const Region &first, const Region &second, const Tolerance &tolerance)
      -> std::expected<Region, Error>
  {
    return Combined(
        {&first, &second}, [](bool inFirst, bool inSecond) -> bool { return inFirst && !inSecond; },
        tolerance);
  }

  auto Intersection(const Region &first, const Region &second, const Tolerance &tolerance)
      -> std::expected<Region, Error>
  {
    return Combined(
        {&first, &second}, [](bool inFirst, bool inSecond) -> bool { return inFirst && inSecond; },
        tolerance);
  }
} // namespace Geoxy
