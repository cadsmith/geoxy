/***************************************************************************************************
 *  @file      source/arrangement.cpp
 *  @brief     Curves arranged once, their gaps bridged, and the loop around a picked point among
 *             them with its islands.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "geoxy/arrangement.hpp"

#include "arranged.hpp"

#include "geoxy/box.hpp"
#include "geoxy/curve.hpp"
#include "geoxy/error.hpp"
#include "geoxy/line.hpp"
#include "geoxy/nearest.hpp"
#include "geoxy/point.hpp"
#include "geoxy/region.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <expected>
#include <memory>
#include <optional>
#include <span>
#include <tuple>
#include <utility>
#include <vector>

namespace Geoxy
{
  /**
   *  What an arrangement keeps, as Boundary reads it.
   */
  struct Arranging
  {
    /**
     *  An arrangement's faces.
     *
     *  @param  arrangement The arrangement.
     *  @return             Its faces; none for one moved from.
     */
    static auto Faces(const Arrangement &arrangement) -> const Arranged *
    {
      return arrangement.m_Faces.get();
    }
  };

  namespace
  {
    /**
     *  The lines that bridge an arrangement's gaps up to a width: between two loose ends, the
     *  nearest pairs first, each end once; then from each loose end left to the nearest point of an
     *  edge that does not reach it.
     *
     *  @param  arranged The arrangement, not traced.
     *  @param  gap      The widest gap.
     *  @return          The lines.
     */
    auto Bridges(const Arranged &arranged, double gap) -> std::vector<Curve>
    {
      std::vector<std::size_t> loose = Frayed(arranged);
      const auto               spot  = [&arranged](std::size_t vertex) -> Point
      {
        return arranged.Points[vertex];
      };
      std::ranges::sort(loose, {},
                        [&spot](std::size_t vertex) -> double { return spot(vertex).X; });
      std::vector<std::tuple<double, std::size_t, std::size_t>> pairs;
      for (std::size_t first = 0; first < loose.size(); ++first)
      {
        for (std::size_t second = first + 1;
             second < loose.size() && spot(loose[second]).X - spot(loose[first]).X <= gap; ++second)
        {
          const double apart = Distance(spot(loose[first]), spot(loose[second]));
          if (apart <= gap)
          {
            pairs.emplace_back(apart, std::min(loose[first], loose[second]),
                               std::max(loose[first], loose[second]));
          }
        }
      }
      std::ranges::sort(pairs);
      std::vector<bool>  bridged(arranged.Points.size(), false);
      std::vector<Curve> bridges;
      for (const auto &[apart, first, second] : pairs)
      {
        if (!bridged[first] && !bridged[second])
        {
          bridged[first]  = true;
          bridged[second] = true;
          bridges.emplace_back(Line{.Start = spot(first), .End = spot(second)});
        }
      }
      for (const std::size_t vertex : loose)
      {
        if (bridged[vertex])
        {
          continue;
        }
        const Point          point = spot(vertex);
        double               best  = gap;
        std::optional<Point> foot;
        for (std::size_t index = 0; index < arranged.Edges.size(); ++index)
        {
          const Arranged::Edge &edge = arranged.Edges[index];
          if (edge.Ends[0] == vertex
              || edge.Ends[1] == vertex
              || !Contains(Padded(arranged.Boxes[index], gap), point))
          {
            continue;
          }
          const Point  near  = Nearest(arranged.Shapes[index], point).At;
          const double apart = Distance(near, point);
          if (apart <= best)
          {
            best = apart;
            foot = near;
          }
        }
        if (foot)
        {
          bridges.emplace_back(Line{.Start = point, .End = *foot});
        }
      }
      return bridges;
    }
  } // namespace

  Arrangement::Arrangement(std::span<const Curve> curves, const Closing &closing)
  {
    Inputs             inputs = {.Runs = {curves}, .Tallies = {Tally{}}, .Spare = true};
    std::vector<Curve> bridges;
    /**
     *  The curves arranged, and arranged again with the bridges across their gaps where there are
     *  any.
     */
    const auto bridged = [&inputs, &bridges, &closing] -> Arranged
    {
      Arranged first = Arrange(inputs, closing.Within);
      if (first.Fault || !(closing.Gap > first.Within.Length))
      {
        return first;
      }
      Outfitted(first);
      bridges = Bridges(first, closing.Gap);
      if (bridges.empty())
      {
        return first;
      }
      inputs.Runs.emplace_back(bridges);
      inputs.Tallies.push_back(Tally{});
      return Arrange(inputs, closing.Within);
    };
    Arranged arranged = bridged();
    if (arranged.Fault)
    {
      m_Faces = std::make_shared<const Arranged>(std::move(arranged));
      return;
    }
    Prune(arranged);
    Trace(arranged);
    /**
     *  Faces no plane's mended as a boolean mends them; the split one, once chosen or not, dropped.
     */
    Arranged chosen = Mended(std::move(arranged), true);
    chosen.Split.reset();
    Outfitted(chosen);
    m_Faces = std::make_shared<const Arranged>(std::move(chosen));
  }

  auto Boundary(const Arrangement &arrangement, Point point) -> std::expected<Region, Error>
  {
    const Arranged *arranged = Arranging::Faces(arrangement);
    if (arranged == nullptr)
    {
      return Region{};
    }
    if (arranged->Fault)
    {
      return std::unexpected(*arranged->Fault);
    }
    if (!std::isfinite(point.X) || !std::isfinite(point.Y))
    {
      return std::unexpected(Error::NotFinite);
    }
    const std::size_t          none = arranged->Samples.size();
    std::optional<std::size_t> face;
    if (!Grazes(*arranged, point))
    {
      face = Enclosing(*arranged, point, none);
    }
    if (!face)
    {
      return Region{};
    }
    Region            region = {.Loops = {Ringed(*arranged, *face)}};
    const std::size_t home   = arranged->Cycles[*face].Component;
    for (std::size_t component = 0; component < none; ++component)
    {
      const Point sample = arranged->Points[arranged->Samples[component]];
      if (component != home
          && Contains(arranged->Cycles[*face].Bounds, sample)
          && Enclosing(*arranged, sample, component) == face)
      {
        region.Loops.push_back(Ringed(*arranged, arranged->Outer[component]));
      }
    }
    return region;
  }
} // namespace Geoxy
