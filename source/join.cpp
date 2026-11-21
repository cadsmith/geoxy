/***************************************************************************************************
 *  @file      source/join.cpp
 *  @brief     Curves end to end, linked into one polyline or one spline.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "geoxy/join.hpp"

#include "chain.hpp"
#include "kernel/piece.hpp"
#include "trig.hpp"

#include "geoxy/angle.hpp"
#include "geoxy/arc.hpp"
#include "geoxy/at.hpp"
#include "geoxy/curve.hpp"
#include "geoxy/line.hpp"
#include "geoxy/point.hpp"
#include "geoxy/polyline.hpp"
#include "geoxy/reverse.hpp"
#include "geoxy/segment.hpp"
#include "geoxy/spline.hpp"
#include "geoxy/tolerance.hpp"
#include "geoxy/vector.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <optional>
#include <span>
#include <variant>
#include <vector>

namespace Geoxy
{
  namespace
  {
    /**
     *  The bulge of an arc, from a quarter of its sweep: its tangent. An angle held in degrees
     *  takes it from its direction, correctly rounded, so that a half circle's is 1 exactly.
     *
     *  @param  quarter A quarter of the sweep, under a quarter turn either way.
     *  @return         Its tangent.
     */
    auto Bulged(Angle quarter) -> double
    {
      if (!quarter.HeldInDegrees())
      {
        return Trig::Tan(quarter.InRadians());
      }
      const Vector way = Direction(quarter);
      return way.Y / way.X;
    }

    /**
     *  Which of a chain's curves are turned, so that each starts where the one before ends.
     *
     *  @param  ends      Each curve's start and end, in order.
     *  @param  tolerance How far apart two ends may be and still touch.
     *  @return           Whether each is turned: the first where its start touches the second and
     *                    its end does not; none where a curve does not touch the one before.
     */
    auto Facing(std::span<const std::array<Point, 2>> ends, const Tolerance &tolerance)
        -> std::optional<std::vector<bool>>
    {
      std::vector<bool> turned(ends.size(), false);
      if (ends.size() > 1)
      {
        const auto touches = [&ends, &tolerance](Point point) -> bool
        {
          return Same(point, ends[1][0], tolerance) || Same(point, ends[1][1], tolerance);
        };
        turned[0] = !touches(ends[0][1]) && touches(ends[0][0]);
      }
      Point reached = turned[0] ? ends[0][0] : ends[0][1];
      for (std::size_t index = 1; index < ends.size(); ++index)
      {
        if (Same(ends[index][0], reached, tolerance))
        {
          reached = ends[index][1];
        }
        else if (Same(ends[index][1], reached, tolerance))
        {
          turned[index] = true;
          reached       = ends[index][0];
        }
        else
        {
          return std::nullopt;
        }
      }
      return turned;
    }

    /**
     *  A line, an arc or a polyline as an open polyline.
     *
     *  @param  curve The curve, a line, an arc or a polyline.
     *  @return       A line as one segment; an arc as one bulge, or as two halves swept a whole
     *                turn; a polyline as it is, a closed one opened with its first vertex again at
     *                its end. None for an arc swept more than a turn or not by a number, and for a
     *                polyline of no segment.
     */
    auto Opened(const Curve &curve) -> std::optional<Polyline>
    {
      if (const Line *line = std::get_if<Line>(&curve))
      {
        return Polyline{.Vertices = {{.Position = line->Start, .Bulge = 0.0},
                                     {.Position = line->End, .Bulge = 0.0}},
                        .Closed   = false};
      }
      if (const Arc *arc = std::get_if<Arc>(&curve))
      {
        const double sweep = arc->Sweep.InRadians();
        if (!(std::abs(sweep) <= 2.0 * Pi))
        {
          return std::nullopt;
        }
        const std::size_t count  = std::abs(sweep) < 2.0 * Pi ? 1 : 2;
        const double      bulge  = Bulged(arc->Sweep / (4.0 * static_cast<double>(count)));
        Polyline          opened = {.Vertices = {}, .Closed = false};
        for (std::size_t index = 0; index <= count; ++index)
        {
          const double along = static_cast<double>(index) / static_cast<double>(count);
          opened.Vertices.push_back(
              {.Position = At(*arc, along), .Bulge = index < count ? bulge : 0.0});
        }
        return opened;
      }
      Polyline opened = std::get<Polyline>(curve);
      if (Count(opened) == 0)
      {
        return std::nullopt;
      }
      if (opened.Closed)
      {
        opened.Vertices.push_back({.Position = opened.Vertices.front().Position, .Bulge = 0.0});
        opened.Closed = false;
      }
      return opened;
    }

    /**
     *  Lines, arcs and polylines linked into one polyline.
     *
     *  @param  curves    The curves.
     *  @param  tolerance How far apart two ends may be and still link.
     *  @return           The polyline, or none.
     */
    auto Polygon(std::span<const Curve> curves, const Tolerance &tolerance) -> std::optional<Curve>
    {
      std::vector<Polyline>             runs;
      std::vector<std::array<Point, 2>> ends;
      for (const Curve &curve : curves)
      {
        const std::optional<Polyline> opened = Opened(curve);
        if (!opened)
        {
          return std::nullopt;
        }
        runs.push_back(*opened);
        ends.push_back({opened->Vertices.front().Position, opened->Vertices.back().Position});
      }
      const std::optional<std::vector<bool>> turned = Facing(ends, tolerance);
      if (!turned)
      {
        return std::nullopt;
      }
      Polyline chain = (*turned)[0] ? Reversed(runs[0]) : runs[0];
      for (std::size_t index = 1; index < runs.size(); ++index)
      {
        const Polyline run          = (*turned)[index] ? Reversed(runs[index]) : runs[index];
        chain.Vertices.back().Bulge = run.Vertices.front().Bulge;
        chain.Vertices.insert(chain.Vertices.end(), run.Vertices.begin() + 1, run.Vertices.end());
      }
      if (chain.Vertices.size() > 2
          && Same(chain.Vertices.back().Position, chain.Vertices.front().Position, tolerance))
      {
        chain.Vertices.pop_back();
        chain.Closed = true;
      }
      return chain;
    }

    /**
     *  Curves linked into one spline.
     *
     *  @param  curves    The curves.
     *  @param  tolerance How far apart two ends may be and still link.
     *  @return           The spline, or none.
     */
    auto Splined(std::span<const Curve> curves, const Tolerance &tolerance) -> std::optional<Curve>
    {
      std::vector<std::vector<Kernel::Piece>> runs;
      std::vector<std::array<Point, 2>>       ends;
      for (const Curve &curve : curves)
      {
        runs.push_back(Sliced(curve));
        if (runs.back().empty())
        {
          return std::nullopt;
        }
        ends.push_back({Kernel::Start(runs.back().front()), Kernel::End(runs.back().back())});
      }
      const std::optional<std::vector<bool>> turned = Facing(ends, tolerance);
      if (!turned)
      {
        return std::nullopt;
      }
      std::vector<Kernel::Piece> pieces;
      for (std::size_t index = 0; index < runs.size(); ++index)
      {
        std::vector<Kernel::Piece> &run = runs[index];
        if ((*turned)[index])
        {
          std::ranges::reverse(run);
          std::ranges::transform(run, run.begin(), Flipped);
        }
        pieces.insert(pieces.end(), run.begin(), run.end());
      }
      Spline spline = Joined(pieces);
      if (Same(spline.Poles.back(), spline.Poles.front(), tolerance))
      {
        spline.Poles.back() = spline.Poles.front();
      }
      return spline;
    }
  } // namespace

  auto Join(std::span<const Curve> curves, const Tolerance &tolerance) -> std::optional<Curve>
  {
    if (curves.empty())
    {
      return std::nullopt;
    }
    const bool polygonal = std::ranges::all_of(curves,
                                               [](const Curve &curve) -> bool
                                               {
                                                 return std::holds_alternative<Line>(curve)
                                                     || std::holds_alternative<Arc>(curve)
                                                     || std::holds_alternative<Polyline>(curve);
                                               });
    return polygonal ? Polygon(curves, tolerance) : Splined(curves, tolerance);
  }
} // namespace Geoxy
