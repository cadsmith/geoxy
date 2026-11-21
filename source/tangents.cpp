/***************************************************************************************************
 *  @file      source/tangents.cpp
 *  @brief     The places of a curve whose tangent passes through a point: arcs in closed form,
 *             every other curve by the kernel's touches on its pieces.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "geoxy/tangents.hpp"

#include "basis.hpp"
#include "elliptic.hpp"
#include "kernel/piece.hpp"
#include "real.hpp"
#include "scratch.hpp"
#include "seek.hpp"
#include "trig.hpp"

#include "geoxy/angle.hpp"
#include "geoxy/arc.hpp"
#include "geoxy/bezier/cubic.hpp"
#include "geoxy/bezier/quadratic.hpp"
#include "geoxy/curve.hpp"
#include "geoxy/ellipse.hpp"
#include "geoxy/line.hpp"
#include "geoxy/place.hpp"
#include "geoxy/point.hpp"
#include "geoxy/polyline.hpp"
#include "geoxy/ray.hpp"
#include "geoxy/segment.hpp"
#include "geoxy/spline.hpp"
#include "geoxy/vector.hpp"
#include "geoxy/xline.hpp"

#include <cstddef>
#include <variant>
#include <vector>

namespace Geoxy
{
  namespace
  {
    /**
     *  The parameters of an arc whose tangent passes through a point.
     *
     *  @param  arc   The arc.
     *  @param  point The point.
     *  @return       The parameters, unordered.
     */
    auto Touching(const Arc &arc, Point point) -> List<double>
    {
      List<double> alongs;
      if (Zero(arc.Radius))
      {
        return alongs;
      }
      const Vector seen  = (point - arc.Center) / arc.Radius;
      const double reach = Length(seen);
      if (!(reach >= 1.0))
      {
        return alongs;
      }
      const double  toward  = Heading(seen).InRadians();
      const double  spread  = Trig::Acos(1.0 / reach);
      const Ellipse ellipse = Elliptic(arc);
      for (const double angle : {toward - spread, toward + spread})
      {
        const List<double> passes = Passes(ellipse, angle);
        alongs.insert(alongs.end(), passes.begin(), passes.end());
      }
      return alongs;
    }

    /**
     *  What a polyline's segment has of the touches: its own, a bulged one's in closed form, and
     *  which way it turns about the point at its ends, for the turns at its vertices.
     *
     *  @param  segment The segment.
     *  @param  point   The point.
     *  @return         The touches' parameters, and the segment's turn about the point at its ends.
     */
    auto Grazed(const Segment &segment, Point point) -> Kernel::Solution
    {
      Kernel::Solution solution = Leaned(segment, point, [](Vector away, Vector tangent) -> double
                                         { return Cross(away, tangent); });
      if (const Arc *arc = std::get_if<Arc>(&segment))
      {
        solution.Roots = Touching(*arc, point);
      }
      return solution;
    }

    /**
     *  The touches of a piece.
     *
     *  @param  point The point.
     *  @return       What solves for the touches on a piece.
     */
    auto Touched(Point point)
    {
      return [point](const Kernel::Piece &piece) -> Kernel::Solution
      {
        return Kernel::Touches(piece, point, Kernel::Budget);
      };
    }
  } // namespace

  auto Tangents(const Line & /*line*/, Point /*point*/) -> std::vector<Place>
  {
    return {};
  }

  auto Tangents(const Ray & /*ray*/, Point /*point*/) -> std::vector<Place>
  {
    return {};
  }

  auto Tangents(const Xline & /*xline*/, Point /*point*/) -> std::vector<Place>
  {
    return {};
  }

  auto Tangents(const Arc &arc, Point point) -> std::vector<Place>
  {
    return Vectored(Placed(arc, Touching(arc, point)));
  }

  auto Tangents(const Ellipse &ellipse, Point point) -> std::vector<Place>
  {
    return Vectored(Placed(ellipse, Rounded(ellipse, Touched(point), Uncapped{})));
  }

  auto Tangents(const Polyline &polyline, Point point) -> std::vector<Place>
  {
    return Vectored(
        Placed(polyline, Segmented(polyline, [point](const Segment &segment) -> Kernel::Solution
                                   { return Grazed(segment, point); })));
  }

  auto Tangents(const Spline &spline, Point point) -> std::vector<Place>
  {
    if (Rational(spline))
    {
      const auto touches = [point](const Kernel::Tread &tread) -> Kernel::Solution
      {
        return Kernel::Touches(tread, point, Kernel::Budget);
      };
      return Vectored(Sought(spline, touches, [&spline, point](std::size_t span) -> Kernel::Solution
                             { return Kernel::Touches(spline, span, point); }));
    }
    return Vectored(Placed(spline, Kernel::Solved(spline, Touched(point))));
  }

  auto Tangents(const Bezier::Quadratic &quadratic, Point point) -> std::vector<Place>
  {
    return Vectored(Placed(quadratic, Touched(point)(Kernel::Polynomial(quadratic)).Roots));
  }

  auto Tangents(const Bezier::Cubic &cubic, Point point) -> std::vector<Place>
  {
    return Vectored(Placed(cubic, Touched(point)(Kernel::Polynomial(cubic)).Roots));
  }

  auto Tangents(const Curve &curve, Point point) -> std::vector<Place>
  {
    return std::visit([point](const auto &shape) -> std::vector<Place>
                      { return Tangents(shape, point); }, curve);
  }
} // namespace Geoxy
