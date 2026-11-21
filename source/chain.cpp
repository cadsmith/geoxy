/***************************************************************************************************
 *  @file      source/chain.cpp
 *  @brief     Curves as runs of rational Bézier pieces, and runs of pieces as one spline.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "chain.hpp"

#include "basis.hpp"
#include "kernel/piece.hpp"
#include "real.hpp"
#include "trig.hpp"
#include "walk.hpp"

#include "geoxy/angle.hpp"
#include "geoxy/arc.hpp"
#include "geoxy/at.hpp"
#include "geoxy/bezier/cubic.hpp"
#include "geoxy/bezier/quadratic.hpp"
#include "geoxy/curve.hpp"
#include "geoxy/domain.hpp"
#include "geoxy/ellipse.hpp"
#include "geoxy/line.hpp"
#include "geoxy/matrix.hpp"
#include "geoxy/point.hpp"
#include "geoxy/polyline.hpp"
#include "geoxy/ray.hpp"
#include "geoxy/segment.hpp"
#include "geoxy/spline.hpp"
#include "geoxy/vector.hpp"
#include "geoxy/xline.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <span>
#include <variant>
#include <vector>

namespace Geoxy
{
  namespace
  {
    /**
     *  The most one rational quadratic piece here sweeps: a quarter turn, in radians.
     */
    constexpr double QuarterTurn = Pi / 2.0;

    /**
     *  A rational quadratic piece of a circle's arc, or of one mapped.
     *
     *  @param  ends  Where it starts and where it ends.
     *  @param  bend  Its middle control point, where its tangents at its ends meet.
     *  @param  sweep Its sweep, in radians, at most a quarter turn either way.
     *  @return       The piece, its middle weight the cosine of half the sweep and its others 1.
     */
    auto Circlet(std::array<Point, 2> ends, Point bend, double sweep) -> Kernel::Piece
    {
      const double  weight = Trig::Cos(sweep / 2.0);
      Kernel::Piece piece  = Kernel::Polynomial(std::array{ends[0], bend, ends[1]});
      piece.Poles[1]       = {.X = bend.X * weight, .Y = bend.Y * weight, .W = weight};
      return piece;
    }

    /**
     *  Where the tangents at the ends of a circle's arc meet.
     *
     *  @param  start Where it starts.
     *  @param  end   Where it ends.
     *  @param  sweep Its sweep, in radians, under a half turn either way.
     *  @return       The chord's middle moved off it by half the chord times the tangent of half
     *                the sweep: to the right of the way it runs for a counterclockwise sweep, the
     *                side the arc bulges to.
     */
    auto Apex(Point start, Point end, double sweep) -> Point
    {
      return Midpoint(start, end) - (Perpendicular(end - start) * (Trig::Tan(sweep / 2.0) / 2.0));
    }

    /**
     *  How many pieces of at most a quarter turn a sweep takes.
     *
     *  @param  sweep The sweep, in radians, at most a turn either way.
     *  @return       The count, 1 at least.
     */
    auto QuarterCount(double sweep) -> std::size_t
    {
      return static_cast<std::size_t>(std::max(1.0, std::ceil(std::abs(sweep) / QuarterTurn)));
    }

    /**
     *  Whether pieces of a quarter turn hold a sweep: at most a turn, and a number.
     *
     *  @param  sweep The sweep.
     *  @return       Whether they do.
     */
    auto Sweepable(Angle sweep) -> bool
    {
      return std::abs(sweep.InRadians()) <= 2.0 * Pi;
    }

    /**
     *  The fraction of a count of equal parts.
     *
     *  @param  index Which part's start.
     *  @param  count How many parts.
     *  @return       index / count, exactly 1 at the count.
     */
    auto Fraction(std::size_t index, std::size_t count) -> double
    {
      return static_cast<double>(index) / static_cast<double>(count);
    }

    auto Sliced(const Line &line) -> std::vector<Kernel::Piece>
    {
      return {Kernel::Straight(line.Start, line.End)};
    }

    auto Sliced(const Ray & /*ray*/) -> std::vector<Kernel::Piece>
    {
      return {};
    }

    auto Sliced(const Xline & /*xline*/) -> std::vector<Kernel::Piece>
    {
      return {};
    }

    auto Sliced(const Arc &arc) -> std::vector<Kernel::Piece>
    {
      if (!Sweepable(arc.Sweep))
      {
        return {};
      }
      const double               sweep = arc.Sweep.InRadians();
      const std::size_t          count = QuarterCount(sweep);
      const double               step  = sweep / static_cast<double>(count);
      std::vector<Kernel::Piece> pieces;
      for (std::size_t index = 0; index < count; ++index)
      {
        const Point start = At(arc, Fraction(index, count));
        const Point end   = At(arc, Fraction(index + 1, count));
        pieces.push_back(Circlet({start, end}, Apex(start, end, step), step));
      }
      return pieces;
    }

    auto Sliced(const Ellipse &ellipse) -> std::vector<Kernel::Piece>
    {
      if (!Sweepable(ellipse.Sweep))
      {
        return {};
      }
      const double               sweep = ellipse.Sweep.InRadians();
      const std::size_t          count = QuarterCount(sweep);
      const double               step  = sweep / static_cast<double>(count);
      std::vector<Kernel::Piece> pieces;
      for (std::size_t index = 0; index < count; ++index)
      {
        const double start    = Fraction(index, count);
        const double end      = Fraction(index + 1, count);
        const Point  leaving  = Point{} + Direction(ellipse.Start + (ellipse.Sweep * start));
        const Point  arriving = Point{} + Direction(ellipse.Start + (ellipse.Sweep * end));
        pieces.push_back(Circlet({At(ellipse, start), At(ellipse, end)},
                                 ellipse.Map * Apex(leaving, arriving, step), step));
      }
      return pieces;
    }

    auto Sliced(const Spline &spline) -> std::vector<Kernel::Piece>
    {
      std::vector<Kernel::Piece> pieces;
      Kernel::Sections(spline, [&pieces](const Kernel::Piece &piece, Range /*span*/) -> void
                       { pieces.push_back(piece); });
      return pieces;
    }

    auto Sliced(const Bezier::Quadratic &quadratic) -> std::vector<Kernel::Piece>
    {
      return {Kernel::Polynomial(quadratic)};
    }

    auto Sliced(const Bezier::Cubic &cubic) -> std::vector<Kernel::Piece>
    {
      return {Kernel::Polynomial(cubic)};
    }

    /**
     *  A piece lifted to a higher degree, the same curve: each round's control points the mix of
     *  two of the last's, in homogeneous coordinates.
     *
     *  @param  piece  The piece.
     *  @param  degree The degree, at least the piece's and at most MaxDegree.
     *  @return        The piece of that degree.
     */
    auto Promoted(Kernel::Piece piece, std::size_t degree) -> Kernel::Piece
    {
      const std::span<Homogeneous> poles{piece.Poles};
      while (piece.Degree < degree)
      {
        const std::size_t from = piece.Degree;
        const auto        next = static_cast<double>(from + 1);
        poles[from + 1]        = poles[from];
        for (std::size_t index = from; index >= 1; --index)
        {
          poles[index] = Mix(poles[index], poles[index - 1], static_cast<double>(index) / next);
        }
        piece.Degree = from + 1;
      }
      return piece;
    }

    /**
     *  A homogeneous point times a factor: the same point, weighed more or less.
     *
     *  @param  point  The point.
     *  @param  factor The factor.
     *  @return        Each coordinate and the weight times it.
     */
    auto Reweighed(const Homogeneous &point, double factor) -> Homogeneous
    {
      return {.X = point.X * factor, .Y = point.Y * factor, .W = point.W * factor};
    }
  } // namespace

  auto Sliced(const Polyline &polyline) -> std::vector<Kernel::Piece>
  {
    std::vector<Kernel::Piece> pieces;
    const std::size_t          segments = Count(polyline);
    for (std::size_t index = 0; index < segments; ++index)
    {
      const auto first = static_cast<double>(index);
      const Step step  = Locate(polyline, first);
      if (Straightness(step))
      {
        pieces.push_back(Kernel::Straight(step.Start, step.End));
        continue;
      }
      const double      sweep = 4.0 * Trig::Atan(step.Bulge);
      const std::size_t count = QuarterCount(sweep);
      const double      part  = sweep / static_cast<double>(count);
      for (std::size_t piece = 0; piece < count; ++piece)
      {
        const Point start = At(polyline, first + Fraction(piece, count));
        const Point end   = At(polyline, first + Fraction(piece + 1, count));
        pieces.push_back(Circlet({start, end}, Apex(start, end, part), part));
      }
    }
    return pieces;
  }

  auto Sliced(const Curve &curve) -> std::vector<Kernel::Piece>
  {
    return std::visit([](const auto &shape) -> std::vector<Kernel::Piece> { return Sliced(shape); },
                      curve);
  }

  auto Flipped(Kernel::Piece piece) -> Kernel::Piece
  {
    std::ranges::reverse(std::span{piece.Poles}.first(piece.Degree + 1));
    return piece;
  }

  auto Joined(std::span<const Kernel::Piece> pieces) -> Spline
  {
    const std::size_t degree = std::ranges::max_element(pieces, {}, &Kernel::Piece::Degree)->Degree;
    const std::size_t count  = pieces.size();
    Spline spline = {.Degree = static_cast<int>(degree), .Knots = {}, .Poles = {}, .Weights = {}};
    spline.Knots.assign(degree + 1, 0.0);
    Homogeneous last = {};
    for (std::size_t index = 0; index < count; ++index)
    {
      const Kernel::Piece                piece = Promoted(pieces[index], degree);
      const std::span<const Homogeneous> poles = Kernel::Poles(piece);
      const double                       scale = index == 0 ? 1.0 : last.W / poles.front().W;
      for (std::size_t at = index == 0 ? 0 : 1; at <= degree; ++at)
      {
        const Homogeneous pole = Reweighed(poles[at], scale);
        spline.Poles.push_back(Projected(pole));
        spline.Weights.push_back(pole.W);
      }
      last = Reweighed(poles.back(), scale);
      spline.Knots.insert(spline.Knots.end(), index + 1 < count ? degree : degree + 1,
                          static_cast<double>(index + 1));
    }
    if (std::ranges::all_of(spline.Weights,
                            [](double weight) -> bool { return Exactly(weight, 1.0); }))
    {
      spline.Weights.clear();
    }
    return spline;
  }
} // namespace Geoxy
