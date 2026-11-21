/***************************************************************************************************
 *  @file      source/offset.cpp
 *  @brief     The curve at a distance from another.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "geoxy/offset.hpp"

#include "basis.hpp"
#include "belt/belt.hpp"
#include "finite.hpp"
#include "kernel/piece.hpp"
#include "parallel/smooth.hpp"
#include "real.hpp"
#include "sane.hpp"
#include "scratch.hpp"
#include "trig.hpp"

#include "geoxy/angle.hpp"
#include "geoxy/arc.hpp"
#include "geoxy/area.hpp"
#include "geoxy/bezier/cubic.hpp"
#include "geoxy/bezier/quadratic.hpp"
#include "geoxy/box.hpp"
#include "geoxy/curve.hpp"
#include "geoxy/domain.hpp"
#include "geoxy/ellipse.hpp"
#include "geoxy/error.hpp"
#include "geoxy/line.hpp"
#include "geoxy/matrix.hpp"
#include "geoxy/point.hpp"
#include "geoxy/polyline.hpp"
#include "geoxy/ray.hpp"
#include "geoxy/segment.hpp"
#include "geoxy/spline.hpp"
#include "geoxy/tolerance.hpp"
#include "geoxy/vector.hpp"
#include "geoxy/xline.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <expected>
#include <iterator>
#include <limits>
#include <span>
#include <type_traits>
#include <utility>
#include <variant>
#include <vector>

namespace Geoxy
{
  namespace
  {
    /**
     *  The displacement square to a way, to its right, as long as a distance.
     *
     *  @param  way      The way, not of no length.
     *  @param  distance The distance, to the right when positive.
     *  @return          The displacement.
     */
    auto Right(Vector way, double distance) -> Vector
    {
      const Vector side = {.X = way.Y, .Y = -way.X};
      /**
       *  A level or upright way's unit normal is exact, and so the displacement: a point moved by
       *  it is the double nearest the true one, rounded once.
       */
      return Zero(way.X) || Zero(way.Y) ? (side / Length(way)) * distance
                                        : side * (distance / Length(way));
    }

    /**
     *  A straight curve's start and way moved sideways, or why not.
     *
     *  @param  start    Its start, or a point it runs through.
     *  @param  way      Its way.
     *  @param  distance The distance.
     *  @return          The start moved.
     */
    auto Sidestepped(Point start, Vector way, double distance) -> std::expected<Point, Error>
    {
      if (!Finite(start) || !Finite(way) || !std::isfinite(distance))
      {
        return std::unexpected(Error::NotFinite);
      }
      if (Zero(Largest(way)))
      {
        return std::unexpected(Error::Degenerate);
      }
      const Point moved = start + Right(way, distance);
      if (!Finite(moved))
      {
        return std::unexpected(Error::NotFinite);
      }
      return moved;
    }

    /**
     *  A tolerance fit to offset within: sane, and no finer than 16 steps of a double at the size
     *  of the numbers, finer than which where a point lies is not known.
     *
     *  @param  tolerance The tolerance asked for.
     *  @param  size      The largest magnitude among the numbers.
     *  @return           The tolerance; NotFinite where the size is past 2^500, where the squares
     *                    of the numbers an offset works with would be past every double.
     */
    auto Suited(const Tolerance &tolerance, double size) -> std::expected<Tolerance, Error>
    {
      if (!(size < 0x1p500))
      {
        return std::unexpected(Error::NotFinite);
      }
      Tolerance fit = Sane(tolerance);
      fit.Length    = std::max(fit.Length, 16.0 * Grain(size));
      return fit;
    }

    /**
     *  The larger of a size and a number's magnitude, past every double once a number is not
     *  finite.
     *
     *  @param  size  The size so far.
     *  @param  value The number.
     *  @return       The size.
     */
    auto Larger(double size, double value) -> double
    {
      return std::isfinite(value) ? std::max(size, std::abs(value))
                                  : std::numeric_limits<double>::infinity();
    }

    /**
     *  Whether a polyline has a vertex within the tolerance of the one before it, or a closed one
     *  its last within the tolerance of its first.
     *
     *  @param  polyline  The polyline.
     *  @param  tolerance The tolerance.
     *  @return           Whether it has.
     */
    auto Doubled(const Polyline &polyline, const Tolerance &tolerance) -> bool
    {
      const std::vector<Vertex> &vertices = polyline.Vertices;
      for (std::size_t index = 1; index < vertices.size(); ++index)
      {
        if (Same(vertices[index - 1].Position, vertices[index].Position, tolerance))
        {
          return true;
        }
      }
      return polyline.Closed
          && vertices.size() > 1
          && Same(vertices.back().Position, vertices.front().Position, tolerance);
    }

    /**
     *  A polyline with each vertex within the tolerance of the one before taken as one with it, and
     *  a closed one's last within the tolerance of its first.
     *
     *  @param  polyline  The polyline.
     *  @param  tolerance The tolerance.
     *  @return           The polyline, no two vertices in a row one point.
     */
    auto Cleaned(const Polyline &polyline, const Tolerance &tolerance) -> Polyline
    {
      Polyline clean = {.Closed = polyline.Closed};
      clean.Vertices.reserve(polyline.Vertices.size());
      for (const Vertex &vertex : polyline.Vertices)
      {
        if (!clean.Vertices.empty()
            && Same(clean.Vertices.back().Position, vertex.Position, tolerance))
        {
          clean.Vertices.back().Bulge = vertex.Bulge;
          continue;
        }
        clean.Vertices.push_back(vertex);
      }
      if (clean.Closed
          && clean.Vertices.size() > 1
          && Same(clean.Vertices.back().Position, clean.Vertices.front().Position, tolerance))
      {
        clean.Vertices.pop_back();
      }
      return clean;
    }

    /**
     *  The most arcs a smooth curve's offset is laid with.
     */
    constexpr std::size_t Budget = 65'536;

    /**
     *  A smooth curve to offset: its runs, and whether its end is its start.
     */
    struct Smooth
    {
      std::span<const Parallel::Run> Runs   = {};    /**< The runs, in order. */
      bool                           Closed = false; /**< Whether the last runs into the
                                                          first.                         */
    };

    /**
     *  A smooth curve laid with arcs as a polyline: their starts and bulges, and for an open curve
     *  the last's end.
     *
     *  @param  layout The arcs.
     *  @param  closed Whether the curve is closed.
     *  @return        The polyline.
     */
    auto Pathed(const Parallel::Layout &layout, bool closed) -> Polyline
    {
      Polyline path = {.Vertices = {}, .Closed = closed};
      path.Vertices.reserve(layout.Pieces.size() + 1);
      for (const Parallel::Moved &piece : layout.Pieces)
      {
        path.Vertices.push_back({.Position = piece.Start, .Bulge = Trig::Tan(piece.Sweep / 4.0)});
      }
      path.Vertices.push_back({.Position = layout.Pieces.back().End});
      if (closed)
      {
        path.Vertices.pop_back();
      }
      return path;
    }

    /**
     *  A smooth curve's offset, laid with arcs and cleaned as a polyline's.
     *
     *  @param  smooth     The curve's runs and whether it closes.
     *  @param  offsetting The distance, the gap kind at its kinks and the accuracy.
     *  @param  tolerance  The tolerance, fit to the curve's size.
     *  @return            The offset's polylines; Spent where the budget ran out.
     */
    auto Smoothed(const Smooth &smooth, const Offsetting &offsetting, const Tolerance &tolerance)
        -> std::expected<std::vector<Polyline>, Error>
    {
      const double           distance = offsetting.Distance;
      const double           accuracy = std::max(offsetting.Approximation.Error, tolerance.Length);
      const Parallel::Layout curve =
          Parallel::Laid(smooth.Runs, {.Distance = 0.0, .Accuracy = accuracy, .Budget = Budget});
      if (curve.Spent)
      {
        return std::unexpected(Error::Spent);
      }
      /**
       *  The curve laid with arcs is a polyline within the accuracy of it, and so is that
       *  polyline's offset of the curve's.
       */
      return Offset(Pathed(curve, smooth.Closed),
                    Offsetting{.Distance = distance, .Gaps = offsetting.Gaps}, tolerance);
    }

    /**
     *  A Bézier curve's offset, as a smooth curve of one run.
     *
     *  @param  piece      The curve, as a piece.
     *  @param  offsetting The distance and the accuracy.
     *  @param  tolerance  The tolerance asked for.
     *  @return            The offset's polylines.
     */
    auto Smoothing(const Kernel::Piece &piece, const Offsetting &offsetting,
                   const Tolerance &tolerance) -> std::expected<std::vector<Polyline>, Error>
    {
      double size = Larger(0.0, offsetting.Distance);
      Box    box;
      for (const Homogeneous &pole : Kernel::Poles(piece))
      {
        const Point point = Projected(pole);
        size              = Larger(Larger(size, point.X), point.Y);
        box               = Grown(box, point);
      }
      const std::expected<Tolerance, Error> fit = Suited(tolerance, size);
      if (!fit)
      {
        return std::unexpected(fit.error());
      }
      if (!(Largest(Size(box)) > fit->Length))
      {
        return std::unexpected(Error::Degenerate);
      }
      const std::array<Parallel::Run, 1> runs = {piece};
      return Smoothed(
          {.Runs = runs, .Closed = Same(Kernel::Start(piece), Kernel::End(piece), *fit)},
          offsetting, *fit);
    }

    /**
     *  A list of polylines as a list of curves.
     *
     *  @param  polylines The polylines, or why there are none.
     *  @return           The curves.
     */
    auto Several(std::expected<std::vector<Polyline>, Error> &&polylines)
        -> std::expected<std::vector<Curve>, Error>
    {
      return std::move(polylines).transform(
          [](std::vector<Polyline> found) -> std::vector<Curve>
          {
            return {std::make_move_iterator(found.begin()), std::make_move_iterator(found.end())};
          });
    }

    /**
     *  One curve as a list of curves.
     *
     *  @param  curve The curve, or why there is none.
     *  @return       The list.
     */
    template <typename CurveT>
    auto Singleton(std::expected<CurveT, Error> &&curve) -> std::expected<std::vector<Curve>, Error>
    {
      return std::move(curve).transform([](CurveT &&found) -> std::vector<Curve>
                                        { return {Curve{std::move(found)}}; });
    }
  } // namespace

  auto Offset(const Line &line, double distance) -> std::expected<Line, Error>
  {
    const Vector way = line.End - line.Start;
    return Sidestepped(line.Start, way, distance)
        .transform([way](Point start) -> Line { return {.Start = start, .End = start + way}; });
  }

  auto Offset(const Ray &ray, double distance) -> std::expected<Ray, Error>
  {
    return Sidestepped(ray.Start, ray.Way, distance)
        .transform([&ray](Point start) -> Ray { return {.Start = start, .Way = ray.Way}; });
  }

  auto Offset(const Xline &xline, double distance) -> std::expected<Xline, Error>
  {
    return Sidestepped(xline.Through, xline.Way, distance)
        .transform([&xline](Point through) -> Xline
                   { return {.Through = through, .Way = xline.Way}; });
  }

  auto Offset(const Arc &arc, double distance) -> std::expected<Arc, Error>
  {
    const double sweep = arc.Sweep.InRadians();
    if (!Finite(arc.Center)
        || !std::isfinite(arc.Radius)
        || !std::isfinite(sweep)
        || !std::isfinite(arc.Start.InRadians())
        || !std::isfinite(distance))
    {
      return std::unexpected(Error::NotFinite);
    }
    if (!(arc.Radius > 0.0) || Zero(sweep))
    {
      return std::unexpected(Error::Degenerate);
    }
    const double radius = arc.Radius + (sweep < 0.0 ? -distance : distance);
    if (!std::isfinite(radius))
    {
      return std::unexpected(Error::NotFinite);
    }
    if (!(radius > 0.0))
    {
      return std::unexpected(Error::TooLarge);
    }
    Arc moved    = arc;
    moved.Radius = radius;
    return moved;
  }

  auto Offset(const Polyline &polyline, const Offsetting &offsetting, const Tolerance &tolerance)
      -> std::expected<std::vector<Polyline>, Error>
  {
    const double distance = offsetting.Distance;
    double       size     = Larger(0.0, distance);
    for (const Vertex &vertex : polyline.Vertices)
    {
      size = Larger(Larger(Larger(size, vertex.Position.X), vertex.Position.Y), 0.0 * vertex.Bulge);
      /**
       *  0 times a bulge is 0, but for one that is not finite: the size is then past every double.
       */
    }
    const std::expected<Tolerance, Error> fit = Suited(tolerance, size);
    if (!fit)
    {
      return std::unexpected(fit.error());
    }
    const bool      doubled = Doubled(polyline, *fit);
    const Polyline  clean   = doubled ? Cleaned(polyline, *fit) : Polyline{};
    const Polyline &source  = doubled ? clean : polyline;
    if (Count(source) == 0)
    {
      return std::unexpected(Error::Degenerate);
    }
    if (Zero(distance))
    {
      return std::vector{source};
    }
    if (!source.Closed)
    {
      return Belt::Flanked(source, distance, offsetting.Gaps);
    }
    /**
     *  To the right of a closed polyline that runs counterclockwise is outside what it fills, and
     *  inside what one that runs clockwise fills. One that encloses no area, as near as doubles
     *  tell, a stretch run there and back or a figure of eight, has no such side: either way its
     *  offset is the edge of all that lies within the distance of it.
     */
    Point least = source.Vertices.front().Position;
    Point most  = least;
    for (const Vertex &vertex : source.Vertices)
    {
      least = {.X = std::min(least.X, vertex.Position.X),
               .Y = std::min(least.Y, vertex.Position.Y)};
      most  = {.X = std::max(most.X, vertex.Position.X), .Y = std::max(most.Y, vertex.Position.Y)};
    }
    /**
     *  No area: less than the area's own rounding, which goes by how far from 0 the polyline lies
     *  and how large it is.
     */
    const double area = Area(source);
    const bool   none = std::abs(area) <= 1e-12 * size * ((most.X - least.X) + (most.Y - least.Y));
    const bool   outward = (distance > 0.0) == (area > 0.0);
    const Side   sides   = outward ? Side::Outward : Side::Inward;
    return Belt::Edged(source, {.Distance = std::abs(distance),
                                .Sides    = none ? Side::Both : sides,
                                .Fill     = Fills::NonZero,
                                .Gaps     = offsetting.Gaps,
                                .Ends     = Cap::Round});
  }

  auto Buffer(const Polyline &polyline, const Buffering &buffering, const Tolerance &tolerance)
      -> std::expected<std::vector<Polyline>, Error>
  {
    const double distance = buffering.Distance;
    double       size     = Larger(0.0, distance);
    for (const Vertex &vertex : polyline.Vertices)
    {
      size = Larger(Larger(Larger(size, vertex.Position.X), vertex.Position.Y), 0.0 * vertex.Bulge);
    }
    const std::expected<Tolerance, Error> fit = Suited(tolerance, size);
    if (!fit)
    {
      return std::unexpected(fit.error());
    }
    const bool      doubled = Doubled(polyline, *fit);
    const Polyline  clean   = doubled ? Cleaned(polyline, *fit) : Polyline{};
    const Polyline &source  = doubled ? clean : polyline;
    if ((!source.Closed && buffering.Sides != Side::Both)
        || Count(source) == 0
        || !(distance > 0.0))
    {
      return std::unexpected(Error::Degenerate);
    }
    return Belt::Edged(source, buffering);
  }

  auto Offset(const Ellipse &ellipse, const Offsetting &offsetting, const Tolerance &tolerance)
      -> std::expected<std::vector<Polyline>, Error>
  {
    const Matrix &map   = ellipse.Map;
    const double  start = ellipse.Start.InRadians();
    const double  sweep = ellipse.Sweep.InRadians();
    const double  size =
        Larger(Larger(Larger(Larger(0.0, offsetting.Distance), 0.0 * (start + sweep)),
                      std::abs(map.Origin.X) + std::abs(map.X.X) + std::abs(map.Y.X)),
               std::abs(map.Origin.Y) + std::abs(map.X.Y) + std::abs(map.Y.Y));
    const std::expected<Tolerance, Error> fit = Suited(tolerance, size);
    if (!fit)
    {
      return std::unexpected(fit.error());
    }
    if (Zero(Cross(map.X, map.Y)) || Zero(sweep))
    {
      return std::unexpected(Error::Degenerate);
    }
    const bool                         whole = std::abs(sweep) >= 2.0 * Pi;
    const std::array<Parallel::Run, 1> runs  = {
        Parallel::Oval{.Map = map, .Start = start, .Sweep = whole ? 2.0 * Pi : sweep}};
    return Smoothed({.Runs = runs, .Closed = whole}, offsetting, *fit);
  }

  auto Offset(const Spline &spline, const Offsetting &offsetting, const Tolerance &tolerance)
      -> std::expected<std::vector<Polyline>, Error>
  {
    double size = Larger(0.0, offsetting.Distance);
    Box    box;
    for (const Point &pole : spline.Poles)
    {
      size = Larger(Larger(size, pole.X), pole.Y);
      box  = Grown(box, pole);
    }
    const std::expected<Tolerance, Error> fit = Suited(tolerance, size);
    if (!fit)
    {
      return std::unexpected(fit.error());
    }
    if (!Valid(spline) || !(Largest(Size(box)) > fit->Length))
    {
      return std::unexpected(Error::Degenerate);
    }
    List<Parallel::Run> runs;
    Kernel::Sections(spline,
                     [&runs](const Kernel::Piece &piece, Range /*span*/) -> void
                     {
                       /**
                        *  A span whose poles are one point is the curve standing still: it has no
                        *  way, and nothing to offset.
                        */
                       if (Largest(Size(Kernel::Hull(piece))) > 0.0)
                       {
                         runs.emplace_back(piece);
                       }
                     });
    const Point first = Kernel::Start(std::get<Kernel::Piece>(runs.front()));
    const Point last  = Kernel::End(std::get<Kernel::Piece>(runs.back()));
    return Smoothed({.Runs = runs, .Closed = Same(first, last, *fit)}, offsetting, *fit);
  }

  auto Offset(const Curve &curve, const Offsetting &offsetting, const Tolerance &tolerance)
      -> std::expected<std::vector<Curve>, Error>
  {
    const double distance = offsetting.Distance;
    return std::visit(
        [distance, &offsetting,
         &tolerance](const auto &shape) -> std::expected<std::vector<Curve>, Error>
        {
          using ShapeT = std::decay_t<decltype(shape)>;
          if constexpr (std::is_same_v<ShapeT, Arc>)
          {
            std::expected<Arc, Error> arc = Offset(shape, distance);
            if (!arc && arc.error() == Error::TooLarge)
            {
              return std::vector<Curve>{};
            }
            return Singleton(std::move(arc));
          }
          else if constexpr (std::is_same_v<ShapeT, Line>
                             || std::is_same_v<ShapeT, Ray>
                             || std::is_same_v<ShapeT, Xline>)
          {
            return Singleton(Offset(shape, distance));
          }
          else if constexpr (std::is_same_v<ShapeT, Bezier::Quadratic>
                             || std::is_same_v<ShapeT, Bezier::Cubic>)
          {
            return Several(Smoothing(Kernel::Polynomial(shape), offsetting, tolerance));
          }
          else
          {
            return Several(Offset(shape, offsetting, tolerance));
          }
        },
        curve);
  }
} // namespace Geoxy
