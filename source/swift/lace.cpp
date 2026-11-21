/***************************************************************************************************
 *  @file      source/swift/lace.cpp
 *  @brief     A spline, an ellipse or a Bézier curve met with a polyline: each straight segment
 *             near it by the kernel of a line and that curve, each bulge by the ring kernel on its
 *             circle, every segment far from it passed, and the crossings told apart and placed as
 *             the thorough way places them.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "crossing.hpp"
#include "scratch.hpp"
#include "swift/lane.hpp"
#include "swift/pieced.hpp"
#include "swift/ring.hpp"
#include "swift/section.hpp"

#include "geoxy/at.hpp"
#include "geoxy/bezier/cubic.hpp"
#include "geoxy/bezier/quadratic.hpp"
#include "geoxy/box.hpp"
#include "geoxy/curve.hpp"
#include "geoxy/domain.hpp"
#include "geoxy/ellipse.hpp"
#include "geoxy/intersect.hpp"
#include "geoxy/line.hpp"
#include "geoxy/matrix.hpp"
#include "geoxy/place.hpp"
#include "geoxy/point.hpp"
#include "geoxy/polyline.hpp"
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

namespace Geoxy::Swift
{
  namespace
  {
    /**
     *  How far apart, in reaches, a segment's box and the curve's box may lie and the two still be
     *  met: as far as the kernel of a line and a spline keeps a span, so that a segment passed lies
     *  farther from the curve than every length its answers rest on.
     */
    constexpr double Afar = 128.0;

    /**
     *  A crossing found on a segment: its point, the double nearest the true one, the curve's place
     *  there where the kernel proved it, the double nearest the true one, and which segment.
     */
    struct Stitch
    {
      Point       At     = {};  /**< The point.                      */
      double      Along  = 0.0; /**< The curve's place, where owned. */
      std::size_t Index  = 0;   /**< The polyline's segment's index. */
      double      Placed = 0.0; /**< The polyline's place, the double nearest the true one,
                                     on a straight segment where owned.                 */
      bool Bent = false;        /**< Whether the segment is a bulge, on which the point is
                                     placed by its angle from the chord.               */
      bool Own = true;          /**< Whether the places are the curves' own, as a spline's
                                     and a Bézier curve's are; an ellipse's crossing is
                                     placed by its point, by the angle on the ellipse and
                                     by the foot on a straight segment.                 */
    };

    /**
     *  A curve and a polyline being met: the two, the tolerance, the curve's box, a spline's
     *  sections where it is cut so or a Bézier curve's piece, and the crossings found.
     */
    struct Lacing
    {
      const Curve             *Shape     = nullptr; /**< The curve.                 */
      const Spline            *Strand    = nullptr; /**< It as a spline, where one. */
      const Polyline          *Chain     = nullptr; /**< The polyline.              */
      Geoxy::Tolerance         Tolerance = {};      /**< The tolerance.             */
      Box                      Hull      = {};      /**< A box the curve keeps in.  */
      std::span<const Section> Sections  = {};      /**< A spline's sections, where it is
                                                         cut so; none where the kernel of a
                                                         line and a spline of any degree
                                                         meets it.                      */
      std::optional<Section> Piece = {};            /**< A Bézier curve's piece. */
      List<Stitch>           Found = {};            /**< The crossings found.    */
    };

    /**
     *  The gap between two boxes, along x or along y, whichever is wider: nothing where they meet.
     *
     *  @param  one   One box.
     *  @param  other The other.
     *  @return       The gap.
     */
    auto Cloven(const Box &one, const Box &other) -> double
    {
      return std::max({one.Min.X - other.Max.X, other.Min.X - one.Max.X, one.Min.Y - other.Max.Y,
                       other.Min.Y - one.Max.Y, 0.0});
    }

    /**
     *  A straight segment met with the curve: as a line whose place runs its index behind the
     *  polyline's, by the kernel of a line and a spline, a Bézier curve or an ellipse.
     *
     *  @param  lacing The spline and the polyline, and where the crossings go.
     *  @param  index  The segment.
     *  @return        Whether it was sure.
     */
    auto Straddled(Lacing &lacing, std::size_t index) -> bool
    {
      const Polyline     &chain = *lacing.Chain;
      const Point         start = chain.Vertices[index].Position;
      const Point         end   = chain.Vertices[(index + 1) % chain.Vertices.size()].Position;
      const Curve         line  = Line{.Start = start, .End = end};
      std::optional<Lane> lane  = Laned(line);
      if (!lane)
      {
        return false;
      }
      lane->Lead = static_cast<double>(index);
      std::optional<Spots> spots;
      if (lacing.Strand != nullptr)
      {
        spots = lacing.Sections.empty() ? Strung(*lane, *lacing.Strand, lacing.Tolerance)
                                        : Strung(*lane, lacing.Sections, lacing.Tolerance);
      }
      else if (lacing.Piece)
      {
        spots = Bernstein(*lane, *lacing.Piece, lacing.Tolerance);
      }
      else
      {
        spots = Oval(*lane, std::get<Ellipse>(*lacing.Shape), lacing.Tolerance);
      }
      if (!spots)
      {
        return false;
      }
      for (const Spot &spot : *spots)
      {
        lacing.Found.push_back({.At     = spot.At,
                                .Along  = spot.Along[1],
                                .Index  = index,
                                .Placed = spot.Along[0],
                                .Bent   = false,
                                .Own    = spot.Own});
      }
      return true;
    }

    /**
     *  A bulged segment met with the curve by the ring kernel, on its circle.
     *
     *  @param  lacing The curve and the polyline, and where the crossings go.
     *  @param  hump   The segment.
     *  @return        Whether it was sure.
     */
    auto Overarched(Lacing &lacing, const Hump &hump) -> bool
    {
      const std::optional<Spots> spots = Humped(hump, *lacing.Shape, lacing.Tolerance);
      if (!spots)
      {
        return false;
      }
      for (const Spot &spot : *spots)
      {
        lacing.Found.push_back({.At     = spot.At,
                                .Along  = spot.Along[1],
                                .Index  = hump.Index,
                                .Placed = 0.0,
                                .Bent   = true,
                                .Own    = spot.Own});
      }
      return true;
    }

    /**
     *  Whether no two crossings are one, as the thorough way merges them: two are one only where
     *  they lie within the tolerance of each other or where, halfway between them along each curve,
     *  each curve is within twice the tolerance of the other. Each pair lies farther apart than the
     *  tolerance, and the curve halfway between them farther than twice it from the polyline, both
     *  with the margin: halfway the plain way, and the short way round where the curve may close on
     *  itself.
     *
     *  @param  lacing The curve, the polyline and the crossings.
     *  @return        Whether none are one.
     */
    auto Unlaced(const Lacing &lacing) -> bool
    {
      const Curve &strand = *lacing.Shape;
      const double length = lacing.Tolerance.Length;
      const double clear  = 2.0 * length * Margin;
      const Range  domain = Domain(strand);
      const double whole  = domain.End - domain.Start;
      const auto   far    = [&lacing, &strand, clear](double along) -> bool
      {
        return Shunned(*lacing.Chain, At(strand, along)) > clear;
      };
      const std::span<const Stitch> found = lacing.Found;
      for (std::size_t one = 0; one < found.size(); ++one)
      {
        for (std::size_t two = one + 1; two < found.size(); ++two)
        {
          const double from   = found[one].Along;
          const double middle = from + ((found[two].Along - from) / 2.0);
          const double round  = middle + (whole / 2.0);
          if (!(Distance(found[one].At, found[two].At) > length * Margin)
              || !far(middle)
              || !far(round > domain.End ? round - whole : round))
          {
            return false;
          }
        }
      }
      return true;
    }

    /**
     *  The crossings made hits, placed as the thorough way places them: it polishes a spline's or a
     *  Bézier curve's crossing on a straight segment by both places, the polyline's the double
     *  nearest its own, which the kernel proved, and works an ellipse's out exactly, placed by the
     *  foot; on a bulge it places the point by its angle, as Located does; on a spline or a Bézier
     *  curve at the place the kernel proved, on an ellipse by the point's angle. A crossing Located
     *  would not place leaves fewer hits than crossings, and the kernel not sure.
     *
     *  @param  lacing The curve, the polyline and the crossings.
     *  @param  side   Which of the two the curve is: 0 for the first.
     *  @param  into   Where the hits go.
     *  @return        Whether every crossing was placed.
     */
    auto Stitches(const Lacing &lacing, std::size_t side, Intersections &into) -> bool
    {
      const Polyline &chain  = *lacing.Chain;
      std::size_t     placed = 0;
      for (const Stitch &stitch : lacing.Found)
      {
        const std::optional<Place> chained =
            stitch.Bent || !stitch.Own
                ? Crossings::Located(chain, static_cast<double>(stitch.Index) + 0.5, stitch.At)
                : std::optional<Place>{
                      Place{.At = At(chain, stitch.Placed), .Along = stitch.Placed}};
        const std::optional<Place> stranded =
            stitch.Own ? std::optional<Place>{Place{.At    = At(*lacing.Shape, stitch.Along),
                                                    .Along = stitch.Along}}
                       : Crossings::Located(*lacing.Shape, 0.5, stitch.At);
        const Place other  = chained.value_or(Place{});
        const Place mine   = stranded.value_or(Place{});
        placed            += static_cast<std::size_t>(chained.has_value())
                           * static_cast<std::size_t>(stranded.has_value());
        into.Hits.push_back({.At     = stitch.At,
                             .First  = side == 0 ? mine : other,
                             .Second = side == 0 ? other : mine,
                             .Kind   = Contact::Cross});
      }
      return placed == lacing.Found.size();
    }

    /**
     *  Whether a curve is met with a polyline leg by leg: a spline, an ellipse, or a quadratic or
     *  cubic Bézier curve.
     *
     *  @param  curve The curve.
     *  @return       Whether it is.
     */
    auto Laceable(const Curve &curve) -> bool
    {
      return std::holds_alternative<Spline>(curve)
          || std::holds_alternative<Ellipse>(curve)
          || std::holds_alternative<Bezier::Quadratic>(curve)
          || std::holds_alternative<Bezier::Cubic>(curve);
    }

    /**
     *  A quadratic or cubic Bézier curve as a piece.
     *
     *  @param  curve The curve.
     *  @return       Its section; nothing for any other curve.
     */
    auto Beziered(const Curve &curve) -> std::optional<Section>
    {
      if (const auto *quadratic = std::get_if<Bezier::Quadratic>(&curve))
      {
        return Sectioned(*quadratic);
      }
      if (const auto *cubic = std::get_if<Bezier::Cubic>(&curve))
      {
        return Sectioned(*cubic);
      }
      return std::nullopt;
    }

    /**
     *  The box of points, each within the range the kernels take.
     *
     *  @param  points The points.
     *  @return        Their box; nothing where there are none, or one lies out of that range.
     */
    auto Penned(std::span<const Point> points) -> std::optional<Box>
    {
      if (points.empty()
          || !std::ranges::all_of(points,
                                  [](Point point) -> bool { return Held({point.X, point.Y}); }))
      {
        return std::nullopt;
      }
      Box hull = {};
      for (const Point point : points)
      {
        hull = Grown(hull, point);
      }
      return hull;
    }

    /**
     *  A box a curve keeps in: a spline's poles', a Bézier curve's control points', a whole
     *  ellipse's, o +- (|X_x| + |Y_x|, |X_y| + |Y_y|), as its point o + X cos + Y sin keeps each
     *  coordinate within that of o's, grown past the roundings of the sums and o's.
     *
     *  @param  curve The curve: a spline, an ellipse or a quadratic or cubic Bézier curve.
     *  @return       The box; nothing where a number lies outside the range the kernels take, or a
     *                spline has no poles.
     */
    auto Fenced(const Curve &curve) -> std::optional<Box>
    {
      if (const auto *spline = std::get_if<Spline>(&curve))
      {
        return Penned(spline->Poles);
      }
      if (const auto *quadratic = std::get_if<Bezier::Quadratic>(&curve))
      {
        return Penned(std::array{quadratic->Start, quadratic->Control, quadratic->End});
      }
      if (const auto *cubic = std::get_if<Bezier::Cubic>(&curve))
      {
        return Penned(std::array{cubic->Start, cubic->Out, cubic->In, cubic->End});
      }
      const Matrix &map = std::get<Ellipse>(curve).Map;
      if (!Held({map.X.X, map.X.Y, map.Y.X, map.Y.Y, map.Origin.X, map.Origin.Y}))
      {
        return std::nullopt;
      }
      const Vector half = {.X = ((std::abs(map.X.X) + std::abs(map.Y.X)) * (1.0 + 0x1p-50))
                              + (0x1p-50 * std::abs(map.Origin.X)),
                           .Y = ((std::abs(map.X.Y) + std::abs(map.Y.Y)) * (1.0 + 0x1p-50))
                              + (0x1p-50 * std::abs(map.Origin.Y))};
      return Box{.Min = map.Origin - half, .Max = map.Origin + half};
    }
  } // namespace

  auto Laced(const Curves &curves, const Tolerance &tolerance, Intersections &into)
      -> std::optional<bool>
  {
    const std::size_t side  = Laceable(curves[0].get()) ? 0 : 1;
    const Curve      &shape = curves.at(side).get();
    const auto       *chain = std::get_if<Polyline>(&curves.at(1 - side).get());
    if (!Laceable(shape) || chain == nullptr)
    {
      return std::nullopt;
    }
    const auto                        *strand = std::get_if<Spline>(&shape);
    const std::optional<List<Section>> sections =
        strand != nullptr ? Sectioned(*strand) : std::nullopt;
    Lacing                   lacing = {.Shape     = &shape,
                                       .Strand    = strand,
                                       .Chain     = chain,
                                       .Tolerance = tolerance,
                                       .Hull      = {},
                                       .Sections  = sections ? std::span<const Section>{*sections}
                                                             : std::span<const Section>{},
                                       .Piece     = Beziered(shape),
                                       .Found     = {}};
    const std::size_t        count  = Count(*chain);
    const std::optional<Box> hull   = Fenced(shape);
    const bool               held =
        std::ranges::all_of(chain->Vertices, [](const Vertex &vertex) -> bool
                            { return Held({vertex.Position.X, vertex.Position.Y, vertex.Bulge}); });
    if (count == 0 || !hull || !held)
    {
      return false;
    }
    lacing.Hull = *hull;
    /**
     *  Every vertex is the end of a segment of some length, as a run of segments of none lies at
     *  one point with the ends of those beside it, where there is one: each end is then weighed
     *  with it.
     */
    const double clear = Afar * Reach(tolerance);
    bool         some  = false;
    for (std::size_t index = 0; index < count; ++index)
    {
      const Vertex &start = chain->Vertices[index];
      const Vertex &end   = chain->Vertices[(index + 1) % chain->Vertices.size()];
      if (start.Position == end.Position)
      {
        continue;
      }
      some = true;
      /**
       *  A straight segment's box is its ends'; a bulge's is passed by the ring kernel, about its
       *  arc's box.
       */
      const bool   straight = !(std::abs(start.Bulge) >= Straight);
      const double bend     = straight ? 0.0 : start.Bulge;
      if (straight
          && Cloven(Grown(Grown(Box{}, start.Position), end.Position), lacing.Hull) > clear)
      {
        continue;
      }
      const Hump hump = {
          .Start = start.Position, .End = end.Position, .Bend = bend, .Index = index};
      if (!(straight ? Straddled(lacing, index) : Overarched(lacing, hump)))
      {
        return false;
      }
    }
    if (!some || !Unlaced(lacing))
    {
      into.Hits.clear();
      return false;
    }
    /**
     *  A crossing Located would not place leaves the kernel unsure, and the hits none.
     */
    const bool placed = Stitches(lacing, side, into);
    into.Hits.resize(into.Hits.size() * static_cast<std::size_t>(placed));
    std::ranges::sort(into.Hits, {}, [](const Hit &hit) -> double { return hit.First.Along; });
    return placed;
  }
} // namespace Geoxy::Swift
