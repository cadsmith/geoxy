/***************************************************************************************************
 *  @file      source/piecewise.cpp
 *  @brief     A polyline met segment by segment: each segment a line or an arc of its own, met
 *             directly with each piece of the other curve whose box comes near it, where the
 *             thorough way, every pair of parts weighed at once, would give the same.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "piecewise.hpp"

#include "bow.hpp"
#include "crossing.hpp"
#include "elliptic.hpp"
#include "hull.hpp"
#include "met/linear.hpp"
#include "met/round.hpp"
#include "met/slight.hpp"
#include "real.hpp"
#include "straight.hpp"
#include "swift/lane.hpp"
#include "swift/met.hpp"
#include "thorough.hpp"
#include "trig.hpp"

#include "geoxy/angle.hpp"
#include "geoxy/arc.hpp"
#include "geoxy/at.hpp"
#include "geoxy/box.hpp"
#include "geoxy/curve.hpp"
#include "geoxy/derivative.hpp"
#include "geoxy/domain.hpp"
#include "geoxy/extents.hpp"
#include "geoxy/intersect.hpp"
#include "geoxy/line.hpp"
#include "geoxy/nearest.hpp"
#include "geoxy/orient.hpp"
#include "geoxy/place.hpp"
#include "geoxy/point.hpp"
#include "geoxy/polyline.hpp"
#include "geoxy/ray.hpp"
#include "geoxy/segment.hpp"
#include "geoxy/tolerance.hpp"
#include "geoxy/vector.hpp"
#include "geoxy/xline.hpp"

#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <memory>
#include <optional>
#include <span>
#include <variant>
#include <vector>

namespace Geoxy
{
  namespace
  {
    /**
     *  How far, in tolerances, the thorough way may polish a crossing to from where it is exactly:
     *  a place it keeps has the two curves within the tolerance, and two curves crossing at a sine
     *  of an eighth or more are that near only within eight tolerances of the crossing. Twice that.
     */
    constexpr double Drift = 16.0;

    /**
     *  How far, in tolerances, everything that ties one pair of pieces' hits to another's is kept
     *  from them: the curves' ends and vertices, other hits, where the curves come near without
     *  crossing. Two pieces whose boxes come this near are met; any farther apart meet in nothing.
     */
    constexpr double Clear = 256.0;

    /**
     *  How far apart, in tolerances, two hits are kept: the thorough way holds two together that it
     *  polished to within the tolerance of each other, and each lies within Drift tolerances of
     *  where it is exactly. Twice that, and the tolerance.
     */
    constexpr double Apart = (2.0 * Drift) + 32.0;

    /**
     *  The most segments of a polyline met with a straight curve by the swift kernel first.
     */
    constexpr std::size_t Few = 8;

    /**
     *  The sine below which a crossing is left to the thorough way, as the direct ways leave it:
     *  under it the curves stay near each other so long that a band about the crossing may form.
     */
    constexpr double Slant = 0.125;

    /**
     *  The most legs a polyline may have for the legs near a box or a point to be found by weighing
     *  each, not by a grid.
     */
    constexpr std::size_t Handful = 24;

    /**
     *  What a leg is: a polyline's straight segment or a line; a ray; an endless line; a polyline's
     *  arc, its exact circle the bulge's; or an arc of its own.
     */
    enum class Form : std::uint8_t
    {
      Straight, /**< A polyline's straight segment, or a line. */
      Ray,      /**< A ray.                                    */
      Endless,  /**< An endless line.                          */
      Bow,      /**< A polyline's arc.                         */
      Round,    /**< An arc.                                   */
    };

    /**
     *  A piece of a curve: a polyline's segment, as a line or an arc; any other curve whole.
     */
    struct Leg
    {
      mutable Curve Shape = Line{};            /**< The piece, as a curve of its own: a
                                                    polyline's arc only once asked for.  */
      Box                  Hull   = {};        /**< A box it cannot leave.          */
      std::array<Point, 2> Ends   = {};        /**< Its ends.                       */
      std::size_t          Count  = 0;         /**< How many ends it has.           */
      double               Offset = 0.0;       /**< The curve's parameter at its 0. */
      double               Speed  = 0.0;       /**< How far it moves for each 1 of its
                                                    parameter.                           */
      double Bulge = 0.0;                      /**< A polyline's arc's bulge. */
      double Size  = 0.0;                      /**< Its numbers' size, as its grain is
                                                    taken from.                          */
      Point           Center = {};             /**< An arc's center, as its shape's. */
      double          Radius = 0.0;            /**< An arc's radius, as its shape's. */
      Vector          Way    = {};             /**< A straight leg's way.            */
      const Polyline *From   = nullptr;        /**< The polyline it is a segment of. */
      std::size_t     Index  = 0;              /**< Which segment.                   */
      Form            Kind   = Form::Straight; /**< What it is.                      */
      mutable bool    Shaped = true;           /**< Whether its shape is made.       */

      /**
       *  Whether it is a polyline's arc, its exact circle the bulge's, not the arc's.
       *
       *  @return Whether it is.
       */
      [[nodiscard]] auto Bulged() const -> bool { return Kind == Form::Bow; }

      /**
       *  Whether it is an arc, a polyline's or one of its own.
       *
       *  @return Whether it is.
       */
      [[nodiscard]] auto Bent() const -> bool { return Kind == Form::Bow || Kind == Form::Round; }

      /**
       *  Whether it is straight: a line, a ray, an endless line or a straight segment.
       *
       *  @return Whether it is.
       */
      [[nodiscard]] auto Lined() const -> bool { return !Bent(); }
    };

    /**
     *  A leg as a curve of its own, a polyline's arc made from its bulge the first time it is asked
     *  for, as Nth makes it. A leg kept by a curve made ready has every arc made already, so that
     *  nothing of it is written here.
     *
     *  @param  leg The leg.
     *  @return     Its shape.
     */
    auto Formed(const Leg &leg) -> const Curve &
    {
      if (!leg.Shaped)
      {
        leg.Shape  = std::get<Arc>(Nth(*leg.From, leg.Index));
        leg.Shaped = true;
      }
      return leg.Shape;
    }

    /**
     *  How far a polyline's arc moves for each 1 of its parameter: its length, as its sweep takes
     *  an arc tangent.
     *
     *  @param  leg The leg, a polyline's arc.
     *  @return     Its speed.
     */
    auto Swiftness(const Leg &leg) -> double
    {
      return leg.Radius * std::abs(4.0 * Trig::Atan(leg.Bulge));
    }

    /**
     *  How far a leg moves for each 1 of its parameter: a polyline's arc's length, worked out only
     *  once asked for, where it is not kept.
     *
     *  @param  leg The leg.
     *  @return     Its speed.
     */
    auto PerUnit(const Leg &leg) -> double
    {
      return leg.Bulged() && !(leg.Speed > 0.0) ? Swiftness(leg) : leg.Speed;
    }

    /**
     *  The size of a point's numbers grown by a reach, as the thorough way's grain weighs it: the
     *  larger of its coordinates' magnitudes, each with the reach added, finite ones only.
     *
     *  @param  point The point.
     *  @param  reach The reach.
     *  @return       The size.
     */
    auto Outsized(Point point, double reach) -> double
    {
      double size = 0.0;
      for (const double value : {std::abs(point.X) + reach, std::abs(point.Y) + reach})
      {
        size = std::isfinite(value) ? std::max(size, value) : size;
      }
      return size;
    }

    /**
     *  A polyline's segment as a leg: its circle worked out from its bulge as Nth's arc has it, its
     *  arc made only once asked for.
     *
     *  @param  polyline The polyline.
     *  @param  index    The segment.
     *  @return          The leg.
     */
    auto Stranded(const Polyline &polyline, std::size_t index) -> Leg
    {
      const std::span<const Vertex> vertices = polyline.Vertices;
      const Point                   start    = vertices[index].Position;
      const Point                   end      = vertices[(index + 1) % vertices.size()].Position;
      const double                  bulge    = vertices[index].Bulge;
      Leg                           leg      = {.Hull   = Hull(polyline, index),
                                                .Ends   = {start, end},
                                                .Count  = 2,
                                                .Offset = static_cast<double>(index),
                                                .Bulge  = bulge,
                                                .From   = &polyline,
                                                .Index  = index};
      if (!(std::abs(bulge) >= Straight) || start == end)
      {
        leg.Shape = Line{.Start = start, .End = end};
        leg.Way   = end - start;
        leg.Speed = Distance(start, end);
        leg.Size  = std::max(Outsized(start, 0.0), Outsized(end, 0.0));
        return leg;
      }
      /**
       *  The center and the radius as Bent works them out.
       */
      const Hoop hoop = Hooped(start, end, bulge);
      leg.Kind        = Form::Bow;
      leg.Shaped      = false;
      leg.Center      = hoop.Center;
      leg.Radius      = hoop.Radius;
      leg.Size        = Outsized(leg.Center, std::abs(leg.Radius));
      return leg;
    }

    /**
     *  A curve that is not a polyline as one leg: a line, a ray, an endless line, or an arc that
     *  neither closes on itself nor turns round more than once.
     *
     *  @param  curve The curve.
     *  @return       Its leg; nothing for any other curve, nor where a number of it is not finite.
     */
    auto OneLeg(const Curve &curve) -> std::optional<Leg>
    {
      Leg leg = {.Shape = curve, .Hull = {}, .Size = 0.0};
      if (const auto *line = std::get_if<Line>(&curve))
      {
        leg.Ends  = {line->Start, line->End};
        leg.Count = 2;
        leg.Speed = Distance(line->Start, line->End);
        leg.Way   = line->End - line->Start;
      }
      else if (const auto *ray = std::get_if<Ray>(&curve))
      {
        leg.Ends  = {ray->Start, ray->Start};
        leg.Count = 1;
        leg.Speed = Length(ray->Way);
        leg.Way   = ray->Way;
        leg.Kind  = Form::Ray;
      }
      else if (const auto *xline = std::get_if<Xline>(&curve))
      {
        leg.Ends  = {xline->Through, xline->Through};
        leg.Speed = Length(xline->Way);
        leg.Way   = xline->Way;
        leg.Kind  = Form::Endless;
      }
      else if (const auto *arc = std::get_if<Arc>(&curve))
      {
        const double sweep = std::abs(arc->Sweep.InRadians());
        leg.Ends           = {At(*arc, 0.0), At(*arc, 1.0)};
        leg.Count          = 2;
        leg.Speed          = arc->Radius * sweep;
        leg.Center         = arc->Center;
        leg.Radius         = arc->Radius;
        leg.Kind           = Form::Round;
        if (!(sweep < Revolution) || leg.Ends[0] == leg.Ends[1])
        {
          return std::nullopt;
        }
      }
      else
      {
        return std::nullopt;
      }
      /**
       *  Boxed and sized once it is a leg: a spline's tight box, worked out for nothing, cost a
       *  quarter of meeting it with a polyline.
       */
      leg.Hull         = Extents(curve);
      leg.Size         = Sized(curve);
      const double sum = leg.Ends[0].X + leg.Ends[0].Y + leg.Ends[1].X + leg.Ends[1].Y + leg.Speed;
      return std::isfinite(sum) && leg.Speed > 0.0 ? std::optional<Leg>{leg} : std::nullopt;
    }

    /**
     *  A curve as legs: a polyline's segments, any other curve whole.
     *
     *  @param  curve The curve.
     *  @param  legs  Where its legs go, emptied first.
     *  @return       Whether it is met piecewise: not a polyline of no segment, with a number not
     *                finite among them, or with a segment so bulged that its box runs past every
     *                double.
     */
    auto Legged(const Curve &curve, std::vector<Leg> &legs) -> bool
    {
      legs.clear();
      if (const auto *polyline = std::get_if<Polyline>(&curve))
      {
        double sum = 0.0;
        for (const Vertex &vertex : polyline->Vertices)
        {
          sum += vertex.Position.X + vertex.Position.Y + vertex.Bulge;
        }
        const std::size_t count = Count(*polyline);
        if (count == 0 || !std::isfinite(sum))
        {
          return false;
        }
        double span = 0.0;
        for (std::size_t index = 0; index < count; ++index)
        {
          legs.push_back(Stranded(*polyline, index));
          const Box &hull  = legs.back().Hull;
          span            += hull.Max.X - hull.Min.X + hull.Max.Y - hull.Min.Y;
        }
        return std::isfinite(span);
      }
      const std::optional<Leg> whole = OneLeg(curve);
      if (whole)
      {
        legs.push_back(*whole);
      }
      return whole.has_value();
    }

    /**
     *  How far a point lies from a leg's carrier or circle at the least, where a few roundings
     *  leave it, less what rounding may take off; nothing past minus infinity for a straight leg of
     *  no length, which has no carrier.
     *
     *  @param  leg   The leg.
     *  @param  point The point.
     *  @return       The distance, at the least.
     */
    auto Aloof(const Leg &leg, Point point) -> double
    {
      if (leg.Bent())
      {
        const double from  = Distance(point, leg.Center);
        const double slack = 0x1p-46 * (from + std::abs(leg.Radius));
        return std::abs(from - leg.Radius) - slack;
      }
      /**
       *  A straight leg's speed is its way's length, worked out once.
       */
      const Vector way   = leg.Way;
      const Vector gap   = point - leg.Ends[0];
      const double speed = leg.Speed;
      if (!(speed > 0.0))
      {
        return -std::numeric_limits<double>::infinity();
      }
      const double slack = 0x1p-50 * (std::abs(gap.X * way.Y) + std::abs(gap.Y * way.X));
      return (std::abs(Cross(gap, way)) - slack) / (speed * (1.0 + 0x1p-48));
    }

    /**
     *  Whether a point lies farther than a reach from a leg: by the distance at the least from its
     *  box, then from its carrier or its circle, where either settles it, by the leg's nearest
     *  point where neither does.
     *
     *  @param  leg   The leg, its shape made where it is weighed whole.
     *  @param  point The point.
     *  @param  reach The reach.
     *  @return       Whether it does.
     */
    auto Farther(const Leg &leg, Point point, double reach) -> bool
    {
      /**
       *  The point is finite, an end of a leg or a point of a curve between two hits, and a length
       *  as Length works it out is never below its larger component: the square root of a rounded
       *  square is the number itself, and rounding keeps the order of a sum.
       */
      const Vector over = Overhang(leg.Hull, point);
      if (Largest(over) > reach)
      {
        return true;
      }
      const double boxed = Length(over);
      return boxed                                           > reach
          || std::max(boxed, Aloof(leg, point))              > reach
          || Distance(point, Nearest(Formed(leg), point).At) > reach;
    }

    /**
     *  A grid of square cells, about as many as a curve's legs, each listing the legs whose boxes
     *  reach into it: the few legs near a point, found without weighing them all.
     */
    struct Stalls
    {
      Point                      Corner  = {};    /**< Where the first cell starts.      */
      double                     Cell    = 1;     /**< A cell's side.                    */
      std::size_t                Columns = 1;     /**< How many cells across.            */
      std::size_t                Rows    = 1;     /**< How many cells up.                */
      std::vector<std::uint32_t> Starts  = {};    /**< Where each cell's list starts.    */
      std::vector<std::uint32_t> Items   = {};    /**< The legs listed, cell after cell. */
      std::vector<std::uint32_t> Next    = {};    /**< Where each cell's list goes on.   */
      bool                       Laid    = false; /**< Whether the legs are filed.       */

      /**
       *  The cells a box reaches into, clamped to the grid.
       *
       *  @param  box The box.
       *  @return     The first and the last column, the first and the last row.
       */
      [[nodiscard]] auto Cells(const Box &box) const -> std::array<std::size_t, 4>
      {
        const auto cell = [this](double from, double value, std::size_t count) -> std::size_t
        {
          const double index = std::floor((value - from) / Cell);
          return static_cast<std::size_t>(std::clamp(index, 0.0, static_cast<double>(count - 1)));
        };
        return {cell(Corner.X, box.Min.X, Columns), cell(Corner.X, box.Max.X, Columns),
                cell(Corner.Y, box.Min.Y, Rows), cell(Corner.Y, box.Max.Y, Rows)};
      }

      /**
       *  The legs listed in the cells a box reaches into, a row of cells at a time, as a span of
       *  the list: a leg in two cells comes twice.
       *
       *  @tparam VisitT What is done with each row's legs.
       *  @param  box    The box.
       *  @param  visit  What is done with a row's legs, a span of their indices.
       */
      template <typename VisitT>
      void Each(const Box &box, const VisitT &visit) const
      {
        const auto [left, right, low, high] = Cells(box);
        for (std::size_t row = low; row <= high; ++row)
        {
          const std::size_t first = (row * Columns) + left;
          const std::size_t last  = (row * Columns) + right;
          visit(std::span{Items}.subspan(Starts[first], Starts[last + 1] - Starts[first]));
        }
      }
    };

    /**
     *  A polyline's legs filed in a grid: cells about as many as the legs and no smaller than a
     *  leg's box is across on average; one cell where every box is one point.
     *
     *  @param legs The legs, each of a finite box.
     *  @param grid The grid, its lists' room kept from call to call.
     */
    void Filed(std::span<const Leg> legs, Stalls &grid)
    {
      Box    bounds;
      double across = 0.0;
      for (const Leg &leg : legs)
      {
        bounds             = Joined(bounds, leg.Hull);
        const Vector size  = Size(leg.Hull);
        across            += std::max(size.X, size.Y);
      }
      const auto   count = static_cast<double>(legs.size());
      const Vector size  = Size(bounds);
      const double side  = std::max(
          {across / count, std::sqrt(size.X * size.Y / count), std::max(size.X, size.Y) / 1024.0});
      grid.Corner  = bounds.Min;
      grid.Cell    = side > 0.0 ? side : 1.0;
      grid.Columns = 1 + static_cast<std::size_t>(size.X / grid.Cell);
      grid.Rows    = 1 + static_cast<std::size_t>(size.Y / grid.Cell);
      grid.Starts.assign((grid.Columns * grid.Rows) + 1, 0);
      const auto each = [&grid, legs](const auto &file) -> void
      {
        for (std::size_t index = 0; index < legs.size(); ++index)
        {
          const auto [left, right, low, high] = grid.Cells(legs[index].Hull);
          for (std::size_t row = low; row <= high; ++row)
          {
            for (std::size_t column = left; column <= right; ++column)
            {
              file((row * grid.Columns) + column, index);
            }
          }
        }
      };
      each([&grid](std::size_t cell, std::size_t /*index*/) -> void { ++grid.Starts[cell + 1]; });
      for (std::size_t index = 1; index < grid.Starts.size(); ++index)
      {
        grid.Starts[index] += grid.Starts[index - 1];
      }
      grid.Items.resize(grid.Starts.back());
      grid.Next.assign(grid.Starts.begin(), grid.Starts.end() - 1);
      each([&grid](std::size_t cell, std::size_t index) -> void
           { grid.Items[grid.Next[cell]++] = static_cast<std::uint32_t>(index); });
      grid.Laid = true;
    }
  } // namespace

  /**
   *  What a polyline or an arc made ready keeps to be met piece by piece, worked out once.
   */
  struct Legwork
  {
    std::vector<Leg> Legs  = {};    /**< Its legs, every arc made and its speed kept.    */
    Stalls           Filed = {};    /**< Their grid, where they are more than a handful. */
    bool             Sound = false; /**< Whether it is met piecewise, as Legged tells.   */
  };

  namespace
  {
    /**
     *  Where a point is asked about, and how near.
     */
    struct Probe
    {
      Point  At    = {};  /**< The point.       */
      double Reach = 0.0; /**< How near counts. */
    };

    /**
     *  A hit of two legs, put on the curves.
     */
    struct Found
    {
      std::array<std::size_t, 2> Legs = {}; /**< Which leg of each curve it is on. */
      Geoxy::Hit                 Hit  = {}; /**< The hit, on the curves.           */
    };

    /**
     *  An end of one leg lying exactly on the other leg: a vertex or an end of one curve on a
     *  segment of the other, or one both share.
     */
    struct Joint
    {
      Point                      At   = {}; /**< The end.                               */
      std::array<std::size_t, 2> Legs = {}; /**< The two legs, the first curve's first. */
    };

    /**
     *  The room a meeting works in, kept from call to call on each thread, so that a meeting takes
     *  no memory once a call has needed as much.
     */
    struct Room
    {
      std::array<std::vector<Leg>, 2>           Legs   = {}; /**< Each curve's legs, not kept. */
      std::array<Stalls, 2>                     Grids  = {}; /**< Their grids, not kept.       */
      std::vector<std::array<std::size_t, 2>>   Pairs  = {}; /**< The pairs of legs to meet.   */
      std::vector<std::size_t>                  Stamps = {}; /**< The leg last weighed.        */
      std::vector<std::size_t>                  Order  = {}; /**< Hits or arcs, in order.      */
      std::vector<Found>                        Hits   = {}; /**< The hits found.              */
      std::vector<std::array<const Found *, 2>> Asked  = {}; /**< Pairs of hits to tell.       */
      std::vector<Joint>                        Joints = {}; /**< Ends lying on the other.     */
      Intersections                             Met    = {}; /**< What two legs meet in.       */
      std::array<std::array<Vector, 2>, 2>      Ways   = {}; /**< Two arcs' tangents at each
                                                                hit, where Arcwards met them.  */
      Polyline Single = {.Vertices = {Vertex{}, Vertex{}}};  /**< A polyline's arc alone. */
    };

    /**
     *  Two curves met piecewise, and what is worked out once about them.
     */
    struct Pairing
    {
      std::array<const Curve *, 2>        Curves    = {};      /**< The two curves.         */
      Geoxy::Tolerance                    Tolerance = {};      /**< The tolerance.          */
      Room                               *Work      = nullptr; /**< Where the work is done. */
      std::array<std::span<const Leg>, 2> Strides   = {};      /**< Each curve's legs.      */
      std::array<const Legwork *, 2>      Kept      = {};      /**< Each one's legwork, where
                                                                    kept.                   */

      /**
       *  A curve's legs.
       *
       *  @param  side Which curve, 0 the first.
       *  @return      Its legs.
       */
      [[nodiscard]] auto Legs(std::size_t side) const -> std::span<const Leg>
      {
        return Strides.at(side);
      }

      /**
       *  Whether a curve is a polyline.
       *
       *  @param  side Which curve, 0 the first.
       *  @return      Whether it is.
       */
      [[nodiscard]] auto Strung(std::size_t side) const -> bool
      {
        return std::holds_alternative<Polyline>(*Curves.at(side));
      }
    };

    /**
     *  A polyline's legs filed: as kept, or the first time they are asked for in a meeting.
     *
     *  @param  pairing The curves being met.
     *  @param  side    Which curve.
     *  @return         Its grid.
     */
    auto Filing(const Pairing &pairing, std::size_t side) -> const Stalls &
    {
      if (const Legwork *kept = pairing.Kept.at(side); kept != nullptr)
      {
        return kept->Filed;
      }
      Stalls &grid = pairing.Work->Grids.at(side);
      if (!grid.Laid)
      {
        Filed(pairing.Legs(side), grid);
      }
      return grid;
    }

    /**
     *  Whether no leg of a curve comes within a reach of a point: each weighed for a curve of a
     *  handful of legs, those its grid files near the point for any other.
     *
     *  @param  pairing The curves being met.
     *  @param  side    Which curve.
     *  @param  probe   The point and the reach.
     *  @return         Whether none does.
     */
    auto Lonely(const Pairing &pairing, std::size_t side, Probe probe) -> bool
    {
      const std::span<const Leg> legs  = pairing.Legs(side);
      const auto                 alone = [&probe](const Leg &leg) -> bool
      {
        return Farther(leg, probe.At, probe.Reach);
      };
      if (legs.size() <= Handful)
      {
        return std::ranges::all_of(legs, alone);
      }
      const Box   reach   = {.Min = {.X = probe.At.X - probe.Reach, .Y = probe.At.Y - probe.Reach},
                             .Max = {.X = probe.At.X + probe.Reach, .Y = probe.At.Y + probe.Reach}};
      std::size_t crowded = 0;
      Filing(pairing, side)
          .Each(reach,
                [&legs, &alone, &crowded](std::span<const std::uint32_t> row) -> void
                {
                  crowded += static_cast<std::size_t>(
                      !std::ranges::all_of(row, [&legs, &alone](std::uint32_t leg) -> bool
                                           { return alone(legs[leg]); }));
                });
      return crowded == 0;
    }

    /**
     *  Whether two arcs' circles are one to within twice the tolerance, as two rounded from one
     *  circle are: the thorough way's pair of their parts then shares a stretch wherever they lie.
     *
     *  @param  pairing The curves being met.
     *  @param  legs    One leg of each.
     *  @return         Whether they are; never for a straight leg.
     */
    auto Concentric(const Pairing &pairing, std::array<const Leg *, 2> legs) -> bool
    {
      return legs[0]->Bent()
          && legs[1]->Bent()
          && Distance(legs[0]->Center, legs[1]->Center)
                     + std::abs(legs[0]->Radius - legs[1]->Radius)
                 <= 2.0 * pairing.Tolerance.Length;
    }

    /**
     *  Whether two legs share a stretch as the thorough way's pair of their parts weighs one: two
     *  straight legs by the thorough way's own rule, unless their stretch the second's projection
     *  covers is certainly too short or too far from it; two arcs where their circles are one; a
     *  straight leg and an arc never do.
     *
     *  @param  pairing The curves being met.
     *  @param  legs    One leg of each.
     *  @return         Whether they do.
     */
    auto Sharing(const Pairing &pairing, std::array<const Leg *, 2> legs) -> bool
    {
      if (legs[0]->Lined() && legs[1]->Lined())
      {
        const auto run = [](const Leg &leg) -> Run
        {
          return {.Start = leg.Ends[0], .Way = leg.Way, .Along = Domain(leg.Shape)};
        };
        const Run  one   = run(*legs[0]);
        const Run  other = run(*legs[1]);
        const bool parallel =
            Orient(Point{}, Point{} + one.Way, Point{} + other.Way) == Orientation::Collinear;
        return (parallel || !Swift::Unaligned(one, other, pairing.Tolerance))
            && Shares({legs[0]->Shape, legs[1]->Shape}, pairing.Tolerance);
      }
      return Concentric(pairing, {legs[0], legs[1]});
    }

    /**
     *  Whether two legs come so near touching somewhere that the direct way and the thorough way,
     *  which work out their points apart, may part on whether they meet: an arc whose circle comes
     *  within a clearance of touching a straight leg's carrier, or of touching another arc's
     *  circle.
     *
     *  @param  legs  One leg of each.
     *  @param  clear How near is near.
     *  @return       Whether they do.
     */
    auto Skimming(std::array<const Leg *, 2> legs, double clear) -> bool
    {
      const Leg &one   = *legs[0];
      const Leg &other = *legs[1];
      if (one.Bent() && other.Bent())
      {
        const double apart = Distance(one.Center, other.Center);
        return std::abs(apart - (one.Radius + other.Radius))         <= clear
            || std::abs(apart - std::abs(one.Radius - other.Radius)) <= clear;
      }
      if (one.Lined() && other.Lined())
      {
        return false;
      }
      /**
       *  The center's distance from the carrier against the radius, both times the way's length, so
       *  that a straight leg of no length, a point, is as near as its circle is to it.
       */
      const Leg   &round    = one.Bent() ? one : other;
      const Leg   &straight = one.Bent() ? other : one;
      const double way      = Length(straight.Way);
      const double cross =
          std::abs(Cross((round.Center - Point{}) - (straight.Ends[0] - Point{}), straight.Way));
      return std::abs(cross - (round.Radius * way)) <= clear * way;
    }

    /**
     *  The gap between two boxes, along x or along y, whichever is wider: nothing where they meet.
     *
     *  @param  one   One box.
     *  @param  other The other.
     *  @return       The gap.
     */
    auto BoxesGap(const Box &one, const Box &other) -> double
    {
      return std::max({one.Min.X - other.Max.X, other.Min.X - one.Max.X, one.Min.Y - other.Max.Y,
                       other.Min.Y - one.Max.Y, 0.0});
    }

    /**
     *  Whether an arc of one curve and an arc of the other are of one circle, near each other or
     *  not: the second curve's arcs taken in order of their radii, those within twice the tolerance
     *  of each of the first's weighed.
     *
     *  @param  pairing The curves being met.
     *  @return         Whether two are.
     */
    auto Cocircular(const Pairing &pairing) -> bool
    {
      const double               twice  = 2.0 * pairing.Tolerance.Length;
      const std::span<const Leg> ones   = pairing.Legs(0);
      const std::span<const Leg> others = pairing.Legs(1);
      std::vector<std::size_t>  &arcs   = pairing.Work->Order;
      if (std::ranges::none_of(ones, &Leg::Bent))
      {
        return false;
      }
      arcs.clear();
      for (std::size_t index = 0; index < others.size(); ++index)
      {
        if (others[index].Bent())
        {
          arcs.push_back(index);
        }
      }
      const auto radius = [&others](std::size_t index) -> double
      {
        return others[index].Radius;
      };
      std::ranges::sort(arcs, {}, radius);
      for (const Leg &leg : ones)
      {
        if (!leg.Bent())
        {
          continue;
        }
        for (auto at = std::ranges::lower_bound(arcs, leg.Radius - twice, {}, radius);
             at != arcs.end() && radius(*at) <= leg.Radius + twice; ++at)
        {
          if (Concentric(pairing, {&leg, &others[*at]}))
          {
            return true;
          }
        }
      }
      return false;
    }

    /**
     *  The pairs of legs to meet: those whose boxes come within the clearance, every two weighed
     *  where a polyline has a handful of legs, found by a grid of its legs where it has more;
     *  nothing where two arcs share a stretch, their circles one, wherever they lie. Two straight
     *  legs farther apart than twice the tolerance never do: a stretch of one whose ends lie within
     *  the tolerance of the other's carrier has the other within it of the stretch, or crossing it.
     *
     *  @param  pairing The curves being met; the pairs go in its room, by index, the first curve's
     *                  leg first.
     *  @return         Whether there are pairs to meet so.
     */
    auto Paired(const Pairing &pairing) -> bool
    {
      if (Cocircular(pairing))
      {
        return false;
      }
      const bool                 flip   = !pairing.Strung(1);
      const auto                 filing = static_cast<std::size_t>(flip ? 0 : 1);
      const std::span<const Leg> filed  = pairing.Legs(filing);
      const std::span<const Leg> asking = pairing.Legs(1 - filing);
      const double               clear  = Clear * pairing.Tolerance.Length;
      auto                      &pairs  = pairing.Work->Pairs;
      pairs.clear();
      const auto pair = [flip](std::size_t ask, std::size_t leg) -> std::array<std::size_t, 2>
      {
        return flip ? std::array{leg, ask} : std::array{ask, leg};
      };
      if (filed.size() <= Handful)
      {
        for (std::size_t ask = 0; ask < asking.size(); ++ask)
        {
          for (std::size_t leg = 0; leg < filed.size(); ++leg)
          {
            if (BoxesGap(asking[ask].Hull, filed[leg].Hull) <= clear)
            {
              pairs.push_back(pair(ask, leg));
            }
          }
        }
        return true;
      }
      const Stalls             &grid   = Filing(pairing, filing);
      std::vector<std::size_t> &stamps = pairing.Work->Stamps;
      stamps.assign(filed.size(), asking.size());
      for (std::size_t ask = 0; ask < asking.size(); ++ask)
      {
        const Box &hull = asking[ask].Hull;
        grid.Each(Padded(hull, clear),
                  [&](std::span<const std::uint32_t> row) -> void
                  {
                    for (const std::uint32_t leg : row)
                    {
                      if (stamps[leg] != ask && BoxesGap(hull, filed[leg].Hull) <= clear)
                      {
                        pairs.push_back(pair(ask, leg));
                      }
                      stamps[leg] = ask;
                    }
                  });
      }
      return true;
    }

    /**
     *  Whether a leg is gentle enough to look for a crossing on as the direct way does: a straight
     *  one always, an arc where its radius is over a clearance.
     *
     *  @param  leg   The leg.
     *  @param  clear The clearance.
     *  @return       Whether it is.
     */
    auto Mellow(const Leg &leg, double clear) -> bool
    {
      return leg.Lined() || leg.Radius > clear;
    }

    /**
     *  Whether a point lies exactly on a straight segment or a line of some length: at an end of
     *  it, or on its carrier, exactly, between its ends.
     *
     *  @param  leg   The leg.
     *  @param  point The point.
     *  @return       Whether it does.
     */
    auto Upon(const Leg &leg, Point point) -> bool
    {
      if (leg.Kind != Form::Straight || leg.Ends[0] == leg.Ends[1])
      {
        return false;
      }
      if (point == leg.Ends[0] || point == leg.Ends[1])
      {
        return true;
      }
      return Orient(leg.Ends[0], leg.Ends[1], point) == Orientation::Collinear
          && Dot(point - leg.Ends[0], leg.Way)       > 0.0
          && Dot(point - leg.Ends[1], leg.Way)       < 0.0;
    }

    /**
     *  Whether every end of each of two legs lies beyond a clearance of the other, or, where both
     *  are straight segments or lines of some length, exactly on it: each such end filed as a joint
     *  of the two, which the thorough way meets with as the hit both its curves pass through.
     *
     *  @param  pairing The curves being met; the joints go in its room.
     *  @param  pair    The legs, by index, the first curve's first.
     *  @param  clear   The clearance.
     *  @return         Whether they do.
     */
    auto Abutted(const Pairing &pairing, std::array<std::size_t, 2> pair, double clear) -> bool
    {
      const std::array<const Leg *, 2> legs = {&pairing.Legs(0)[pair[0]],
                                               &pairing.Legs(1)[pair[1]]};
      for (const std::size_t side : {0UZ, 1UZ})
      {
        const Leg &ended   = *legs.at(side);
        const Leg &against = *legs.at(1 - side);
        for (const Point end : std::span{ended.Ends}.first(ended.Count))
        {
          if (Farther(against, end, clear))
          {
            continue;
          }
          if (!Upon(ended, end) || !Upon(against, end))
          {
            return false;
          }
          pairing.Work->Joints.push_back({.At = end, .Legs = pair});
        }
      }
      return true;
    }

    /**
     *  The tolerance two legs are held to: the tolerance, or their own grain where an arc's center
     *  and radius are coarser than the polyline's vertices, as Usable takes it from their shapes.
     *
     *  @param  legs      One leg of each.
     *  @param  tolerance The tolerance.
     *  @return           The tolerance they are held to.
     */
    auto Tolerated(std::array<const Leg *, 2> legs, Tolerance tolerance) -> Tolerance
    {
      const double size = std::max(legs[0]->Size, legs[1]->Size);
      tolerance.Length  = std::max(tolerance.Length, 16.0 * Above(size));
      return tolerance;
    }

    /**
     *  A crossing put on the curves made a hit of two legs, in the meeting's room.
     *
     *  @param pairing The curves being met.
     *  @param pair    The legs, by index, the first curve's first.
     *  @param crossed The crossing, put on the curves.
     */
    void Put(const Pairing &pairing, std::array<std::size_t, 2> pair,
             const Crossings::Crossed &crossed)
    {
      pairing.Work->Hits.push_back(
          {.Legs = pair,
           .Hit  = {.At     = crossed.At,
                    .First  = {.At = crossed.On[0], .Along = crossed.Along[0]},
                    .Second = {.At = crossed.On[1], .Along = crossed.Along[1]},
                    .Kind   = Contact::Cross}});
    }

    /**
     *  Whether two ways cross at a sine above a least one.
     *
     *  @param  ways  The two ways.
     *  @param  least The least sine.
     *  @return       Whether they do.
     */
    auto Steeply(std::array<Vector, 2> ways, double least) -> bool
    {
      return std::abs(Cross(ways[0], ways[1])) > least * Length(ways[0]) * Length(ways[1]);
    }

    /**
     *  A straight leg and a polyline's arc met by the swift kernel of a line and a polyline, the
     *  arc a polyline of its own: each crossing worked out from the bulge's exact circle in
     *  double-double and proved the double nearest the true one, as the exact way rounds it, and
     *  placed on the curves as the exact way places it.
     *
     *  @param  pairing The curves being met; the hits go in its room.
     *  @param  pair    The legs, by index, the first curve's first, one straight and one an arc.
     *  @param  steep   The least sine.
     *  @return         Whether they were met so: nothing where the kernel is not sure or a crossing
     *                  lies on an axis, false where a crossing is shallower than the least sine or
     *                  placed by none.
     */
    auto Bellied(const Pairing &pairing, std::array<std::size_t, 2> pair, double steep)
        -> std::optional<bool>
    {
      const Leg &one      = pairing.Legs(0)[pair[0]];
      const Leg &other    = pairing.Legs(1)[pair[1]];
      const bool flip     = one.Bulged();
      const Leg &bent     = flip ? one : other;
      const Leg &straight = flip ? other : one;
      Polyline  &single   = pairing.Work->Single;
      single.Vertices[0]  = {.Position = bent.Ends[0], .Bulge = bent.Bulge};
      single.Vertices[1]  = {.Position = bent.Ends[1], .Bulge = 0.0};
      const std::optional<Swift::Spots> spots =
          Swift::Laned(straight.Shape)
              .and_then([&single, &pairing](const Swift::Lane &lane) -> std::optional<Swift::Spots>
                        { return Swift::Bulged(lane, single, pairing.Tolerance); });
      /**
       *  A crossing on an axis is left to the direct way, which works it out again by the exact
       *  way, whose nothing has a sign of its own.
       */
      const auto axes = [](const Swift::Spots &found) -> std::size_t
      {
        std::size_t count = 0;
        for (const Swift::Spot &spot : found)
        {
          count +=
              static_cast<std::size_t>(Zero(spot.At.X)) + static_cast<std::size_t>(Zero(spot.At.Y));
        }
        return count;
      };
      if (spots.transform(axes).value_or(1) > 0)
      {
        return std::nullopt;
      }
      const auto placed = [&pairing, pair, flip, steep](const Swift::Spots &found) -> bool
      {
        const Leg  &curved = pairing.Legs(flip ? 0 : 1)[pair.at(flip ? 0 : 1)];
        const Leg  &level  = pairing.Legs(flip ? 1 : 0)[pair.at(flip ? 1 : 0)];
        std::size_t lost   = 0;
        for (const Swift::Spot &spot : found)
        {
          const double                along = level.Offset + spot.Along[0];
          const double                round = curved.Offset + 0.5;
          const std::array<double, 2> places =
              flip ? std::array{round, along} : std::array{along, round};
          if (!Steeply({level.Way, Perpendicular(spot.At - curved.Center)}, steep))
          {
            return false;
          }
          const std::optional<Crossings::Crossed> exact =
              Crossings::Placed(pairing.Curves, places, spot.At);
          lost += static_cast<std::size_t>(!exact.has_value());
          Put(pairing, pair, exact.value_or(Crossings::Crossed{}));
        }
        return lost == 0;
      };
      return spots.transform(placed);
    }

    /**
     *  Two arcs, a polyline's or ones of their own, met as the direct way meets two arcs, its steps
     *  taken here: Straightly takes no arc, Piecewise no curve but a polyline, Concentrically none
     *  but arcs of one circle, which Sharing leaves no pair of, and Swiftless asks no swift kernel;
     *  so Roughly's Circularly gives its crossings, each made exact as Exactly makes it, by Aiming
     *  as by Crossing with the same roots, put in order along the first arc, and made an end where
     *  Circularly keeps an end at it. Each hit's tangents are Derivative's, from the directions
     *  Aiming kept: the turn At would take, once. Hits on a whole circle's seam, ordered by Seamed,
     *  are left to the direct way.
     *
     *  @param  pairing The curves being met; the hits and tangents go in its room.
     *  @param  legs    The two arcs, the first curve's first.
     *  @param  held    The tolerance they are held to.
     *  @return         Whether they were met so; where not, the direct way turns to Banded, and the
     *                  caller to the direct way.
     */
    auto Arcwards(const Pairing &pairing, std::array<const Leg *, 2> legs, const Tolerance &held)
        -> bool
    {
      const Curve                     &first  = Formed(*legs[0]);
      const Curve                     &second = Formed(*legs[1]);
      const Arc                       &round  = std::get<Arc>(first);
      const Arc                       &across = std::get<Arc>(second);
      const std::optional<Met::Nearby> near   = Met::Circularly(round, across, held);
      if (!near)
      {
        return false;
      }
      /**
       *  Anything the direct way's steps do otherwise is counted, not branched on, and leaves the
       *  pair to the direct way itself: a whole circle's seam, a crossing not made exact or placed
       *  off its arc's sweep, a touch, an end at a crossing.
       */
      Intersections         &met = pairing.Work->Met;
      const Crossings::Roots roots =
          near->Held().empty() ? Crossings::Roots{} : Crossings::Quickly({&first, &second});
      auto        other = static_cast<std::size_t>(!(std::abs(round.Sweep.InRadians()) < 2.0 * Pi));
      std::size_t index = 0;
      for (const Approach &approach : near->Held())
      {
        const std::optional<Crossings::Bearings> aimed =
            Crossings::Aiming({&first, &second}, approach, roots);
        const Crossings::Bearings found = aimed.value_or(Crossings::Bearings{});
        other += static_cast<std::size_t>(!aimed.has_value())
               + static_cast<std::size_t>(approach.Touch)
               + static_cast<std::size_t>(near->Ending(found.At, held.Length))
               + static_cast<std::size_t>(!found.Ways[0].has_value())
               + static_cast<std::size_t>(!found.Ways[1].has_value());
        met.Hits.push_back({.At     = found.At,
                            .First  = {.At = {}, .Along = found.Along[0]},
                            .Second = {.At = {}, .Along = found.Along[1]},
                            .Kind   = Contact::Cross});
        const auto tangent = [&found](const Arc &arc, std::size_t side) -> Vector
        {
          return Perpendicular(found.Ways.at(side).value_or(Vector{}))
               * (arc.Radius * arc.Sweep.InRadians());
        };
        pairing.Work->Ways.at(index++) = {tangent(round, 0), tangent(across, 1)};
      }
      /**
       *  In order along the first arc, as Exactly sorts them: two at the most, swapped where the
       *  second comes first.
       */
      if (met.Hits.size() == 2 && met.Hits[1].First.Along < met.Hits[0].First.Along)
      {
        std::swap(met.Hits[0], met.Hits[1]);
        std::swap(pairing.Work->Ways[0], pairing.Work->Ways[1]);
      }
      return other == 0;
    }

    /**
     *  Two legs met directly, their hits put on the curves: nothing where the direct ways leave the
     *  pair to the thorough way, where a hit is no plain crossing at the least sine or more, or
     *  where no exact way places a crossing on the curves themselves.
     *
     *  @param  pairing The curves being met; the hits go in its room.
     *  @param  pair    The legs, by index, the first curve's first.
     *  @param  held    The tolerance they are held to.
     *  @return         Whether they were met so.
     */
    auto Directly(const Pairing &pairing, std::array<std::size_t, 2> pair, const Tolerance &held)
        -> bool
    {
      const Leg     &one    = pairing.Legs(0)[pair[0]];
      const Leg     &other  = pairing.Legs(1)[pair[1]];
      const Curve   &first  = *pairing.Curves[0];
      const Curve   &second = *pairing.Curves[1];
      Intersections &met    = pairing.Work->Met;
      met.Hits.clear();
      met.Overlaps.clear();
      met.Spent = false;
      /**
       *  A pair that shares no stretch meets in no stretch, and spends no budget: a line's and an
       *  arc's parts are solved in closed form.
       */
      const bool arced = one.Bent() && other.Bent() && Arcwards(pairing, {&one, &other}, held);
      if (!arced)
      {
        met.Hits.clear();
      }
      if (!arced && !Direct({Formed(one), Formed(other)}, held, met))
      {
        return false;
      }
      const double steep = std::max(Slant, 2.0 * Trig::Sin(held.Radians));
      std::size_t  lost  = 0;
      for (std::size_t index = 0; index < met.Hits.size(); ++index)
      {
        const Hit   &hit = met.Hits[index];
        const Vector way =
            arced ? pairing.Work->Ways.at(index)[0] : Derivative(one.Shape, hit.First.Along);
        const Vector across =
            arced ? pairing.Work->Ways.at(index)[1] : Derivative(other.Shape, hit.Second.Along);
        if (hit.Kind != Contact::Cross || !Steeply({way, across}, steep))
        {
          return false;
        }
        const std::array<double, 2> places = {one.Offset + hit.First.Along,
                                              other.Offset + hit.Second.Along};
        /**
         *  Each crossing ends as the thorough way ends it. Two straight pieces, or a straight one
         *  and an arc, cross at the double nearest the carriers' exact crossing, which the direct
         *  way found already: it is placed on each curve as the exact way places it. On a
         *  polyline's arc, which the direct way took as an arc rounded from its bulge, the crossing
         *  is worked out again on the curves themselves, the bulge by its exact circle, and
         *  polished in double-double where that finds none; and so is one on an axis, whose nothing
         *  the exact way gives a sign of its own. A crossing placed by none leaves the pair to the
         *  thorough way, which keeps its own polish there.
         */
        const bool                              axis  = Zero(hit.At.X) || Zero(hit.At.Y);
        const bool                              again = one.Bulged() || other.Bulged() || axis;
        const std::optional<Crossings::Crossed> exact =
            (again ? Crossings::Crossing({&first, &second}, places, hit.At)
                   : Crossings::Placed({&first, &second}, places, hit.At))
                .or_else([&first, &second, &places, &hit] -> std::optional<Crossings::Crossed>
                         { return Crossings::Polished({&first, &second}, places, hit.At); });
        lost += static_cast<std::size_t>(!exact.has_value());
        Put(pairing, pair, exact.value_or(Crossings::Crossed{}));
      }
      return lost == 0;
    }

    /**
     *  Whether a straight leg certainly does not cross a polyline's arc whose circle its carrier
     *  cuts: each of the two points they meet at, worked out in doubles, lies clearly past an end
     *  of the leg or clearly on the side of the chord the arc is not on, by many times what
     *  rounding may move it.
     *
     *  @param  straight The straight leg.
     *  @param  bent     The polyline's arc.
     *  @return          Whether it certainly does not; never where the carrier misses the circle.
     */
    auto Aside(const Leg &straight, const Leg &bent) -> bool
    {
      const Vector way    = straight.Way;
      const double square = Dot(way, way);
      const Vector toward = bent.Center - straight.Ends[0];
      const double off    = Cross(toward, way);
      const double whole  = bent.Radius * bent.Radius * square;
      /**
       *  The carrier cuts the circle deeper than the clearance, as nothing grazes: the root is of a
       *  number well above nothing, its rounding small beside it.
       */
      const double root = std::sqrt(std::max(whole - (off * off), 0.0));
      const double foot = Dot(toward, way) / square;
      const double half = root / square;
      const double size = Length(toward) + bent.Radius;
      const double slack =
          0x1p-40 * ((size * std::sqrt(square)) + ((whole + (off * off)) / root)) / square;
      const Range  span  = Domain(straight.Shape);
      const Vector chord = bent.Ends[1] - bent.Ends[0];
      const double turn  = std::copysign(1.0, bent.Bulge);
      /**
       *  A point clearly past the leg's stretch, or clearly on the side of the chord the arc is not
       *  on: a bulge above zero lies right of its chord.
       */
      const auto apart = [&](double sign) -> bool
      {
        const double along  = foot + (sign * half);
        const Point  point  = straight.Ends[0] + (way * along);
        const double side   = Cross(chord, point - bent.Ends[0]) * turn;
        const double margin = 0x1p-36 * Length(chord) * (size + (slack * std::sqrt(square)));
        const double past   = std::max(span.Start - slack - along, along - span.End - slack);
        return std::max(past, side - margin) > 0.0;
      };
      return apart(-1.0) && apart(1.0);
    }

    /**
     *  Whether a point lies clearly on the side of an arc's chord the arc is not on, by more than a
     *  reach: a bulge or an arc swept counterclockwise lies right of its chord, as Swept tells.
     *
     *  @param  leg   The arc, a polyline's or one of its own, swept less than a whole turn.
     *  @param  point The point.
     *  @param  reach How far past the chord's line.
     *  @return       Whether it does.
     */
    auto Unswept(const Leg &leg, Point point, double reach) -> bool
    {
      const Vector chord = leg.Ends[1] - leg.Ends[0];
      const double turn  = leg.Bulged() ? leg.Bulge : std::get<Arc>(leg.Shape).Sweep.InRadians();
      return Cross(chord, point - leg.Ends[0]) * std::copysign(1.0, turn) > reach * Length(chord);
    }

    /**
     *  Whether two arcs whose circles cross meet in nothing, as Met::Circularly finds on the arcs
     *  the direct way takes them as, its center and radius each leg's: their circles cross at a
     *  sine above Steep and do not kiss, worked out as Circularly works them out; each crossing,
     *  worked out as it works them out, lies clearly off one arc, on the side of its chord it is
     *  not on, farther from the chord's line than the tolerance, so that it is neither swept nor
     *  near an end; and every end of each lies clearly farther than four tolerances from the
     *  other's circle, so that Settle keeps none. Clearly is by a slack of 2^-30 of their size,
     *  many times what Circularly's cuts, the library's sine and cosine of the arcs' angles, may
     *  lie from the ends here. Circularly then gives nothing near, and the direct way nothing.
     *
     *  @param  legs   One leg of each, both arcs swept less than a whole turn, of radii over the
     *                 clearance, their circles crossing and not skimming.
     *  @param  length The tolerance the pair is held to.
     *  @return        Whether they meet in nothing so.
     */
    auto Estranged(std::array<const Leg *, 2> legs, double length) -> bool
    {
      /**
       *  The circles do not kiss: Skimming tells so at the clearance, by the same distance and
       *  radii, far past Circularly's four tolerances.
       */
      const Leg   &one   = *legs[0];
      const Leg   &other = *legs[1];
      const double clear = 4.0 * length;
      const double inner = one.Radius;
      const double outer = other.Radius;
      const double apart = Distance(one.Center, other.Center);
      const double cosine =
          ((inner * inner) + (outer * outer) - (apart * apart)) / (2.0 * inner * outer);
      if (!((1.0 - (cosine * cosine)) > Met::Steep * Met::Steep))
      {
        return false;
      }
      const Leg   &small  = inner <= outer ? one : other;
      const Leg   &large  = inner <= outer ? other : one;
      const double less   = small.Radius;
      const double more   = large.Radius;
      const Vector toward = (large.Center - small.Center) / apart;
      const double along  = ((apart * apart) + (less * less) - (more * more)) / (2.0 * apart);
      const double height = std::sqrt(std::max(0.0, (less * less) - (along * along)));
      const Point  base   = small.Center + (toward * along);
      const Vector rise   = Perpendicular(toward) * height;
      const double slack =
          0x1p-30 * std::max(Outsized(one.Center, inner), Outsized(other.Center, outer));
      const auto far = [clear, slack](const Leg &ended, const Leg &against) -> bool
      {
        return std::ranges::all_of(
            ended.Ends, [&against, clear, slack](Point end) -> bool
            { return std::abs(Distance(end, against.Center) - against.Radius) > clear + slack; });
      };
      const auto off = [&one, &other, length, slack](Point point) -> bool
      {
        return Unswept(one, point, length + slack) || Unswept(other, point, length + slack);
      };
      return far(one, other) && far(other, one) && off(base + rise) && off(base - rise);
    }

    /**
     *  Whether two legs certainly do not cross, where every end of each lies beyond the clearance
     *  of the other, no circle grazes the other's carrier or circle and no stretch is shared: two
     *  straight segments where both ends of one lie on one side of the other's carrier, exactly; a
     *  straight leg whose carrier misses an arc's circle, or cuts a polyline's arc's circle only
     *  off the leg or off the arc; two arcs whose circles lie apart, or one within the other, or
     *  cross only off the arcs, as Estranged tells. Two such legs lie farther apart than the
     *  clearance, as their ends do, and meet in nothing.
     *
     *  @param  legs   One leg of each.
     *  @param  length The tolerance the pair is held to.
     *  @return        Whether they certainly do not.
     */
    auto Missed(std::array<const Leg *, 2> legs, double length) -> bool
    {
      const Leg &one   = *legs[0];
      const Leg &other = *legs[1];
      if (one.Bent() && other.Bent())
      {
        const double apart = Distance(one.Center, other.Center);
        return apart > one.Radius + other.Radius
            || apart < std::abs(one.Radius - other.Radius)
            || Estranged(legs, length);
      }
      if (one.Bent() != other.Bent())
      {
        const Leg &round    = one.Bent() ? one : other;
        const Leg &straight = one.Bent() ? other : one;
        return std::abs(Cross(round.Center - straight.Ends[0], straight.Way))
                 > round.Radius * Length(straight.Way)
            || (round.Bulged() && Aside(straight, round));
      }
      if (one.Kind != Form::Straight || other.Kind != Form::Straight)
      {
        return false;
      }
      const auto aside = [](const Leg &carrier, const Leg &ended) -> bool
      {
        const Orientation start = Orient(carrier.Ends[0], carrier.Ends[1], ended.Ends[0]);
        return start != Orientation::Collinear
            && start == Orient(carrier.Ends[0], carrier.Ends[1], ended.Ends[1]);
      };
      return aside(one, other) || aside(other, one);
    }

    /**
     *  Which side of a polyline's arc's chord a point lies on, clearly: by more than a reach of the
     *  chord's length past its line, the side the arc is on or the other.
     *
     *  @param  leg   The polyline's arc.
     *  @param  point The point.
     *  @param  reach How far past the chord's line.
     *  @return       1 on the arc's side, -1 on the other, 0 where neither clearly.
     */
    auto Flank(const Leg &leg, Point point, double reach) -> int
    {
      const Vector chord = leg.Ends[1] - leg.Ends[0];
      const double turn  = leg.Bulged() ? leg.Bulge : std::get<Arc>(leg.Shape).Sweep.InRadians();
      const double off   = Cross(chord, point - leg.Ends[0]) * std::copysign(1.0, turn);
      const double edge  = reach * Length(chord);
      return static_cast<int>(off < -edge) - static_cast<int>(off > edge);
    }

    /**
     *  The most a polyline's arc may be bulged for it to be met from its circle alone: its sweep,
     *  four arc tangents of its bulge, then lies clear of a whole turn.
     */
    constexpr double Tame = 0x1p20;

    /**
     *  Whether an arc, a polyline's or one of its own, sweeps clearly less than a whole turn: a
     *  polyline's arc bulged less than Tame, an arc swept short of a turn by 2^-20 of one.
     *
     *  @param  leg The arc.
     *  @return     Whether it does.
     */
    auto Tamed(const Leg &leg) -> bool
    {
      return leg.Bulged() ? std::abs(leg.Bulge) < Tame
                          : std::abs(std::get<Arc>(leg.Shape).Sweep.InRadians())
                                < Revolution * (1.0 - 0x1p-20);
    }

    /**
     *  A polyline's arc and another's, or an arc of its own, met from their circles alone, where
     *  the direct way's steps for two arcs (Arcwards, then Directly) come clearly to the same hits:
     *  each swept clearly short of a whole turn, as Tamed tells; their circles crossing, as
     *  Circularly crosses them, at a sine clearly above Steep, Slant and twice the tolerance's
     *  angle; every end of each clearly farther than four tolerances from the other's circle, so
     *  that Circularly keeps none, by a slack of 2^-30 of their size as Estranged takes it; and
     *  each crossing, worked out as Circularly works it out, clearly on one side of each arc's
     *  chord or the other by the tolerance and that slack, so that it lies on the arc or off it as
     *  Circularly sweeps it, and far from its ends. A crossing on both is put on the curves as
     *  Directly puts it, polished on the exact circles, which Crossing leaves to Polished, from as
     *  near as Aiming's point, on the segments the hit's places name; no arc is made, nor its
     *  circles crossed again, its roots found in double-double or its crossings aimed at. A
     *  crossing not polished so leaves the pair to the thorough way, as Directly leaves it. An arc
     *  of its own and a bulge are polished by Newton's method alone, whose last bit may hang on
     *  where it starts at a coordinate under 2^-8 of the size about it: such a crossing is left to
     *  the direct way's steps, as Retold leaves two bulges' to Bowed.
     *
     *  @param  pairing The curves being met; the hits go in its room.
     *  @param  pair    The legs, by index, the first curve's first, two arcs, a polyline's and
     *                  another's or one of its own, of radii over the clearance, whose circles
     *                  cross, as Missed and Skimming leave them, no nearer touching than it.
     *  @param  held    The tolerance they are held to.
     *  @return         Whether they were met so, as Directly would tell; nothing where it is not
     *                  clear.
     */
    auto Bowwise(const Pairing &pairing, std::array<std::size_t, 2> pair, const Tolerance &held)
        -> std::optional<bool>
    {
      const Leg   &one   = pairing.Legs(0)[pair[0]];
      const Leg   &other = pairing.Legs(1)[pair[1]];
      const double inner = one.Radius;
      const double outer = other.Radius;
      const double apart = Distance(one.Center, other.Center);
      /**
       *  The sine as Circularly weighs it, clear of what the directions Aiming keeps may make of
       *  it: its cosine's rounding grows with the ratio of the radii.
       */
      const double cosine =
          ((inner * inner) + (outer * outer) - (apart * apart)) / (2.0 * inner * outer);
      const double steep = std::max({Met::Steep, Slant, 2.0 * Trig::Sin(held.Radians)});
      const double ratio = std::max(inner, outer) / std::min(inner, outer);
      const double clear = 4.0 * held.Length;
      const double size  = std::max(Outsized(one.Center, inner), Outsized(other.Center, outer));
      const double slack = 0x1p-30 * size;
      const auto   far   = [clear, slack](const Leg &ended, const Leg &against) -> bool
      {
        return std::ranges::all_of(
            ended.Ends, [&against, clear, slack](Point end) -> bool
            { return std::abs(Distance(end, against.Center) - against.Radius) > clear + slack; });
      };
      std::size_t unclear =
          static_cast<std::size_t>(!Tamed(one))
          + static_cast<std::size_t>(!Tamed(other))
          + static_cast<std::size_t>(
              !((1.0 - (cosine * cosine)) > (steep * steep) + 0x1p-24 + (0x1p-40 * ratio)))
          + static_cast<std::size_t>(!far(one, other))
          + static_cast<std::size_t>(!far(other, one));
      const Leg   &small  = inner <= outer ? one : other;
      const Leg   &large  = inner <= outer ? other : one;
      const double less   = small.Radius;
      const double more   = large.Radius;
      const Vector toward = (large.Center - small.Center) / apart;
      const double along  = ((apart * apart) + (less * less) - (more * more)) / (2.0 * apart);
      const double height = std::sqrt(std::max(0.0, (less * less) - (along * along)));
      const Point  base   = small.Center + (toward * along);
      const Vector rise   = Perpendicular(toward) * height;
      const double reach  = held.Length + slack;
      const bool   bows   = one.Bulged() && other.Bulged();
      std::array<Point, 2> points = {};
      std::size_t          count  = 0;
      for (const Point point : {base + rise, base - rise})
      {
        const double least   = std::min(std::abs(point.X), std::abs(point.Y));
        const double around  = size + std::abs(point.X) + std::abs(point.Y);
        unclear             += static_cast<std::size_t>(!bows && !(least >= 0x1p-8 * around));
        const int onOne      = Flank(one, point, reach);
        const int onOther    = Flank(other, point, reach);
        unclear += static_cast<std::size_t>(onOne == 0) + static_cast<std::size_t>(onOther == 0);
        points.at(count)  = point;
        count            += static_cast<std::size_t>(std::min(onOne, onOther) > 0);
      }
      if (unclear > 0)
      {
        return std::nullopt;
      }
      std::size_t                       lost  = 0;
      std::array<Crossings::Crossed, 2> found = {};
      for (std::size_t index = 0; index < count; ++index)
      {
        const std::optional<Crossings::Crossed> exact = Crossings::Polished(
            pairing.Curves, {one.Offset + 0.5, other.Offset + 0.5}, points.at(index));
        lost            += static_cast<std::size_t>(!exact.has_value());
        found.at(index)  = exact.value_or(Crossings::Crossed{});
      }
      for (const Crossings::Crossed &crossed : std::span{found}.first(count))
      {
        Put(pairing, pair, crossed);
      }
      return lost == 0;
    }

    /**
     *  Two legs met, their hits put on the curves: nothing where the thorough way may answer
     *  otherwise, an arc too tight, an end within the clearance of the other leg but not exactly on
     *  it, a stretch shared or a touch grazed; two straight legs an end of one lies on by the joint
     *  filed, a straight leg and a polyline's arc by the swift kernel, any other pair directly.
     *
     *  @param  pairing The curves being met; the hits go in its room.
     *  @param  pair    The legs, by index, the first curve's first.
     *  @return         Whether they were met so.
     */
    auto Intersected(const Pairing &pairing, std::array<std::size_t, 2> pair) -> bool
    {
      const Leg &one   = pairing.Legs(0)[pair[0]];
      const Leg &other = pairing.Legs(1)[pair[1]];
      /**
       *  The pair is held to the tolerance, or to its own grain where an arc's center and radius
       *  are coarser than the polyline's vertices: the tolerance its own thorough way, and so the
       *  direct way, holds it to, and the clearance kept that many times over.
       */
      const Tolerance   held   = Tolerated({&one, &other}, pairing.Tolerance);
      const double      clear  = Clear * held.Length;
      const std::size_t joints = pairing.Work->Joints.size();
      const bool        clean  = Mellow(one, clear)
                              && Mellow(other, clear)
                              && !Skimming({&one, &other}, clear)
                              && Abutted(pairing, pair, clear)
                              && !Sharing(pairing, {&one, &other});
      if (!clean)
      {
        return false;
      }
      /**
       *  Two straight legs, not along each other, meet at most once: where an end of one lies on
       *  the other, they meet there and nowhere else.
       */
      if (pairing.Work->Joints.size() > joints || Missed({&one, &other}, held.Length))
      {
        return true;
      }
      if (one.Bent() && other.Bent())
      {
        if (const std::optional<bool> told = Bowwise(pairing, pair, held))
        {
          return *told;
        }
      }
      const bool bowed =
          one.Bulged() != other.Bulged() && one.Kind != Form::Round && other.Kind != Form::Round;
      /**
       *  Where the swift kernel is not sure, as where an end of the straight leg lies near the
       *  arc's circle away from the arc, the pair is met directly.
       */
      const std::optional<bool> swift =
          bowed ? Bellied(pairing, pair, std::max(Slant, 2.0 * Trig::Sin(held.Radians)))
                : std::nullopt;
      return swift.has_value() ? *swift : Directly(pairing, pair, held);
    }

    /**
     *  How a curve passes a joint: inside one leg, through the vertex two legs share, or at its
     *  end, the legs it passes by.
     */
    struct Passage
    {
      std::array<std::size_t, 2> Legs   = {};    /**< Its legs there.               */
      std::size_t                Count  = 0;     /**< How many: one or two.         */
      bool                       Ending = false; /**< Whether the joint is its end. */
    };

    /**
     *  How a curve passes a joint, from the legs it was met with there: inside the one leg; at the
     *  vertex both of the legs it lies between share; or at the curve's end, a line's or an open
     *  polyline's. Nothing where those legs are not all and only the ones it passes by, nor at a
     *  closed polyline's seam, where the thorough way may place it at either end of the domain.
     *
     *  @param  pairing The curves being met.
     *  @param  side    Which curve.
     *  @param  joints  The joints at one point, each with the two legs met there.
     *  @return         How it passes there.
     */
    auto Passed(const Pairing &pairing, std::size_t side, std::span<const Joint> joints)
        -> std::optional<Passage>
    {
      const std::span<const Leg> legs  = pairing.Legs(side);
      const Point                point = joints.front().At;
      std::array<std::size_t, 2> met   = {};
      std::size_t                count = 0;
      for (const Joint &joint : joints)
      {
        const std::size_t leg = joint.Legs.at(side);
        if (std::ranges::contains(std::span{met}.first(count), leg))
        {
          continue;
        }
        if (count == met.size())
        {
          return std::nullopt;
        }
        met.at(count++) = leg;
      }
      const std::span<const std::size_t> taken = std::span{met}.first(count);
      const auto                         ended = std::ranges::find_if(
          taken, [&legs, point](std::size_t leg) -> bool
          { return point == legs[leg].Ends[0] || point == legs[leg].Ends[1]; });
      if (ended == taken.end())
      {
        return count == 1 ? std::optional<Passage>{{.Legs = met, .Count = 1}} : std::nullopt;
      }
      const auto *polyline = std::get_if<Polyline>(pairing.Curves.at(side));
      if (polyline == nullptr)
      {
        return Passage{.Legs = met, .Count = 1, .Ending = true};
      }
      const std::size_t total  = legs.size();
      const std::size_t vertex = point == legs[*ended].Ends[0] ? *ended : *ended + 1;
      const bool        seam   = vertex == 0 || vertex == total;
      if (polyline->Closed && seam)
      {
        return std::nullopt;
      }
      const std::size_t last    = vertex == 0 ? 0 : total - 1;
      const Passage     passage = seam ? Passage{.Legs = {last, last}, .Count = 1, .Ending = true}
                                       : Passage{.Legs = {vertex - 1, vertex}, .Count = 2};
      return std::ranges::is_permutation(taken, std::span{passage.Legs}.first(passage.Count))
               ? std::optional<Passage>{passage}
               : std::nullopt;
    }

    /**
     *  Whether two crossings are one to the bit: their point, and each curve's place and point.
     *
     *  @param  one   One.
     *  @param  other The other.
     *  @return       Whether they are.
     */
    auto Twin(const Crossings::Crossed &one, const Crossings::Crossed &other) -> bool
    {
      const auto bits = [](const Crossings::Crossed &crossed) -> std::array<std::uint64_t, 8>
      {
        return {std::bit_cast<std::uint64_t>(crossed.At.X),
                std::bit_cast<std::uint64_t>(crossed.At.Y),
                std::bit_cast<std::uint64_t>(crossed.Along[0]),
                std::bit_cast<std::uint64_t>(crossed.Along[1]),
                std::bit_cast<std::uint64_t>(crossed.On[0].X),
                std::bit_cast<std::uint64_t>(crossed.On[0].Y),
                std::bit_cast<std::uint64_t>(crossed.On[1].X),
                std::bit_cast<std::uint64_t>(crossed.On[1].Y)};
      };
      return bits(one) == bits(other);
    }

    /**
     *  The hit where both curves pass through a joint, as the thorough way finds it, in the
     *  meeting's room. Its candidates there are polished to the joint or kept within the tolerance
     *  of it, each on one of the legs either curve passes by, or at the vertex between two; their
     *  run is merged into one hit, an end where the joint is an end of either curve, else a
     *  crossing where every way in or out of the joint of one curve passes every way of the other
     *  to one side, steeper than the least sine, so that a candidate wherever placed is judged one;
     *  and the hit is put exactly where the legs' carriers cross, the joint itself, placed on the
     *  leg a candidate's place names: taken only where every leg either curve passes by places it
     *  to the bit alike.
     *
     *  @param  pairing The curves being met; the hit goes in its room.
     *  @param  joints  The joints at one point, each with the two legs met there.
     *  @return         Whether it was found so: never where a curve passes the joint otherwise than
     *                  at a vertex, inside a leg or at its end, nor where it may be a touch.
     */
    auto AtJoint(const Pairing &pairing, std::span<const Joint> joints) -> bool
    {
      const std::optional<Passage> one   = Passed(pairing, 0, joints);
      const std::optional<Passage> other = one ? Passed(pairing, 1, joints) : std::nullopt;
      if (!other)
      {
        return false;
      }
      const Point  point  = joints.front().At;
      const bool   ending = one->Ending || other->Ending;
      const double steep  = std::max(Slant, 2.0 * Trig::Sin(pairing.Tolerance.Radians));
      const std::span<const std::size_t> ones   = std::span{one->Legs}.first(one->Count);
      const std::span<const std::size_t> others = std::span{other->Legs}.first(other->Count);
      /**
       *  Every leg either curve passes by places the joint alike: where none does, or two place it
       *  apart, the joint is left, though the carriers of two straight legs at a sine of an eighth
       *  or more cross at it exactly, wherever it is placed.
       */
      int                passes = 0;
      std::size_t        faults = 0;
      std::size_t        lost   = 0;
      Crossings::Crossed made   = {};
      for (const std::size_t first : ones)
      {
        for (const std::size_t second : others)
        {
          const Leg &leg     = pairing.Legs(0)[first];
          const Leg &across  = pairing.Legs(1)[second];
          const int  passing = Cross(across.Way, leg.Way) > 0.0 ? 1 : -1;
          const std::optional<Crossings::Crossed> exact =
              Crossings::Crossing(pairing.Curves, {leg.Offset + 0.5, across.Offset + 0.5}, point);
          const Crossings::Crossed crossed  = exact.value_or(Crossings::Crossed{});
          faults                           += Steeply({leg.Way, across.Way}, steep) ? 0U : 1U;
          faults                           += !ending && passes != 0 && passing != passes ? 1 : 0;
          lost                             += static_cast<std::size_t>(!exact.has_value());
          lost   += static_cast<std::size_t>(passes != 0 && !Twin(made, crossed));
          passes  = passing;
          made    = crossed;
        }
      }
      if (faults + lost > 0)
      {
        return false;
      }
      const auto slowest = [&pairing](std::size_t                  side,
                                      std::span<const std::size_t> legs) -> std::size_t
      {
        return *std::ranges::min_element(legs, {}, [&pairing, side](std::size_t leg) -> double
                                         { return PerUnit(pairing.Legs(side)[leg]); });
      };
      pairing.Work->Hits.push_back({.Legs = {slowest(0, ones), slowest(1, others)},
                                    .Hit  = {.At     = made.At,
                                             .First  = {.At = made.On[0], .Along = made.Along[0]},
                                             .Second = {.At = made.On[1], .Along = made.Along[1]},
                                             .Kind   = ending ? Contact::End : Contact::Cross}});
      return true;
    }

    /**
     *  The hits at every joint, the joints at one point taken together.
     *
     *  @param  pairing The curves being met; the joints are in its room, and the hits go there.
     *  @return         Whether every one was found so.
     */
    auto Jointed(const Pairing &pairing) -> bool
    {
      std::vector<Joint> &joints = pairing.Work->Joints;
      std::ranges::sort(
          joints, [](const Joint &one, const Joint &other) -> bool
          { return one.At.X < other.At.X || (!(other.At.X < one.At.X) && one.At.Y < other.At.Y); });
      for (auto from = joints.begin(); from != joints.end();)
      {
        const auto until = std::find_if(from, joints.end(), [from](const Joint &joint) -> bool
                                        { return joint.At != from->At; });
        if (!AtJoint(pairing, {from, until}))
        {
          return false;
        }
        from = until;
      }
      return true;
    }

    /**
     *  Whether no point of a curve comes within a reach of a point: weighed whole for a curve that
     *  is one leg, by its legs near the point for a polyline.
     *
     *  @param  pairing The curves being met, whose grids are filed when first asked for.
     *  @param  side    Which curve.
     *  @param  probe   The point and the reach.
     *  @return         Whether none does.
     */
    auto Alone(const Pairing &pairing, std::size_t side, Probe probe) -> bool
    {
      if (!pairing.Strung(side))
      {
        return Farther(pairing.Legs(side).front(), probe.At, probe.Reach);
      }
      return Lonely(pairing, side, probe);
    }

    /**
     *  How the thorough way would look between two hits to tell them apart: on which curve, and how
     *  far between them along it.
     */
    struct Look
    {
      std::size_t Side = 0;   /**< The curve looked along, 0 the first. */
      double      Part = 0.5; /**< How far from the first hit's place.  */
    };

    /**
     *  Whether the thorough way certainly tells two hits apart by a look between them: the place it
     *  looks at, part of the way from one hit's place on a curve to the other's, lies farther from
     *  the other curve than twice the tolerance it holds hits together by, however far it polished
     *  the hits from where they are exactly. Each place it polished lies within Drift tolerances of
     *  the exact one, so the place looked at moves by at most that over the slower leg's speed, as
     *  fast as the fastest leg it may move along; on a closed polyline the short way round is
     *  taken, and nothing is certain where the two ways are nearly as long.
     *
     *  @param  pairing The curves being met.
     *  @param  pair    The hits, the earlier along the first curve first.
     *  @param  look    Where to look.
     *  @return         Whether it does.
     */
    auto Discerned(const Pairing &pairing, std::array<const Found *, 2> pair, Look look) -> bool
    {
      const std::size_t          side     = look.Side;
      const Curve               &curve    = *pairing.Curves.at(side);
      const std::span<const Leg> legs     = pairing.Legs(side);
      const double               length   = pairing.Tolerance.Length;
      const auto                *polyline = std::get_if<Polyline>(&curve);
      const auto                 along    = [side](const Found &found) -> double
      {
        return side == 0 ? found.Hit.First.Along : found.Hit.Second.Along;
      };
      const double slowest =
          std::min(PerUnit(legs[pair[0]->Legs.at(side)]), PerUnit(legs[pair[1]->Legs.at(side)]));
      const double drift  = Drift * length / slowest;
      const Range  domain = Domain(curve);
      const double whole  = domain.End - domain.Start;
      const bool   closed = polyline != nullptr && polyline->Closed;
      double       step   = along(*pair[1]) - along(*pair[0]);
      if (closed && std::abs(std::abs(step) - (whole / 2)) <= 2.0 * drift)
      {
        return false;
      }
      if (closed && std::abs(step) > whole / 2)
      {
        step -= std::copysign(whole, step);
      }
      double middle = along(*pair[0]) + (step * look.Part);
      middle        = middle > domain.End ? middle - whole : middle;
      middle        = middle < domain.Start ? middle + whole : middle;
      /**
       *  A leg with a crossing on it runs longer than twice the clearance, the drift is then less
       *  than a thirty-second of one, and the place looked at moves along two legs at the most.
       */
      double fastest = PerUnit(legs.front());
      if (polyline != nullptr)
      {
        const auto count = static_cast<double>(legs.size());
        for (const double index : {std::floor(middle - drift), std::floor(middle + drift)})
        {
          const double wrapped = closed ? index - (count * std::floor(index / count))
                                        : std::clamp(index, 0.0, count - 1.0);
          fastest = std::max(fastest, PerUnit(legs[static_cast<std::size_t>(wrapped)]));
        }
      }
      const double reach = (drift * fastest * (1.0 + 0x1p-20)) + (3.0 * length);
      return reach <= Clear * length
          && Alone(pairing, 1 - side, {.At = At(curve, middle), .Reach = reach});
    }

    /**
     *  Whether the thorough way certainly tells two hits apart: by a look at a quarter, half or
     *  three quarters of the way between them along either curve, any one of which settles it, a
     *  curve that is one leg asked first, as it costs least.
     *
     *  @param  pairing The curves being met.
     *  @param  pair    The hits, the earlier along the first curve first.
     *  @return         Whether it does.
     */
    auto Separable(const Pairing &pairing, std::array<const Found *, 2> pair) -> bool
    {
      const std::size_t cheap = pairing.Strung(1) ? 1 : 0;
      for (const double part : {0.5, 0.25, 0.75})
      {
        for (const std::size_t side : {cheap, 1 - cheap})
        {
          if (Discerned(pairing, pair, {.Side = side, .Part = part}))
          {
            return true;
          }
        }
      }
      return false;
    }

    /**
     *  Whether the thorough way would keep every hit apart from every other, as each is here: every
     *  two more than Apart tolerances apart, so that none it polished lies within the tolerance of
     *  another; and each told apart from the hit before it along the first curve, and from the
     *  first hit where the first curve closes on itself, as the thorough way looks to join them.
     *
     *  @param  pairing The curves being met, the hits in its room in order along the first curve.
     *  @return         Whether it would.
     */
    auto Distinct(const Pairing &pairing) -> bool
    {
      const std::span<const Found> found = pairing.Work->Hits;
      const double                 clear = Apart * pairing.Tolerance.Length;
      std::vector<std::size_t>    &order = pairing.Work->Order;
      order.clear();
      for (std::size_t index = 0; index < found.size(); ++index)
      {
        order.push_back(index);
      }
      const auto spot = [found](std::size_t index) -> Point
      {
        return found[index].Hit.At;
      };
      std::ranges::sort(order, {}, [&spot](std::size_t index) -> double { return spot(index).X; });
      for (std::size_t one = 0; one < order.size(); ++one)
      {
        for (std::size_t other = one + 1;
             other < order.size() && spot(order[other]).X - spot(order[one]).X <= clear; ++other)
        {
          if (Distance(spot(order[one]), spot(order[other])) <= clear)
          {
            return false;
          }
        }
      }
      const auto *first = std::get_if<Polyline>(pairing.Curves[0]);
      auto       &asked = pairing.Work->Asked;
      asked.clear();
      for (std::size_t index = 1; index < found.size(); ++index)
      {
        asked.push_back({&found[index - 1], &found[index]});
      }
      for (std::size_t index = 2; first != nullptr && first->Closed && index < found.size();
           ++index)
      {
        asked.push_back({found.data(), &found[index]});
      }
      return std::ranges::all_of(asked, [&pairing](const std::array<const Found *, 2> &pair) -> bool
                                 { return Separable(pairing, pair); });
    }

    /**
     *  Whether a straight curve across a polyline of a few segments is met by the swift kernel: it
     *  meets the two whole as the thorough way does where it is sure, at a fraction of the cost of
     *  meeting the legs; past a few its weighing of every two hits costs more than it saves. Asked
     *  before the legs are laid: where they cannot be, the thorough way's direct ways ask the swift
     *  kernel the same next.
     *
     *  @param  curves    The two curves.
     *  @param  tolerance The tolerance.
     *  @param  into      Where the hits go, where it is sure.
     *  @return           Whether it was.
     */
    auto Quick(const Stance &curves, const Tolerance &tolerance, Intersections &into) -> bool
    {
      const Curve &first  = *curves.Shapes[0];
      const Curve &second = *curves.Shapes[1];
      const auto   lined  = [](const Curve &curve) -> bool
      {
        return std::holds_alternative<Line>(curve)
            || std::holds_alternative<Ray>(curve)
            || std::holds_alternative<Xline>(curve);
      };
      /**
       *  One of the two is a polyline, and a straight curve is none: the other is.
       */
      const Curve &strung = lined(first) ? second : first;
      return (lined(first) || lined(second))
          && Count(std::get<Polyline>(strung)) <= Few
          && Swift::Met({first, second}, tolerance, into);
    }
  } // namespace

  auto Legworked(const Curve &curve) -> std::shared_ptr<const Legwork>
  {
    if (!std::holds_alternative<Polyline>(curve) && !std::holds_alternative<Arc>(curve))
    {
      return nullptr;
    }
    const std::shared_ptr<Legwork> work = std::make_shared<Legwork>();
    work->Sound                         = Legged(curve, work->Legs);
    for (Leg &leg : work->Legs)
    {
      if (leg.Bulged())
      {
        Formed(leg);
        leg.Speed = Swiftness(leg);
      }
    }
    if (work->Sound && work->Legs.size() > Handful)
    {
      Filed(work->Legs, work->Filed);
    }
    return work;
  }

  auto Piecewise(const Stance &curves, const Tolerance &tolerance, Intersections &into) -> bool
  {
    const Curve &first  = *curves.Shapes[0];
    const Curve &second = *curves.Shapes[1];
    if (!std::holds_alternative<Polyline>(first) && !std::holds_alternative<Polyline>(second))
    {
      return false;
    }
    /**
     *  A slight arc's center lies too far for a leg's circle: left to the thorough way.
     */
    if (Met::Remotely(first, second))
    {
      return false;
    }
    if (Quick(curves, tolerance, into))
    {
      return true;
    }
    thread_local Room room;
    for (Stalls &grid : room.Grids)
    {
      grid.Laid = false;
    }
    room.Hits.clear();
    room.Joints.clear();
    Pairing pairing = {.Curves = {&first, &second}, .Tolerance = tolerance, .Work = &room};
    for (const std::size_t side : {0UZ, 1UZ})
    {
      const Legwork *kept = curves.Legs.at(side);
      if (kept != nullptr ? !kept->Sound : !Legged(*curves.Shapes.at(side), room.Legs.at(side)))
      {
        return false;
      }
      pairing.Strides.at(side) = kept != nullptr ? std::span<const Leg>{kept->Legs}
                                                 : std::span<const Leg>{room.Legs.at(side)};
      pairing.Kept.at(side)    = kept;
    }
    if (!Paired(pairing))
    {
      return false;
    }
    for (const std::array<std::size_t, 2> &pair : room.Pairs)
    {
      if (!Intersected(pairing, pair))
      {
        return false;
      }
    }
    if (!Jointed(pairing))
    {
      return false;
    }
    std::ranges::sort(room.Hits, {},
                      [](const Found &one) -> double { return one.Hit.First.Along; });
    if (!Distinct(pairing))
    {
      return false;
    }
    for (const Found &one : room.Hits)
    {
      into.Hits.push_back(one.Hit);
    }
    return true;
  }
} // namespace Geoxy
