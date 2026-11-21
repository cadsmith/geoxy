/***************************************************************************************************
 *  @file      source/swift/arcade.hpp
 *  @brief     What the kernels of a quadratic and a cubic Bézier curve share: a carrier's equation
 *             along a piece, its roots honed, and the hits gathered.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#ifndef GEOXY_SWIFT_ARCADE_HPP
#define GEOXY_SWIFT_ARCADE_HPP

#include "scratch.hpp"
#include "sure.hpp"
#include "swift/circuit.hpp"
#include "swift/gather.hpp"
#include "swift/lane.hpp"
#include "swift/piece.hpp"
#include "swift/pieced.hpp"
#include "swift/roots.hpp"
#include "swift/section.hpp"
#include "swift/terms.hpp"
#include "twofold.hpp"

#include "geoxy/at.hpp"
#include "geoxy/box.hpp"
#include "geoxy/curve.hpp"
#include "geoxy/domain.hpp"
#include "geoxy/intersect.hpp"
#include "geoxy/place.hpp"
#include "geoxy/point.hpp"
#include "geoxy/spline.hpp"
#include "geoxy/tolerance.hpp"
#include "geoxy/vector.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <optional>
#include <span>
#include <utility>

namespace Geoxy::Swift
{
  /**
   *  A Bézier curve's equation at a grid place of a span, in double-double: the span's true point,
   *  as Exactly or, where rational, Homogeneously gives it, and the parabola's or the cubic's
   *  equation there, as Keystone works it out.
   *
   *  @tparam CarrierT The Bézier curve's kind: a parabola's or a cubic's.
   *  @param  carrier  The Bézier curve.
   *  @param  piece    The span.
   *  @param  along    The place, on the grid.
   *  @return          H there, and the point's numerator's coordinates and its weight.
   */
  template <typename CarrierT>
  auto Spandrel(const CarrierT &carrier, const Piece<Ample> &piece, double along)
      -> std::pair<Number, std::array<Number, 3>>
  {
    std::array<Number, 3> point = {};
    if (piece.Source->Rated)
    {
      point = Homogeneously(*piece.Source, Sure::Exact(along));
    }
    else
    {
      const std::array<Number, 2> exact = Exactly(*piece.Source, along);
      point                             = {exact[0], exact[1], Sure::Exact(1.0)};
    }
    return {Keystone(carrier, point), point};
  }

  /**
   *  A Bézier curve and a span being met, and what is worked out once about them.
   *
   *  @tparam CarrierT The Bézier curve's kind: a parabola's or a cubic's.
   */
  template <typename CarrierT>
  struct Arcade
  {
    const CarrierT               *Carrier   = nullptr; /**< The Bézier curve.        */
    const Piece<Ample>           *Part      = nullptr; /**< The span it is met with. */
    Geoxy::Swift::Equation<Ample> Equation  = {};      /**< Its equation on it.      */
    Geoxy::Tolerance              Tolerance = {};      /**< The tolerance.           */
  };

  /**
   *  A root of a Bézier curve's equation on a span, polished once in double-double, and the
   *  crossing's point.
   */
  struct Honed
  {
    Number                Place    = {};  /**< The place on the span, with its bound. */
    std::array<Number, 2> Crossing = {};  /**< The point, with its bounds.            */
    double                Near     = 0.0; /**< How far from the crossing a point of the
                                               curve may lie where it comes within the
                                               reach of the span about it.             */
  };

  /**
   *  A root of a Bézier curve's equation on a span polished and its point worked out as a ring's on
   *  a span: the root t* lies within R = (|H(t0)| + E) / g of the grid place t0, g the least |H'|
   *  on the stretch, and the step t1 = t0 - H0 / H'(t0) within [E + 2u |H0| + F (slack + |H''| R) /
   *  g] / |H'(t0)| of it, F the most |H(t0)| may be. Where the span comes within the reach of the
   *  curve, |H| is no more than T, so within T / g of a root on its stretch, and its point within
   *  speed T / g of that root's: a point of the curve that near lies within near = reach + speed T
   *  / g of the crossing.
   *
   *  @tparam CarrierT The Bézier curve's kind: a parabola's or a cubic's.
   *  @param  arcade   The curve and the span.
   *  @param  root     The root's stretch.
   *  @return          The place, the point and how near.
   */
  template <typename CarrierT>
  auto Honing(const Arcade<CarrierT> &arcade, const Root &root) -> Honed
  {
    const Piece<Ample>    &piece    = *arcade.Part;
    const Equation<Ample> &equation = arcade.Equation;
    const auto [along, slope]       = Burnished(equation, root);
    const auto [value, point]       = Spandrel(*arcade.Carrier, piece, along);
    const double slack              = Sloping(equation.Degree, equation.First, equation.Largest);
    const double most               = (Sure::Size(value) + value.Bound) * Grow;
    const double far                = most / root.Slope * Grow;
    const double step               = (value.Value.High + value.Value.Low) / slope;
    const double off                = ((((value.Bound + (2.0 * Unit * Sure::Size(value)))
                                         + (most * (slack + (equation.Turn * far)) / root.Slope))
                                        / std::abs(slope))
                                       + (2.0 * Unit * std::abs(step)))
                                    * Grow;
    return {.Place    = {.Value = Twofold::Add(Twofold::Lifted(along), Twofold::Lifted(-step)),
                         .Bound = off},
            .Crossing = Moved(piece, point, {along, step, off}),
            .Near = (Reach(arcade.Tolerance) + (piece.Speed * equation.Reach / root.Slope * Grow))
                  * Grow};
  }

  /**
   *  A crossing of a Bézier curve and a span: its point, the double nearest the true one, and its
   *  places, each the double nearest the true one.
   */
  struct Ribbed
  {
    Point  At    = {};  /**< The point.                   */
    double Along = 0.0; /**< The spline's place.          */
    double Bent  = 0.0; /**< The Bézier curve's place, s. */
  };

  /**
   *  The crossings a Bézier curve's kernel gathers, and the tolerance it meets within.
   */
  struct Springs
  {
    Geoxy::Tolerance Tolerance = {}; /**< The tolerance.       */
    List<Ribbed>     Found     = {}; /**< The crossings found. */
  };

  /**
   *  A Bézier curve met with a spline span by span: a span whose poles' box lies farther than twice
   *  the reach from the Bézier curve's control points' box is far from the curve, which keeps to
   *  their hull; each span's roots isolated as a ring's are, every one decided; each span crossed
   *  must be one about whose crossings the thorough way holds its candidates to the tolerance
   *  alone, a step of the curve's place moving it no more than its stride.
   *
   *  @tparam CarrierT The Bézier curve's kind: a parabola's or a cubic's.
   *  @param  carrier  The Bézier curve.
   *  @param  spline   The spline.
   *  @param  springs  The tolerance, and where the crossings go.
   *  @return          Whether it was sure.
   */
  template <typename CarrierT>
  auto Coursing(const CarrierT &carrier, const Spline &spline, Springs &springs) -> bool
  {
    const Box                       &hull  = carrier.Hull;
    const double                     reach = Reach(springs.Tolerance);
    const double                     clear = (2.0 * reach)
                                           + (0x1p-50
                                              * (std::abs(hull.Min.X)
                                                 + std::abs(hull.Min.Y)
                                                 + std::abs(hull.Max.X)
                                                 + std::abs(hull.Max.Y)));
    const Vector                     pad   = {.X = clear, .Y = clear};
    const std::optional<List<Broad>> sections =
        Weighed(spline, {.Min = hull.Min - pad, .Max = hull.Max + pad});
    if (!sections)
    {
      return false;
    }
    for (const Broad &section : *sections)
    {
      std::optional<Piece<Ample>> piece = Ready(section);
      if (!piece)
      {
        return false;
      }
      Shifting(section, carrier.Hub, *piece);
      const std::optional<Equation<Ample>>     equation = Springing(carrier, *piece, reach);
      const std::optional<Roots<Ranks<Ample>>> roots =
          equation ? Isolates(*equation) : std::nullopt;
      if (!roots)
      {
        return false;
      }
      const Arcade<CarrierT> arcade = {.Carrier   = &carrier,
                                       .Part      = &*piece,
                                       .Equation  = *equation,
                                       .Tolerance = springs.Tolerance};
      const std::size_t      before = springs.Found.size();
      if (!std::ranges::all_of(std::span{roots->At}.first(roots->Count),
                               [&arcade, &springs](const Root &root) -> bool
                               { return Sprung(arcade, root, springs.Found); })
          || (springs.Found.size() > before && !Fine(section, carrier.Stride, springs.Tolerance)))
      {
        return false;
      }
    }
    return true;
  }

  /**
   *  Whether no two hits of a Bézier curve and a spline are one, as the thorough way merges them:
   *  each pair lies farther apart than the tolerance, and the spline halfway between them, the
   *  plain way and the short way round, farther than twice it from the parabola or the cubic, both
   *  with the margin, as RingApart weighs a ring's.
   *
   *  @tparam CarrierT The Bézier curve's kind: a parabola's or a cubic's.
   *  @param  carrier  The Bézier curve.
   *  @param  pairing  The two curves, which is the spline, and the tolerance.
   *  @param  hits     The hits.
   *  @return          Whether none are one.
   */
  template <typename CarrierT>
  auto Unparted(const CarrierT &carrier, const RingPair &pairing, std::span<const Hit> hits) -> bool
  {
    const Curve &shape  = pairing.Pair->at(pairing.Side);
    const double length = pairing.Tolerance.Length;
    const double clear  = 2.0 * length * Margin;
    const Range  domain = Domain(shape);
    const double whole  = domain.End - domain.Start;
    const auto   place  = [&pairing](const Hit &hit) -> double
    {
      return pairing.Side == 0 ? hit.First.Along : hit.Second.Along;
    };
    const auto far = [&carrier, &shape, clear](double along) -> bool
    {
      return Clearing(carrier, At(shape, along), clear);
    };
    for (std::size_t one = 0; one < hits.size(); ++one)
    {
      for (std::size_t two = one + 1; two < hits.size(); ++two)
      {
        const double from   = place(hits[one]);
        const double middle = from + ((place(hits[two]) - from) / 2.0);
        const double round  = middle + (whole / 2.0);
        if (!(Distance(hits[one].At, hits[two].At) > length * Margin)
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
   *  A Bézier curve met with a spline by its parabola's or its cubic's kernel: the spline's spans
   *  coursed, each crossing made a hit at the places proved, each the curve's own and the double
   *  nearest the true one, and no two hits one.
   *
   *  @pre             The pairing's other curve is a spline.
   *  @tparam CarrierT The Bézier curve's kind: a parabola's or a cubic's.
   *  @param  carrier  The Bézier curve.
   *  @param  pairing  The two curves, which is the spline, and the tolerance.
   *  @param  into     Where the hits go, empty on entry, and left so where it is not sure.
   *  @return          Whether it was sure.
   */
  template <typename CarrierT>
  auto Trussed(const CarrierT &carrier, const RingPair &pairing, Intersections &into) -> bool
  {
    const Curves &curves  = *pairing.Pair;
    Springs       springs = {.Tolerance = pairing.Tolerance, .Found = {}};
    if (!Coursing(carrier, std::get<Spline>(curves.at(pairing.Side).get()), springs))
    {
      return false;
    }
    for (const Ribbed &ribbed : springs.Found)
    {
      const std::array<double, 2> places = pairing.Side == 1
                                             ? std::array{ribbed.Bent, ribbed.Along}
                                             : std::array{ribbed.Along, ribbed.Bent};
      into.Hits.push_back({.At     = ribbed.At,
                           .First  = {.At = At(curves[0], places[0]), .Along = places[0]},
                           .Second = {.At = At(curves[1], places[1]), .Along = places[1]},
                           .Kind   = Contact::Cross});
    }
    if (!Unparted(carrier, pairing, into.Hits))
    {
      into.Hits.clear();
      return false;
    }
    std::ranges::sort(into.Hits, {}, [](const Hit &hit) -> double { return hit.First.Along; });
    return true;
  }
} // namespace Geoxy::Swift

#endif /* GEOXY_SWIFT_ARCADE_HPP */
