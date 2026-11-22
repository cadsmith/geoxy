/***************************************************************************************************
 *  @file      source/swift/gather.cpp
 *  @brief     The hits of a ring and a curve gathered piece by piece, and told apart.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "swift/gather.hpp"

#include "crossing.hpp"
#include "scratch.hpp"
#include "swift/circuit.hpp"
#include "swift/lane.hpp"
#include "swift/piece.hpp"
#include "swift/pieced.hpp"
#include "swift/rim.hpp"
#include "swift/ring.hpp"
#include "swift/section.hpp"
#include "swift/terms.hpp"
#include "trig.hpp"

#include "geoxy/arc.hpp"
#include "geoxy/at.hpp"
#include "geoxy/bezier/cubic.hpp"
#include "geoxy/bezier/quadratic.hpp"
#include "geoxy/box.hpp"
#include "geoxy/curve.hpp"
#include "geoxy/domain.hpp"
#include "geoxy/ellipse.hpp"
#include "geoxy/intersect.hpp"
#include "geoxy/matrix.hpp"
#include "geoxy/place.hpp"
#include "geoxy/point.hpp"
#include "geoxy/spline.hpp"
#include "geoxy/tolerance.hpp"
#include "geoxy/vector.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <numbers>
#include <optional>
#include <span>
#include <utility>
#include <variant>

namespace Geoxy::Swift
{
  namespace
  {
    /**
     *  Whether a quarter of an arc or an ellipse lies outside its sweep, clear of it by a margin
     *  past what rounding the angles may do: then no point of the curve is on it, and it is not
     *  met. Its angles run from q pi / 2 to (q + 1) pi / 2, the sweep's over |sweep| from its least
     *  angle, both taken round to a turn from the quarter's start.
     *
     *  @param  ring    The arc or the ellipse.
     *  @param  quarter Which quarter.
     *  @return         Whether it lies outside.
     */
    auto Beyond(const Ring &ring, std::size_t quarter) -> bool
    {
      constexpr double Whole = 2.0 * std::numbers::pi;
      constexpr double Verge = 0x1p-20;
      const double     least = std::min(ring.Start, ring.Start + ring.Sweep);
      const double     shift = least - (static_cast<double>(quarter) * std::numbers::pi / 2.0);
      const double     from  = shift - (Whole * std::floor(shift / Whole));
      return !ring.Whole
          && from > (std::numbers::pi / 2.0) + Verge
          && from + std::abs(ring.Sweep) < Whole - Verge;
    }

    /**
     *  A Bézier curve as a polynomial piece.
     *
     *  @param  curve The curve.
     *  @return       Its section; nothing for any other curve.
     */
    auto Bezier(const Curve &curve) -> std::optional<Section>
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
     *  What the hits of a ring and a piece's curve are told apart by.
     */
    struct Parting
    {
      const Ring  *Carrier = nullptr; /**< The ring whose equation.                       */
      const Curve *Shape   = nullptr; /**< The piece's curve.                             */
      std::size_t  Side    = 0;       /**< Which of a hit's places is on it: 0 the first. */
      double       Length  = 0.0;     /**< The tolerance's length.                        */
    };

    /**
     *  Whether no two hits are one, as the thorough way merges them: two are one only where they
     *  lie within the tolerance of each other or where, halfway between them along each curve, each
     *  curve is within twice the tolerance of the other. Each pair lies farther apart than the
     *  tolerance, and the piece's curve halfway between them farther than twice it from the ring's
     *  carrier, both with the margin: halfway the plain way, and the short way round where the
     *  curve may close on itself.
     *
     *  @param  parting The ring, the piece's curve, and the tolerance's length.
     *  @param  hits    The hits.
     *  @return         Whether none are one.
     */
    auto RingApart(const Parting &parting, std::span<const Hit> hits) -> bool
    {
      const double clear = 2.0 * parting.Length * Margin;
      const auto   place = [&parting](const Hit &hit) -> double
      {
        return parting.Side == 0 ? hit.First.Along : hit.Second.Along;
      };
      const auto far = [&parting, clear](double along) -> bool
      {
        return RingGap(*parting.Carrier, At(*parting.Shape, along)) > clear;
      };
      const Range  domain = Domain(*parting.Shape);
      const double whole  = domain.End - domain.Start;
      for (std::size_t one = 0; one < hits.size(); ++one)
      {
        for (std::size_t two = one + 1; two < hits.size(); ++two)
        {
          const double from   = place(hits[one]);
          const double middle = from + ((place(hits[two]) - from) / 2.0);
          const double round  = middle + (whole / 2.0);
          if (!(Distance(hits[one].At, hits[two].At) > parting.Length * Margin)
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
     *  A crossing made a hit, placed as the thorough way places it: on an arc or an ellipse by its
     *  point's angle, on a Bézier curve at the place the kernel proved.
     *
     *  @param  curves The two curves.
     *  @param  side   Which is the piece's curve: 0 for the first.
     *  @param  found  The crossing.
     *  @return        The hit; nothing where the point is not placed on a curve.
     */
    auto Hitting(const Curves &curves, std::size_t side, const Found &found) -> std::optional<Hit>
    {
      const auto place = [&curves, side, &found](std::size_t which) -> std::optional<Place>
      {
        const Curve &curve = curves.at(which);
        return which == side && found.Own
                 ? std::optional<Place>{Place{.At = At(curve, found.Along), .Along = found.Along}}
                 : Crossings::Located(curve, 0.5, found.At);
      };
      return place(0).and_then(
          [&place, &found](const Place &first) -> std::optional<Hit>
          {
            return place(1).transform(
                [&first, &found](const Place &second) -> Hit
                {
                  return {.At = found.At, .First = first, .Second = second, .Kind = Contact::Cross};
                });
          });
    }

    /**
     *  The crossings a ring's kernel gathers, and the tolerance it meets within.
     */
    struct Gathering
    {
      Geoxy::Tolerance          Tolerance = {}; /**< The tolerance.       */
      List<Geoxy::Swift::Found> Found     = {}; /**< The crossings found. */
    };

    /**
     *  A ring met with a polynomial piece, its crossings gathered.
     *
     *  @tparam CountT    The piece's room.
     *  @param  carrier   The ring.
     *  @param  section   The piece.
     *  @param  gathering The tolerance, and where the crossings go.
     *  @return           Whether it was sure.
     */
    template <std::size_t CountT>
    auto Laid(const Ring &carrier, const Sectional<CountT> &section, Gathering &gathering) -> bool
    {
      const std::optional<Piece<CountT>> piece = Plain(section, carrier);
      return piece
          && Encountered(Scene<CountT>{.Carrier   = &carrier,
                                       .Part      = &*piece,
                                       .Tolerance = gathering.Tolerance},
                         gathering.Found);
    }

    /**
     *  A box holding a bulge's arc: every point of it lies within |c| / 2 max(1, |b|) of the
     *  chord's middle, c the chord, as the arc's middle lies |b| |c| / 2 from it, and past a half
     *  turn the circle's far side r + |c| |1 - b^2| / (4 |b|) = |c| |b| / 2. Grown past what
     *  rounding takes.
     *
     *  @param  hump The bulge.
     *  @param  pad  How far past the arc it reaches.
     *  @return      The box.
     */
    auto Boxing(const Hump &hump, double pad) -> Box
    {
      const Vector chord  = hump.End - hump.Start;
      const Point  middle = hump.Start + (chord * 0.5);
      const double half =
          (Length(chord) * 0.5 * std::max(1.0, std::abs(hump.Bend)) * (1.0 + 0x1p-40))
          + (0x1p-40 * (std::abs(middle.X) + std::abs(middle.Y)))
          + pad;
      return {.Min = middle - Vector{.X = half, .Y = half},
              .Max = middle + Vector{.X = half, .Y = half}};
    }

    /**
     *  A ring met with a spline span by span: a span whose poles' box lies farther than twice the
     *  reach from the ring's is far from it; each span crossed must be one about whose crossings
     *  the thorough way holds its candidates to the tolerance alone, a step of the ring's place,
     *  from 0 to 1, moving it by no more than |M|_F times its sweep times a step of the doubles at
     *  1.
     *
     *  @param  carrier   The ring.
     *  @param  spline    The spline.
     *  @param  gathering The tolerance, and where the crossings go.
     *  @return           Whether it was sure.
     */
    auto Threaded(const Ring &carrier, const Spline &spline, Gathering &gathering) -> bool
    {
      const Matrix &map    = carrier.Map;
      const Vector  extent = {.X = Trig::Hypot(map.X.X, map.Y.X),
                              .Y = Trig::Hypot(map.X.Y, map.Y.Y)};
      const double  clear =
          (2.0 * Reach(gathering.Tolerance))
          + (0x1p-50 * (std::abs(map.Origin.X) + std::abs(map.Origin.Y) + extent.X + extent.Y))
          + carrier.Drift;
      const Vector reach = extent + Vector{.X = clear, .Y = clear};
      /**
       *  A bulge's spans are met about its arc's box alone: one farther from it lies farther from
       *  the arc, and its crossings with the rest of the circle are not the polyline's.
       */
      const Box    near   = carrier.Chord != nullptr
                              ? Boxing(*carrier.Chord, clear)
                              : Box{.Min = map.Origin - reach, .Max = map.Origin + reach};
      const double stride = carrier.Chord != nullptr
                              ? carrier.Stride
                              : carrier.Most * std::abs(carrier.Sweep) * 0x1p-52 * Grow;
      const std::optional<List<Broad>> sections = Weighed(spline, near);
      if (!sections)
      {
        return false;
      }
      for (const Broad &section : *sections)
      {
        const std::size_t before = gathering.Found.size();
        if (!Laid(carrier, section, gathering)
            || (gathering.Found.size() > before && !Fine(section, stride, gathering.Tolerance)))
        {
          return false;
        }
      }
      return true;
    }

    /**
     *  Where a bulge's circle meets pieces of the other curve: its arc's box, grown by twice the
     *  reach and the roundings of its center and the box, as Threaded keeps a spline's spans. A
     *  piece whose box lies outside it lies farther than the reach from the arc, and its crossings
     *  with the rest of the circle are not the polyline's.
     *
     *  @param  carrier   The ring, a bulge's circle.
     *  @param  tolerance The tolerance.
     *  @return           The box.
     */
    auto Bulging(const Ring &carrier, const Tolerance &tolerance) -> Box
    {
      const double clear = (2.0 * Reach(tolerance))
                         + (0x1p-50
                            * (std::abs(carrier.Map.Origin.X)
                               + std::abs(carrier.Map.Origin.Y)
                               + (2.0 * carrier.High)))
                         + carrier.Drift;
      return Boxing(*carrier.Chord, clear);
    }

    /**
     *  Whether two boxes lie apart, neither reaching the other along x or along y.
     *
     *  @param  one   One box.
     *  @param  other The other.
     *  @return       Whether they do.
     */
    auto Disjoint(const Box &one, const Box &other) -> bool
    {
      return one.Min.X   > other.Max.X
          || other.Min.X > one.Max.X
          || one.Min.Y   > other.Max.Y
          || other.Min.Y > one.Max.Y;
    }

    /**
     *  A box a quarter of an arc or an ellipse keeps in: its three control points', as a rational
     *  quadratic piece of weights above zero lies in their hull. Each is o' + M c, c's coordinates
     *  0 or 1 turned, worked out in two sums within u of the sizes of their terms, |X_x| + |Y_x| no
     *  more than twice |M|_F: the box grown by 2^-48 (|o'_x| + |o'_y| + |M|_F) holds the true
     *  points.
     *
     *  @param  bent    The arc or the ellipse.
     *  @param  quarter Which quarter.
     *  @return         The box.
     */
    auto Kerbed(const Ring &bent, std::size_t quarter) -> Box
    {
      constexpr std::array<std::array<double, 2>, 3> Corners = {
          {{1.0, 0.0}, {1.0, 1.0}, {0.0, 1.0}}};
      const Matrix &map  = bent.Map;
      Box           hull = {};
      for (const std::array<double, 2> &corner : Corners)
      {
        const auto turned = Quartered(corner, quarter);
        hull              = Grown(hull, map.Origin + (map.X * turned[0]) + (map.Y * turned[1]));
      }
      const double pad = 0x1p-48 * (std::abs(map.Origin.X) + std::abs(map.Origin.Y) + bent.Most);
      return {.Min = hull.Min - Vector{.X = pad, .Y = pad},
              .Max = hull.Max + Vector{.X = pad, .Y = pad}};
    }
  } // namespace

  auto Conic(const Curve &curve) -> bool
  {
    return std::holds_alternative<Arc>(curve) || std::holds_alternative<Ellipse>(curve);
  }

  auto Polynomial(const Curve &curve) -> bool
  {
    return std::holds_alternative<Bezier::Quadratic>(curve)
        || std::holds_alternative<Bezier::Cubic>(curve)
        || std::holds_alternative<Spline>(curve);
  }

  auto Gathered(const Ring &carrier, const Curve &other, const Tolerance &tolerance)
      -> std::optional<List<Found>>
  {
    Gathering                gathering = {.Tolerance = tolerance, .Found = {}};
    const std::optional<Box> near =
        carrier.Chord != nullptr ? std::optional{Bulging(carrier, tolerance)} : std::nullopt;
    if (const std::optional<Section> section = Bezier(other))
    {
      if (near && Disjoint(section->Hull, *near))
      {
        return std::move(gathering.Found);
      }
      return Laid(carrier, *section, gathering) ? std::optional{std::move(gathering.Found)}
                                                : std::nullopt;
    }
    if (const auto *spline = std::get_if<Spline>(&other))
    {
      return Threaded(carrier, *spline, gathering) ? std::optional{std::move(gathering.Found)}
                                                   : std::nullopt;
    }
    const std::optional<Ring> bent = Rimmed(other);
    if (!bent)
    {
      return std::nullopt;
    }
    const Bend bend = Bendings(*bent, carrier);
    for (std::size_t quarter = 0; quarter < 4; ++quarter)
    {
      if (Beyond(*bent, quarter) || (near && Disjoint(Kerbed(*bent, quarter), *near)))
      {
        continue;
      }
      const Piece<Thin> piece = QuarterPiece(bend, quarter);
      if (!Encountered(Scene<Thin>{.Carrier = &carrier, .Part = &piece, .Tolerance = tolerance},
                       gathering.Found))
      {
        return std::nullopt;
      }
    }
    return std::move(gathering.Found);
  }

  auto Rung(const Ring &carrier, const RingPair &pairing, Intersections &into) -> bool
  {
    const Curves                    &curves = *pairing.Pair;
    const std::optional<List<Found>> found =
        Gathered(carrier, curves.at(pairing.Side), pairing.Tolerance);
    const std::span<const Found> crossings =
        found ? std::span<const Found>{*found} : std::span<const Found>{};
    const std::size_t wanted = found ? found->size() : std::numeric_limits<std::size_t>::max();
    std::size_t       placed = 0;
    for (const Found &crossing : crossings)
    {
      const std::optional<Hit> hit  = Hitting(curves, pairing.Side, crossing);
      placed                       += static_cast<std::size_t>(hit.has_value());
      into.Hits.push_back(hit.value_or(Hit{}));
    }
    return placed == wanted
        && RingApart({.Carrier = &carrier,
                      .Shape   = &curves.at(pairing.Side).get(),
                      .Side    = pairing.Side,
                      .Length  = pairing.Tolerance.Length},
                     into.Hits);
  }
} // namespace Geoxy::Swift
