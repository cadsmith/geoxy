/***************************************************************************************************
 *  @file      source/inertia.cpp
 *  @brief     Second moments of what a curve encloses: the polygon of its pieces' chords summed in
 *             closed form, each piece with its chord by Green's theorem in the chord's frame.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "geoxy/inertia.hpp"

#include "legendre.hpp"
#include "real.hpp"
#include "trig.hpp"

#include "geoxy/accuracy.hpp"
#include "geoxy/arc.hpp"
#include "geoxy/at.hpp"
#include "geoxy/curve.hpp"
#include "geoxy/derivative.hpp"
#include "geoxy/domain.hpp"
#include "geoxy/ellipse.hpp"
#include "geoxy/matrix.hpp"
#include "geoxy/point.hpp"
#include "geoxy/polyline.hpp"
#include "geoxy/ray.hpp"
#include "geoxy/segment.hpp"
#include "geoxy/spline.hpp"
#include "geoxy/transform.hpp"
#include "geoxy/vector.hpp"
#include "geoxy/xline.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <functional>
#include <iterator>
#include <numbers>
#include <optional>
#include <span>
#include <type_traits>
#include <variant>
#include <vector>

namespace Geoxy
{
  namespace
  {
    /**
     *  The integrals over an area: of 1, x, y, x^2, y^2 and x y.
     */
    struct Sums
    {
      double Area = 0.0; /**< Of 1.   */
      double X    = 0.0; /**< Of x.   */
      double Y    = 0.0; /**< Of y.   */
      double XX   = 0.0; /**< Of x^2. */
      double YY   = 0.0; /**< Of y^2. */
      double XY   = 0.0; /**< Of x y. */
    };

    /**
     *  Two sums added.
     *
     *  @param  lhs One.
     *  @param  rhs The other.
     *  @return     Their sum.
     */
    auto operator+(const Sums &lhs, const Sums &rhs) -> Sums
    {
      return {.Area = lhs.Area + rhs.Area,
              .X    = lhs.X + rhs.X,
              .Y    = lhs.Y + rhs.Y,
              .XX   = lhs.XX + rhs.XX,
              .YY   = lhs.YY + rhs.YY,
              .XY   = lhs.XY + rhs.XY};
    }

    /**
     *  Sums times a number.
     *
     *  @param  sums   The sums.
     *  @param  factor The number.
     *  @return        Each times it.
     */
    auto operator*(const Sums &sums, double factor) -> Sums
    {
      return {.Area = sums.Area * factor,
              .X    = sums.X * factor,
              .Y    = sums.Y * factor,
              .XX   = sums.XX * factor,
              .YY   = sums.YY * factor,
              .XY   = sums.XY * factor};
    }

    /**
     *  A frame: an origin and the unit way of its x axis, its y axis a quarter turn on.
     */
    struct Frame
    {
      Point  Origin = {};                   /**< The origin. */
      Vector Way    = {.X = 1.0, .Y = 0.0}; /**< The x axis. */
    };

    /**
     *  A piece of a curve to integrate: the curve, moved near the origin, and the stretch of its
     *  parameter; or a polyline's bulge, by its half sweep a and its signed radius r, its points r
     *  (sin t, cos a - cos t) for t from -a to a, in the frame of its chord.
     */
    struct Swath
    {
      const Curve *Shape  = nullptr; /**< The curve; none for a bulge.         */
      Geoxy::Range Range  = {};      /**< The stretch: a bulge's from -a to a. */
      double       Radius = 0.0;     /**< A bulge's signed radius.             */
      double       Half   = 0.0;     /**< A bulge's half sweep.                */
      Frame        Chord  = {};      /**< A bulge's chord: its middle and way. */
    };

    /**
     *  What Green's theorem gives a piece at a parameter in a frame whose x axis is its chord: the
     *  forms -f(x, y) x' whose integrals along the piece are those over the area between it and its
     *  chord, each f with a factor of y, so that the chord adds nothing.
     *
     *  @param  stretch The piece.
     *  @param  frame   The chord's frame.
     *  @param  along   The parameter.
     *  @return         The forms.
     */
    auto Forms(const Swath &stretch, const Frame &frame, double along) -> Sums
    {
      const double radius = stretch.Radius;
      const double half   = stretch.Half;
      double       x      = radius * Trig::Sin(along);
      double y = -2.0 * radius * Trig::Sin((half + along) / 2.0) * Trig::Sin((half - along) / 2.0);
      double run = -radius * Trig::Cos(along);
      if (stretch.Shape != nullptr)
      {
        const Vector offset = At(*stretch.Shape, along) - frame.Origin;
        const Vector across = Perpendicular(frame.Way);
        x                   = Dot(offset, frame.Way);
        y                   = Dot(offset, across);
        run                 = -Dot(Derivative(*stretch.Shape, along), frame.Way);
      }
      return {.Area = y * run,
              .X    = x * y * run,
              .Y    = y * y * run / 2.0,
              .XX   = x * x * y * run,
              .YY   = y * y * y * run / 3.0,
              .XY   = x * y * y * run / 2.0};
    }

    /**
     *  The forms integrated over a stretch by a 16-point Gauss-Legendre rule, and the integral of
     *  their magnitudes.
     *
     *  @param  stretch The piece's stretch.
     *  @param  frame   The chord's frame.
     *  @return         The integrals, and the magnitudes'.
     */
    auto Gaussian(const Swath &stretch, const Frame &frame) -> std::array<Sums, 2>
    {
      const double         middle = (stretch.Range.Start + stretch.Range.End) / 2.0;
      const double         half   = (stretch.Range.End - stretch.Range.Start) / 2.0;
      const Legendre::Rule rule   = Legendre::Ruled(8);
      Sums                 sum    = {};
      Sums                 size   = {};
      for (std::size_t index = 0; index < rule.Nodes.size(); ++index)
      {
        for (const double node : {-rule.Nodes[index], rule.Nodes[index]})
        {
          const Sums   forms = Forms(stretch, frame, middle + (half * node));
          const double share = rule.Weights[index] * half;
          sum                = sum + (forms * share);
          size               = size
                             + (Sums{.Area = std::abs(forms.Area),
                                     .X    = std::abs(forms.X),
                                     .Y    = std::abs(forms.Y),
                                     .XX   = std::abs(forms.XX),
                                     .YY   = std::abs(forms.YY),
                                     .XY   = std::abs(forms.XY)}
                                * std::abs(share));
        }
      }
      return {sum, size};
    }

    /**
     *  How near two rules must agree, for each integral: 2^-44 of its magnitude's, or the accuracy,
     *  times the chord's length to the power the integral's units are past an area.
     *
     *  @param  size     The magnitudes' integrals.
     *  @param  accuracy The accuracy.
     *  @param  chord    The chord's length.
     *  @return          The slack of each.
     */
    auto Agreement(const Sums &size, Accuracy accuracy, double chord) -> Sums
    {
      constexpr double Relative = 0x1p-44;
      const double     error    = std::max(accuracy.Error, 0.0);
      return {.Area = std::max(size.Area * Relative, error),
              .X    = std::max(size.X * Relative, error * chord),
              .Y    = std::max(size.Y * Relative, error * chord),
              .XX   = std::max(size.XX * Relative, error * chord * chord),
              .YY   = std::max(size.YY * Relative, error * chord * chord),
              .XY   = std::max(size.XY * Relative, error * chord * chord)};
    }

    /**
     *  Whether two rules agree within the slack.
     *
     *  @param  whole  The rule over the stretch.
     *  @param  halves The rules over its halves, summed.
     *  @param  slack  The slack.
     *  @return        Whether every integral does.
     */
    auto Agreed(const Sums &whole, const Sums &halves, const Sums &slack) -> bool
    {
      const std::array<double, 6> gaps = {
          std::abs(whole.Area - halves.Area), std::abs(whole.X - halves.X),
          std::abs(whole.Y - halves.Y),       std::abs(whole.XX - halves.XX),
          std::abs(whole.YY - halves.YY),     std::abs(whole.XY - halves.XY)};
      const std::array<double, 6> slacks = {slack.Area, slack.X,  slack.Y,
                                            slack.XX,   slack.YY, slack.XY};
      return std::ranges::equal(gaps, slacks, std::less_equal<>{});
    }

    /**
     *  How the halving goes: the frame, the accuracy and the chord's length, and how many halvings
     *  are left.
     */
    struct Halving
    {
      Geoxy::Frame    Frame    = {};  /**< The chord's frame.          */
      Geoxy::Accuracy Accuracy = {};  /**< The accuracy.               */
      double          Chord    = 0.0; /**< The chord's length.         */
      int             Left     = 0;   /**< How many halvings are left. */
    };

    /**
     *  The forms integrated over a stretch, halved until a rule and the two over its halves agree.
     *
     *  @param  stretch The stretch.
     *  @param  whole   The rule over it.
     *  @param  halving The frame, the accuracy, and how deep it may go.
     *  @return         The integrals.
     */
    auto Integrated(const Swath &stretch, const Sums &whole, const Halving &halving) -> Sums
    {
      const double middle = (stretch.Range.Start + stretch.Range.End) / 2.0;
      Swath        first  = stretch;
      Swath        second = stretch;
      first.Range.End     = middle;
      second.Range.Start  = middle;

      const std::array<Sums, 2> before = Gaussian(first, halving.Frame);
      const std::array<Sums, 2> after  = Gaussian(second, halving.Frame);
      const Sums                halves = before[0] + after[0];
      const Sums slack = Agreement(before[1] + after[1], halving.Accuracy, halving.Chord);
      if (halving.Left == 0 || Agreed(whole, halves, slack))
      {
        return halves;
      }
      const Halving deeper = {.Frame    = halving.Frame,
                              .Accuracy = halving.Accuracy,
                              .Chord    = halving.Chord,
                              .Left     = halving.Left - 1};
      return Integrated(first, before[0], deeper) + Integrated(second, after[0], deeper);
    }

    /**
     *  Sums in a frame moved to sums about a point, along the plane's own axes.
     *
     *  @param  sums   The sums in the frame.
     *  @param  frame  The frame.
     *  @param  origin The point.
     *  @return        The sums about it.
     */
    auto Shifted(const Sums &sums, const Frame &frame, Point origin) -> Sums
    {
      const Vector shift    = frame.Origin - origin;
      const Vector way      = frame.Way;
      const Vector across   = Perpendicular(way);
      const double x        = (sums.X * way.X) + (sums.Y * across.X);
      const double y        = (sums.X * way.Y) + (sums.Y * across.Y);
      const double squaredX = (way.X * way.X * sums.XX)
                            + (2.0 * way.X * across.X * sums.XY)
                            + (across.X * across.X * sums.YY);
      const double squaredY = (way.Y * way.Y * sums.XX)
                            + (2.0 * way.Y * across.Y * sums.XY)
                            + (across.Y * across.Y * sums.YY);
      const double crossed  = (way.X * way.Y * sums.XX)
                            + (((way.X * across.Y) + (across.X * way.Y)) * sums.XY)
                            + (across.X * across.Y * sums.YY);
      return {.Area = sums.Area,
              .X    = (sums.Area * shift.X) + x,
              .Y    = (sums.Area * shift.Y) + y,
              .XX   = (sums.Area * shift.X * shift.X) + (2.0 * shift.X * x) + squaredX,
              .YY   = (sums.Area * shift.Y * shift.Y) + (2.0 * shift.Y * y) + squaredY,
              .XY   = (sums.Area * shift.X * shift.Y) + (shift.X * y) + (shift.Y * x) + crossed};
    }

    /**
     *  What a piece encloses with its chord, about a point.
     *
     *  @param  stretch  The piece.
     *  @param  origin   The point.
     *  @param  accuracy The accuracy.
     *  @return          The sums.
     */
    auto Segmental(const Swath &stretch, Point origin, Accuracy accuracy) -> Sums
    {
      Frame  frame  = stretch.Chord;
      double length = 2.0 * std::abs(stretch.Radius * Trig::Sin(stretch.Half));
      if (stretch.Shape != nullptr)
      {
        const Point  start = At(*stretch.Shape, stretch.Range.Start);
        const Vector chord = At(*stretch.Shape, stretch.Range.End) - start;
        length             = Length(chord);
        frame              = {.Origin = start,
                              .Way    = Zero(length) ? Vector{.X = 1.0, .Y = 0.0} : chord / length};
      }
      const Sums whole = Gaussian(stretch, frame)[0];
      const Sums sums  = Integrated(
          stretch, whole, {.Frame = frame, .Accuracy = accuracy, .Chord = length, .Left = 16});
      return Shifted(sums, frame, origin);
    }

    /**
     *  What a polygon encloses, about its first corner, by its sums in closed form.
     *
     *  @param  corners The corners, in order, the last joined to the first.
     *  @return         The sums about the first.
     */
    auto Cornerwise(std::span<const Point> corners) -> Sums
    {
      Sums sums = {};
      for (std::size_t index = 1; index + 1 < corners.size(); ++index)
      {
        const Vector one   = corners[index] - corners[0];
        const Vector two   = corners[index + 1] - corners[0];
        const double cross = Cross(one, two);
        sums = sums
             + Sums{.Area = cross / 2.0,
                    .X    = cross * (one.X + two.X) / 6.0,
                    .Y    = cross * (one.Y + two.Y) / 6.0,
                    .XX   = cross * ((one.X * one.X) + (one.X * two.X) + (two.X * two.X)) / 12.0,
                    .YY   = cross * ((one.Y * one.Y) + (one.Y * two.Y) + (two.Y * two.Y)) / 12.0,
                    .XY   = cross
                          * ((2.0 * one.X * one.Y)
                             + (one.X * two.Y)
                             + (two.X * one.Y)
                             + (2.0 * two.X * two.Y))
                          / 24.0};
      }
      return sums;
    }

    /**
     *  The pieces of a curve to integrate, each smooth and short enough for the rules: a polyline's
     *  bulges, a spline's spans, a Bézier curve whole, an arc or an ellipse by quarter turns; the
     *  polygon of their ends; all about the curve's first point, which they are moved from.
     */
    struct Cuttings
    {
      Point              Origin  = {}; /**< The curve's first point.           */
      std::vector<Curve> Shapes  = {}; /**< The curve, moved to start at 0, 0. */
      std::vector<Swath> Pieces  = {}; /**< The pieces.                        */
      std::vector<Point> Corners = {}; /**< The polygon of the pieces' ends.   */
    };

    /**
     *  A polyline cut into its bulges and the polygon of its vertices, about its first vertex.
     *
     *  @param polyline The polyline.
     *  @param cut      The cut, filled.
     */
    void Bulges(const Polyline &polyline, Cuttings &cut)
    {
      const std::size_t count = Count(polyline);
      cut.Origin              = At(polyline, 0.0);
      for (std::size_t index = 0; index < count; ++index)
      {
        const Vertex &vertex = polyline.Vertices[index];
        const Point   end    = At(polyline, static_cast<double>(index + 1));
        const Vector  chord  = end - vertex.Position;
        cut.Corners.push_back(Point{} + (vertex.Position - cut.Origin));
        if (std::holds_alternative<Arc>(Nth(polyline, index)))
        {
          const double length = Length(chord);
          const double half   = 2.0 * Trig::Atan(vertex.Bulge);
          cut.Pieces.push_back(
              {.Shape  = nullptr,
               .Range  = {.Start = -half, .End = half},
               .Radius = length / (2.0 * Trig::Sin(half)),
               .Half   = half,
               .Chord  = {.Origin = Point{} + (vertex.Position - cut.Origin) + (chord / 2.0),
                          .Way    = chord / length}});
        }
      }
      cut.Corners.push_back(Point{} + (At(polyline, static_cast<double>(count)) - cut.Origin));
    }

    /**
     *  How many pieces a curve is integrated in: an arc's or an ellipse's quarter turns, one for a
     *  Bézier curve.
     *
     *  @tparam ShapeT The curve's kind.
     *  @param  shape  The curve.
     *  @return        The count.
     */
    template <typename ShapeT>
    auto Quarters(const ShapeT &shape) -> std::size_t
    {
      if constexpr (std::is_same_v<ShapeT, Arc> || std::is_same_v<ShapeT, Ellipse>)
      {
        return static_cast<std::size_t>(
            std::max(1.0, std::ceil(std::abs(shape.Sweep.InRadians()) * 2.0 / std::numbers::pi)));
      }
      else
      {
        return 1;
      }
    }

    /**
     *  A curve other than a polyline cut into pieces, moved to start at the origin: a spline at its
     *  knots, an arc or an ellipse at its quarter turns, a Bézier curve whole.
     *
     *  @tparam ShapeT The curve's kind.
     *  @param  shape  The curve.
     *  @param  curve  The curve, as any curve.
     *  @param  cut    The cut, filled.
     */
    template <typename ShapeT>
    void Sown(const ShapeT &shape, const Curve &curve, Cuttings &cut)
    {
      const Range domain = Domain(shape);
      cut.Origin         = At(shape, domain.Start);
      cut.Shapes.push_back(Transformed(curve, Translation(Point{} - cut.Origin)));
      const Curve        &moved = cut.Shapes.front();
      std::vector<double> cuts;
      if constexpr (std::is_same_v<ShapeT, Spline>)
      {
        std::ranges::copy_if(shape.Knots, std::back_inserter(cuts), [domain](double knot) -> bool
                             { return domain.Start <= knot && knot <= domain.End; });
        cuts.erase(std::ranges::unique(cuts).begin(), cuts.end());
      }
      else
      {
        const std::size_t count = Quarters(shape);
        for (std::size_t index = 0; index <= count; ++index)
        {
          cuts.push_back(static_cast<double>(index) / static_cast<double>(count));
        }
      }
      for (std::size_t index = 0; index + 1 < cuts.size(); ++index)
      {
        cut.Corners.push_back(At(moved, cuts[index]));
        cut.Pieces.push_back(
            {.Shape = &moved, .Range = {.Start = cuts[index], .End = cuts[index + 1]}});
      }
      cut.Corners.push_back(At(moved, cuts.back()));
    }

    /**
     *  A curve cut into pieces, moved to start at the origin.
     *
     *  @param  curve The curve, finite.
     *  @return       The pieces.
     */
    auto Diced(const Curve &curve) -> Cuttings
    {
      Cuttings cut = {};
      std::visit(
          [&cut, &curve](const auto &shape) -> void
          {
            using ShapeT = std::decay_t<decltype(shape)>;
            if constexpr (std::is_same_v<ShapeT, Polyline>)
            {
              Bulges(shape, cut);
            }
            else
            {
              Sown(shape, curve, cut);
            }
          },
          curve);
      return cut;
    }

    /**
     *  Whether a curve can enclose an area: it ends, and a spline is valid.
     *
     *  @param  curve The curve.
     *  @return       Whether it can.
     */
    auto Enclosable(const Curve &curve) -> bool
    {
      return !std::holds_alternative<Ray>(curve)
          && !std::holds_alternative<Xline>(curve)
          && (!std::holds_alternative<Spline>(curve) || Valid(std::get<Spline>(curve)));
    }
  } // namespace

  auto Inertia(const Curve &curve, Accuracy accuracy) -> std::optional<Moments>
  {
    if (!Enclosable(curve))
    {
      return std::nullopt;
    }
    const Cuttings cut  = Diced(curve);
    Sums           sums = Cornerwise(cut.Corners);
    for (const Swath &piece : cut.Pieces)
    {
      sums = sums + Segmental(piece, {}, accuracy);
    }
    if (Zero(sums.Area) || !std::isfinite(sums.Area))
    {
      return std::nullopt;
    }
    const double x = sums.X / sums.Area;
    const double y = sums.Y / sums.Area;
    return Moments{.Area     = sums.Area,
                   .Centroid = cut.Origin + Vector{.X = x, .Y = y},
                   .AboutX   = sums.YY - (sums.Area * y * y),
                   .AboutY   = sums.XX - (sums.Area * x * x),
                   .Product  = sums.XY - (sums.Area * x * y)};
  }

  auto Axes(const Moments &moments) -> Principal
  {
    const double sign    = moments.Area < 0.0 ? -1.0 : 1.0;
    const double alongX  = moments.AboutY * sign;
    const double alongY  = moments.AboutX * sign;
    const double product = moments.Product * sign;
    const double mean    = (alongX + alongY) / 2.0;
    const double spread  = Trig::Hypot((alongX - alongY) / 2.0, product);
    const double angle   = Trig::Atan2(2.0 * product, alongX - alongY) / 2.0;
    return {.Least = mean - spread, .Most = mean + spread, .Axis = Trig::Way(angle)};
  }
} // namespace Geoxy
