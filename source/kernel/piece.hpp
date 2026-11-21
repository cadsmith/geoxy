/***************************************************************************************************
 *  @file      source/kernel/piece.hpp
 *  @brief     A rational Bézier piece: what every curved question is answered on.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#ifndef GEOXY_KERNEL_PIECE_HPP
#define GEOXY_KERNEL_PIECE_HPP

#include "basis.hpp"
#include "scratch.hpp"

#include "geoxy/bezier/cubic.hpp"
#include "geoxy/bezier/quadratic.hpp"
#include "geoxy/box.hpp"
#include "geoxy/domain.hpp"
#include "geoxy/matrix.hpp"
#include "geoxy/point.hpp"
#include "geoxy/spline.hpp"
#include "geoxy/vector.hpp"

#include <array>
#include <cstddef>
#include <cstdint>
#include <span>
#include <vector>

namespace Geoxy::Kernel
{
  /**
   *  A rational Bézier piece, its control points homogeneous, their weights positive.
   *
   *  It lies in the hull of its control points, so its box is within theirs, and halving it is one
   *  de Casteljau: every curved question is a recursion over pieces with its own test to prune by.
   *  Held whole on the stack, so no question allocates.
   */
  struct Piece
  {
    std::size_t                            Degree = 0; /**< The degree, at most MaxDegree. */
    std::array<Homogeneous, MaxDegree + 1> Poles =
        {}; /**< The control points, the first Degree + 1. */
  };

  /**
   *  A piece's control points.
   *
   *  @param  piece The piece.
   *  @return       Its first Degree + 1 control points.
   */
  inline auto Poles(const Piece &piece) -> std::span<const Homogeneous>
  {
    return std::span{piece.Poles}.first(piece.Degree + 1);
  }

  /**
   *  Where a piece starts.
   *
   *  @param  piece The piece.
   *  @return       Its first control point, projected.
   */
  inline auto Start(const Piece &piece) -> Point
  {
    return Projected(Poles(piece).front());
  }

  /**
   *  Where a piece ends.
   *
   *  @param  piece The piece.
   *  @return       Its last control point, projected.
   */
  inline auto End(const Piece &piece) -> Point
  {
    return Projected(Poles(piece).back());
  }

  /**
   *  A piece cut in two at a parameter, by de Casteljau's algorithm.
   *
   *  @param  piece The piece.
   *  @param  along Where, from 0 to 1.
   *  @return       The part before and the part after.
   */
  auto Split(const Piece &piece, double along) -> std::array<Piece, 2>;

  /**
   *  A point of a piece and how fast it moves there.
   */
  struct Sample
  {
    Point  At  = {}; /**< The point.                       */
    Vector Way = {}; /**< Its derivative by the parameter. */
  };

  /**
   *  A piece's point and derivative at a parameter, by one pass of de Casteljau's algorithm: the
   *  point from the last two control points but one round, the derivative from their difference
   *  times the degree, by the quotient rule where the weights differ.
   *
   *  @param  piece The piece, of degree 1 at least.
   *  @param  along Where, from 0 to 1.
   *  @return       The point and the derivative.
   */
  auto Sampled(const Piece &piece, double along) -> Sample;

  /**
   *  A piece's point at a parameter, by de Casteljau's algorithm, and nothing of how fast it moves:
   *  a piece cut from beside a heavy weight has weights whose square is nothing to a double, which
   *  the derivative would divide by.
   *
   *  @param  piece The piece.
   *  @param  along Where, from 0 to 1.
   *  @return       The point.
   */
  auto Reached(const Piece &piece, double along) -> Point;

  /**
   *  The box of a piece's control points: a box that holds the piece.
   *
   *  @param  piece The piece.
   *  @return       The box.
   */
  auto Hull(const Piece &piece) -> Box;

  /**
   *  Whether a piece is within a distance of the chord between its ends: every control point is,
   *  and the piece lies in their hull.
   *
   *  @param  piece     The piece.
   *  @param  tolerance The distance.
   *  @return           Whether it is.
   */
  auto Flat(const Piece &piece, double tolerance) -> bool;

  /**
   *  How many pieces a question may halve one into: past it, it answers from the pieces it has.
   */
  inline constexpr std::size_t Budget = 1U << 16U;

  /**
   *  A box grown to hold a piece, within a few steps of a double.
   *
   *  A piece whose control points reach past the box is halved, the point between the halves added,
   *  until what a piece's control points add is within four steps of a double of their coordinates;
   *  those are then added whole. Past the budget, every piece left adds its control points whole,
   *  so the box still holds the piece, a little larger.
   *
   *  @param  piece  The piece, its ends already in the box.
   *  @param  box    The box.
   *  @param  budget How many halvings it may take.
   *  @return        The box grown.
   */
  auto Bound(const Piece &piece, Box box, std::size_t budget) -> Box;

  /**
   *  An axis of the plane.
   */
  enum class Axis : std::uint8_t
  {
    X, /**< The x axis. */
    Y, /**< The y axis. */
  };

  /**
   *  A line across an axis: the points whose coordinate along it is a value.
   */
  struct Level
  {
    Axis   Along = Axis::X; /**< The axis.                */
    double Value = 0.0;     /**< The coordinate along it. */
  };

  /**
   *  A coordinate along an axis.
   *
   *  @tparam PairT What has an x and a y: a point, a vector or a homogeneous point.
   *  @param  pair  It.
   *  @param  axis  The axis.
   *  @return       Its x or its y.
   */
  template <typename PairT>
  auto Coordinate(const PairT &pair, Axis axis) -> double
  {
    return axis == Axis::X ? pair.X : pair.Y;
  }

  /**
   *  The most coefficients of a polynomial the kernel solves: a foot's on a piece of the greatest
   *  degree.
   */
  inline constexpr std::size_t Terms = 3 * static_cast<std::size_t>(MaxDegree);

  /**
   *  The roots of a polynomial in Bernstein form over 0 to 1.
   *
   *  Halved while its coefficients change sign more than once, it holds no root where they keep
   *  their sign and one where they change it once, which bisection then finds to the last bit.
   *
   *  @param  values Its coefficients, Terms at most, one at least.
   *  @param  budget How many halvings it may take; once spent, a part whose ends differ in sign
   *                 gives one root, and any other none.
   *  @return        Where it changes sign and where it is zero at an end, in order; a root where it
   *                 keeps its sign may or may not be among them.
   */
  auto Roots(std::span<const double> values, std::size_t budget) -> List<double>;

  /**
   *  What solving on a piece found: its roots, and the function's values at its ends, whose signs
   *  tell a root that rounding put just past an end from one that is not there.
   */
  struct Solution
  {
    List<double> Roots = {};  /**< The roots, from 0 to 1, in order. */
    double       Start = 0.0; /**< The function's value at 0.        */
    double       End   = 0.0; /**< Its value at 1.                   */
  };

  /**
   *  A place on a chain of pieces: which piece, and its parameter there.
   */
  struct Spot
  {
    std::size_t Index = 0;   /**< Which piece.         */
    double      Along = 0.0; /**< Its parameter there. */
  };

  /**
   *  The roots of a chain of pieces, and each seam the function changes sign across.
   *
   *  A root on a seam leaves the function zero there in truth, but each piece works it out on its
   *  own, and rounding may leave it a little off zero and of opposite signs on the two sides:
   *  neither piece then finds it, and the seam has it. Where the chain has a corner, the seam has
   *  the turn of the function there. A piece with no roots and both end values zero, one of no
   *  length or one the function is zero all along, is looked across: the pieces before and after it
   *  meet as at a seam, at the start of the one after.
   *
   *  @param  solutions What each piece found, in the chain's order.
   *  @param  closed    Whether the last piece's end is the first's start.
   *  @return           The roots and seams, a seam as the start of the piece after it, in the
   *                    chain's order but for a closed chain's seam at its start, which comes last.
   */
  auto Chained(std::span<const Solution> solutions, bool closed) -> List<Spot>;

  /**
   *  Where a piece meets a level.
   *
   *  The piece's coordinate less the level, times its weight, is a polynomial in Bernstein form.
   *
   *  @param  piece  The piece.
   *  @param  level  The level.
   *  @param  budget How many halvings the roots may take.
   *  @return        Where it crosses the level and where an end of it lies on it, in order; a touch
   *                 within may or may not be among them.
   */
  auto Roots(const Piece &piece, Level level, std::size_t budget) -> Solution;

  /**
   *  Where a piece's tangent is square to the way to a point: the feet of the perpendiculars from
   *  the point.
   *
   *  The piece less the point, dotted with its derivative, is a rational function; its numerator, W
   *  (U U' + V V') - W' (U^2 + V^2) for U, V and W the piece's homogeneous coordinates less the
   *  point's times the weight, is a polynomial of degree 3 n - 1, whose roots are found in
   *  Bernstein form; with one positive weight throughout, U U' + V V' alone, of degree 2 n - 1.
   *
   *  @param  piece  The piece.
   *  @param  point  The point.
   *  @param  budget How many halvings the roots may take.
   *  @return        The parameters, from 0 to 1, in order, and the numerator's values at the ends,
   *                 positive where the piece runs away from the point; none, and both values 0,
   *                 when every coefficient is within its rounding of 0, every point of the piece as
   *                 far from the point.
   */
  auto Perpendiculars(const Piece &piece, Point point, std::size_t budget) -> Solution;

  /**
   *  The place of a piece nearest a point: the nearest of its feet and its ends.
   *
   *  @param  piece The piece.
   *  @param  point The point.
   *  @return       The parameter there, from 0 to 1.
   */
  auto Closest(const Piece &piece, Point point) -> double;

  /**
   *  Where a piece's tangent passes through a point: the places it touches the lines from the
   *  point.
   *
   *  The piece less the point, crossed with its derivative, is a rational function whose numerator,
   *  the weight dropping out, is U V' - V U' for U and V the piece's homogeneous coordinates less
   *  the point's times the weight: a polynomial of degree 2 n - 1, whose roots are found in
   *  Bernstein form.
   *
   *  @param  piece  The piece.
   *  @param  point  The point.
   *  @param  budget How many halvings the roots may take.
   *  @return        The parameters, from 0 to 1, in order, and the numerator's values at the ends,
   *                 positive where the piece turns counterclockwise about the point; none, and both
   *                 values 0, when every coefficient is within its rounding of 0, the piece
   *                 straight through the point. Where the piece passes through the point, the place
   *                 there is a double root, and may or may not be among them.
   */
  auto Touches(const Piece &piece, Point point, std::size_t budget) -> Solution;

  /**
   *  The numerator of a piece's derivative: W X' - W' X and W Y' - W' Y for X, Y and W the piece's
   *  homogeneous coordinates, polynomials of degree 2 n - 2 in Bernstein form. Whatever point the
   *  piece is seen from, they are the same.
   */
  struct Hodograph
  {
    std::size_t               Degree = 0;  /**< Their degree, 2 n - 2.         */
    std::array<double, Terms> Across = {}; /**< Along x, the first Degree + 1. */
    std::array<double, Terms> Rising = {}; /**< Along y, the first Degree + 1. */
  };

  /**
   *  A piece with the numerator of its derivative carried beside it.
   *
   *  A part of a piece far smaller than its coordinates has control points a step of a double or
   *  less apart, and which way it runs cannot be read from them. Cut from the whole piece's
   *  numerator, the part's own keeps every digit of it.
   */
  struct Tread
  {
    Piece     Shape = {}; /**< The piece.                                             */
    Hodograph Speed = {}; /**< Its derivative's numerator, times any positive number. */
  };

  /**
   *  The numerator of a piece's derivative, worked out pair by pair of its control points, so that
   *  a pair of light weights keeps its digits beside a pair of heavy ones.
   *
   *  @pre          The degree is 1 at least.
   *  @param  piece The piece.
   *  @return       The numerator.
   */
  auto Derived(const Piece &piece) -> Hodograph;

  /**
   *  The numerator of the derivative of a spline's span, worked out from the spline's own poles.
   *
   *  A span's Bézier control points are sums of its poles, and beside a heavy pole the sum rounds
   *  the light ones away: two control points then differ by nothing, and the span's piece does not
   *  tell which way it runs there. Each pair of control points is instead every pair of poles, by
   *  how much of the two the control points hold, and the poles' differences are exact.
   *
   *  @pre           The spline is valid and the span is not empty.
   *  @param  spline The spline.
   *  @param  span   The index of the knot that starts the span.
   *  @return        The numerator over the span, times a positive number.
   */
  auto Derived(const Spline &spline, std::size_t span) -> Hodograph;

  /**
   *  Where a spline's span has its tangent pass through a point, worked out from the spline's own
   *  poles: each pair of them the two legs from the point crossed, weighed and carried to the
   *  span's control points as Derived carries theirs. Of degree 2 n - 2.
   *
   *  @pre           The spline is valid and the span is not empty.
   *  @param  spline The spline.
   *  @param  span   The index of the knot that starts the span.
   *  @param  point  The point.
   *  @return        As Touches of a piece gives them, the span's parameter from 0 to 1.
   */
  auto Touches(const Spline &spline, std::size_t span, Point point) -> Solution;

  /**
   *  The numerator of a piece's derivative over each half of the piece, times a positive number.
   *
   *  @param  speed The numerator over the whole.
   *  @return       It over the first half and over the second.
   */
  auto Halved(const Hodograph &speed) -> std::array<Hodograph, 2>;

  /**
   *  Where a piece's tangent is square to the way to a point, its derivative's numerator given: the
   *  roots of U times its x part and V times its y part, of degree 3 n - 2.
   *
   *  @param  tread  The piece and its derivative's numerator.
   *  @param  point  The point.
   *  @param  budget How many halvings the roots may take.
   *  @return        As Perpendiculars of a piece gives them.
   */
  auto Perpendiculars(const Tread &tread, Point point, std::size_t budget) -> Solution;

  /**
   *  Where a piece's tangent passes through a point, its derivative's numerator given: the roots of
   *  U times its y part less V times its x part, of degree 3 n - 2.
   *
   *  @param  tread  The piece and its derivative's numerator.
   *  @param  point  The point.
   *  @param  budget How many halvings the roots may take.
   *  @return        As Touches of a piece gives them.
   */
  auto Touches(const Tread &tread, Point point, std::size_t budget) -> Solution;

  /**
   *  A line through a point along a way: the carrier of a straight run, zero where the cross
   *  product of a point less its point with its way is.
   */
  struct Carrier
  {
    Point  Through = {}; /**< A point on it. */
    Vector Way     = {}; /**< Its direction. */
  };

  /**
   *  What solving a curve's implicit equation on a piece found: where it is zero, with its values
   *  at the piece's ends for the seams of a chain, and where it turns, every tangency among them.
   */
  struct Zeros
  {
    Solution     Crossing = {};    /**< Where the equation is zero.               */
    List<double> Turns    = {};    /**< Where its numerator's derivative is zero. */
    bool         On       = false; /**< Whether all the piece lies on the curve.  */
  };

  /**
   *  Where a piece meets a line: the roots of the cross product of the piece less the line's point,
   *  times the weight, with its way, a polynomial of the piece's degree.
   *
   *  @param  piece   The piece.
   *  @param  carrier The line.
   *  @param  budget  How many halvings the roots may take.
   *  @return         The zeros and turns, from 0 to 1; none, and On, when the piece lies on the
   *                  line, to within rounding.
   */
  auto Zeroed(const Piece &piece, const Carrier &carrier, std::size_t budget) -> Zeros;

  /**
   *  Where a piece meets an ellipse: the roots of U^2 + V^2 - W^2, U and V the piece's homogeneous
   *  coordinates mapped so the ellipse is the unit circle, a polynomial of twice the piece's
   *  degree.
   *
   *  @param  piece  The piece.
   *  @param  unit   The map that takes the ellipse to the unit circle.
   *  @param  budget How many halvings the roots may take.
   *  @return        The zeros and turns, from 0 to 1; none, and On, when the piece lies on the
   *                 ellipse, to within rounding.
   */
  auto Zeroed(const Piece &piece, const Matrix &unit, std::size_t budget) -> Zeros;

  /**
   *  What a stretch of a curve encloses with a point: the signed area of the fan from the point and
   *  its first moment about the point, as Green's theorem gives them from the curve alone.
   */
  struct Enclosure
  {
    double Area   = 0.0; /**< The signed area, counterclockwise positive.    */
    Vector Moment = {};  /**< The first moment, the area times its centroid. */
  };

  /**
   *  What a polynomial piece encloses with a point.
   *
   *  By Green's theorem the area is a half of the integral of the piece less the point crossed with
   *  its derivative, and the moment a third of the piece less the point times that; for a
   *  polynomial piece both are polynomials in Bernstein form, whose integral over 0 to 1 is the
   *  mean of their coefficients: exact but for rounding.
   *
   *  @pre          The piece's weights are all one value.
   *  @param  piece The piece.
   *  @param  point The point, best near the piece so that its coordinates are small.
   *  @return       The area and the moment, taken about the point.
   */
  auto Enclosed(const Piece &piece, Point point) -> Enclosure;

  /**
   *  A quarter of an ellipse as a rational quadratic piece.
   *
   *  The quarter of the unit circle between index and index + 1 quarter turns has the corners of
   *  the unit square about it for control points, weighted 1, the square root of a half and 1; the
   *  ellipse's map takes them to its own.
   *
   *  @param  map   The ellipse's map of the unit circle.
   *  @param  index Which quarter, from 0 to 3.
   *  @return       The piece.
   */
  auto Quarter(const Matrix &map, std::size_t index) -> Piece;

  /**
   *  The angle on the unit circle of a quarter piece's parameter.
   *
   *  The piece runs over its quarter as tan(angle / 2) = u sin(q / 2) / (1 - u + u cos(q / 2)), q a
   *  quarter turn, its angle measured from the quarter's start.
   *
   *  @param  index Which quarter.
   *  @param  along The piece's parameter.
   *  @return       The angle, in radians; at 1, the next quarter's start to the bit.
   */
  auto Quartered(std::size_t index, double along) -> double;

  /**
   *  A spline's parameter at a parameter of one of its pieces.
   *
   *  @param  span  The span the piece is of.
   *  @param  along The piece's parameter.
   *  @return       The spline's, the span's end exactly at 1.
   */
  auto Across(Range span, double along) -> double;

  /**
   *  A line as a piece of degree 1.
   *
   *  @param  start Where it starts.
   *  @param  end   Where it ends.
   *  @return       The piece.
   */
  auto Straight(Point start, Point end) -> Piece;

  /**
   *  A polynomial Bézier curve's control points as a piece, every weight 1.
   *
   *  @param  points The control points, at most MaxDegree + 1.
   *  @return        The piece.
   */
  auto Polynomial(std::span<const Point> points) -> Piece;

  /**
   *  A quadratic Bézier curve as a piece.
   *
   *  @param  quadratic The curve.
   *  @return           The piece, of degree 2.
   */
  auto Polynomial(const Bezier::Quadratic &quadratic) -> Piece;

  /**
   *  A cubic Bézier curve as a piece.
   *
   *  @param  cubic The curve.
   *  @return       The piece, of degree 3.
   */
  auto Polynomial(const Bezier::Cubic &cubic) -> Piece;

  /**
   *  A spline's span as a piece: its Bézier control points, each the span's blossom at its ends,
   *  the start repeated as often as the end is not.
   *
   *  @pre           The spline is valid and the span is not empty.
   *  @param  spline The spline.
   *  @param  span   The index of the knot that starts the span.
   *  @return        The piece.
   */
  auto Section(const Spline &spline, std::size_t span) -> Piece;

  /**
   *  Each of a spline's pieces, one per span that is not empty, in order, and the stretch of the
   *  spline's parameter it runs over. A spline that is not valid has none: it is a curve of
   *  nothing, as a polyline with no vertices is.
   *
   *  @tparam VisitT What takes each piece and its stretch.
   *  @param  spline The spline.
   *  @param  visit  What takes them.
   */
  template <typename VisitT>
  void Sections(const Spline &spline, const VisitT &visit)
  {
    if (!Valid(spline))
    {
      return;
    }
    const auto degree = static_cast<std::size_t>(spline.Degree);
    for (std::size_t span = degree; span < spline.Poles.size(); ++span)
    {
      if (spline.Knots[span] < spline.Knots[span + 1])
      {
        visit(Section(spline, span),
              Range{.Start = spline.Knots[span], .End = spline.Knots[span + 1]});
      }
    }
  }

  /**
   *  What a solver finds on each of a spline's pieces, chained, in the spline's parameter: nothing
   *  on a spline that is not valid.
   *
   *  @tparam SolveT What solves on a piece, giving a Solution.
   *  @param  spline The spline.
   *  @param  solve  The solver.
   *  @return        The roots and seams, in order.
   */
  template <typename SolveT>
  auto Solved(const Spline &spline, const SolveT &solve) -> List<double>
  {
    List<Solution> solutions;
    List<Range>    spans;
    Sections(spline,
             [&solutions, &spans, &solve](const Piece &piece, Range span) -> void
             {
               solutions.push_back(solve(piece));
               spans.push_back(span);
             });
    List<double> alongs;
    for (const Spot spot : Chained(solutions, false))
    {
      alongs.push_back(Across(spans[spot.Index], spot.Along));
    }
    return alongs;
  }

  /**
   *  What a solver finds on each of a spline's spans that is not empty, chained, in the spline's
   *  parameter: nothing on a spline that is not valid.
   *
   *  @tparam SolveT What solves on a span, told by the index of the knot that starts it, giving a
   *                 Solution.
   *  @param  spline The spline.
   *  @param  solve  The solver.
   *  @return        The roots and seams, in order.
   */
  template <typename SolveT>
  auto Stepped(const Spline &spline, const SolveT &solve) -> List<double>
  {
    if (!Valid(spline))
    {
      return {};
    }
    List<Solution> solutions;
    List<Range>    spans;
    const auto     degree = static_cast<std::size_t>(spline.Degree);
    for (std::size_t span = degree; span < spline.Poles.size(); ++span)
    {
      if (spline.Knots[span] < spline.Knots[span + 1])
      {
        solutions.push_back(solve(span));
        spans.push_back({.Start = spline.Knots[span], .End = spline.Knots[span + 1]});
      }
    }
    List<double> alongs;
    for (const Spot spot : Chained(solutions, false))
    {
      alongs.push_back(Across(spans[spot.Index], spot.Along));
    }
    return alongs;
  }
} // namespace Geoxy::Kernel

#endif /* GEOXY_KERNEL_PIECE_HPP */
