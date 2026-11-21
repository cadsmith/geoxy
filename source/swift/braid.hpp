/***************************************************************************************************
 *  @file      source/swift/braid.hpp
 *  @brief     A cubic Bézier curve met with a curve: its implicit equation, a Bézout determinant,
 *             along each piece of the other.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#ifndef GEOXY_SWIFT_BRAID_HPP
#define GEOXY_SWIFT_BRAID_HPP

#include "scratch.hpp"
#include "swift/arcade.hpp"
#include "swift/circuit.hpp"
#include "swift/piece.hpp"
#include "swift/roots.hpp"
#include "swift/terms.hpp"

#include "geoxy/box.hpp"
#include "geoxy/intersect.hpp"
#include "geoxy/point.hpp"
#include "geoxy/tolerance.hpp"
#include "geoxy/vector.hpp"

#include <array>
#include <cstddef>
#include <optional>

namespace Geoxy::Swift
{
  /**
   *  How many entries a cubic's Bézout matrix holds: a, b, c, d, e and g of the symmetric [[a, b,
   *  c], [b, d, e], [c, e, g]].
   */
  inline constexpr std::size_t Fibres = 6;

  /**
   *  A cubic Bézier curve's implicit cubic, its Bézout form about the first control point P0
   *  (Sederberg's): with g = p - P0, a_i = P_i - P0, A_ij = a_i x a_j and l_ij = (P_i - p) x (P_j
   *  - p) = A_ij - (P_i - P_j) x g, every point of the curve's cubic, at every place s, has f = det
   *    [[a, b, c], [b, d, e], [c, e, g]] nothing for a = 3 l32, b = 3 l31, c = l30, d = l30 + 9
   *    l21, e = 3 l20 and g = 3 l10, each entry a linear form K - S x g (tests/proof/braid.py's
   *    implicit). There (s^2, s (1 - s), (1 - s)^2) is the matrix's kernel; the column of c, e and
   *    g, three forms that vanish at P0, is killed by lambda = (3 A21, A13, A32), so that N = 3 A21
   *    a + A13 b + A32 c and D = 3 A21 b + A13 d + A32 e, linear forms too, give the place of every
   *    point of it but its double point: s (D - N) = D (inverse).
   */
  struct Braid
  {
    Point                                     Hub   = {}; /**< P0.             */
    std::array<std::array<Number, 2>, Fibres> Sides = {}; /**< Each entry's S. */
    std::array<Number, Fibres>                Areas = {}; /**< Each entry's K. */
    std::array<Vector, Fibres>                Ways  = {}; /**< Each S rounded. */
    std::array<double, Fibres>                Slips = {}; /**< How far a rounded S's
                                                               coordinate may be off.           */
    std::array<double, Fibres>           Long   = {};     /**< Each |S|, at the most. */
    std::array<double, Fibres>           Levels = {};     /**< Each K rounded.        */
    std::array<double, Fibres>           Drifts = {};     /**< How far it may be off. */
    std::array<std::array<Number, 2>, 2> Turns  = {};     /**< N's and D's S.         */
    std::array<Number, 2>                Bases  = {};     /**< N's and D's K.         */
    double                               Stride = 0.0;    /**< What the thorough way
                                                               charges a step of its place:
                                                               Strict's.                        */
    Box Hull = {};                                        /**< Its control points' box. */
  };

  /**
   *  The cubic's equation on a spline's span, H(t) = W(t)^3 f(N(t) / W(t)) = a (d g - e^2) - b (b g
   *  - c e) + c (b e - c d) of the entries on the span, of degree 3 n, held by its Bernstein
   *    coefficients: each minor's as Twilled has them, and h_k the sum over i + j = k of C(n, i)
   *    (a_i m0_j - b_i m1_j + c_i m2_j) over C(3 n, k). Each term three products, a difference, a
   *    sum and the binomial weight's product, the running sum of no more than n + 1 terms and the
   *    quotient are within 16 u of the terms' sizes for every degree to Tightest (braid.py's
   *    summed, term and wide), past what the entries and the minors may be off times each other
   *    (ring.py's product). Every binomial weight is whole and exact.
   *
   *  A point p within r of the cubic, at q, has |f(p)| = |f(p) - f(q)| no more than U1 r + U2 r^2 +
   *  U3 r^3 as Teased gives them from the entries' most sizes on the span (braid.py's monomial):
   *  where |H| is above T = reach W_max^3 (U1 + U2 reach + U3 reach^2), the span lies farther than
   *  the reach from the cubic. A T or a bound past the doubles, which only weights far apart could
   *  make, leaves no stretch far and no root's: nothing is isolated.
   *
   *  @param  braid The cubic.
   *  @param  piece The span, moved to P0.
   *  @param  reach The reach.
   *  @return       The equation; nothing where an entry may be past Thickest.
   */
  auto Springing(const Braid &braid, const Piece<Ample> &piece, double reach)
      -> std::optional<Equation<Ample>>;

  /**
   *  The cubic's equation at a homogeneous point (N, W), W^3 f(N / W), in double-double.
   *
   *  @param  braid The cubic.
   *  @param  point The point's numerator's coordinates and its weight.
   *  @return       H there.
   */
  auto Keystone(const Braid &braid, const std::array<Number, 3> &point) -> Number;

  /**
   *  One of the cubic's equation's roots on a span: polished and its point worked out as Honing
   *  works them out; then kept where it lies on the Bézier curve, decided clear of where the curves
   *  may come within the reach about it, steep, and its point and places settled.
   *
   *  Every point q of the cubic within near of the crossing p but its double point has its place
   *  D(q) / (D(q) - N(q)) (braid.py's inverse), and each of D and D - N, affine, moves by no more
   *  than its side's length times near from its value at p: so a quotient whose bounds hold those
   *  moves holds every such place. Where it lies farther than its bound past 0 or 1 every point of
   *  the cubic about the crossing is off the curve, and where it lies that far inside every one is
   *  on it, the curve's ends farther than the reach from the span; where D - N may be nothing, the
   *  double point may be near, and it is not sure. Only a crossing kept need be steep: the sine
   *  there is |H'| / (W^3 |grad f| |p'|) (braid.py's sine), |grad f| no more than Pitched gives at
   *  the entries' sizes at the crossing.
   *
   *  @param  arcade The cubic and the span.
   *  @param  root   The root's stretch.
   *  @param  found  Where the crossing goes, where it is one.
   *  @return        Whether it was decided.
   */
  auto Sprung(const Arcade<Braid> &arcade, const Root &root, List<Ribbed> &found) -> bool;

  /**
   *  Whether a point is farther than a length c from a cubic: |f(p)| above U1 c + U2 c^2 + U3 c^3
   *  as Teased gives them from the entries' sizes at the point, each entry and f in double-double
   *  from the exact gap (braid.py's monomial).
   *
   *  @param  braid The cubic.
   *  @param  point The point.
   *  @param  clear The length, c.
   *  @return       Whether it is.
   */
  auto Clearing(const Braid &braid, Point point, double clear) -> bool;

  /**
   *  A cubic Bézier curve met with a rational spline or a polynomial one past a cubic, of degree to
   *  five, by the ring kernel on the curve's implicit cubic: the determinant of its Bézout matrix,
   *  whose entries are linear forms, a cubic form in a span's homogeneous point, of three times the
   *  span's degree on it, its roots isolated in doubles and each polished once in double-double.
   *  The curve's place at a point of its cubic is a quotient of two linear forms in the point, so a
   *  crossing is kept where its place lies within the curve, clear of its ends by the margin about
   *  the crossing. It is sure only where the thorough way would find these crossings and nothing
   *  else: wherever the curves come within eight tolerances of each other, a steep crossing there,
   *  far from either curve's ends and from any other. Each crossing is the double nearest the true
   *  point, both its places the doubles nearest the true ones.
   *
   *  @param  curves    The two curves.
   *  @param  tolerance The tolerance, no finer than the curves' grain.
   *  @param  into      Where the hits go, empty on entry, and left so where it is not sure.
   *  @return           Whether it was sure; nothing where the curves are not such a pair, a spline
   *                    past degree five aside, of which it is never sure.
   */
  auto Braided(const Curves &curves, const Tolerance &tolerance, Intersections &into)
      -> std::optional<bool>;
} // namespace Geoxy::Swift

#endif /* GEOXY_SWIFT_BRAID_HPP */
