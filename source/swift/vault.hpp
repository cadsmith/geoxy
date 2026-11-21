/***************************************************************************************************
 *  @file      source/swift/vault.hpp
 *  @brief     A quadratic Bézier curve met with a curve: its implicit equation along each piece of
 *             the other.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#ifndef GEOXY_SWIFT_VAULT_HPP
#define GEOXY_SWIFT_VAULT_HPP

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
   *  A quadratic Bézier curve's parabola as the kernel takes it: the points p where Q(p) = L1^2 - 4
   *  L0 L2 is nothing, L_k(p) the barycentric coordinates of p in the control triangle P0, P1, P2
   *  times twice its signed area D = (P1 - P0) x (P2 - P1). About the middle control point, g = p -
   *  P1, they are L0 = e0 x g, L2 = e2 x g and L1 = D - e1 x g for the sides e0 = P2 - P1, e1 = P2
   *  - P0 and e2 = P1 - P0; on the curve at s they are (1 - s)^2 D, 2 s (1 - s) D and s^2 D, so
   *    that s = 1/2 + (e2 - e0) x g / (2 D), affine in the point (tests/proof/arch.py's implicit
   *    and place). Q's gradient there is 2 D times B'(s) turned a quarter, and its second
   *    derivative along a unit direction no more than 2 |e1|^2 + 8 |e0| |e2| (gradient, curving).
   *    Two points of the parabola at s and s' lie |s' - s| |B'((s + s') / 2)| apart, and |B'| is
   *    nowhere less than 2 |D| / |e2 - e0| (chord, least).
   */
  struct Vault
  {
    Point                                Hub   = {}; /**< P1, where L0 and L2 cross. */
    std::array<std::array<Number, 2>, 3> Sides = {}; /**< e0, e1 and e2, exactly.    */
    std::array<Vector, 3>                Ways  = {}; /**< Each side rounded.         */
    std::array<double, 3>                Slips = {}; /**< How far a rounded side's
                                                          coordinate may be off.           */
    std::array<double, 3> Long = {};                 /**< Each side's length, at the most. */
    Number                Area = {};                 /**< D, with its bound.               */
    double                Bend = 0.0;                /**< 2 |e1|^2 + 8 |e0| |e2|, at the
                                                          most.                            */
    double Fast   = 0.0;                             /**< |B'| from 0 to 1, at the most. */
    double Slow   = 0.0;                             /**< |B'| anywhere, at the least.   */
    double Stride = 0.0;                             /**< What the thorough way charges a
                                                          step of its place: Strict's.     */
    Box Hull = {};                                   /**< Its control points' box. */
  };

  /**
   *  The parabola's equation on a spline's span, H(t) = W(t)^2 Q(N(t) / W(t)) = L1^2 - 4 L0 L2 for
   *  the span's numerator N and weight W, each L_k a polynomial of the span's degree n in t whose
   *  Bernstein coefficient at control point i is L_k at w_i (P_i - P1) and w_i: e0 x G_i, e2 x G_i
   *  and D w_i - e1 x G_i for G_i the moved control point. H is held by its coefficients of degree
   *  2 n, c_k the sum over i + j = k of C(n, i) C(n, j) (l1_i l1_j - 4 l0_i l2_j), over C(2 n, k):
   *  each term two products and a difference, the binomial weight's product, the running sum and
   *  the quotient, which ring.py's term, running and finished (its rho^2 nothing) hold within 16 u
   *  of the terms' sizes for every degree to Broadest (its wide), past what the l_k may be off.
   *  Each l_k is two products and a difference from the rounded side and the moved point, within 4
   *  u of the products' sizes, past the side's slip times the point and the side times the point's
   *  slack; l1's weighted area adds its own rounding, D's bound times the weight and D times the
   *  heft (arch.py's line).
   *
   *  A point p within r of the parabola has |Q(p)| no more than |grad Q(p)| r + (|e1|^2 + 4 |e0|
   *  |e2|) r^2 (arch.py's reach), and |grad Q(p)| = |2 L1 grad L1 - 4 (L0 grad L2 + L2 grad L0)| no
   *  more than U = 2 |L1| |e1| + 4 |L0| |e2| + 4 |L2| |e0|, the most of which on the span is no
   *  more than on its control points, each |L_k| at a point of their hull no more than at its
   *  corners: where |H| is above T = reach W_max^2 (U + Bend reach / 2), the span lies farther than
   *  the reach from the parabola.
   *
   *  @param  vault The parabola.
   *  @param  piece The span, moved to P1.
   *  @param  reach The reach.
   *  @return       The equation.
   */
  auto Springing(const Vault &vault, const Piece<Ample> &piece, double reach) -> Equation<Ample>;

  /**
   *  The parabola's equation at a homogeneous point (N, W), L1^2 - 4 L0 L2 of g = N - P1 W and W,
   *  in double-double.
   *
   *  @param  vault The parabola.
   *  @param  point The point's numerator's coordinates and its weight.
   *  @return       H there.
   */
  auto Keystone(const Vault &vault, const std::array<Number, 3> &point) -> Number;

  /**
   *  One of the parabola's equation's roots on a span: polished and its point worked out as Honing
   *  works them out; then kept where it lies on the Bézier curve, decided clear of where the curves
   *  may come within the reach about it, steep, and its point and places settled.
   *
   *  A point of the parabola within near of the crossing has its place within near / |B'|_min of
   *  the crossing's (arch.py's chord and least). The crossing's place is s = 1/2 + (e2 - e0) x (X
   *  - P1) / (2 D) of its point, in double-double, the point's bound carried: where s lies farther
   *    than that margin past 0 or 1 every point of the parabola about the crossing is off the
   *    curve, and where it lies that far inside every one is on it, the curve's ends farther than
   *    the reach from the span. Only a crossing kept need be steep: the sine there is |H'| / (W^2 2
   *    |D| |B'| |p'|) (arch.py's sine), |B'| no more than Fast.
   *
   *  @param  arcade The parabola and the span.
   *  @param  root   The root's stretch.
   *  @param  found  Where the crossing goes, where it is one.
   *  @return        Whether it was decided.
   */
  auto Sprung(const Arcade<Vault> &arcade, const Root &root, List<Ribbed> &found) -> bool;

  /**
   *  Whether a point is farther than a length c from a parabola: |Q(p)| above U(p) c + (|e1|^2 + 4
   *  |e0| |e2|) c^2, U(p) = 2 |L1| |e1| + 4 |L0| |e2| + 4 |L2| |e0| at the point, each L_k in
   *  double-double from the exact gap (arch.py's reach).
   *
   *  @param  vault The parabola.
   *  @param  point The point.
   *  @param  clear The length, c.
   *  @return       Whether it is.
   */
  auto Clearing(const Vault &vault, Point point, double clear) -> bool;

  /**
   *  A quadratic Bézier curve met with a rational spline or a polynomial one past a cubic, of
   *  degree to eight, by the ring kernel on the curve's parabola: L1^2 - 4 L0 L2 = 0 for the
   *  barycentric coordinates L_k of its control triangle, a quadratic form in a span's homogeneous
   *  point, of twice the span's degree on it, its roots isolated in doubles and each polished once
   *  in double-double. The curve's place at a point of the parabola is affine in the point, so a
   *  crossing is kept where its place lies within the curve, clear of its ends by the margin about
   *  the crossing. It is sure only where the thorough way would find these crossings and nothing
   *  else: wherever the curves come within eight tolerances of each other, a steep crossing there,
   *  far from either curve's ends and from any other. Each crossing is the double nearest the true
   *  point, both its places the doubles nearest the true ones.
   *
   *  @param  curves    The two curves.
   *  @param  tolerance The tolerance, no finer than the curves' grain.
   *  @param  into      Where the hits go, empty on entry, and left so where it is not sure.
   *  @return           Whether it was sure; nothing where the curves are not such a pair.
   */
  auto Vaulted(const Curves &curves, const Tolerance &tolerance, Intersections &into)
      -> std::optional<bool>;
} // namespace Geoxy::Swift

#endif /* GEOXY_SWIFT_VAULT_HPP */
