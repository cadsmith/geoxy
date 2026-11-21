/***************************************************************************************************
 *  @file      source/corner/surd.hpp
 *  @brief     Numbers with a square root inside a square root, held exactly: the signs they take
 *             against a double, and the double nearest one, ties to even.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#ifndef GEOXY_CORNER_SURD_HPP
#define GEOXY_CORNER_SURD_HPP

#include "expansion.hpp"

/**
 *  What a fillet and a chamfer work out: where a circle of a radius touches two lines or circles,
 *  and the curves cut there.
 */
namespace Geoxy::Corner
{
  /**
   *  A number held exactly, as an expansion, with the arithmetic of one: a sum, a difference and a
   *  product of two never round.
   */
  struct Big
  {
    Expansion::Number Value = {}; /**< The number. */
  };

  /**
   *  The sum of two numbers, exactly.
   *
   *  @param  lhs One.
   *  @param  rhs The other.
   *  @return     The sum.
   */
  inline auto operator+(const Big &lhs, const Big &rhs) -> Big
  {
    return {.Value = Expansion::Sum(lhs.Value, rhs.Value)};
  }

  /**
   *  The difference of two numbers, exactly.
   *
   *  @param  lhs The number.
   *  @param  rhs The number taken from it.
   *  @return     The difference.
   */
  inline auto operator-(const Big &lhs, const Big &rhs) -> Big
  {
    return {.Value = Expansion::Difference(lhs.Value, rhs.Value)};
  }

  /**
   *  The product of two numbers, exactly.
   *
   *  @param  lhs One.
   *  @param  rhs The other.
   *  @return     The product.
   */
  inline auto operator*(const Big &lhs, const Big &rhs) -> Big
  {
    return {.Value = Expansion::Product(lhs.Value, rhs.Value)};
  }

  /**
   *  A number with its sign changed, exactly.
   *
   *  @param  number The number.
   *  @return        Its negative.
   */
  inline auto operator-(const Big &number) -> Big
  {
    return {.Value = Expansion::Scaled(number.Value, -1.0)};
  }

  /**
   *  A coordinate a fillet or a chamfer works out, as the numbers it is made of: (Whole + Inner
   *  sqrt(Root) + Outer sqrt(Base + Bend sqrt(Root))) / Below. A point where a line meets another
   *  takes Whole and Below alone; one where a circle touches a line and a circle, all of them.
   *
   *  @tparam NumberT The numbers: doubles to choose among the answers, sure double-doubles to
   *                  certify one, exact ones where no bound settles it.
   */
  template <typename NumberT>
  struct Form
  {
    NumberT Whole = {}; /**< What is added to the roots' shares.            */
    NumberT Inner = {}; /**< What the root of Root is multiplied by.        */
    NumberT Root  = {}; /**< Under the inner root, not below nothing.       */
    NumberT Outer = {}; /**< What the outer root is multiplied by.          */
    NumberT Base  = {}; /**< Under the outer root, less Bend sqrt(Root).    */
    NumberT Bend  = {}; /**< What sqrt(Root) is multiplied by in the outer. */
    NumberT Below = {}; /**< What all is divided by, not nothing.           */
  };

  /**
   *  The sign of a + b sqrt(p), exactly.
   *
   *  @param  whole a.
   *  @param  share b.
   *  @param  root  p, not below nothing.
   *  @return       1, -1 or 0.
   */
  auto Signed(const Big &whole, const Big &share, const Big &root) -> int;

  /**
   *  The sign of a form less a number, exactly: of its numerator less the number times its
   *  denominator, times the denominator's sign.
   *
   *  @param  form The form, its roots' arguments not below nothing.
   *  @param  mark The number.
   *  @return      1 where the form is above it, -1 where below, 0 where equal.
   */
  auto Signed(const Form<Big> &form, const Expansion::Number &mark) -> int;

  /**
   *  The double nearest a form's value, the even one where two are as near.
   *
   *  The doubles about a guess are searched for the least whose halfway point to the next one up
   *  lies at or above the value, each step one exact comparison: within twice a bound of the guess,
   *  and two places either side, where the value lies within the bound; across every finite double
   *  where it does not, or where the guess is not finite or past 2^1000. The value's products with
   *  halfway points must stay within the doubles, as they do for the forms of curves scaled to near
   *  1.
   *
   *  @param  form  The form.
   *  @param  guess A double near the value.
   *  @param  bound How far from the guess the value may be.
   *  @return       The double nearest.
   */
  auto Pinpointed(const Form<Big> &form, double guess, double bound) -> double;
} // namespace Geoxy::Corner

#endif /* GEOXY_CORNER_SURD_HPP */
