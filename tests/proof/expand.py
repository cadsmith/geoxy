"""Bounds by expansion, for steps whose error Z3 does not settle in time: a quotient, a square root.

A value worked out in floating point is written exactly in the relative errors of its roundings:
each rounding of an exact z clear of underflow gives z (1 + e), |e| <= u, as rounding.v's
rounded_relative proves, and a product that underflows gives z (1 + e) + h a, |a| <= 1. Divided
by the exact result's own scale, the computed value less the exact one is a polynomial in those
errors and in the inputs' relative parts, each in a box. Its magnitude is at most the sum, over
its terms, of each coefficient's magnitude times the bounds of its variables: the triangle
inequality, term by term. The coefficients are exact rationals, so nothing here rounds.

A bound so proved holds for every value the errors may take, as Z3's would. A mutant is refuted
by one point of the boxes, worked in exact rationals, where the error is past the bound it asks:
so the bound the proof gives is not loose past what the code needs.
"""

import fractions
import itertools
import random

Fraction = fractions.Fraction
U = Fraction(1, 2 ** 53)


class Poly:
    """A polynomial with exact rational coefficients, each term a product of named variables."""

    def __init__(self, terms=None):
        self.terms = {key: value for key, value in (terms or {}).items() if value != 0}

    @classmethod
    def var(cls, name):
        return cls({((name, 1),): Fraction(1)})

    @classmethod
    def const(cls, value):
        return cls({(): Fraction(value)})

    @staticmethod
    def lift(other):
        return other if isinstance(other, Poly) else Poly.const(other)

    def __add__(self, other):
        other = Poly.lift(other)
        terms = dict(self.terms)
        for key, value in other.terms.items():
            terms[key] = terms.get(key, Fraction(0)) + value
        return Poly(terms)

    __radd__ = __add__

    def __neg__(self):
        return Poly({key: -value for key, value in self.terms.items()})

    def __sub__(self, other):
        return self + (-Poly.lift(other))

    def __rsub__(self, other):
        return Poly.lift(other) - self

    def __mul__(self, other):
        other = Poly.lift(other)
        terms = {}
        for (left, first), (right, second) in itertools.product(self.terms.items(),
                                                                 other.terms.items()):
            powers = dict(left)
            for name, power in right:
                powers[name] = powers.get(name, 0) + power
            key = tuple(sorted(powers.items()))
            terms[key] = terms.get(key, Fraction(0)) + first * second
        return Poly(terms)

    __rmul__ = __mul__

    def bound(self, bounds):
        """The most the polynomial's magnitude may be, each variable within its bound."""
        total = Fraction(0)
        for key, value in self.terms.items():
            term = abs(value)
            for name, power in key:
                term *= bounds[name] ** power
            total += term
        return total

    def value(self, values):
        """The polynomial at a point, exactly."""
        total = Fraction(0)
        for key, coefficient in self.terms.items():
            term = coefficient
            for name, power in key:
                term *= values[name] ** power
            total += term
        return total

    def parts(self, names):
        """The polynomial as a part free of the named variables, and each named variable's
        coefficient; no term may hold two of them, nor one twice."""
        free = {}
        shares = {name: {} for name in names}
        for key, value in self.terms.items():
            held = [(name, power) for name, power in key if name in names]
            if not held:
                free[key] = value
                continue
            if len(held) > 1 or held[0][1] > 1:
                raise ValueError("a term holds two of the named variables")
            rest = tuple(pair for pair in key if pair[0] != held[0][0])
            shares[held[0][0]][rest] = value
        return Poly(free), {name: Poly(terms) for name, terms in shares.items()}


def rounding(name):
    """The factor 1 + e of a rounding, e within u."""
    return 1 + Poly.var(name)


def errors(*names):
    """Bounds of u for the named relative errors."""
    return {name: U for name in names}


def verdict(name, held, told):
    """A check's result, as lemma.py's run prints it."""

    def check():
        return held, told

    check.__name__ = name
    return check


def shares(poly, names, bounds):
    """The most each named variable's terms may be, per unit of it, and the most the rest may be:
    a value homogeneous in its scales, each term holding one scale or one underflow, is at most
    the sum of each share times its variable's magnitude, and the rest."""
    free, parts = poly.parts(names)
    return {name: part.bound(bounds) for name, part in parts.items()}, free.bound(bounds)


def peak(poly, bounds, fixed, tries=256):
    """The largest magnitude the polynomial takes at some corners of its boxes, each variable at
    plus or minus its bound, the fixed ones as given: a point every bound must cover, and so the
    refutation of a bound below it."""
    chooser = random.Random(1)
    names = [name for name in bounds if name not in fixed]
    best = Fraction(0)
    for _ in range(tries):
        values = {name: chooser.choice((-1, 1)) * bounds[name] for name in names}
        values.update(fixed)
        best = max(best, abs(poly.value(values)))
    return best
