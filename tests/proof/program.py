"""The floating point operations of the exact predicates' proved functions, as the code performs
them: each function a graph of operations on its inputs, written once.

translation.py checks every graph against the compiler's LLVM IR of the code, operation for
operation and constant for constant. lemma.py checks every rounding a Z3 lemma states against a
graph: each lemma binds its variables to the graph's values, and each value it rounds must be
the graph's operation on the graph's operands. So the proofs are about the code the compiler
compiles, not about a copy of it.

The error-free transformations are the code's, two.hpp's, operation for operation, as
rounding.v proves them: TwoSum (two_sum), TwoDiff with TwoDiffTail (two_diff_tail), TwoProduct
with a fused multiply-add (two_product_fma). orient.cpp's CrossSign is rounding.v's cross_sign,
its Safe safe_product.
"""

import struct

from lemma import (CIRCLE_BOUND_A, CIRCLE_BOUND_B, CIRCLE_BOUND_C, CIRCLE_BOUND_COMPENSATED,
                   CIRCLE_BOUND_UNDERFLOW, NORMAL, ORIENT_BOUND_A, RESULT_BOUND)

COMMUTATIVE = {"add", "mul"}


class Value:
    """A value of a graph: an input, a constant, or an operation on values."""

    def __init__(self, op, *args, name=None, constant=None):
        self.op = op
        self.args = args
        self.name = name
        self.constant = constant
        self._key = None

    def key(self):
        """A canonical text of the value: see signed."""
        if self._key is None:
            sign, text = signed(self)
            self._key = ("-" if sign < 0 else "") + text
        return self._key

    def __eq__(self, other):
        return isinstance(other, Value) and self.key() == other.key()

    def __hash__(self):
        return hash(self.key())


def signed(value):
    """The value as a sign and a canonical text of a computation whose result is the same real
    number up to that sign. Rounding to nearest is odd, rnd(-z) = -rnd(z) (rounding.v's
    symmetric), so a negation moves out of any rounded sum, difference, product or fused
    multiply-add: -a flips a's sign, a product takes its factors' signs, a sum sorts its terms
    and factors out the first one's sign, a fused multiply-add its product's. Commutative
    operands are sorted, and a doubling, 2 a or a + a, written one way. Equal computations get
    one text, whatever negations and orders the code or the compiler wrote them with; only the
    sign of a zero may differ, which no proof and no comparison of the predicates tells apart."""
    op = value.op
    if op == "input":
        return 1, value.name
    if op == "constant":
        sign = -1 if struct.pack(">d", value.constant)[0] & 0x80 else 1
        return sign, "c" + struct.pack(">d", abs(value.constant)).hex()
    args = [signed(arg) for arg in value.args]
    if op == "neg":
        return -args[0][0], args[0][1]
    if op == "abs":
        return 1, f"abs({args[0][1]})"
    if op == "mul":
        (s1, a), (s2, b) = args
        if a == "c" + struct.pack(">d", 2.0).hex():
            return s1 * s2, f"twice({b})"
        if b == "c" + struct.pack(">d", 2.0).hex():
            return s1 * s2, f"twice({a})"
        return s1 * s2, f"mul({','.join(sorted((a, b)))})"
    if op in ("add", "sub"):
        (s1, a), (s2, b) = args
        if op == "sub":
            s2 = -s2
        terms = sorted(((a, s1), (b, s2)))
        (first, sign), (second, other) = terms
        relative = "+" if sign * other > 0 else "-"
        if first == second and relative == "+":
            return sign, f"twice({first})"
        return sign, f"sum({first},{relative}{second})"
    if op == "fma":
        (s1, a), (s2, b), (s3, c) = args
        sign = s1 * s2
        relative = "+" if sign * s3 > 0 else "-"
        return sign, f"fma({','.join(sorted((a, b)))},{relative}{c})"
    texts = [("-" if s < 0 else "") + t for s, t in args]
    return 1, f"{op}({','.join(texts)})"


def given(name):
    return Value("input", name=name)


def constant(value):
    return Value("constant", constant=value)


def add(a, b):
    return Value("add", a, b)


def sub(a, b):
    return Value("sub", a, b)


def mul(a, b):
    return Value("mul", a, b)


def neg(a):
    return Value("neg", a)


def magnitude(a):
    return Value("abs", a)


def fma(a, b, c):
    return Value("fma", a, b, c)


def greater(a, b):
    return Value("gt", a, b)


def two_sum(first, second):
    """two.hpp's TwoSum: the sum and its roundoff, rounding.v's two_sum."""
    total = add(first, second)
    virtual_second = sub(total, first)
    virtual_first = sub(total, virtual_second)
    second_round = sub(second, virtual_second)
    first_round = sub(first, virtual_first)
    return total, add(first_round, second_round)


def two_diff(first, second):
    """two.hpp's TwoDiff, by TwoDiffTail: rounding.v's two_diff_tail."""
    difference = sub(first, second)
    virtual_second = sub(first, difference)
    virtual_first = add(difference, virtual_second)
    second_round = sub(virtual_second, second)
    first_round = sub(first, virtual_first)
    return difference, add(first_round, second_round)


def two_product(first, second):
    """two.hpp's TwoProduct where products fuse: rounding.v's two_product_fma."""
    product = mul(first, second)
    return product, fma(first, second, neg(product))


CORNERS = ("a", "b", "c")
TURNS = (("a", "b", "c"), ("b", "c", "a"), ("c", "a", "b"))


def orient_first():
    """orient.hpp's Approximate and Decided, as Orient runs them on (a, b, c)."""
    ax, ay, bx, by, cx, cy = (given(name) for name in ("ax", "ay", "bx", "by", "cx", "cy"))
    v = {"dx1": sub(ax, cx), "dy2": sub(by, cy), "dy1": sub(ay, cy), "dx2": sub(bx, cx)}
    v["left"] = mul(v["dx1"], v["dy2"])
    v["right"] = mul(v["dy1"], v["dx2"])
    v["det"] = sub(v["left"], v["right"])
    v["sum"] = add(v["left"], v["right"])
    v["w"] = add(magnitude(v["sum"]), constant(NORMAL))
    v["bound"] = mul(constant(ORIENT_BOUND_A), v["w"])
    v["decided"] = greater(magnitude(v["det"]), v["bound"])
    return v


SAFE_LOW = 2.0 ** -968
SAFE_HIGH = 2.0 ** 1020


def orient_adapt():
    """orient.cpp's OrientAdapt to CrossSign, on (a, b, c): the differences from c and their
    roundoffs, Safe's two products, and CrossSign's two fused ones."""
    ax, ay, bx, by, cx, cy = (given(name) for name in ("ax", "ay", "bx", "by", "cx", "cy"))
    v = {}
    v["fx"], v["fx tail"] = two_diff(ax, cx)
    v["fy"], v["fy tail"] = two_diff(ay, cy)
    v["sx"], v["sx tail"] = two_diff(bx, cx)
    v["sy"], v["sy tail"] = two_diff(by, cy)
    v["safe first"] = magnitude(mul(v["fx"], v["sy"]))
    v["safe second"] = magnitude(mul(v["fy"], v["sx"]))
    v["product"] = mul(v["fy"], v["sx"])
    v["roundoff"] = fma(v["fy"], v["sx"], neg(v["product"]))
    v["rest"] = fma(v["fx"], v["sy"], neg(v["product"]))
    return v


def differences():
    """Each corner's difference from the point d, rounded: x and y of a, b and c."""
    point = {"x": given("dx"), "y": given("dy")}
    return {f"{corner}{axis}": sub(given(f"{corner}{axis}"), point[axis])
            for corner in CORNERS for axis in "xy"}


def circle_first():
    """circle.hpp's Lift and InCircle's test, on the circle (a, b, c) and the point d."""
    v = differences()
    for own, nxt, before in TURNS:
        v[f"{own}xx"] = mul(v[f"{own}x"], v[f"{own}x"])
        v[f"{own}yy"] = mul(v[f"{own}y"], v[f"{own}y"])
        v[f"{own}lift"] = add(v[f"{own}xx"], v[f"{own}yy"])
        v[f"{own}p1"] = mul(v[f"{nxt}x"], v[f"{before}y"])
        v[f"{own}p2"] = mul(v[f"{before}x"], v[f"{nxt}y"])
        v[f"{own}m"] = sub(v[f"{own}p1"], v[f"{own}p2"])
        v[f"{own}t"] = mul(v[f"{own}lift"], v[f"{own}m"])
        v[f"{own}pm"] = add(magnitude(v[f"{own}p1"]), magnitude(v[f"{own}p2"]))
        v[f"{own}q"] = mul(v[f"{own}pm"], v[f"{own}lift"])
    v["s"] = add(v["at"], v["bt"])
    v["det"] = add(v["s"], v["ct"])
    v["ps"] = add(v["aq"], v["bq"])
    v["perm"] = add(v["ps"], v["cq"])
    v["s1"] = add(v["alift"], v["blift"])
    v["lifts"] = add(v["s1"], v["clift"])
    v["w"] = add(v["lifts"], constant(1.0))
    v["k1"] = mul(constant(CIRCLE_BOUND_A), v["perm"])
    v["k2"] = mul(constant(CIRCLE_BOUND_UNDERFLOW), v["w"])
    v["bound"] = add(v["k1"], v["k2"])
    v["decided"] = greater(magnitude(v["det"]), v["bound"])
    return v


def compensated():
    """circle.cpp's CircleAdapt from its differences to its third stage: Compensated, the two
    thresholds it is decided against, Correction and the third stage's bound. The permanent is
    an input, the first stage's."""
    v = differences()
    point = {"x": given("dx"), "y": given("dy")}
    permanent = given("permanent")
    for corner in CORNERS:
        for axis in "xy":
            _, v[f"{corner}{axis}tail"] = two_diff(given(f"{corner}{axis}"), point[axis])
    for own, nxt, before in TURNS:
        v[f"{own}xB"], v[f"{own}xS"] = two_product(v[f"{own}x"], v[f"{own}x"])
        v[f"{own}yB"], v[f"{own}yS"] = two_product(v[f"{own}y"], v[f"{own}y"])
        v[f"{own}hB"], v[f"{own}hS"] = two_sum(v[f"{own}xB"], v[f"{own}yB"])
        v[f"{own}ls"] = add(v[f"{own}xS"], v[f"{own}yS"])
        v[f"{own}Lo"] = add(v[f"{own}ls"], v[f"{own}hS"])
        v[f"{own}f1"], v[f"{own}e1"] = two_product(v[f"{nxt}x"], v[f"{before}y"])
        v[f"{own}f2"], v[f"{own}e2"] = two_product(v[f"{nxt}y"], v[f"{before}x"])
        v[f"{own}mB"], v[f"{own}mS"] = two_diff(v[f"{own}f1"], v[f"{own}f2"])
        v[f"{own}ms"] = sub(v[f"{own}e1"], v[f"{own}e2"])
        v[f"{own}Mlo"] = add(v[f"{own}ms"], v[f"{own}mS"])
        v[f"{own}p"], v[f"{own}pe"] = two_product(v[f"{own}hB"], v[f"{own}mB"])
        v[f"{own}c1"] = mul(v[f"{own}hB"], v[f"{own}Mlo"])
        v[f"{own}c2"] = mul(v[f"{own}Lo"], v[f"{own}mB"])
        v[f"{own}c3"] = add(v[f"{own}c1"], v[f"{own}c2"])
        v[f"{own}lo"] = add(v[f"{own}pe"], v[f"{own}c3"])
    v["o"], v["os"] = two_sum(v["ap"], v["bp"])
    v["t"], v["ts"] = two_sum(v["o"], v["cp"])
    v["a1"] = add(v["alo"], v["blo"])
    v["a2"] = add(v["a1"], v["clo"])
    v["a3"] = add(v["os"], v["ts"])
    v["low"] = add(v["a2"], v["a3"])
    v["E"] = add(v["t"], v["low"])
    v["exact bound"] = mul(constant(CIRCLE_BOUND_COMPENSATED), permanent)
    v["second bound"] = mul(constant(CIRCLE_BOUND_B), permanent)
    for own, nxt, before in TURNS:
        a = (v[f"{own}x"], v[f"{own}y"])
        b = (v[f"{nxt}x"], v[f"{nxt}y"])
        c = (v[f"{before}x"], v[f"{before}y"])
        at = (v[f"{own}xtail"], v[f"{own}ytail"])
        bt = (v[f"{nxt}xtail"], v[f"{nxt}ytail"])
        ct = (v[f"{before}xtail"], v[f"{before}ytail"])
        v[f"{own}Lh"] = add(mul(a[0], a[0]), mul(a[1], a[1]))
        v[f"{own}q"] = sub(add(mul(b[0], ct[1]), mul(c[1], bt[0])),
                           add(mul(b[1], ct[0]), mul(c[0], bt[1])))
        v[f"{own}f"] = mul(v[f"{own}Lh"], v[f"{own}q"])
        v[f"{own}w"] = add(mul(a[0], at[0]), mul(a[1], at[1]))
        v[f"{own}m"] = sub(mul(b[0], c[1]), mul(b[1], c[0]))
        v[f"{own}g"] = mul(mul(constant(2.0), v[f"{own}w"]), v[f"{own}m"])
        v[f"{own}corr"] = add(v[f"{own}f"], v[f"{own}g"])
    v["cs"] = add(v["acorr"], v["bcorr"])
    v["Corr"] = add(v["cs"], v["ccorr"])
    v["C"] = add(v["E"], v["Corr"])
    v["third k"] = mul(constant(CIRCLE_BOUND_C + CIRCLE_BOUND_COMPENSATED), permanent)
    v["third r"] = mul(constant(RESULT_BOUND), magnitude(v["E"]))
    v["third bound"] = add(v["third k"], v["third r"])
    return v


def rotated(value, mapping):
    """A value with its inputs renamed: a corner's computation turned into another's."""
    if value.op == "input":
        return given(mapping.get(value.name, value.name))
    if value.op == "constant":
        return value
    return Value(value.op, *(rotated(arg, mapping) for arg in value.args))


def corners_alike(values, name):
    """Whether a corner's value of a graph is the next corner's turned: so a lemma proved of one
    corner holds of all three."""
    turn = {"ax": "bx", "ay": "by", "bx": "cx", "by": "cy", "cx": "ax", "cy": "ay"}
    return (rotated(values[f"a{name}"], turn) == values[f"b{name}"]
            and rotated(values[f"b{name}"], turn) == values[f"c{name}"])
