"""An exact oracle for where two curves meet, sharing nothing with geoxy.

Every number is a Fraction. A line or a conic has a rational implicit equation; the other curve,
a polynomial or a conic's rational parametrization, put into it gives a polynomial whose real
roots Sturm's theorem counts and isolates exactly, a repeated one a tangency, an identically zero
one a shared carrier. Two polynomial pieces meet where the resultant of their coordinates'
differences, a determinant of polynomials worked out by Bareiss's fraction-free elimination,
vanishes. Roots are narrowed to 2^-140; a place within 2^-100 of a boundary is on it, as only a
case built to be there is.
"""

from dataclasses import dataclass, field
from fractions import Fraction as F

NARROW = F(1, 2 ** 140)
INF = F(2 ** 4000)
ON = F(1, 2 ** 100)


# Polynomials over the rationals, their coefficients from the constant up.

def trim(p):
    """The polynomial with its trailing zero coefficients dropped."""
    p = list(p)
    while p and p[-1] == 0:
        p.pop()
    return p


def add(p, q):
    """The sum of two polynomials."""
    n = max(len(p), len(q))
    return trim([(p[i] if i < len(p) else 0) + (q[i] if i < len(q) else 0) for i in range(n)])


def sub(p, q):
    """The difference of two polynomials, p minus q."""
    return add(p, [-c for c in q])


def mul(p, q):
    """The product of two polynomials."""
    if not p or not q:
        return []
    out = [F(0)] * (len(p) + len(q) - 1)
    for i, a in enumerate(p):
        if a:
            for j, b in enumerate(q):
                out[i + j] += a * b
    return trim(out)


def scale(p, c):
    """The polynomial with every coefficient multiplied by c."""
    return trim([a * c for a in p])


def deriv(p):
    """The derivative of a polynomial."""
    return trim([i * p[i] for i in range(1, len(p))])


def value(p, x):
    """The polynomial's value at x, by Horner's rule."""
    out = F(0)
    for c in reversed(p):
        out = out * x + c
    return out


def divmod_(p, q):
    """The quotient and remainder of dividing polynomial p by q."""
    p, q = trim(p), trim(q)
    out = [F(0)] * max(len(p) - len(q) + 1, 1)
    while len(p) >= len(q) and p:
        k = len(p) - len(q)
        c = p[-1] / q[-1]
        out[k] = c
        p = sub(p, [F(0)] * k + scale(q, c))
    return trim(out), p


def gcd(p, q):
    """The monic greatest common divisor of two polynomials, or [] where both are zero."""
    p, q = trim(p), trim(q)
    while q:
        p, q = q, divmod_(p, q)[1]
    return scale(p, 1 / p[-1]) if p else []


def power(p, n):
    """The polynomial raised to the power n."""
    out = [F(1)]
    for _ in range(n):
        out = mul(out, p)
    return out


def sturm(p):
    """The Sturm chain of a polynomial."""
    chain = [p, deriv(p)]
    while chain[-1]:
        rest = divmod_(chain[-2], chain[-1])[1]
        if not rest:
            break
        chain.append(scale(rest, -1))
    return chain


def changes(chain, x):
    """The number of sign changes along the chain's values at x, zeros skipped."""
    signs = [s for s in ((value(p, x) > 0) - (value(p, x) < 0) for p in chain) if s]
    return sum(1 for a, b in zip(signs, signs[1:]) if a != b)


def bound(p):
    """A bound that every real root of the polynomial lies within."""
    return 1 + max(abs(c / p[-1]) for c in p[:-1]) if len(p) > 1 else F(1)


def roots(p, lo, hi):
    """The real roots of a polynomial on [lo, hi], each narrowed to a rational within 2^-140,
    and whether it is a repeated root."""
    p = trim(p)
    if len(p) <= 1:
        return []
    repeated = gcd(p, deriv(p))
    free = divmod_(p, repeated)[0] if len(repeated) > 1 else p
    lo = -bound(free) if lo is None else F(lo)
    hi = bound(free) if hi is None else F(hi)
    simple = isolated(free, lo, hi)
    doubled = isolated(repeated, lo, hi) if len(repeated) > 1 else []
    return [(r, any(abs(r - d) <= ON for d in doubled)) for r in simple]


def isolated(free, lo, hi):
    """The roots of a square-free polynomial on [lo, hi], by Sturm's theorem and bisection."""
    chain = sturm(free)
    found = [lo] if value(free, lo) == 0 else []
    pending = [(lo, hi)]
    while pending:
        a, b = pending.pop()
        count = changes(chain, a) - changes(chain, b)
        if count == 0:
            continue
        if count == 1 and value(free, a) != 0:
            found.append(b if value(free, b) == 0 else narrowed(free, a, b))
            continue
        middle = (a + b) / 2
        pending.append((a, middle))
        pending.append((middle, b))
    return sorted(set(found))


def simplest(a, b):
    """The fraction of the least denominator from a to b, both taken in, by their continued
    fractions."""
    whole = a.numerator // a.denominator
    if whole == a:
        return a
    if whole < b.numerator // b.denominator:
        return F(whole + 1)
    return whole + 1 / simplest(1 / (b - whole), 1 / (a - whole))


def narrowed(p, a, b):
    """The root of p in (a, b), found by bisection down to 2^-140: the root itself where it is
    rational, as the fraction of the least denominator left in so narrow a stretch then is, and
    is tried; so that a meeting at a rational place is written exactly, a coordinate that is
    nothing as nothing, and not as what is left of 2^-140."""
    fa = value(p, a)
    while b - a > NARROW:
        m = (a + b) / 2
        fm = value(p, m)
        if fm == 0:
            return m
        if (fm > 0) == (fa > 0):
            a, fa = m, fm
        else:
            b = m
    plain = simplest(a, b)
    return plain if value(p, plain) == 0 else (a + b) / 2


# Points and the curves' pieces.

def cross(a, b):
    """The cross product of two vectors."""
    return a[0] * b[1] - a[1] * b[0]


def dot(a, b):
    """The dot product of two vectors."""
    return a[0] * b[0] + a[1] * b[1]


def minus(a, b):
    """The vector from b to a."""
    return (a[0] - b[0], a[1] - b[1])


def near(a, b):
    """Whether two points are within 2^-100 of each other on both axes."""
    return abs(a[0] - b[0]) <= ON and abs(a[1] - b[1]) <= ON


@dataclass
class Piece:
    """A part of a curve: polynomial coordinates and weight over a parameter's stretch, or a
    conic's rational parametrization with its arc. Ends are where the whole curve ends."""
    x: list
    y: list
    w: list
    lo: object
    hi: object
    implicit: dict = None
    arc: tuple = None
    unit: tuple = None
    straight: tuple = None
    infinity: tuple = None
    handed: int = 1


def point_at(piece, t):
    """The point a piece has at parameter t."""
    w = value(piece.w, t)
    return (value(piece.x, t) / w, value(piece.y, t) / w)


def line(start, way, lo, hi):
    """A piece along start + way * t over [lo, hi], a single point if way is zero."""
    if way[0] == 0 and way[1] == 0:
        return Piece(x=[start[0]] if start[0] else [], y=[start[1]] if start[1] else [], w=[F(1)],
                     lo=F(0), hi=F(0))
    p = Piece(x=trim([start[0], way[0]]), y=trim([start[1], way[1]]), w=[F(1)], lo=lo, hi=hi)
    p.implicit = {(0, 0): -cross(start, way), (1, 0): way[1], (0, 1): -way[0]}
    p.straight = (start, way, lo, hi)
    return p


def bezier(points, weights=None):
    """The piece for a Bézier curve with these control points, over [0, 1]; rational where
    weights are given, each positive: x and y the weighted sums, w the weights' sum."""
    n = len(points) - 1
    weights = [F(1)] * (n + 1) if weights is None else [F(w) for w in weights]
    xs, ys, ws = [], [], []
    for axis, out in ((0, xs), (1, ys), (2, ws)):
        poly = []
        for i, pt in enumerate(points):
            basis = scale(mul(power([F(0), F(1)], i), power([F(1), F(-1)], n - i)), F(comb(n, i)))
            value_ = weights[i] if axis == 2 else weights[i] * pt[axis]
            poly = add(poly, scale(basis, value_))
        out.extend(poly or [F(0)])
    return Piece(x=xs, y=ys, w=trim(ws), lo=F(0), hi=F(1))


def level(piece, pt):
    """The equation along one axis that a polynomial piece meets a point's coordinate on, times
    the weight: x - px w, or y - py w where x - px w vanishes everywhere."""
    across_x = sub(piece.x, scale(piece.w, pt[0]))
    return across_x if across_x else sub(piece.y, scale(piece.w, pt[1]))


def homogeneous(a_axis, a_w, b):
    """The coefficients in a's parameter of a's coordinate times b's weight less b's coordinate
    times a's weight, each a polynomial in b's parameter."""
    b_axis, b_w = b
    count = max(len(a_axis), len(a_w))
    out = []
    for k in range(count):
        ak = a_axis[k] if k < len(a_axis) else F(0)
        wk = a_w[k] if k < len(a_w) else F(0)
        out.append(sub(scale(b_w, ak), scale(b_axis, wk)))
    while len(out) > 1 and not out[-1]:
        out.pop()
    return out


def comb(n, k):
    """The binomial coefficient n choose k."""
    out = 1
    for i in range(k):
        out = out * (n - i) // (i + 1)
    return out


def conic(origin, xaxis, yaxis, start, end, counterclockwise, full):
    """An ellipse arc: origin + xaxis cos + yaxis sin, from the unit vector start to end, both
    rational, the way given, or the whole of it. Parametrized by u through a rational turn that
    keeps the point at u = infinity off both ends."""
    for c, s in ((F(3, 5), F(4, 5)), (F(5, 13), F(12, 13)), (F(8, 17), F(15, 17)),
                 (F(20, 29), F(21, 29))):
        far = (-c, -s)
        if full or (not near(far, start) and not near(far, end)):
            break
    one_minus = [F(1), F(0), F(-1)]
    two_u = [F(0), F(2)]
    cn = sub(scale(one_minus, c), scale(two_u, s))
    sn = add(scale(one_minus, s), scale(two_u, c))
    den = [F(1), F(0), F(1)]
    x = add(add(scale(den, origin[0]), scale(cn, xaxis[0])), scale(sn, yaxis[0]))
    y = add(add(scale(den, origin[1]), scale(cn, xaxis[1])), scale(sn, yaxis[1]))
    piece = Piece(x=x, y=y, w=den, lo=None, hi=None)
    det = cross(xaxis, yaxis)
    inv = ((yaxis[1] / det, -xaxis[1] / det), (-yaxis[0] / det, xaxis[0] / det))
    # unit = inv applied to (P - origin): u = inv[0] * dx + inv[1] * dy, columns.
    a, b = inv[0][0], inv[1][0]
    cc, d = inv[0][1], inv[1][1]
    ox, oy = origin
    ux = {(1, 0): a, (0, 1): b, (0, 0): -(a * ox + b * oy)}
    uy = {(1, 0): cc, (0, 1): d, (0, 0): -(cc * ox + d * oy)}
    piece.implicit = add2(add2(sq2(ux), sq2(uy)), {(0, 0): F(-1)})
    piece.unit = (origin, inv)
    piece.arc = (start, end, counterclockwise, full)
    piece.handed = 1 if det > 0 else -1
    piece.infinity = (origin[0] + xaxis[0] * far[0] + yaxis[0] * far[1],
                      origin[1] + xaxis[1] * far[0] + yaxis[1] * far[1])
    return piece


def add2(p, q):
    """The sum of two polynomials in two variables, as dictionaries of exponent pairs."""
    out = dict(p)
    for k, v in q.items():
        out[k] = out.get(k, 0) + v
    return {k: v for k, v in out.items() if v}


def sq2(p):
    """The square of a polynomial in two variables."""
    out = {}
    for (i, j), a in p.items():
        for (k, m), b in p.items():
            out[(i + k, j + m)] = out.get((i + k, j + m), 0) + a * b
    return {k: v for k, v in out.items() if v}


def unit_of(piece, pt):
    """The point of the unit circle that a point on the conic's plane maps to."""
    (ox, oy), inv = piece.unit
    dx, dy = pt[0] - ox, pt[1] - oy
    return (inv[0][0] * dx + inv[1][0] * dy, inv[0][1] * dx + inv[1][1] * dy)


def half(a, v):
    """Which half turn from a the vector v lies in: 0 for the first, 1 for the second."""
    c = cross(a, v)
    return 0 if c > ON or (abs(c) <= ON and dot(a, v) > 0) else 1


def turn_before(a, v, w):
    """Whether, turning counterclockwise from a, v comes no later than w."""
    hv, hw = half(a, v), half(a, w)
    if hv != hw:
        return hv < hw
    return cross(v, w) >= -ON


def on_arc(piece, pt):
    """Whether a point of a conic's carrier lies on the arc."""
    if piece.arc is None:
        return True
    start, end, counterclockwise, full = piece.arc
    if full:
        return True
    q = unit_of(piece, pt)
    if near(q, start) or near(q, end):
        return True
    if counterclockwise:
        return turn_before(start, q, end)
    return turn_before(end, q, start)


def on_straight(piece, pt):
    """Whether a point of a line's carrier lies in the line's stretch."""
    if piece.straight is None:
        return True
    start, way, lo, hi = piece.straight
    s = dot(minus(pt, start), way) / dot(way, way)
    return (lo is None or s >= lo - ON) and (hi is None or s <= hi + ON)


def compose(implicit, piece):
    """An implicit equation of degree d along a piece, times its weight to the d."""
    degree = max(i + j for i, j in implicit)
    out = []
    for (i, j), c in implicit.items():
        term = mul(mul(power(piece.x, i), power(piece.y, j)), power(piece.w, degree - i - j))
        out = add(out, scale(term, c))
    return out


@dataclass
class Found:
    """What a meeting of two pieces finds: each hit with whether it is a tangency and the sine of
    the angle between the curves there, None where no angle is known."""
    hits: list = field(default_factory=list)      # (point, tangent)
    shared: list = field(default_factory=list)    # (point, point)
    sines: list = field(default_factory=list)     # one for each hit
    places: list = field(default_factory=list)    # each hit's (place, way) on a, then on b

    def hit(self, point, tangent, sine=None, places=None):
        """Records a hit, with where it is on each piece and the way each runs there."""
        self.hits.append((point, tangent))
        self.sines.append(sine)
        self.places.append(places)

    def swapped(self):
        """The same meeting with the pieces the other way round."""
        return Found(self.hits, self.shared, self.sines,
                     [None if p is None else (p[1], p[0]) for p in self.places])


def place(piece, pt, t=None):
    """Where a point is on a piece, and the way the piece runs there, exactly: a conic's point
    by its unit vector and its way by the sense its arc runs; a line's by its parameter; a
    polynomial piece's by the parameter given or found."""
    if piece.arc is not None:
        gx, gy = gradient(piece.implicit, pt)
        sense = piece.handed * (1 if piece.arc[2] else -1)
        return unit_of(piece, pt), (-gy * sense, gx * sense)
    if piece.straight is not None:
        start, way, _, _ = piece.straight
        return dot(minus(pt, start), way) / dot(way, way), way
    t = placed(piece, pt) if t is None else t
    return t, (None if t is None else tangent(piece, t))


def tangent(piece, t):
    """A piece's direction at a parameter, exactly: the derivative of x/w and y/w, times w^2."""
    w, dw = value(piece.w, t), value(deriv(piece.w), t)
    return (value(deriv(piece.x), t) * w - value(piece.x, t) * dw,
            value(deriv(piece.y), t) * w - value(piece.y, t) * dw)


def gradient(implicit, pt):
    """The gradient of an implicit equation at a point, exactly."""
    gx = sum(c * i * pt[0] ** (i - 1) * pt[1] ** j for (i, j), c in implicit.items() if i)
    gy = sum(c * j * pt[0] ** i * pt[1] ** (j - 1) for (i, j), c in implicit.items() if j)
    return (gx, gy)


def sine(u, v):
    """The sine of the angle between two directions, as a float; 0 where either vanishes."""
    lengths = float(dot(u, u)) * float(dot(v, v))
    return 0.0 if lengths == 0 else min(1.0, abs(float(cross(u, v))) / lengths ** 0.5)


def across(a, pt):
    """A piece's direction at a point on it: along its carrier's level set where it has an
    implicit equation, else its derivative at the parameter that reaches the point."""
    if a.implicit is not None:
        gx, gy = gradient(a.implicit, pt)
        return (-gy, gx)
    t = placed(a, pt)
    return None if t is None else tangent(a, t)


def implicit_meet(a, b):
    """Where a piece with an implicit equation meets another, by the other's parametrization."""
    found = Found()
    g = compose(a.implicit, b)
    if not g:
        pairs, touched = carrier(a, b)
        found.shared.extend(pairs)
        for p in touched:
            found.hit(p, False)
        return found
    for t, repeated in roots(g, b.lo, b.hi):
        pt = point_at(b, t)
        if on_arc(b, pt) and on_arc(a, pt) and on_straight(a, pt):
            way = across(b, pt) if b.implicit is not None else tangent(b, t)
            found.hit(pt, repeated, sine(across(a, pt), way), (place(a, pt), place(b, pt, t)))
    if b.infinity is not None:
        pt = b.infinity
        if (value_implicit(a.implicit, pt) == 0 and on_arc(b, pt) and on_arc(a, pt)
                and on_straight(a, pt)):
            found.hit(pt, False, sine(across(a, pt), across(b, pt)), (place(a, pt), place(b, pt)))
    return found


def value_implicit(implicit, pt):
    """The value of an implicit equation at a point."""
    return sum(c * pt[0] ** i * pt[1] ** j for (i, j), c in implicit.items())


def ends(piece):
    """Where a piece starts and ends, or None for a whole conic or an endless line."""
    if piece.arc is not None:
        start, end, _, full = piece.arc
        if full:
            return []
        (ox, oy), inv = piece.unit
        det = inv[0][0] * inv[1][1] - inv[1][0] * inv[0][1]
        # back from the unit circle: P = origin + M q, M the inverse of inv
        m = ((inv[1][1] / det, -inv[0][1] / det), (-inv[1][0] / det, inv[0][0] / det))
        out = []
        for q in (start, end):
            out.append((ox + m[0][0] * q[0] + m[1][0] * q[1], oy + m[0][1] * q[0] + m[1][1] * q[1]))
        return out
    out = []
    for t in (piece.lo, piece.hi):
        if t is not None:
            out.append(point_at(piece, t))
    return out


def carrier(a, b):
    """Where b lies on a's carrier and on a too: the stretches shared, as pairs of ends (None for
    an end at infinity), and the points they only touch at."""
    if a.straight is not None:
        return straight_carrier(a, b)
    if a.arc is not None and b.arc is not None:
        return conic_carrier(a, b)
    return [], []


def straight_carrier(a, b):
    """The stretches and touching points where b lies on a straight piece a."""
    start, way, lo, hi = a.straight
    square = dot(way, way)

    def at(q):
        """How far along the line's way a point lies, in units of way."""
        return dot(minus(q, start), way) / square

    if b.straight is not None:
        other, bway, blo, bhi = b.straight
        base, rate = at(other), dot(bway, way) / square
        values = []
        for t, end in ((blo, -1), (bhi, 1)):
            if t is None:
                values.append(end * (1 if rate > 0 else -1) * INF)
            else:
                values.append(base + rate * t)
        low_b, high_b = min(values), max(values)
    else:
        ends_b = ends(b)
        values = sorted(at(q) for q in ends_b)
        low_b, high_b = values[0], values[-1]
    low = max(-INF if lo is None else lo, low_b)
    high = min(INF if hi is None else hi, high_b)

    def place(s):
        """The point at distance s along the line, or None at infinity."""
        if abs(s) == INF:
            return None
        return (start[0] + way[0] * s, start[1] + way[1] * s)

    if low < high:
        return [(place(low), place(high))], []
    if low == high:
        return [], [place(low)]
    return [], []


def conic_carrier(a, b):
    """The stretches and touching points where two conic pieces share a carrier."""
    if a.arc[3] and b.arc[3]:
        start = ends_of_full(a)
        return [(start, start)], []
    points = []
    for p in ends(a) + ends(b):
        if on_arc(a, p) and on_arc(b, p) and not any(near(p, q) for q in points):
            points.append(p)
    if not points:
        return [], []
    keyed = sorted(points, key=lambda p: angle_key(a.arc[0], unit_of(a, p)))
    pairs = []
    for p, q in zip(keyed, keyed[1:] + keyed[:1]):
        if near(p, q):
            continue
        middle = arc_middle(a, p, q)
        if on_arc(a, middle) and on_arc(b, middle):
            pairs.append((p, q))
    touched = [p for p in keyed if not any(near(p, e) for pair in pairs for e in pair)]
    return pairs, touched


def ends_of_full(piece):
    """The point of a whole conic where its unit circle's start direction lands."""
    (ox, oy), inv = piece.unit
    det = inv[0][0] * inv[1][1] - inv[1][0] * inv[0][1]
    m = ((inv[1][1] / det, -inv[0][1] / det), (-inv[1][0] / det, inv[0][0] / det))
    q = piece.arc[0]
    return (ox + m[0][0] * q[0] + m[1][0] * q[1], oy + m[0][1] * q[0] + m[1][1] * q[1])


def angle_key(start, v):
    """A sort key for a direction going counterclockwise from start."""
    return (half(start, v), Slope(start, v))


class Slope:  # pylint: disable=too-few-public-methods
    """Orders directions within a half turn from a start by cross products, exactly."""

    def __init__(self, start, v):
        """Keeps the direction and the start it is measured from."""
        self.start, self.v = start, v

    def __lt__(self, other):
        """Whether this direction comes before the other, by their cross product."""
        return cross(self.v, other.v) > 0


def mid(p, q):
    """The midpoint of two points."""
    return ((p[0] + q[0]) / 2, (p[1] + q[1]) / 2)


def arc_middle(piece, p, q):
    """A point of a conic's carrier between two of its points, going counterclockwise in the
    unit circle's frame from p to q: a direction is all on_arc reads."""
    a, b = unit_of(piece, p), unit_of(piece, q)
    if abs(a[0] + b[0]) <= ON and abs(a[1] + b[1]) <= ON:
        s = (-a[1], a[0])
    elif cross(a, b) >= 0:
        s = (a[0] + b[0], a[1] + b[1])
    else:
        s = (-(a[0] + b[0]), -(a[1] + b[1]))
    (ox, oy), inv = piece.unit
    det = inv[0][0] * inv[1][1] - inv[1][0] * inv[0][1]
    m = ((inv[1][1] / det, -inv[0][1] / det), (-inv[1][0] / det, inv[0][0] / det))
    return (ox + m[0][0] * s[0] + m[1][0] * s[1], oy + m[0][1] * s[0] + m[1][1] * s[1])


def on_all(piece, pt):
    """Whether a point of any piece's carrier lies on the piece."""
    if piece.arc is not None:
        return on_arc(piece, pt)
    if piece.straight is not None:
        return on_straight(piece, pt)
    return matches(piece, pt)


def sylvester_det(p, q):
    """The resultant in s of two polynomials whose coefficients are polynomials in t, by
    Bareiss's fraction-free elimination over polynomials in t."""
    m, n = len(p) - 1, len(q) - 1
    size = m + n
    rows = []
    for i in range(n):
        rows.append([[] for _ in range(i)] + [list(c) for c in reversed(p)]
                    + [[] for _ in range(size - m - 1 - i)])
    for i in range(m):
        rows.append([[] for _ in range(i)] + [list(c) for c in reversed(q)]
                    + [[] for _ in range(size - n - 1 - i)])
    sign = 1
    previous = [F(1)]
    for k in range(size - 1):
        if not trim(rows[k][k]):
            swap = next((r for r in range(k + 1, size) if trim(rows[r][k])), None)
            if swap is None:
                return []
            rows[k], rows[swap] = rows[swap], rows[k]
            sign = -sign
        for i in range(k + 1, size):
            for j in range(k + 1, size):
                top = sub(mul(rows[i][j], rows[k][k]), mul(rows[i][k], rows[k][j]))
                rows[i][j] = divmod_(top, previous)[0]
            rows[i][k] = []
        previous = rows[k][k]
    return scale(rows[size - 1][size - 1], sign)


def is_point(piece):
    """Whether a polynomial piece is a single point."""
    start = point_at(piece, F(piece.lo))
    return not sub(piece.x, scale(piece.w, start[0])) and not sub(piece.y, scale(piece.w, start[1]))


def bezier_meet(a, b):
    """Where two polynomial pieces meet: the roots in b's parameter of the resultant in a's of
    their coordinates' differences, each matched to a's parameter."""
    found = Found()
    if is_point(a) or is_point(b):
        point, other = (point_at(a, F(0)), b) if is_point(a) else (point_at(b, F(0)), a)
        if matches(other, point):
            found.hit(point, False)
        return found
    p = homogeneous(a.x, a.w, (b.x, b.w))
    q = homogeneous(a.y, a.w, (b.y, b.w))
    if len(p) == 1:
        p, q = q, p
    if len(q) == 1:
        # a runs along one axis only: b meets it where b's coordinate on that axis is a's.
        g = q[0]
        for t, repeated in roots(g, b.lo, b.hi) if g else []:
            pt = point_at(b, t)
            ways_on(found, (a, b), (pt, t, repeated))
        if not g:
            found.shared.extend(shared_polynomial(a, b))
        return found
    resultant = sylvester_det(p, q)
    if not resultant:
        found.shared.extend(shared_polynomial(a, b))
        return found
    for t, repeated in roots(resultant, b.lo, b.hi):
        pt = point_at(b, t)
        ways_on(found, (a, b), (pt, t, repeated))
    return found


def ways_on(found, pieces, where):
    """Records a hit of b's at a point and parameter on every place of a's there: two where a
    passes through the point twice, as at a loop's double point, each then a crossing of its
    own however the resultant's root repeats."""
    a, b = pieces
    pt, t, repeated = where
    on = places_on(a, pt)
    for s in on:
        found.hit(pt, repeated and len(on) == 1, angle(a, b, (pt, t)),
                  ((s, tangent(a, s)), (t, tangent(b, t))))


def places_on(piece, pt):
    """Every parameter at which a polynomial piece passes within 2^-80 of a point."""
    if is_point(piece):
        return [F(0)] if close(point_at(piece, F(0)), pt) else []
    lo, hi = F(piece.lo) - ON, F(piece.hi) + ON
    return [t for t, _ in roots(level(piece, pt), lo, hi) if close(point_at(piece, t), pt)]


def angle(a, b, where):
    """The sine of the angle two polynomial pieces meet at, b's point and parameter given."""
    way = across(a, where[0])
    return None if way is None else sine(way, tangent(b, where[1]))


def placed(piece, pt):
    """The parameter at which a polynomial piece passes within 2^-80 of a point, in its
    stretch, or None."""
    found = places_on(piece, pt)
    return found[0] if found else None


def matches(piece, pt):
    """Whether a polynomial piece passes within 2^-80 of a point, at a parameter in its
    stretch."""
    return placed(piece, pt) is not None


def close(a, b):
    """Whether two points are within 2^-80 of each other on both axes."""
    return abs(a[0] - b[0]) <= F(1, 2 ** 80) and abs(a[1] - b[1]) <= F(1, 2 ** 80)


def shared_polynomial(a, b):
    """The stretch two polynomial pieces share, as a pair of ends, if they share one."""
    points = []
    for p in ends(a):
        if matches(b, p):
            points.append(p)
    for p in ends(b):
        if matches(a, p) and not any(near(p, q) for q in points):
            points.append(p)
    if len(points) == 2:
        return [(points[0], points[1])]
    return []


def reach2(piece, q):
    """The square of the distance from a point to a piece, exactly to the 2^-140 its roots are
    narrowed to: the least over its ends and the places where (P - q).P' vanishes, P = (x, y)/w,
    whose numerator is (X - qx W)(X' W - X W') + (Y - qy W)(Y' W - Y W'); for a conic over its
    whole parametrization and its point at infinity, kept to its arc."""
    if piece.straight is not None:
        start, way, lo, hi = piece.straight
        s = dot(minus(q, start), way) / dot(way, way)
        s = max(s, lo) if lo is not None else s
        s = min(s, hi) if hi is not None else s
        foot = (start[0] + way[0] * s, start[1] + way[1] * s)
        return dot(minus(q, foot), minus(q, foot))
    x, y, w = piece.x, piece.y, piece.w
    dx = sub(mul(deriv(x), w), mul(x, deriv(w)))
    dy = sub(mul(deriv(y), w), mul(y, deriv(w)))
    turn = add(mul(sub(x, scale(w, q[0])), dx), mul(sub(y, scale(w, q[1])), dy))
    places = [point_at(piece, t) for t, _ in roots(turn, piece.lo, piece.hi)] if turn else []
    places += ends(piece)
    if piece.infinity is not None:
        places.append(piece.infinity)
    if piece.arc is not None:
        places = [p for p in places if on_arc(piece, p)]
    return min(dot(minus(q, p), minus(q, p)) for p in places)


def meet_pieces(a, b):
    """Where two pieces meet, as hits and shared stretches."""
    if a.implicit is not None:
        return implicit_meet(a, b)
    if b.implicit is not None:
        return implicit_meet(b, a).swapped()
    return bezier_meet(a, b)
