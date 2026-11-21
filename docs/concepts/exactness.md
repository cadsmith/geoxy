# Exactness

Floating point arithmetic rounds every step. A geometry library that ignores this gives answers that are nearly right most of the time, and wrong in ways that are hard to see: an intersection missed because two lines were almost parallel, a point reported on the wrong side of a line. geoxy is written so that, wherever it can be done, the answer is the one exact arithmetic would give, rounded once.

## What "exact" means here {#meaning}

A double cannot hold most numbers. The best a function returning doubles can do is return the double nearest the true answer. That is what geoxy means by exact: the true result, rounded once at the end, not the accumulated rounding of every step along the way.

Where this holds, two things follow. The answer does not depend on how the input happens to be written: a line given from left to right or right to left meets a circle in the same points. And decisions are right: a point is reported on the side of a line it is truly on.

## Where it holds {#where}

- **Decisions.** `Orient`, the side of a line a point is on, is always right, whatever the coordinates. So are the decisions built on it.
- **Lines and arcs.** Where two of them meet, the length of an arc, the area of a polyline, the extents of an arc, the arc of a fillet between them: each is the nearest double to the true value, or within a step or two of it.
- **Offsets of polylines.** Each vertex of an offset is the nearest double to its true place in the great majority of cases, and within a few steps otherwise.

## Where it is approximated {#approximate}

Some results have no closed form: the length of an ellipse or a spline, the offset of a spline, where two splines meet. geoxy computes these to a few steps of a double where it can, and to the `Accuracy` you pass where the answer is a chain of chords or arcs. Such calls take an `Accuracy` in their signature; that is how to tell them from the exact ones.

## What the tolerance is not {#tolerance}

A tolerance does not make a computation looser. It answers one question, whether two things that are close count as one, and it is asked only where that question arises: two curves that graze, two ends that should join, a loop that should close. The points themselves are computed from the curves' own numbers.

## Far from the origin {#far}

A double has about 16 digits. At coordinates in the millions, as survey drawings have, 7 of them are spent before the decimal point. geoxy's exact results stay exact there; what changes is the size of a step, and `Tolerance::Of` takes it into account. Nothing in the library assumes a drawing sits near the origin or is of unit size.

## The same answer on every run {#same}

geoxy is built so that each floating point operation is rounded on its own (`-ffp-contract=off`), so a compiler does not fuse or reorder its arithmetic. What does not depend on the system's math library, which is everything decided between lines, gives the same bits wherever it is built. Sines, cosines and arc tangents come from the system's library by default; the `GEOXY_ENABLE_TRIGO` option builds geoxy's own instead.
