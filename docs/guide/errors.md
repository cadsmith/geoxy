# Errors

geoxy never throws. A call that may have no answer says so in its return type, and a caller reads it before using the answer.

## Three ways a call answers {#ways}

- **A plain value**, where there is always an answer: `Length` of an arc, `Nearest` on a curve, `Extents`.
- **A `std::optional`**, where there may be none and the reason is obvious: `Unit` of a vector of no length, `Tangent` at a place where the curve stands still, `Join` for curves that do not touch.
- **A `std::expected`**, where the reason matters: `Offset`, the region operations and `Fit` give an `Error`; `Fillet` and `Chamfer` give a `Misfit`.

## Error {#error}

| Value | Means |
| :--- | :--- |
| `Degenerate` | What was given cannot have an answer, such as a loop that does not close |
| `NotFinite` | A number given is infinite or not a number |
| `TooLarge` | A size asked for is past what the shapes allow |
| `NoConvergence` | An approximation did not reach its accuracy |
| `Spent` | A work budget ran out, so part of the answer may be missing |

## Misfit {#misfit}

`Fillet` and `Chamfer` say why a corner could not be made: `NotFinite`, `Degenerate`, `Parallel` (the two curves never form a corner), `TooLarge` (the radius or distance does not fit) or `Apart` (the curves do not reach each other).

```cpp [main.cpp]
#include "geoxy/geoxy.hpp"

#include <print>

auto main() -> int
{
  const Geoxy::Tolerance tolerance = {.Length = 1e-9, .Radians = 1e-9};
  const Geoxy::Polyline  dot       = {.Vertices = {{.Position = {.X = 1, .Y = 1}}}, .Closed = false};

  const auto offset = Geoxy::Offset(dot, {.Distance = 2}, tolerance);
  if (!offset)
  {
    std::println("no offset: {}", offset.error() == Geoxy::Error::Degenerate ? "degenerate" : "other");
  }

  const auto unit = Geoxy::Unit(Geoxy::Vector{.X = 0, .Y = 0});
  std::println("a vector of no length has {} direction", unit ? "a" : "no");
}
```

```text [Output]
no offset: degenerate
a vector of no length has no direction
```

## Spent {#spent}

A few operations bound their own work so that no input can hang them. `Intersect` sets `Intersections::Spent` when it stopped before it was sure it had every hit; the hits it did find are still right. It is never silent: check the flag where a missing hit would matter.

## Numbers that are not finite {#finite}

Passing an infinity or a NaN is not undefined behavior. A measure of such a curve comes back as a NaN or as its documented fallback, a decision as "no", and an operation with an `Error` as `NotFinite`.
