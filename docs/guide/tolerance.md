# Tolerance and accuracy

Two small types tell an operation how to answer. They are passed to every call that needs one; geoxy keeps no global setting.

## Tolerance: when two things are one {#tolerance}

A `Tolerance` holds a `Length`, the distance within which two points count as the same point, and `Radians`, the angle within which two directions count as the same. Operations that decide something take it: `Intersect` to tell a touch from a near miss, `Join` to tell whether two ends join, the region operations to tell whether a loop closes.

`Tolerance::Of(extents)` gives one fitted to a model: 1e-10 of the model's larger side, the ratio AutoCAD's default has to a unit drawing, and never less than 16 steps of a double at the farthest corner, so that it still means something at survey coordinates.

```cpp [main.cpp]
#include "geoxy/geoxy.hpp"

#include <print>

auto main() -> int
{
  const Geoxy::Box       plan = {.Min = {.X = 0, .Y = 0}, .Max = {.X = 50, .Y = 30}};
  const Geoxy::Tolerance near = Geoxy::Tolerance::Of(plan);
  std::println("a 50 unit plan:      {}", near.Length);

  const Geoxy::Box       survey = {.Min = {.X = 5'000'000, .Y = 3'000'000},
                                   .Max = {.X = 5'000'050, .Y = 3'000'030}};
  const Geoxy::Tolerance far    = Geoxy::Tolerance::Of(survey);
  std::println("the same plan, far:  {}", far.Length);
}
```

```text [Output]
a 50 unit plan:      5e-09
the same plan, far:  1.4901161193847656e-08
```

A program usually has several: the model's, the pick aperture under the cursor, a command's own fuzz. That is why it is an argument and not a setting.

## Accuracy: how near an approximation comes {#accuracy}

An `Accuracy` holds an `Error`, the most an approximated answer may be off, in drawing units. Only operations that approximate take one: the chords `Flatten` lays along a curve, the length of a spline, the offset of an ellipse. Finer costs more work, within a budget that bounds it whatever is asked. An error at or below 0 asks for as near as can be.

What is worked out exactly takes none: the length of a line or an arc, the area of a polyline, where two lines meet.

```cpp [main.cpp]
#include "geoxy/geoxy.hpp"

#include <print>

auto main() -> int
{
  const Geoxy::Curve circle = Geoxy::Circle({.X = 0, .Y = 0}, 100);
  for (const double error : {1.0, 0.01, 0.0001})
  {
    const auto points = Geoxy::Flatten(circle, Geoxy::Accuracy{.Error = error});
    std::println("within {:>6}: {} points", error, points.size());
  }
}
```

```text [Output]
within      1: 24 points
within   0.01: 224 points
within  1e-04: 2223 points
```

## Which is which {#which}

A tolerance answers "are these the same?". An accuracy answers "how good must this be?". Passing a loose tolerance does not make anything faster, and passing a fine accuracy does not make two near points merge.
