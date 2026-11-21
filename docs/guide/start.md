# Getting started

geoxy is a C++23 library of 2D geometry for CAD: points, vectors and curves, and what drawing, editing and snapping need of them. It finds where curves meet, the place on a curve nearest a point, lengths and areas, offsets, fillets and chamfers, the regions curves bound and their unions, and it indexes curves and points so that a large drawing is asked quickly. Its answers are worked out to the last digit a double holds wherever that can be done, and it never throws.

## Add it to a project {#add}

geoxy is a CMake project that depends only on the standard library and needs C++23. With its source in your project, add it and link the target `geoxy::geoxy`:

```cmake [CMakeLists.txt]
add_subdirectory(geoxy)
target_link_libraries(app PRIVATE geoxy::geoxy)
```

[Installing geoxy](/guide/install/) gives the compilers it needs and the other ways to add it.

## A first program {#first}

`geoxy/geoxy.hpp` includes the whole interface, in the namespace `Geoxy`. This program meets a line with a circle and prints where they cross.

```cpp [main.cpp]
#include "geoxy/geoxy.hpp"

#include <print>

auto main() -> int
{
  const Geoxy::Curve     line      = Geoxy::Line{.Start = {.X = -10, .Y = 3}, .End = {.X = 10, .Y = 3}};
  const Geoxy::Curve     circle    = Geoxy::Circle({.X = 0, .Y = 0}, 5);
  const Geoxy::Tolerance tolerance = {.Length = 1e-9, .Radians = 1e-9};

  const Geoxy::Intersections found = Geoxy::Intersect(line, circle, tolerance);
  for (const Geoxy::Hit &hit : found.Hits)
  {
    std::println("({}, {}) at {} along the line", hit.At.X, hit.At.Y, hit.First.Along);
  }
}
```

```text [Output]
(-4, 3) at 0.3 along the line
(4, 3) at 0.7 along the line
```

## How the interface reads {#shape}

Three things hold everywhere:

- **Curves are plain values.** A `Line` is two points, an `Arc` a center, a radius and two angles. `Curve` is a `std::variant` of the nine kinds, and every operation takes either one kind or a `Curve`. [Curves](/guide/curves/) lists them.
- **Operations are free functions.** `Intersect`, `Nearest`, `Length`, `Offset`, `Fillet`: a verb, the curves, and what the verb needs. Nothing is a method of a curve.
- **Nothing throws.** An operation that can fail returns a `std::expected` or a `std::optional`, and says why. [Errors](/guide/errors/) lists the reasons.

Two small types say how to answer. A `Tolerance` says when two points count as one; an `Accuracy` says how near an approximation must come. [Tolerance and accuracy](/guide/tolerance/) explains both.

## Where to go next {#next}

The guide:

- [Installing geoxy](/guide/install/): compilers, CMake and the ways to link it.
- [Curves](/guide/curves/): the nine kinds, their parameters, and places on them.
- [Tolerance and accuracy](/guide/tolerance/): what each call asks for and why.
- [Errors](/guide/errors/): how a call says it has no answer.

The tutorials, one task each:

- [Where curves meet](/tutorials/intersecting/) and [Snapping](/tutorials/snapping/): intersections, nearest points, tangents.
- [Measuring](/tutorials/measuring/): length, area, extents, centroids.
- [Editing curves](/tutorials/editing/): trim, extend, break, join, fillet, chamfer.
- [Offsetting](/tutorials/offsetting/): one side of a curve, or a shape grown and shrunk.
- [Regions](/tutorials/regions/): the area a pick lies in, and unions, differences and intersections.
- [Moving and scaling](/tutorials/transforming/): matrices applied to curves.
- [Indexing a drawing](/tutorials/indexing/): picking and snapping among many curves.
- [Drawing curves](/tutorials/displaying/): chords for the screen, triangles to fill, dashes.

The concepts explain the model behind the code: [Exactness](/concepts/exactness/), [How it is tested](/concepts/quality/) and [Limits](/concepts/limits/).
