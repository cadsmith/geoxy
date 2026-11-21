# Where curves meet

`Intersect` finds every point two curves share, and every stretch they run along together. It takes any two curves of any kinds.

## Hits {#hits}

An `Intersections` holds the `Hits`, in order along the first curve. Each hit has the point, the place on each curve, and how they meet there: `Cross` where they pass through each other, `Touch` where they are tangent, `End` where the meeting is at an end of one of them.

```cpp [main.cpp]
#include "geoxy/geoxy.hpp"

#include <print>

auto main() -> int
{
  const Geoxy::Tolerance tolerance = {.Length = 1e-9, .Radians = 1e-9};
  const Geoxy::Curve     circle    = Geoxy::Circle({.X = 0, .Y = 0}, 5);
  const Geoxy::Curve     chord     = Geoxy::Line{.Start = {.X = -10, .Y = 3}, .End = {.X = 10, .Y = 3}};
  const Geoxy::Curve     tangent   = Geoxy::Line{.Start = {.X = -10, .Y = 5}, .End = {.X = 10, .Y = 5}};
  const Geoxy::Curve     miss      = Geoxy::Line{.Start = {.X = -10, .Y = 7}, .End = {.X = 10, .Y = 7}};

  const auto show = [&](const char *name, const Geoxy::Curve &line)
  {
    const Geoxy::Intersections found = Geoxy::Intersect(line, circle, tolerance);
    std::println("{}: {} hits", name, found.Hits.size());
    for (const Geoxy::Hit &hit : found.Hits)
    {
      const char *kind = hit.Kind == Geoxy::Contact::Cross ? "cross" : hit.Kind == Geoxy::Contact::Touch ? "touch" : "end";
      std::println("  ({}, {}) {}", hit.At.X, hit.At.Y, kind);
    }
  };
  show("chord", chord);
  show("tangent", tangent);
  show("miss", miss);
}
```

```text [Output]
chord: 2 hits
  (-4, 3) cross
  (4, 3) cross
tangent: 1 hits
  (0, 5) touch
miss: 0 hits
```

## Overlaps {#overlaps}

Where two curves lie on each other along a stretch, there is no single point to give. `Intersections::Overlaps` holds each shared stretch as a range of each curve's parameter, and whether the second runs the same way as the first.

```cpp [main.cpp]
#include "geoxy/geoxy.hpp"

#include <print>

auto main() -> int
{
  const Geoxy::Tolerance tolerance = {.Length = 1e-9, .Radians = 1e-9};
  const Geoxy::Curve     first     = Geoxy::Line{.Start = {.X = 0, .Y = 0}, .End = {.X = 10, .Y = 0}};
  const Geoxy::Curve     second    = Geoxy::Line{.Start = {.X = 15, .Y = 0}, .End = {.X = 5, .Y = 0}};

  const Geoxy::Intersections found = Geoxy::Intersect(first, second, tolerance);
  for (const Geoxy::Overlap &overlap : found.Overlaps)
  {
    std::println("first from {} to {}, second from {} to {}, {}", overlap.First.Start, overlap.First.End,
                 overlap.Second.Start, overlap.Second.End, overlap.Along ? "the same way" : "the other way");
  }
}
```

```text [Output]
first from 0.5 to 1, second from 0.5 to 1, the other way
```

## A curve with itself {#self}

`SelfIntersect` finds where one curve crosses itself: a polyline that loops over its own path, a spline with a loop.

## Many curves at once {#many}

Intersecting every pair of a drawing's curves one by one is slow. A `Grid` prepares the curves once and intersects only those whose boxes come near; [Indexing a drawing](/tutorials/indexing/) shows it.

## What the tolerance does {#tolerance}

Two curves that pass within `tolerance.Length` of each other without crossing meet in a `Touch`. Everything else is decided from the curves' own numbers: between lines and arcs each hit is the double nearest the true point, whatever the tolerance.
