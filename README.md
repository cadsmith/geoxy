# geoxy

[![C++23](https://img.shields.io/badge/C%2B%2B-23-blue)](docs/cpp.md)
[![License: MIT](https://img.shields.io/badge/license-MIT-blue)](LICENSE)

The geometry kernel of a 2D CAD application: points, vectors and curves, and
what drawing, editing and snapping need of them. Where curves meet, the place
nearest a point, lengths and areas, offsets, fillets and chamfers, regions
and their unions, and indexes that keep a large drawing quick to ask. It
sits beside [dxfio](https://github.com/cadsmith/dxfio) and dwgio, which read
and write the drawings.

Answers are worked out to the last digit a double holds wherever that can be
done, nothing throws, and it depends on the standard library alone.

```cpp
#include "geoxy/geoxy.hpp"

const Geoxy::Curve     line      = Geoxy::Line{.Start = {.X = -10, .Y = 3}, .End = {.X = 10, .Y = 3}};
const Geoxy::Curve     circle    = Geoxy::Circle({.X = 0, .Y = 0}, 5);
const Geoxy::Tolerance tolerance = {.Length = 1e-9, .Radians = 1e-9};

for (const Geoxy::Hit &hit : Geoxy::Intersect(line, circle, tolerance).Hits)
{
  std::println("({}, {})", hit.At.X, hit.At.Y);   // (-4, 3) and (4, 3)
}
```

## Documentation

[Getting started](docs/guide/start.md) is the place to begin. The guide, the
tutorials and the concepts are in `docs/`, one page a file:

| | |
| :--- | :--- |
| Guide | [Getting started](docs/guide/start.md), [Installing](docs/guide/install.md), [Curves](docs/guide/curves.md), [Tolerance and accuracy](docs/guide/tolerance.md), [Errors](docs/guide/errors.md) |
| Tutorials | [Where curves meet](docs/tutorials/intersecting.md), [Snapping](docs/tutorials/snapping.md), [Measuring](docs/tutorials/measuring.md), [Editing curves](docs/tutorials/editing.md), [Offsetting](docs/tutorials/offsetting.md), [Regions](docs/tutorials/regions.md), [Moving and scaling](docs/tutorials/transforming.md), [Indexing a drawing](docs/tutorials/indexing.md), [Drawing curves](docs/tutorials/displaying.md) |
| Concepts | [Exactness](docs/concepts/exactness.md), [How it is tested](docs/concepts/quality.md), [Limits](docs/concepts/limits.md) |

Every program in those pages is built and run, and held to the output
printed under it:

```sh
python3 -I tools/docs/examples.py
```

## Build

C++23, CMake 3.28 or later and Ninja. Every step goes through a preset.

```sh
cmake --preset debug
cmake --build --preset debug
ctest --preset debug
```

The other presets, the targets and the coverage gate are in
[docs/cmake.md](docs/cmake.md).

## Use

As a subdirectory or an installed package, the target is the same:

```cmake
find_package(geoxy CONFIG REQUIRED)
target_link_libraries(app PRIVATE geoxy::geoxy)
```

## Layout

| Path | Holds |
| :--- | :--- |
| `include/geoxy/` | The public headers, each listed in `include/CMakeLists.txt` |
| `source/` | The library's sources and private headers |
| `cmake/` | The build's own modules: how geoxy is built, the developer's setup and the install rules |
| `vendor/` | Third party code, fetched at a pinned release for the tests and the benchmarks alone: doctest, Shewchuk's predicates as the tests' oracle, Google Benchmark |
| `tests/geoxy/` | The doctest suite, one file per unit, mirroring `source/` |
| `tests/suites/` | Cases of curves that meet, each answered exactly by `tools/suites`, and the program that scores them |
| `tests/proof/` | Proofs of the bounds the exact arithmetic rests on, checked as a test |
| `tests/fuzz/` | The libFuzzer targets and the inputs they found |
| `tests/bench/` | `geoxy_bench`, on Google Benchmark |
| `tests/api/` | The public interface's snapshot and the headers-alone check |
| `tests/package/` | A program built against the installed package |
| `tests/support/` | What the tests share: curves as text, to write a case down and read it back |
| `tools/suites/`, `tools/reference/` | The exact answers the suites and the tests are held to, worked out in Python |
| `tools/gate/` | Every check a change must pass, run with one command |
| `tools/format/` | The formatter: clang-format, then every comment framed and aligned |
| `tools/tidy/` | The script CI uses to clang-tidy only what a change reaches |
| `tools/docs/` | The check that builds and runs every program in the docs |
| `docs/` | The user's guide, tutorials and concepts; and the C++, CMake and commit guidelines |
| `site/` | The documentation website, VitePress over `docs/` |

## Public interface

`tests/api/snapshot.txt` holds the public headers with their comments
removed, and the `geoxy_api` test fails when they differ from it, so that a
change to the interface is always a deliberate one. When it is, take the new
snapshot:

```sh
cmake --build --preset debug --target geoxy_api
```

## License

[MIT](LICENSE)
