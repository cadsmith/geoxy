# CMake Guidelines

## Rules

- No `#` comments are allowed in any CMake file.
- Lines must not exceed 80 columns.
- The library target is `geoxy`, its alias `geoxy::geoxy`; other targets use
  the `geoxy_` prefix.
- Libraries must be placed in `lib/` and executables in `bin/`.
- No flag is set for a directory: `add_compile_options`,
  `add_compile_definitions` and `add_link_options` are not used. What
  every target of ours is built with is on `geoxy_build`, which each links.
- Every option starts `GEOXY_`, is declared in the top level
  `CMakeLists.txt`, and one that turns something on is `GEOXY_ENABLE_`.
- Nothing geoxy needs only to develop itself reaches a project that uses
  it: tests, vendored code, warnings, hardening and the output folders are
  for geoxy as the top level project alone.
- Configuration, building, and testing must use CMake presets.
- Every public header must be listed in `include/CMakeLists.txt`; configuring
  fails on one that is not.
- Every test file lives under `tests/geoxy/` and is listed in its
  `CMakeLists.txt`; CTest runs each of its test cases on its own, as
  `geoxy/<file>/<test case>`.

## Layout

| File | Holds |
| :--- | :--- |
| `CMakeLists.txt` | The project, every option, and the order of the rest |
| `cmake/build.cmake` | `geoxy_build`: C++23 and each product and sum rounded on its own, always; the warnings, hardening, coverage, sanitizers and fuzzing where geoxy is the top level project |
| `cmake/developer.cmake` | What geoxy as the top level project sets up: the output folders, ccache, CTest and `geoxy_docs` |
| `cmake/install.cmake`, `cmake/config.cmake.in` | The install rules and the package, where `GEOXY_INSTALL` is on |
| `cmake/cases.cmake`, `cmake/stamp.cmake` | Test cases listed for CTest; work kept until what it reads changes |
| `source/`, `include/` | The library's sources, and its public headers as a file set |
| `tests/` | Every test target, each linking `geoxy_internal`: the library, `geoxy_build` and the private headers |
| `vendor/` | Third party code, fetched for the tests and the benchmarks alone |

## Options

| Option | Default | Does |
| :--- | :--- | :--- |
| `GEOXY_INSTALL` | on at the top level | Adds the install rules and the package |
| `GEOXY_ENABLE_TRIGO` | off | Takes sines, arc tangents and lengths from `source/trigo` |
| `GEOXY_ENABLE_HARDENING` | on | Checks the standard library's preconditions |
| `GEOXY_ENABLE_CCACHE` | on | Compiles through ccache where it is found |
| `GEOXY_ENABLE_COVERAGE` | off | Instruments for LLVM coverage |
| `GEOXY_ENABLE_SANITIZE` | off | Instruments for ASan and UBSan |
| `GEOXY_ENABLE_FUZZ` | off | Builds the libFuzzer targets; needs LLVM's clang |
| `GEOXY_ENABLE_TIDY` | off | Runs clang-tidy on the library as it builds |
| `GEOXY_ENABLE_BENCH` | off | Adds `geoxy_bench` |

Only the first two reach a project that uses geoxy; the rest are read
where geoxy is the top level project.

## Presets

| Name | Type | Purpose |
| :--- | :--- | :--- |
| `debug` | `Debug` | Development with symbols and assertions; adds `geoxy_bench`, unbuilt, so that clangd sees it |
| `release` | `Release` | Optimized build |
| `coverage` | `Debug` at `-O1` | Instrumented build for LLVM coverage |
| `sanitize` | `Debug` at `-O1` | ASan and UBSan instrumented build, float division by zero included |
| `tidy` | `Debug` | The library built through clang-tidy |
| `bench` | `Release` | `geoxy_bench` on Google Benchmark, `GEOXY_ENABLE_HARDENING` off |
| `fuzz` | `RelWithDebInfo` | libFuzzer targets with ASan and UBSan; configured with LLVM's clang, as below |
| `unity` | `Debug` | The library and the tests built as a unity build, eight files compiled as one, so that no two files' helpers may share a name |

Every build when geoxy is the top level project defines
`_GLIBCXX_ASSERTIONS` and a libc++ hardening mode, so the standard library
checks its preconditions: the fast mode in `Release` and `MinSizeRel`,
bounds for a few percent, and the debug mode in every other build, which
checks each sort's ordering too. `GEOXY_ENABLE_HARDENING` off, as `bench`
turns it, defines neither, to time what a caller runs. Warnings are
errors, `-Wfloat-equal` and `-Wdouble-promotion` among them: floating
point values compare through a tolerance, never with `==`.

Every build rounds each product and sum on its own: `-ffp-contract=off`
with GCC and Clang, `/fp:precise` with MSVC, and no `-ffast-math`. Clang
fuses `a * b - c * d` by default on arm64, which breaks the exact
predicates and makes arm64 and x86-64 disagree. `tests/geoxy/bits.cpp` holds
a checksum of a million seeded operations that must be one number on every
platform; it changes when contraction comes back. A caller's build compiles
what a public header computes inline with the caller's flags, so a public
header keeps its floating point arithmetic out of line.

`GEOXY_ENABLE_TRIGO`, off by default, takes every sine, cosine, tangent, arc
function, length, logarithm, power and cube root the library calls outside
its correctly rounded `Direction` and `Heading` from `source/trigo`, geoxy's
own port of trigo, instead of the C library: within 0.56 steps of a double
for the sines, degrees reduced exactly, and the same bits on every machine,
arm64 and x86-64, with a fused multiply-add in its hardware or without,
where the C libraries differ. `tests/geoxy/trigo.cpp` holds a checksum of
its functions that must be one number on every platform. It stays off while
it costs speed: intersections of arcs are up to 9% faster with it, but an ellipse's
offset, flattening and extents and a polyline's flattening are 8 to 10%
slower, their chains of tangents and arc tangents slower than Apple's.

## Targets

| Target | Alias | Type | Description |
| :--- | :--- | :--- | :--- |
| `geoxy` | `geoxy::geoxy` | Static, or shared with `BUILD_SHARED_LIBS` | The geometry kernel |
| `geoxy_build` | - | Interface | What every target of ours is built with; never installed |
| `geoxy_internal` | - | Interface, tests | The library with `geoxy_build` and its private headers, for the tests and the benchmarks |
| `geoxy_tests` | - | Executable | Doctest test suite, doctest and the standard headers precompiled; CTest runs it once per test case, as `geoxy/<file>/<test case>` |
| `geoxy_tests_cases` | - | Custom | Lists the test cases of `geoxy_tests` for CTest, each time it links |
| `geoxy_suites`, `geoxy_scoreboard` | - | Executable | The cases of `tests/suites` run against their exact answers, and the count of those answered to the bit, held to the count kept |
| `geoxy_support` | - | Static, tests | What the tests share: curves written as text and read back |
| `geoxy_doctest` | - | Object | Doctest's runner, `tests/main.cpp`, compiled once for `geoxy_tests` and `geoxy_suites` |
| `geoxy_headers` | - | Object | Each public header compiled on its own, twice |
| `geoxy_api` | - | Custom, test | Takes the public interface's snapshot; the test compares it |
| `coverage` | - | Custom | Runs the test cases eight at a time and gates lines, functions, regions and branches at 100% |
| `geoxy_proofs` | - | Custom, test | Brings each proof of `tests/proof` up to date: run again only when what it proves from changes, once for every preset |
| `geoxy_predicates` | - | Static, tests and bench | Shewchuk's `predicates.c`, fetched: the oracle whose signs `source/exact` is checked against, and the reference its benchmarks run beside |
| `geoxy_bench` | - | Executable, on request | Google Benchmark of the library: the exact predicates, `Intersect` per pair kind, as given and prepared, and `Nearest`, `Length`, `Area`, `Extents` and `Flatten` per curve kind |
| `geoxy_fuzz_orient`, `geoxy_fuzz_circle`, `geoxy_fuzz_batch`, `geoxy_fuzz_two` | - | Executable, `fuzz` | libFuzzer targets of `source/exact` against `geoxy_predicates` |
| `geoxy_docs` | - | Custom, with Doxygen | The API's documentation in `build/docs`, failing on any warning |

## Install

`cmake --install` puts the library, its headers and a package under a
prefix, where `find_package(geoxy CONFIG)` finds the `geoxy::geoxy` target
with the same usage requirements the subdirectory gives. The library
links the platform's threads, `Threads::Threads`, for the pools a grid
is intersected on; the package finds them for the program that uses it.
`tests/package` is a program built that way, by CI against a fresh
install; given `-DGEOXY_PACKAGE_SOURCE=<geoxy>` it adds geoxy as a
subdirectory instead, the other way a project uses it. The package is
compatible within one minor version until 1.0.

## Commands

```sh
cmake --preset debug
cmake --build --preset debug
ctest --preset debug
cmake --workflow --preset debug
ctest --preset debug -R geoxy/version
./build/debug/bin/geoxy_tests -tc="version*"

cmake --preset coverage
cmake --build --preset coverage --target coverage

cmake --preset sanitize
cmake --build --preset sanitize
ctest --preset sanitize

cmake --preset tidy
cmake --build --preset tidy

cmake --build --preset debug --target geoxy_api

cmake --preset unity
cmake --build --preset unity
cmake --preset unity -B build/one -DCMAKE_UNITY_BUILD_BATCH_SIZE=0
cmake --build build/one --target geoxy geoxy_tests

cmake --preset bench
cmake --build --preset bench
./build/bench/bin/geoxy_bench
./build/bench/bin/geoxy_bench --benchmark_filter='ByIntersect|ByPrepared'

cmake --preset fuzz -DCMAKE_CXX_COMPILER=/opt/homebrew/opt/llvm/bin/clang++ \
  -DCMAKE_C_COMPILER=/opt/homebrew/opt/llvm/bin/clang
cmake --build --preset fuzz
mkdir -p build/fuzz/corpus/orient
./build/fuzz/bin/geoxy_fuzz_orient -max_total_time=600 -timeout=2 \
  -artifact_prefix=tests/fuzz/found/orient/ \
  build/fuzz/corpus/orient tests/fuzz/corpus/orient

python3 -I tools/gate/run.py source/intersect.cpp tests/geoxy/intersect.cpp
```

The last runs every check a change must pass, side by side: the debug,
release and sanitize presets built and tested, coverage gated, the library
linted by the `tidy` preset, `geoxy_bench` built by the `bench` preset, the
library and the tests built as a unity build by the `unity` preset,
configured first where it is not, the documentation and the interface's
snapshot built, and the files named checked for their formatting and for what
an editor running clangd shows. Under twenty seconds where a file or two
changed.

A unity build compiles a file's anonymous helpers with the next files', so
two helpers of one name and kind in one namespace, a helper that shadows or
hides another file's name, or a using-directive at namespace scope breaks it:
every helper of the library has a name of its own, and so does every test's
that another test or a header of the library shared. The `unity` preset
compiles eight files as one, which a change builds again in a few seconds;
`CMAKE_UNITY_BUILD_BATCH_SIZE` at 0 compiles each target as one file, the
strict check, which fails only on the names `abut.cpp`, `graze.cpp`,
`prepared.cpp` and their tests still share.

Every build is incremental: Ninja compiles again what a change reaches and
nothing else, and so with the checks that are steps of a build, the `tidy`
preset's clang-tidy and the proofs of `tests/proof` among them. Where
`ccache` is installed and geoxy is the top level project, every compile goes
through it, so a file compiled before with the same flags, as after a change
taken back, is not compiled again; `GEOXY_ENABLE_CCACHE` off, or a
`CMAKE_CXX_COMPILER_LAUNCHER` of one's own, says otherwise. ccache 4.8 or
later is told that the precompiled header's macros and times are not its
input, so it caches what is compiled with it too.

`geoxy_tests` precompiles doctest and the standard headers the tests include
most, and only those: precompiling geoxy's own headers too puts the whole
interface in scope of every test, where a test includes a part of it, and
`Length` turns ambiguous in `tests/geoxy/vector.cpp`. Doctest's runner,
`tests/main.cpp`, is compiled once, as `geoxy_doctest`, for `geoxy_tests` and
`geoxy_suites` both.

CTest runs each test case of `geoxy_tests` as a test of its own, named
`geoxy/<file>/<test case>`, so that the costliest file is no longer one long
process: `geoxy_tests_cases` lists them as `geoxy_tests` links, by
`cmake/cases.cmake`, and the test folder includes that list. A test preset
runs eight tests at a time; `--parallel` on the command line says otherwise.
The longest tests start first: `suites` and `geoxy_proofs` say they are costly,
and CTest learns what the others cost from the runs before. The `coverage`
target runs the test cases the same way, eight at a time, each merging what it
covered into one of eight raw profiles.

A randomized test runs all its rounds in a build that is optimized, as every
preset's but `debug` is, and one in sixteen of them in `debug`, whose rounds
cost seventeen times as much: `tests/geoxy/draws.hpp` counts them. A test that
bounds how long a call takes counts the processor time its process spends,
`std::clock`, after a first call, and not the wall clock's: gates running side
by side load the machine past its cores, and a process pays once for the
memory it touches first.

The checks of `tools/suites`, `suites_oracle` and `suites_forged`, are
Python and never touch the C++ build; `suites_forged` runs on macOS alone,
where the cases were forged: on Linux the forge makes cases that differ,
for a reason not yet found. Each runs through
`cmake/stamp.cmake`, which stamps a check that passes in `build/stamps`, the
folder every preset shares (`GEOXY_STAMPS`), by a hash of what it reads:
`tools/suites/*.py`, its cases, and Python's path and date. A check whose
stamp is there passes at once; so each runs once for every preset, and again
only when what it reads changes. A preset that comes to a check while another
runs it waits for it.

Each fuzz target decodes its bytes into the inputs the exact predicates find
hardest (`tests/fuzz/draw.hpp`): points on a line or a circle nudged by a few
steps of a double, coordinates differing in their low bits, mixed magnitudes,
survey-scale offsets, and both sides of Shewchuk's limits. Every sign must be
the sign of Shewchuk's `predicates.c`, on the places scaled by a power of two
into his limits where they lie outside, and keep the laws of its determinant;
a NaN or infinite coordinate must give zero. A file in
`tests/fuzz/corpus/<target>` or `tests/fuzz/found/<target>` runs again in
`geoxy/exact/fuzz`: an input the fuzzer finds lands in `found/<target>`, and
stays there once the predicate is fixed.

`tests/proof` checks the exact predicates' arithmetic by machine, down to
the hardware's rounding and up to the compiler's code:

- `rounding.v` proves in Rocq, with Flocq, every fact of IEEE 754 binary64
  the other proofs use: the relative and the underflowing error of a
  rounding, its monotony and symmetry, `two.hpp`'s TwoSum, TwoDiff, fused
  and split TwoProduct exact, each in the code's own order, the orientation
  test's one fused comparison, and IEEE 754's operations that rounding.
- `integers.v` proves in Rocq the whole numbers `grid.cpp` and `whole.cpp`
  sum the last stages in, step for step in C++'s unsigned arithmetic: the
  carries, the products of words, the signed sums and the grids' two's
  complement, a double's bits read as Flocq reads them, and each sign the
  sign of the determinant of the doubles. `integers.py` runs the code
  itself, built by Clang and by GCC, with 128 bit products and without,
  against exact rational arithmetic on thousands of hard cases, exact zeros
  among them.
- `rocq.py` builds Flocq, fetched by `vendor/flocq`, then both proofs, and
  fails on any admitted step or any axiom but Rocq's own for the reals.
- `circle.py` and `first.py` prove with Z3, from those facts, the circle
  test's compensated second stage and both first stages, through underflow,
  overflow and NaN, down to the signs their bounds decide. Each run also
  takes a guard out or shrinks a bound, and Z3 must find the counterexample.
- `program.py` writes each proved function once as a graph of its
  operations. Every rounding a lemma states must be the graph's, and
  `translation.py` reads the compiler's LLVM IR of the code back into
  graphs that must equal them, constant for constant: the first stages,
  inline and in the batch loop, the circle test's later stages, and the
  orientation test's fused comparison with its guards.
- `convert.py` checks the two-double constants `Angle` converts degrees
  and radians with: that one multiplication gives the nearest double,
  always.
- `machine.py` disassembles what Apple's and LLVM's Clang emit for arm64
  and x86-64, and Homebrew's GCC for the machine it runs on, and finds no
  fused multiply-add the source did not ask for and no x87 instruction,
  whatever the build's flags, since
  `source/exact/strict.hpp` turns fusing off and stops fast math.

Each is a step of the build with a stamp of its own under `tests/proof` in
the build folder, made again when a file of `tests/proof` or `source/exact`,
a header they read, `source/angle.cpp`, `source/degrees.hpp` or the tool that
checks it changes, and not otherwise. A proof that passes is stamped in
`build/stamps` too, as the Python checks are, by a hash of those files and of
the tools' paths and dates, so a fresh tree proves once for every preset, not
once a preset. Where a proof runs: `rocq` wherever `rocq` is found,
`circle` and `first` wherever `z3` is, `convert` everywhere, `integers`
with Clang, `translation` with Clang on arm64, whose fused multiply-add
its graphs are written for, and `machine` with Clang on macOS. CTest runs one
test for them all, `geoxy_proofs`, which builds the target of that name: it
proves what is out of date, a few minutes in a fresh tree, and nothing where
nothing they prove from has changed. Deleting `build/stamps` and building the
target `clean` take the stamps away, to prove everything again, as after a
compiler changes behind a path and a date that do not.
