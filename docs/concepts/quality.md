# How it is tested

A geometry library is trusted for its answers, so most of the work in geoxy is checking them.

## Answers checked against the truth {#truth}

A test that compares a result with a number typed by hand only checks that the code still does what it did. geoxy's tests compare against answers computed another way, to many more digits than a double holds: exact rational arithmetic, arithmetic at hundreds of bits, and for offsets a separate program that works out the true offset to 120 digits. A result is then judged by how many steps of a double it is from the truth.

## Every line run {#coverage}

Every change must leave the test suite running 100% of the library's lines, branches and regions. A branch no test can reach is removed, not excused.

## The gate {#gate}

No change is kept unless all of this passes:

- the tests, in debug, release and under the address and undefined behavior sanitizers, with division by zero in floating point trapped;
- full coverage, as above;
- clang-tidy with warnings as errors, including a limit on how complicated one function may be;
- a unity build, the formatter, and a check that every public header compiles on its own;
- the public interface compared with a recorded snapshot, so that it never changes by accident.

## Random shapes {#random}

Beyond fixed cases, shapes are generated at random, at small, ordinary and very large coordinates, and every answer is compared with the truth or with asking every curve one by one. The indexes are tested this way to the bit: a `Tree` must return exactly what a loop over all the curves returns.
