# Installing geoxy

geoxy builds with CMake 3.28 or later and a C++23 compiler, and depends on nothing but the standard library and the system's threads.

## As a subdirectory {#subdirectory}

With geoxy's source inside your project:

```cmake [CMakeLists.txt]
add_subdirectory(geoxy)
target_link_libraries(app PRIVATE geoxy::geoxy)
```

Only the library is built this way. The tests, the benchmarks, their dependencies, geoxy's own warning flags and its install rules are for geoxy as the top level project alone; set `GEOXY_INSTALL` on to have its install rules in your project too.

## With FetchContent {#fetch}

```cmake [CMakeLists.txt]
include(FetchContent)
FetchContent_Declare(geoxy
  GIT_REPOSITORY https://github.com/cadsmith/geoxy.git
  GIT_TAG        main
)
FetchContent_MakeAvailable(geoxy)
target_link_libraries(app PRIVATE geoxy::geoxy)
```

Pin `GIT_TAG` to a commit for a build that does not change under you.

## As an installed package {#package}

Build and install geoxy once, then find it:

```sh [Shell]
cmake --preset release
cmake --build --preset release
cmake --install build/release --prefix /usr/local
```

```cmake [CMakeLists.txt]
find_package(geoxy CONFIG REQUIRED)
target_link_libraries(app PRIVATE geoxy::geoxy)
```

## Static or shared {#linkage}

geoxy follows `BUILD_SHARED_LIBS`: a static library by default, a shared one with `-DBUILD_SHARED_LIBS=ON`. The target carries what the headers need either way, so nothing changes in your code.

## One header or several {#headers}

`geoxy/geoxy.hpp` includes everything. Each part also stands alone, and a file that only meets curves can include `geoxy/intersect.hpp` and nothing else; every public header compiles on its own.

## Floating point {#floating}

geoxy is built with `-ffp-contract=off`, so that each product and sum is rounded on its own and the same input gives the same answer on every machine. Your own code needs no such flag to call it. Do not build geoxy with `-ffast-math`: its exact arithmetic relies on IEEE 754 as written.
