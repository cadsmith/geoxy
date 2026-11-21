cmake_minimum_required(VERSION 3.28)

project(geoxy
  VERSION 0.1.0
  DESCRIPTION "Exact and fast 2D geometry for CAD"
  HOMEPAGE_URL "https://github.com/cadsmith/geoxy"
  LANGUAGES CXX
)

option(GEOXY_INSTALL "Add geoxy's install rules" ${PROJECT_IS_TOP_LEVEL})
option(GEOXY_ENABLE_TRIGO
  "Take sines, arc tangents and lengths from source/trigo" OFF)
option(GEOXY_ENABLE_HARDENING
  "Check the standard library's preconditions" ON)
option(GEOXY_ENABLE_CCACHE "Compile through ccache where it is found" ON)
option(GEOXY_ENABLE_COVERAGE "Instrument for LLVM coverage" OFF)
option(GEOXY_ENABLE_SANITIZE "Instrument for ASan and UBSan" OFF)
option(GEOXY_ENABLE_FUZZ "Build the libFuzzer targets" OFF)
option(GEOXY_ENABLE_TIDY "Run clang-tidy on the library as it builds" OFF)
option(GEOXY_ENABLE_BENCH "Add geoxy_bench, on Google Benchmark" OFF)

set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_SCAN_FOR_MODULES OFF)

if (PROJECT_IS_TOP_LEVEL)
  include(cmake/developer.cmake)
endif()

include(cmake/build.cmake)

add_subdirectory(vendor SYSTEM)
add_subdirectory(source)
add_subdirectory(include)

if (GEOXY_INSTALL)
  include(cmake/install.cmake)
endif()

if (PROJECT_IS_TOP_LEVEL AND BUILD_TESTING)
  add_subdirectory(tests)
endif()
