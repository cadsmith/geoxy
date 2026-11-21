set(BENCHMARK_ENABLE_TESTING OFF CACHE BOOL "")
set(BENCHMARK_ENABLE_GTEST_TESTS OFF CACHE BOOL "")
set(BENCHMARK_ENABLE_INSTALL OFF CACHE BOOL "")
set(BENCHMARK_ENABLE_WERROR OFF CACHE BOOL "")
set(BENCHMARK_INSTALL_DOCS OFF CACHE BOOL "")

FetchContent_Declare(
  benchmark
  GIT_REPOSITORY https://github.com/google/benchmark.git
  GIT_TAG        v1.9.4
  GIT_SHALLOW    TRUE
  GIT_PROGRESS   TRUE
  SYSTEM
  FIND_PACKAGE_ARGS 1.9.4 GLOBAL
)

FetchContent_MakeAvailable(benchmark)

if (benchmark_SOURCE_DIR)
  set_target_properties(benchmark benchmark_main PROPERTIES
    COMPILE_WARNING_AS_ERROR OFF
  )
  if (NOT MSVC)
    target_compile_options(benchmark PRIVATE -w)
    target_compile_options(benchmark_main PRIVATE -w)
  endif()
endif()
