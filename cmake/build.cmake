add_library(geoxy_build INTERFACE)

target_compile_features(geoxy_build INTERFACE cxx_std_23)

if (MSVC)
  target_compile_options(geoxy_build INTERFACE /utf-8 /fp:strict)
else()
  target_compile_options(geoxy_build INTERFACE -ffp-contract=off)
endif()

if (NOT PROJECT_IS_TOP_LEVEL)
  return()
endif()

if (MSVC)
  target_compile_options(geoxy_build INTERFACE /W4)
else()
  target_compile_options(geoxy_build
    INTERFACE
      -Wall -Wextra -Wpedantic -Wshadow
      -Wconversion -Wsign-conversion -Wold-style-cast
      -Wnon-virtual-dtor -Woverloaded-virtual
      -Wfloat-equal -Wdouble-promotion
      $<$<CXX_COMPILER_ID:Clang,AppleClang>:-Wnull-dereference>
      $<$<CXX_COMPILER_ID:GNU>:-Wno-maybe-uninitialized>
  )
endif()

if (GEOXY_ENABLE_HARDENING)
  set(hardening
    $<IF:$<CONFIG:Release,MinSizeRel>,FAST,DEBUG>)
  target_compile_definitions(geoxy_build
    INTERFACE
      _GLIBCXX_ASSERTIONS
      _LIBCPP_HARDENING_MODE=_LIBCPP_HARDENING_MODE_${hardening}
  )
endif()

if (GEOXY_ENABLE_COVERAGE)
  if (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    target_compile_options(geoxy_build
      INTERFACE
        -fprofile-instr-generate
        -fcoverage-mapping
        -fprofile-update=atomic
    )
    target_link_options(geoxy_build INTERFACE -fprofile-instr-generate)
  elseif (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    target_compile_options(geoxy_build
      INTERFACE --coverage -fprofile-update=atomic)
    target_link_options(geoxy_build INTERFACE --coverage)
  endif()
endif()

if (GEOXY_ENABLE_SANITIZE)
  set(sanitizers address,undefined,float-divide-by-zero)
  if (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    string(APPEND sanitizers ,implicit-conversion,local-bounds,nullability)
    target_compile_options(geoxy_build
      INTERFACE
        -fsanitize-ignorelist=${PROJECT_SOURCE_DIR}/tests/ignorelist.txt
    )
  endif()
  target_compile_options(geoxy_build
    INTERFACE
      -fsanitize=${sanitizers}
      -fno-sanitize-recover=all
      -fno-omit-frame-pointer
  )
  target_link_options(geoxy_build INTERFACE -fsanitize=${sanitizers})
endif()

if (GEOXY_ENABLE_FUZZ)
  if (NOT CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    message(FATAL_ERROR
      "GEOXY_ENABLE_FUZZ needs LLVM's clang, which has libFuzzer: "
      "configure with -DCMAKE_CXX_COMPILER=<llvm>/bin/clang++ "
      "and -DCMAKE_C_COMPILER=<llvm>/bin/clang")
  endif()
  target_compile_options(geoxy_build INTERFACE -fsanitize=fuzzer-no-link)
endif()
