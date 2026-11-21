enable_language(C)

FetchContent_Declare(
  predicates
  URL https://www.cs.cmu.edu/afs/cs/project/quake/public/code/predicates.c
  URL_HASH
    SHA256=f8662c3f407d1c1c5dcd4dd49ea8b8ddd801a71827d65734206ddc3741792029
  DOWNLOAD_NO_EXTRACT TRUE
)

FetchContent_MakeAvailable(predicates)

add_library(geoxy_predicates STATIC ${predicates_SOURCE_DIR}/predicates.c)

set_target_properties(geoxy_predicates PROPERTIES
  C_STANDARD 11
  C_EXTENSIONS ON
  COMPILE_WARNING_AS_ERROR OFF
)

if (MSVC)
  set(shim ${CMAKE_CURRENT_BINARY_DIR}/shim)
  file(WRITE ${shim}/sys/time.h "")
  target_include_directories(geoxy_predicates PRIVATE ${shim})
  target_compile_definitions(geoxy_predicates PRIVATE random=rand)
  target_compile_options(geoxy_predicates PRIVATE /fp:precise /w)
else()
  target_compile_options(geoxy_predicates PRIVATE -ffp-contract=off -w)
endif()
