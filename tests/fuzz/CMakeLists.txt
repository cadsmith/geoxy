foreach(target batch circle orient two)
  add_executable(geoxy_fuzz_${target}
    ${target}.cpp
    draw.cpp
    oracle.cpp
  )

  target_link_libraries(geoxy_fuzz_${target}
    PRIVATE
      geoxy_internal
      geoxy_predicates
  )

  target_link_options(geoxy_fuzz_${target} PRIVATE -fsanitize=fuzzer)
endforeach()
