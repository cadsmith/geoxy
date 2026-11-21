cmake_minimum_required(VERSION 3.28)
project(geoxy_package LANGUAGES CXX)

set(GEOXY_PACKAGE_SOURCE "" CACHE PATH
  "geoxy's source, to add as a subdirectory; empty to find it installed")

if (GEOXY_PACKAGE_SOURCE)
  add_subdirectory(${GEOXY_PACKAGE_SOURCE} geoxy)
else()
  find_package(geoxy 0.1 CONFIG REQUIRED)
endif()

add_executable(geoxy_package main.cpp)

target_link_libraries(geoxy_package PRIVATE geoxy::geoxy)
