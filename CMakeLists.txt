cmake_minimum_required(VERSION 2.8)

project(mfg)

set(CMAKE_CXX_STANDARD 14)

set(HEADER_FILES ${CMAKE_CURRENT_SOURCE_DIR})

set(SOURCE_FILES solver1.cpp solver1.cpp malgo_tests.cpp mfg_tests.cpp catch.hpp inc/solver1.h inc/utils.h inc/malgo.h)

include_directories(inc)

add_executable(mfg ${SOURCE_FILES} ${HEADER_FILES})
