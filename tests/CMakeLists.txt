set(CMAKE_CXX_STANDARD 14)

set(CPP_SOURCES malgo_tests.cpp mfg_tests.cpp catch.hpp)
add_executable(run_tests ${CPP_SOURCES})

target_link_libraries(run_tests mfg)