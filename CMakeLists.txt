cmake_minimum_required(VERSION 3.1)
project(matslise)

if (NOT DEFINED EMSCRIPTEN)
    add_subdirectory(cmake/pybind11)
endif ()
add_subdirectory(main)
add_subdirectory(src)
add_subdirectory(test)

enable_testing()
add_test(NAME MatsliseTesting COMMAND matslise_test)