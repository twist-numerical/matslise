cmake_minimum_required(VERSION 3.5)
project(schrodinger)

find_library(matslise REQUIRED)
set(CMAKE_CXX_STANDARD 11)

add_subdirectory (matslise)

add_executable(schrodinger main.cpp)
target_link_libraries(schrodinger LINK_PUBLIC matslise)


# add_subdirectory(swig/python_wrapper)

# add_subdirectory(python)
