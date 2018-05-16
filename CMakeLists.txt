cmake_minimum_required(VERSION 3.5)
project(schrodinger)

set(CMAKE_CXX_FLAGS "-O3 -Wall -Wextra -fPIC")
set(CMAKE_CXX_STANDARD 11)

add_subdirectory (matslise)

add_executable(schrodinger main.cpp)
target_link_libraries(schrodinger LINK_PUBLIC matslise)


add_subdirectory(swig/python_wrapper)