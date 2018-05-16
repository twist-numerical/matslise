cmake_minimum_required(VERSION 3.9)
project(schrodinger)

set(CMAKE_CXX_STANDARD 11)

add_executable(schrodinger main.cpp matslise/matslise.cpp matslise/matslise.h matslise/legendre.cpp matslise/legendre.h matslise/Array2D.h)