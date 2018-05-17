cmake_minimum_required(VERSION 3.5)
project(pyslise)

add_subdirectory(pybind11)

pybind11_add_module(pyslise SHARED pyslise.cpp)
target_link_libraries(pyslise PRIVATE matslise)