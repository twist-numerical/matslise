cmake_minimum_required(VERSION 3.9)
project(matslise)

add_subdirectory(matslise)
add_subdirectory(test)

enable_testing()
add_test(NAME MatsliseTesting COMMAND matslise_test)