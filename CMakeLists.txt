cmake_minimum_required(VERSION 3.5)
project(matslise)

add_subdirectory(main)
add_subdirectory(src)
add_subdirectory(test)

enable_testing ()
add_test (NAME matslise_test COMMAND Test)