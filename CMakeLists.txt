cmake_minimum_required(VERSION 3.5)
project(schrodinger)

set(CMAKE_CXX_STANDARD 11)

find_package(Eigen3 REQUIRED NO_MODULE)
add_subdirectory (matslise)

add_executable(schrodinger main.cpp matslise/matslise.cpp matslise/matscs.cpp)
target_link_libraries(schrodinger Eigen3::Eigen)


# add_subdirectory(swig/python_wrapper)

# add_subdirectory(python)
