cmake_minimum_required(VERSION 3.5)
project(pyslise)

set(CMAKE_CXX_FLAGS "-O2 ${CMAKE_CXX_FLAGS}")

find_package(Eigen3 REQUIRED NO_MODULE)

add_subdirectory(pybind11)

pybind11_add_module(pyslise SHARED pyslise.cpp)
if (MKL_FOUND)
    target_link_libraries(pyslise PUBLIC
            m
            Threads::Threads
            ${CMAKE_DL_LIBS}
            ${BLAS_LIBRARIES})
    add_definitions(-DEIGEN_USE_BLAS -DEIGEN_USE_MKL -DMKL_DIRECT_CALL)
endif ()
target_link_libraries(pyslise PRIVATE matslise Eigen3::Eigen)