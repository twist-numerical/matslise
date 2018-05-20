cmake_minimum_required(VERSION 3.5)
project(matslise)

set(CMAKE_BUILD_TYPE Release)
set(CMAKE_CXX_FLAGS "-O2 -Wall -Wextra -fPIC")
set(CMAKE_INSTALL_PREFIX /usr/)

find_package(Eigen3 REQUIRED NO_MODULE)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
add_library(matslise SHARED matslise.cpp legendre.h Array2D.h matscs.cpp)
target_link_libraries(matslise Eigen3::Eigen)

install(TARGETS matslise DESTINATION lib)
install(FILES matslise.h Array2D.h matscs.h DESTINATION include)