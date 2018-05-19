cmake_minimum_required(VERSION 3.5)
project(matslise)

set(CMAKE_BUILD_TYPE Release)
set(CMAKE_CXX_FLAGS "-Wall -Wextra -fPIC")
set(CMAKE_INSTALL_PREFIX /usr/)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
add_library(matslise SHARED matslise.cpp legendre.cpp Array2D.h)


install(TARGETS matslise DESTINATION lib)
install(FILES matslise.h Array2D.h DESTINATION include)