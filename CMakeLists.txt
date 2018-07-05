cmake_minimum_required(VERSION 3.5)
project(matslise)

set(CMAKE_CXX_FLAGS "-O2 -Wall -Wextra -fPIC -fno-var-tracking-assignments")
set(CMAKE_CXX_STANDARD 11)

cmake_minimum_required(VERSION 3.5)

set(CMAKE_INSTALL_PREFIX /usr/)

find_package(Eigen3 REQUIRED NO_MODULE)

SET(MATSLISE_SRC src/util/calculateEta.cpp
        src/matslise/matslise.cpp src/matslise/halfrange.cpp src/matslise/sector.cpp
        src/matscs/matscs_formulas.cpp src/matscs/matscs.cpp src/matscs/sector.cpp
        src/se2d/se2d.cpp src/Matrix2D.h)

add_library(matslise SHARED ${MATSLISE_SRC})
target_link_libraries(matslise Eigen3::Eigen)

install(TARGETS matslise DESTINATION lib)
install(FILES src/matslise.h src/Array2D.h src/matscs.h src/se2d.h src/Matrix2D.h DESTINATION include/matslise/)


add_executable(schrodinger src/main.cpp ${MATSLISE_SRC})
target_link_libraries(schrodinger Eigen3::Eigen)

# add_subdirectory(swig/python_wrapper)

# add_subdirectory(python)
