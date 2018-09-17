cmake_minimum_required(VERSION 3.5)
project(matslise)

set(CMAKE_CXX_FLAGS "-O2 -Wall -Wextra -fPIC")
set(CMAKE_CXX_STANDARD 11)

cmake_minimum_required(VERSION 3.5)

set(CMAKE_INSTALL_PREFIX /usr/)

INCLUDE_DIRECTORIES("/home/toon/Programmeren/cpp/eigen")

SET(MATSLISE_SRC src/util/calculateEta.cpp src/util/lobatto.cpp
        src/matslise/matslise.cpp src/matslise/halfrange.cpp src/matslise/sector.cpp src/matslise/matslise_formulas.cpp src/matslise/eigenfunction_calculator.cpp
        src/matscs/matscs.cpp src/matscs/sector.cpp src/matscs/matscs_formulas.cpp src/matscs/eigenfunction_calculator.cpp
        src/se2d/se2d.cpp src/se2d/sector.cpp)

if (DEFINED EMSCRIPTEN)
    add_executable(matslise src/js/main.cpp ${MATSLISE_SRC})

    set_target_properties(matslise
            PROPERTIES SUFFIX ".js"
            LINK_FLAGS " --bind -s  NO_EXIT_RUNTIME=1 -O2")
    em_link_js_library(matslise ${libraryJsFiles})
else ()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fno-var-tracking-assignments")
    add_library(matslise SHARED ${MATSLISE_SRC} src/matslise/eigenfunction_calculator.cpp)
    # target_link_libraries(matslise Eigen3::Eigen)

    install(TARGETS matslise DESTINATION lib)
    install(FILES src/matslise.h src/Array2D.h src/matscs.h src/se2d.h src/Matrix2D.h src/Evaluator.h DESTINATION include/matslise/)

    add_executable(schrodinger src/main.cpp ${MATSLISE_SRC})
    # target_link_libraries(schrodinger Eigen3::Eigen)
endif (DEFINED EMSCRIPTEN)
