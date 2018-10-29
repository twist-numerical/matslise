set(CMAKE_CXX_FLAGS "-O2 -Wall -Wextra -fPIC")
set(CMAKE_CXX_STANDARD 11)

find_package(Eigen3 3.3 REQUIRED NO_MODULE)

SET(MATSLISE_SRC util/calculateEta.cpp util/lobatto.cpp
        matslise/matslise.cpp matslise/halfrange.cpp matslise/sector.cpp matslise/matslise_formulas.cpp matslise/eigenfunction_calculator.cpp
        matscs/matscs.cpp matscs/sector.cpp matscs/matscs_formulas.cpp matscs/eigenfunction_calculator.cpp
        se2d/se2d.cpp se2d/sector.cpp)

if (DEFINED EMSCRIPTEN)
    add_executable(matslise js/main.cpp ${MATSLISE_SRC})
    target_link_libraries(matslise Eigen3::Eigen)

    set_target_properties(matslise
            PROPERTIES SUFFIX ".js"
            LINK_FLAGS " --bind -s  NO_EXIT_RUNTIME=1 -O2 -s TOTAL_MEMORY=128mb")
    em_link_js_library(matslise ${libraryJsFiles})
else ()
    add_library(matslise SHARED ${MATSLISE_SRC})
    target_link_libraries(matslise Eigen3::Eigen)

    install(TARGETS matslise LIBRARY DESTINATION lib)
    install(FILES matslise.h Array2D.h matscs.h se2d.h Matrix2D.h Evaluator.h DESTINATION include/matslise/)


endif (DEFINED EMSCRIPTEN)
