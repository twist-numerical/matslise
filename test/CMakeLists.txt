find_package(Eigen3 3.3 REQUIRED NO_MODULE)

if (NOT (DEFINED EMSCRIPTEN))
    add_executable(matslise_test EXCLUDE_FROM_ALL
            main.cpp
            matslise/prufer.cpp matslise/mathieu.cpp matslise/simple.cpp
            matscs/propagatePsi.cpp
            se2d/henon.cpp se2d/ixaru.cpp se2d/zero.cpp se2d/two_sector.cpp
            util/lobatto.cpp)
    target_link_libraries(matslise_test matslise Eigen3::Eigen)
endif ()
