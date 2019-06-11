find_package(Eigen3 3.3 REQUIRED NO_MODULE)
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/../cmake)

set(CMAKE_CXX_FLAGS "-O2 -Wall")

if (NOT (DEFINED EMSCRIPTEN))
    set(libraries matslise Eigen3::Eigen)
    find_package(Boost)
    if (Boost_FOUND)
        add_definitions(-DBOOST)
        include_directories(${Boost_INCLUDE_DIRS})


        find_package(Quadmath)
        if (Quadmath_FOUND)
            list(APPEND libraries Quadmath::quadmath)
        endif ()
    endif ()

    add_executable(matslise_test EXCLUDE_FROM_ALL
            matslise/prufer.cpp matslise/mathieu.cpp matslise/simple.cpp matslise/coffey_evans.cpp
            matslise/hr_mathieu.cpp
            matscs/propagatePsi.cpp
            se2d/henon.cpp se2d/ixaru.cpp se2d/zero.cpp se2d/quartic.cpp
            util/lobatto.cpp
            main.cpp)
    target_link_libraries(matslise_test ${libraries})
endif ()
