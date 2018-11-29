find_package(Eigen3 3.3 REQUIRED NO_MODULE)

if (NOT (DEFINED EMSCRIPTEN))
    add_executable(test EXCLUDE_FROM_ALL main.cpp matslise/prufer.cpp matslise/mathieu.cpp matscs/propagatePsi.cpp)
    target_link_libraries(test matslise Eigen3::Eigen)
endif ()
