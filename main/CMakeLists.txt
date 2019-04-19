if (NOT (DEFINED EMSCRIPTEN))
    set(CMAKE_CXX_FLAGS "-O2 -Wall")
    find_package(Eigen3 3.3 REQUIRED NO_MODULE)

    add_executable(schrodinger EXCLUDE_FROM_ALL main.cpp)
    target_link_libraries(schrodinger matslise Eigen3::Eigen)
endif ()
