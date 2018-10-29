
find_package(Eigen3 3.3 REQUIRED NO_MODULE)

if (NOT (DEFINED EMSCRIPTEN))
    add_executable(schrodinger EXCLUDE_FROM_ALL main.cpp)
    target_link_libraries(schrodinger matslise Eigen3::Eigen)
endif ()
