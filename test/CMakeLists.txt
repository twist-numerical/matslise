set(CMAKE_CXX_STANDARD 17)
find_package(Eigen3 3.3 REQUIRED NO_MODULE)
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/../cmake)

cmake_policy(SET CMP0077 NEW)

message("Matslise test")
if (MSVC)
    add_compile_options(
            /wd4661 # no suitable definition provided for explicit template instantiation request
            /wd4820 #'bytes' bytes padding added after construct 'member_name'
            /wd4464  # relative include path contains '..'
            /DNOMINMAX
            /D_USE_MATH_DEFINES
            /bigobj)
else ()
    add_compile_options(-Wall -Wextra -Wno-register -Wno-undefined-var-template)

    OPTION(OPTIMIZE "Enable extra optimization flags (-O3 -march=native)" OFF)
    IF (OPTIMIZE)
        message("-- Using -O3 -march=native")
        add_compile_options(-O3 -march=native)
    ENDIF (OPTIMIZE)
endif ()

if (DEFINED EMSCRIPTEN)
    add_custom_target(matslise_test_copy_files ALL
            COMMAND ${CMAKE_COMMAND} -E remove_directory ${CMAKE_BINARY_DIR}/test-js
            COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_SOURCE_DIR}/test/js ${CMAKE_BINARY_DIR}/test-js
            DEPENDS matslise)
    option(NODE_BIN "Binary directory of node installation. Make sure 'ava' is installed." "")
    add_custom_target(matslise_test
            COMMAND ${CMAKE_COMMAND} -E env PATH=$ENV{PATH}:${NODE_BIN} npm link ava
            COMMAND ${CMAKE_COMMAND} -E env PATH=$ENV{PATH}:${NODE_BIN} ava --verbose -T 5m
            DEPENDS matslise_test_copy_files
            WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/test-js)
else ()
    Include(FetchContent)

    FetchContent_Declare(
            Catch2
            GIT_REPOSITORY https://github.com/catchorg/Catch2.git
            GIT_TAG        v3.3.2 # or a later release
    )

    FetchContent_MakeAvailable(Catch2)

    set(MATSLISE_TEST_SRC
            util/legendre.cpp util/rectangle.cpp util/calculateEta.cpp util/polynomial.cpp
            matslise/prufer.cpp matslise/mathieu.cpp matslise/simple.cpp matslise/coffey_evans.cpp
            matslise/hr_mathieu.cpp matslise/ixaru.cpp matslise/marletta.cpp matslise/periodic.cpp
            matslise/sturm_liouville.cpp
            liouville/liouville.cpp liouville/klotter.cpp)

    add_executable(matslise_test EXCLUDE_FROM_ALL ${MATSLISE_TEST_SRC})

    target_link_libraries(matslise_test PRIVATE matslise Eigen3::Eigen Catch2::Catch2WithMain)
    add_dependencies(matslise_test matslise)
endif ()
