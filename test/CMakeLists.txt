set(CMAKE_CXX_STANDARD 17)
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/../cmake)

cmake_policy(SET CMP0077 NEW)

message("Matslise test")

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
            GIT_TAG        v3.3.2
    )

    FetchContent_MakeAvailable(Catch2)

    set(MATSLISE_TEST_SRC
            util/legendre.cpp util/rectangle.cpp util/calculateEta.cpp util/polynomial.cpp
            matslise/prufer.cpp matslise/mathieu.cpp matslise/simple.cpp matslise/coffey_evans.cpp
            matslise/hr_mathieu.cpp matslise/ixaru.cpp matslise/marletta.cpp matslise/periodic.cpp
            matslise/sturm_liouville.cpp matslise/harmonic.cpp
            liouville/liouville.cpp liouville/klotter.cpp)

    add_executable(matslise_test EXCLUDE_FROM_ALL ${MATSLISE_TEST_SRC})

    target_link_libraries(matslise_test PRIVATE matslise Catch2::Catch2WithMain)
    add_dependencies(matslise_test matslise)
endif ()
