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
    set(libraries matslise Eigen3::Eigen)

    OPTION(MATSLISE_LONG_DOUBLE "Compile also for the 'long double' type." OFF)
    IF (MATSLISE_LONG_DOUBLE)
        ADD_DEFINITIONS(-DMATSLISE_LONG_DOUBLE)
    ENDIF (MATSLISE_LONG_DOUBLE)

    OPTION(MATSLISE_QUADMATH "Compile also for quadruple precision (_float128). Boost and a supporting compiler are required." OFF)
    IF (MATSLISE_QUADMATH)
        find_package(Boost REQUIRED)
        if (Boost_FOUND)
            add_definitions(-DBOOST)
            include_directories(${Boost_INCLUDE_DIRS})

            find_package(Quadmath REQUIRED)
            if (Quadmath_FOUND)
                add_definitions(-DMATSLISE_float128)
                list(APPEND libraries Quadmath::quadmath)
            endif ()
        endif ()
    ENDIF (MATSLISE_QUADMATH)

    option(MATSLISE_MATSCS "Build matscs for coupled systems" ON)

    set(MATSLISE_TEST_SRC
            util/legendre.cpp util/rectangle.cpp util/heap.cpp util/calculateEta.cpp
            matslise/prufer.cpp matslise/mathieu.cpp matslise/simple.cpp matslise/coffey_evans.cpp
            matslise/hr_mathieu.cpp matslise/ixaru.cpp matslise/marletta.cpp matslise/periodic.cpp)

    if (MATSLISE_MATSCS)
        set(MATSLISE_TEST_SRC ${MATSLISE_TEST_SRC}
                matscs/propagatePsi.cpp matscs/propagateTheta.cpp)
    endif ()

    add_executable(matslise_test EXCLUDE_FROM_ALL ${MATSLISE_TEST_SRC} main.cpp)

    find_package(OpenMP)
    if (OpenMP_CXX_FOUND)
        target_link_libraries(matslise_test PUBLIC OpenMP::OpenMP_CXX)
    endif ()

    OPTION(MATSLISE_USE_MKL "Use the Intel's MKL if available" OFF)
    if (USE_MKL)
        find_package(MKL)
        if (MKL_FOUND)
            set(BLA_INCLUDE_DIR "${MKL_INCLUDE_DIR}")
            set(BLAS_LIBRARIES "${MKL_LIBRARIES}")
            set(BLA_DEFINITIONS USE_MKL)
            target_compile_definitions(matslise_test PUBLIC ${BLA_DEFINITIONS})
            target_include_directories(matslise_test PUBLIC ${BLA_INCLUDE_DIR})

            set(THREADS_PREFER_PTHREAD_FLAG ON)
            find_package(Threads REQUIRED)
            target_link_libraries(matslise_test PUBLIC
                    m
                    Threads::Threads
                    ${CMAKE_DL_LIBS}
                    ${BLAS_LIBRARIES})
            add_definitions(-DEIGEN_USE_BLAS -DEIGEN_USE_MKL -DMKL_DIRECT_CALL)
        endif ()
    endif ()

    target_link_libraries(matslise_test PRIVATE ${libraries})

    add_dependencies(matslise_test matslise)
endif ()
