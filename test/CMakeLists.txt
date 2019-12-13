find_package(Eigen3 3.3 REQUIRED NO_MODULE)
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/../cmake)

if(MSVC)
    add_compile_options(
	/wd4661 # no suitable definition provided for explicit template instantiation request
	/wd4820 #'bytes' bytes padding added after construct 'member_name'
	/wd4464  # relative include path contains '..'
	/DNOMINMAX
	/D_USE_MATH_DEFINES
	/bigobj)
else()
    add_compile_options(-Wall -Wextra -Werror)
endif()

if (NOT (DEFINED EMSCRIPTEN))
    set(libraries matslise Eigen3::Eigen)

    OPTION(LONG_DOUBLE "Option description" OFF)
    IF (LONG_DOUBLE)
        ADD_DEFINITIONS(-DMATSLISE_long_double)
    ENDIF (LONG_DOUBLE)

    OPTION(QUADMATH "Option description" OFF)
    IF (QUADMATH)
        find_package(Boost REQUIRED)
        if (Boost_FOUND)
            add_definitions(-DBOOST)
            include_directories(${Boost_INCLUDE_DIRS})

            find_package(Quadmath)
            if (Quadmath_FOUND)
                add_definitions(-DMATSLISE_float128)
                list(APPEND libraries Quadmath::quadmath)
            endif ()
        endif ()
    ENDIF (QUADMATH)

    add_executable(matslise_test EXCLUDE_FROM_ALL
            matslise/prufer.cpp matslise/mathieu.cpp matslise/simple.cpp matslise/coffey_evans.cpp
            matslise/hr_mathieu.cpp
            matscs/propagatePsi.cpp matscs/propagateTheta.cpp
            se2d/henon.cpp se2d/ixaru.cpp se2d/zero.cpp se2d/quartic.cpp
            util/lobatto.cpp
            main.cpp)

    find_package(OpenMP)
    if(OpenMP_CXX_FOUND)
        target_link_libraries(matslise_test PUBLIC OpenMP::OpenMP_CXX)
    endif()

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

    target_link_libraries(matslise_test PRIVATE ${libraries})


    configure_file(benchmark.py benchmark.py)
    file(TO_NATIVE_PATH "${CMAKE_SOURCE_DIR}/benchmark.txt" benchmark_results)
    add_custom_target(benchmark_tests
            COMMAND ${CMAKE_COMMAND} -E env MATSLISE_TEST_EXECUTABLE="$<TARGET_FILE:matslise_test>" MATSLISE_BENCHMARK_RESULTS="${benchmark_results}" ${PYTHON_EXECUTABLE} benchmark.py
            WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
    add_dependencies(matslise_test matslise)
endif ()
