set(CMAKE_CXX_STANDARD 17)
set(PYSLISE_VERSION "3.1.1")
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/../cmake)

message("Matslise")
if (MSVC)
    add_compile_options(
            /wd4267 # possible loss of data size_t -> int
            /wd4244 # possible loss of data size_t -> int
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

set(MATSLISE_SRC
        matslise/matslise.cpp matslise/halfrange.cpp matslise/abstract_matslise.cpp matslise/sector.cpp
        util/quadrature.cpp util/calculateEta.cpp util/sectorbuilder.impl.h)

option(MATSLISE_MATSCS "Build matscs for coupled systems" ON)
option(MATSLISE_2D "Build matslise for 2D equations" ON)
option(MATSLISE_3D "Build matslise for 3D equations" OFF)

if (MATSLISE_3D)
    set(MATSLISE_2D ON)
endif ()
if (MATSLISE_2D)
    set(MATSLISE_MATSCS ON)
endif ()

if (MATSLISE_MATSCS)
    set(MATSLISE_SRC ${MATSLISE_SRC}
            matscs/matscs.cpp matscs/sector.cpp matscs/matscs_formulas.cpp)
endif ()

if (MATSLISE_2D)
    set(MATSLISE_SRC ${MATSLISE_SRC}
            matslise2d/halfrange.cpp matslise2d/eigenfunctions.cpp
            matslise2d/matslise2d.cpp matslise2d/sector.cpp matslise2d/basisQuadrature.cpp
            matsliseNd/sector.cpp)
endif ()

if (MATSLISE_3D)
    set(MATSLISE_SRC ${MATSLISE_SRC}
            matslise3d/matslise3d.cpp matslise3d/sector.cpp)
endif ()

if (DEFINED EMSCRIPTEN)
    add_executable(matslise js/module.cpp js/matslise.cpp js/matslise2d.cpp ${MATSLISE_SRC})
    find_package(Eigen3 3.3 REQUIRED NO_MODULE)
    target_link_libraries(matslise Eigen3::Eigen)

    configure_file(js/matslise.d.ts matslise.d.ts COPYONLY)
    set_target_properties(matslise
            PROPERTIES SUFFIX ".js"
            LINK_FLAGS " --bind -s NO_EXIT_RUNTIME=1 -O3 -s ASSERTIONS=1 -s DISABLE_EXCEPTION_CATCHING=0 -s SAFE_HEAP=1 -s TOTAL_MEMORY=128mb -s MODULARIZE=1 -s USE_ES6_IMPORT_META=0")
    em_link_js_library(matslise ${libraryJsFiles})
else ()
    set(libraries Eigen3::Eigen)

    add_library(matslise ${MATSLISE_SRC})
    install(TARGETS matslise
            ARCHIVE DESTINATION lib
            LIBRARY DESTINATION lib)
    install(FILES matslise.h matscs.h matslise2d.h matslise3d.h matsliseNd.h formula_constants.h
            DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/matslise)
    install(FILES util/y.h util/eigen.h util/sectorbuilder.h util/rectangle.h util/scoped_timer.h
            DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/matslise/util)
    install(FILES matslise2d/basisQuadrature.h
            DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/matslise/matslise2d)

    OPTION(PROFILE "Compile with profiling information." OFF)
    if (PROFILE)
        ADD_DEFINITIONS(-DMATSLISE_PROFILING)
    ENDIF (PROFILE)

    OPTION(PARALLEL "Compile with parallel algorithms." OFF)
    IF (PARALLEL)
        if (PROFILE)
            message("only one of 'parallel' or 'profile' can be on, ignoring 'parallel'")
        else ()
            ADD_DEFINITIONS(-DMATSLISE_parallel -DEIGEN_DONT_PARALLELIZE)
            find_package(OpenMP REQUIRED)
            target_link_libraries(matslise PUBLIC OpenMP::OpenMP_CXX)
        ENDIF ()
    ENDIF (PARALLEL)

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

    OPTION(USE_MKL "Use the Intel's MKL if available" OFF)
    if (USE_MKL)
        find_package(MKL)
        if (MKL_FOUND)
            set(BLA_INCLUDE_DIR "${MKL_INCLUDE_DIR}")
            set(BLAS_LIBRARIES "${MKL_LIBRARIES}")
            set(BLA_DEFINITIONS USE_MKL)
            foreach (target matslise)
                target_compile_definitions(${target} PUBLIC ${BLA_DEFINITIONS})
                target_include_directories(${target} PUBLIC ${BLA_INCLUDE_DIR})

                set(THREADS_PREFER_PTHREAD_FLAG ON)
                find_package(Threads REQUIRED)
                target_link_libraries(${target} PUBLIC
                        m
                        Threads::Threads
                        ${CMAKE_DL_LIBS}
                        ${CMAKE_DL_LIBS}
                        ${BLAS_LIBRARIES})
            endforeach ()
            add_definitions(-DEIGEN_USE_BLAS -DEIGEN_USE_MKL -DMKL_DIRECT_CALL)
        endif ()
    endif ()

    find_package(Eigen3 3.3 REQUIRED NO_MODULE)
    message("-- Eigen3 version ${Eigen3_VERSION}")
    target_link_libraries(matslise PUBLIC ${libraries})

    OPTION(WITH_PYTHON "Also build the python bindings" ON)
    if (WITH_PYTHON)
        if (EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/../cmake/pybind11/CMakeLists.txt)
            add_subdirectory(../cmake/pybind11 pybind11)
            message("-- python: ${PYTHON_EXECUTABLE}")
            pybind11_add_module(pyslise MODULE ${MATSLISE_SRC}
                    pyslise/module.cpp pyslise/pyslise.cpp pyslise/pyscs.cpp pyslise/pyslise2d.cpp pyslise/pyslise3d.cpp)
            target_link_libraries(pyslise PRIVATE ${libraries})
            IF (PARALLEL)
                target_link_libraries(pyslise PUBLIC OpenMP::OpenMP_CXX)
            ENDIF (PARALLEL)

            configure_file(pyslise/setup.py setup.py NEWLINE_STYLE LF)
            configure_file(pyslise/description.md description.md NEWLINE_STYLE LF)
            add_custom_target(build_wheel
                    COMMAND ${CMAKE_COMMAND} -E env PYSLISE_LIBRARY="$<TARGET_FILE:pyslise>"
                    ${PYTHON_EXECUTABLE} setup.py bdist_wheel
                    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
            add_dependencies(build_wheel pyslise)

            add_custom_target(pyslise_install
                    # COMMAND ${PYTHON_EXECUTABLE} -m pip uninstall -y --verbose pyslise
                    COMMAND ${CMAKE_COMMAND} -E env PYSLISE_LIBRARY="$<TARGET_FILE:pyslise>"
                    ${PYTHON_EXECUTABLE} -m pip install --upgrade --ignore-installed .
                    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
            add_dependencies(pyslise_install pyslise)
        else ()
            message("To enable python bindings, make sure pybind11 is available in /cmake/pybind11. Or disable the python bindings with -DWITH_PYTHON=OFF")
        endif ()
    endif ()
endif (DEFINED EMSCRIPTEN)
