set(CMAKE_CXX_STANDARD 17)
set(MATSLISE_VERSION "3.2.1")
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/../cmake)

cmake_policy(SET CMP0077 NEW)

message("Matslise")
if (MSVC)
    add_compile_options(
            /wd4267 # possible loss of data size_t -> int
            /wd4244 # possible loss of data size_t -> int
            /wd4661 # no suitable definition provided for explicit template instantiation request
            /wd4820 #'bytes' bytes padding added after construct 'member_name'
            /wd4464 # relative include path contains '..'
            /EHsc
    )
else ()
    add_compile_options(-Wall -Wextra -Wno-register -Wno-undefined-var-template)
endif ()

Include(FetchContent)
FetchContent_Declare(
        eigen
        GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
        GIT_TAG 3.4.0
)
FetchContent_MakeAvailable(eigen)

set(MATSLISE_SRC
        matslise/matslise.cpp matslise/halfrange.cpp matslise/abstract_matslise.cpp matslise/periodic.cpp
        matslise/matslise_formulas.cpp matslise/sector.cpp liouville/liouville.cpp
        liouville/sturm_liouville.cpp util/sectorbuilder.cpp)

set(MATSLISE_MATSCS_SRC
        matscs/matscs.cpp matscs/sector.cpp matscs/matscs_formulas.cpp)
set(PYSLISE_SRC pyslise/module.cpp pyslise/pyslise.cpp pyslise/sturm_liouville.cpp)

set(MATSLISE_FULL_SRC ${MATSLISE_SRC} ${MATSLISE_MATSCS_SRC})

set(defines "")

option(MATSLISE_MATSCS "Build matscs for coupled systems" ON)

if (MATSLISE_MATSCS)
    set(MATSLISE_SRC ${MATSLISE_SRC} ${MATSLISE_MATSCS_SRC})
    set(PYSLISE_SRC ${PYSLISE_SRC} pyslise/pyscs.cpp)
    list(APPEND defines WITH_MATSCS)
endif ()

list(APPEND defines MATSLISE_VERSION="${MATSLISE_VERSION}")

if (DEFINED EMSCRIPTEN)
    add_executable(matslise js/module.cpp js/matslise.cpp js/sturm_liouville.cpp ${MATSLISE_FULL_SRC})
    target_link_libraries(matslise Eigen3::Eigen)

    configure_file(js/matslise.d.ts matslise.d.ts COPYONLY)
    set_target_properties(matslise
            PROPERTIES SUFFIX ".js"
            LINK_FLAGS " --bind -s NO_EXIT_RUNTIME=1 -O3 -s ASSERTIONS=1 -s DISABLE_EXCEPTION_CATCHING=0 -s SAFE_HEAP=1 -s TOTAL_MEMORY=128mb -s MODULARIZE=1 -s USE_ES6_IMPORT_META=0")
    em_link_js_library(matslise ${libraryJsFiles})
else ()
    set(libraries Eigen3::Eigen)

    option(MATSLISE_STATIC "Build matslise as a static library" OFF)
    if (MATSLISE_STATIC)
        add_library(matslise STATIC ${MATSLISE_SRC})
    else ()
        add_library(matslise ${MATSLISE_SRC})
    endif ()

    set(MATSLISE_COMPILE_OPTIONS)

    IF(MSVC)
        set(MATSLISE_COMPILE_OPTIONS ${MATSLISE_COMPILE_OPTIONS} /DNOMINMAX /D_USE_MATH_DEFINES /bigobj)
    endif ()

    OPTION(MATSLISE_OPTIMIZE "Enable extra optimization flags (-O3 -march=native)" OFF)
    if (MATSLISE_OPTIMIZE AND NOT MSVC)
        message("-- Using -O3 -march=native")
        set(MATSLISE_COMPILE_OPTIONS ${MATSLISE_COMPILE_OPTIONS} -O3 -march=native)
    endif ()
    target_compile_options(matslise PUBLIC ${MATSLISE_COMPILE_OPTIONS})

    install(TARGETS matslise
            ARCHIVE DESTINATION lib
            LIBRARY DESTINATION lib)
    install(FILES matslise.h matscs.h matsliseNd.h formula_constants.h
            DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/matslise)
    install(FILES util/y.h util/eigen.h util/sectorbuilder.h util/rectangle.h util/scoped_timer.h
            DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/matslise/util)

    OPTION(MATSLISE_PROFILING "Compile with profiling information." OFF)
    if (MATSLISE_PROFILING)
        list(APPEND defines MATSLISE_PROFILING)
    ENDIF (MATSLISE_PROFILING)

    OPTION(MATSLISE_LONG_DOUBLE "Option description" OFF)
    IF (MATSLISE_LONG_DOUBLE)
        list(APPEND defines MATSLISE_LONG_DOUBLE)
    ENDIF (MATSLISE_LONG_DOUBLE)

    OPTION(MATSLISE_QUADMATH "Option description" OFF)
    IF (MATSLISE_QUADMATH)
        find_package(Boost REQUIRED)
        if (Boost_FOUND)
            list(APPEND defines BOOST)
            include_directories(${Boost_INCLUDE_DIRS})

            find_package(Quadmath REQUIRED)
            if (Quadmath_FOUND)
                list(APPEND defines MATSLISE_QUADMATH)
                list(APPEND libraries Quadmath::quadmath)
            endif ()
        endif ()
    ENDIF (MATSLISE_QUADMATH)

    OPTION(MATSLISE_USE_MKL "Use the Intel's MKL if available" OFF)
    if (MATSLISE_USE_MKL)
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
                        ${BLAS_LIBRARIES})
            endforeach ()
            list(APPEND defines EIGEN_USE_BLAS EIGEN_USE_MKL MKL_DIRECT_CALL)
        endif ()
    endif ()

    message("-- Eigen3 version ${Eigen3_VERSION}")
    target_link_libraries(matslise PUBLIC ${libraries})
    target_compile_definitions(matslise PUBLIC ${defines})

    OPTION(MATSLISE_PYTHON "Also build the python bindings" ON)
    if (MATSLISE_PYTHON)
        if (EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/../cmake/pybind11/CMakeLists.txt)
            add_subdirectory(../cmake/pybind11 pybind11)
            message("-- python: ${PYTHON_EXECUTABLE}")
            pybind11_add_module(pyslise MODULE ${MATSLISE_SRC} ${PYSLISE_SRC})
            target_link_libraries(pyslise PRIVATE ${libraries})
            target_compile_definitions(pyslise PUBLIC ${defines})
            target_compile_options(pyslise PUBLIC ${MATSLISE_COMPILE_OPTIONS})

            configure_file(pyslise/setup.py setup.py NEWLINE_STYLE LF)
            configure_file(pyslise/description.md description.md NEWLINE_STYLE LF)
            add_custom_target(build_wheel
                    COMMAND ${CMAKE_COMMAND} -E env PYSLISE_LIBRARY="$<TARGET_FILE:pyslise>"
                    ${PYTHON_EXECUTABLE} setup.py bdist_wheel
                    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
            add_dependencies(build_wheel pyslise)

            add_custom_target(pyslise_install
                    COMMAND ${CMAKE_COMMAND} -E env PYSLISE_LIBRARY="$<TARGET_FILE:pyslise>"
                    ${PYTHON_EXECUTABLE} -m pip install --upgrade --force-reinstall .
                    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
            add_dependencies(pyslise_install pyslise)

            add_custom_target(pyslise_install_test
                    COMMAND ${PYTHON_EXECUTABLE} test/pyslise/test_docs.py
                    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
            add_dependencies(pyslise_install_test pyslise_install)
        else ()
            message("To enable python bindings, make sure pybind11 is available in /cmake/pybind11. Or disable the python bindings with -DMATSLISE_PYTHON=OFF")
        endif ()
    endif ()
endif (DEFINED EMSCRIPTEN)
