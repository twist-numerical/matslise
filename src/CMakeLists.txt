set(CMAKE_CXX_STANDARD 14)
set(PYSLISE_VERSION 0.0.4)
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/../cmake)

message("Matslise")
if(MSVC)
    add_compile_options(
    /wd4661 # no suitable definition provided for explicit template instantiation request
    /wd4820 #'bytes' bytes padding added after construct 'member_name'
    /wd4464  # relative include path contains '..'
    /DNOMINMAX
    /D_USE_MATH_DEFINES
    /bigobj)
else()
    add_compile_options(-Wall -Wextra)
    OPTION(OPTIMIZE "Enable extra optimization flags (-O3 -march=native)" OFF)
    IF (OPTIMIZE)
        message("-- Using -O3 -march=native")
        add_compile_options(-O3 -march=native)
    ENDIF (OPTIMIZE)
endif()

SET(MATSLISE_SRC
        matslise/matslise.cpp matslise/halfrange.cpp matslise/sector.cpp
        matscs/matscs.cpp matscs/sector.cpp
        se2d/halfrange.cpp se2d/se_eigenvalues.cpp se2d/se_eigenfunctions.cpp se2d/se2d.cpp se2d/sector.cpp schrodinger.h
        util/legendre.cpp)

if (DEFINED EMSCRIPTEN)
    add_executable(matslise js/main.cpp ${MATSLISE_SRC})
    find_package(Eigen3 3.3 REQUIRED NO_MODULE)
    target_link_libraries(matslise Eigen3::Eigen)

    set_target_properties(matslise
            PROPERTIES SUFFIX ".js"
            LINK_FLAGS " --bind -s  NO_EXIT_RUNTIME=1 -O2 -s TOTAL_MEMORY=128mb -s EXPORT_ES6=1 -s MODULARIZE=1 -s USE_ES6_IMPORT_META=0")
    em_link_js_library(matslise ${libraryJsFiles})
else ()
    set(targets matslise pyslise)
    set(libraries Eigen3::Eigen)

    add_subdirectory(../cmake/pybind11 pybind11)
    pybind11_add_module(pyslise MODULE ${MATSLISE_SRC}
            pyslise/module.cpp pyslise/pyslise.cpp pyslise/pyscs.cpp pyslise/pyse2d.cpp)

    add_library(matslise ${MATSLISE_SRC})

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
            foreach (target ${targets})
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
            add_definitions(-DEIGEN_USE_BLAS -DEIGEN_USE_MKL -DMKL_DIRECT_CALL)
        endif ()
    endif ()

    find_package(Eigen3 3.3 REQUIRED NO_MODULE)
    message("-- Eigen3 version ${Eigen3_VERSION}")
    target_link_libraries(matslise PUBLIC ${libraries})
    target_link_libraries(pyslise PRIVATE ${libraries})

    message("-- python: ${PYTHON_EXECUTABLE}")
    message("-- auditwheel: ${AUDITWHEEL_repair_plat}")
    configure_file(pyslise/setup.py setup.py)
    configure_file(pyslise/description.md description.md)
    file(TO_NATIVE_PATH "${CMAKE_SOURCE_DIR}/wheelhouse" WHEEL_COPY_TO)
    add_custom_target(build_wheel
            COMMAND ${CMAKE_COMMAND} -E env PYSLISE_LIBRARY="$<TARGET_FILE:pyslise>" PYSLISE_MOVE="${WHEEL_COPY_TO}"
            ${PYTHON_EXECUTABLE} setup.py bdist_wheel
            WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
    add_dependencies(build_wheel pyslise)

    add_custom_target(pyslise_install
            COMMAND ${CMAKE_COMMAND} -E env PYSLISE_LIBRARY="$<TARGET_FILE:pyslise>"
            ${PYTHON_EXECUTABLE} setup.py bdist_wheel
            COMMAND ${CMAKE_COMMAND} -E chdir ${CMAKE_BINARY_DIR}
            ${PYTHON_EXECUTABLE} -m pip uninstall --verbose pyslise -y
            COMMAND ${CMAKE_COMMAND} -E chdir ${CMAKE_CURRENT_BINARY_DIR}/dist
            ${PYTHON_EXECUTABLE} -m pip install --upgrade --ignore-installed ${CMAKE_CURRENT_BINARY_DIR}/dist/*.whl
            WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
    add_dependencies(pyslise_install pyslise)

    install(FILES matslise.h Array2D.h matscs.h se2d.h schrodinger.h DESTINATION include/matslise/)
    install(FILES util/lobatto.h util/y.h util/eigen.h util/SectorBuilder.h DESTINATION include/matslise/util)

endif (DEFINED EMSCRIPTEN)
