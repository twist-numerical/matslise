set(CMAKE_CXX_STANDARD 11)
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/../cmake)

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang" OR "${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2 -ansi -Wall -Wextra -Wconversion -fPIC")
    message("-- added GCC/Clang flags")
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /O2 /DNOMINMAX /D_USE_MATH_DEFINES")
    message("-- added MSVC flags")
endif ()

SET(MATSLISE_SRC
        matslise/matslise.cpp matslise/halfrange.cpp matslise/sector.cpp
        matscs/matscs.cpp matscs/sector.cpp
        se2d/se_eigenvalues.cpp se2d/se_eigenfunctions.cpp se2d/se2d.cpp se2d/sector.cpp schrodinger.h
        util/lobatto.cpp util/legendre.cpp)

if (DEFINED EMSCRIPTEN)
    add_executable(matslise js/main.cpp ${MATSLISE_SRC})
    find_package(Eigen3 3.3 REQUIRED NO_MODULE)
    target_link_libraries(matslise Eigen3::Eigen)

    set_target_properties(matslise
            PROPERTIES SUFFIX ".js"
            LINK_FLAGS " --bind -s  NO_EXIT_RUNTIME=1 -O2 -s TOTAL_MEMORY=128mb -s BINARYEN_TRAP_MODE='clamp'")
    em_link_js_library(matslise ${libraryJsFiles})
else ()
    set(targets matslise pyslise)
    set(libraries Eigen3::Eigen)
    add_subdirectory(python/pybind11)
    pybind11_add_module(pyslise MODULE ${MATSLISE_SRC} python/pyslise.cpp)
    add_library(matslise SHARED ${MATSLISE_SRC})

    find_package(Boost)
    if (Boost_FOUND)
        add_definitions(-DBOOST)
        include_directories(${Boost_INCLUDE_DIRS})


        find_package(Quadmath)
        if (Quadmath_FOUND)
            list(APPEND libraries Quadmath::quadmath)
        endif ()
    endif ()

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

    find_package(Eigen3 3.3 REQUIRED NO_MODULE)
    target_link_libraries(matslise PUBLIC ${libraries})
    target_link_libraries(pyslise PRIVATE ${libraries})

    install(FILES matslise.h Array2D.h matscs.h se2d.h schrodinger.h DESTINATION include/matslise/)
    install(FILES util/lobatto.h util/y.h util/eigen.h util/SectorBuilder.h DESTINATION include/matslise/util)

endif (DEFINED EMSCRIPTEN)
