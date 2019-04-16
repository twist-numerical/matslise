set(CMAKE_CXX_FLAGS "-O2 -Wall -fPIC")
set(CMAKE_CXX_STANDARD 11)


SET(MATSLISE_SRC util/calculateEta.cpp util/lobatto.cpp util/y.cpp
        matslise/matslise.cpp matslise/halfrange.cpp matslise/sector.cpp matslise/matslise_formulas.cpp matslise/eigenfunction_calculator.cpp
        matscs/matscs.cpp matscs/sector.cpp matscs/matscs_formulas.cpp matscs/eigenfunction_calculator.cpp
        se2d/se_eigenvalues.cpp se2d/se_eigenfunctions.cpp se2d/se2d.cpp se2d/sector.cpp schrodinger.h util/y.cpp util/y.h util/eigen.h)

if (DEFINED EMSCRIPTEN)
    add_executable(matslise js/main.cpp ${MATSLISE_SRC})
    find_package(Eigen3 3.3 REQUIRED NO_MODULE)
    target_link_libraries(matslise Eigen3::Eigen)

    set_target_properties(matslise
            PROPERTIES SUFFIX ".js"
            LINK_FLAGS " --bind -s  NO_EXIT_RUNTIME=1 -O2 -s TOTAL_MEMORY=128mb -s BINARYEN_TRAP_MODE='clamp'")
    em_link_js_library(matslise ${libraryJsFiles})
else ()
    # set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fno-var-tracking-assignments")
    # set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsave-optimization-record")
    add_library(matslise SHARED ${MATSLISE_SRC})

    list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
    find_package(MKL)
    if (MKL_FOUND)
        set(BLA_INCLUDE_DIR "${MKL_INCLUDE_DIR}")
        set(BLAS_LIBRARIES "${MKL_LIBRARIES}")
        set(BLA_DEFINITIONS USE_MKL)
        message(STATUS "MKL found")
        target_compile_definitions(matslise PUBLIC ${BLA_DEFINITIONS})
        target_include_directories(matslise PUBLIC ${BLA_INCLUDE_DIR})

        set(THREADS_PREFER_PTHREAD_FLAG ON)
        find_package(Threads REQUIRED)
        target_link_libraries(matslise PUBLIC
                m
                Threads::Threads
                ${CMAKE_DL_LIBS}
                ${BLAS_LIBRARIES})
        add_definitions(-DEIGEN_USE_BLAS -DEIGEN_USE_MKL -DMKL_DIRECT_CALL)
    endif ()

    find_package(Eigen3 3.3 REQUIRED NO_MODULE)
    target_link_libraries(matslise PUBLIC Eigen3::Eigen)
    if (WIN32)
        install(TARGETS matslise RUNTIME DESTINATION lib)
    else ()
        install(TARGETS matslise LIBRARY DESTINATION lib)
    endif ()
    install(FILES matslise.h Array2D.h matscs.h se2d.h Evaluator.h schrodinger.h DESTINATION include/matslise/)
    install(FILES util/lobatto.h util/y.h util/eigen.h DESTINATION include/matslise/util)


endif (DEFINED EMSCRIPTEN)
