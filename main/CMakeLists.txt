if (NOT (DEFINED EMSCRIPTEN))
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

    find_package(Eigen3 3.3 REQUIRED NO_MODULE)

    add_executable(schrodinger EXCLUDE_FROM_ALL main.cpp)
    target_link_libraries(schrodinger PRIVATE matslise Eigen3::Eigen)
endif ()
