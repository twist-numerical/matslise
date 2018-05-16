find_package(SWIG REQUIRED)
include(${SWIG_USE_FILE})

# Include python
find_package(PythonLibs 2.7 REQUIRED)
include_directories(${PYTHON_INCLUDE_PATH})


set(CMAKE_SWIG_FLAGS "")
set_source_files_properties(pyslise.i PROPERTIES CPLUSPLUS ON)
include_directories(../..)

# Add swig module
swig_add_module(pyslise python pyslise.i)
swig_link_libraries(pyslise matslise ${PYTHON_LIBRARIES})

# Files to install with Python
set(PYTHON_INSTALL_FILES
        ${CMAKE_CURRENT_BINARY_DIR}/pyslise.py
        ${CMAKE_CURRENT_BINARY_DIR}/_pyslise.so)

# Configure setup.py.in and copy to output directory
set(SETUP_PY_IN setup.py.in)
set(SETUP_PY_OUT ${CMAKE_CURRENT_BINARY_DIR}/setup.py)
configure_file(${SETUP_PY_IN} ${SETUP_PY_OUT})

# Declare install target for python
#install(TARGETS swig_example
#        COMMAND "${PYTHON_EXECUTABLE} setup.py.in"
#        COMPONENT swig-python)

# Install target to call setup.py.in
add_custom_target(install-python
        DEPENDS _pyslise
        #COMMAND /usr/lib/sagemath/local/bin/python2.7 ${SETUP_PY_OUT} install)
        COMMAND /home/toon/miniconda3/bin/python ${SETUP_PY_OUT} install)