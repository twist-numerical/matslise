# Matslise

This [1] is the C++ version of Matslise [2].

Matslise is a collection of algorithms to solve the one and two dimensional time-independent Schrödinger equations. These algorithms are based upon constant perturbation methods to efficiently solve these eigenvalue problems.

The code (and name) is based on Matslise [1][2]. This is a feature-rich matlab library for solving the one dimensional time independent Schrödinger equation.

To solve the two dimensional problem an algorithm is developed on the basis of a method proposed by Ixaru [3].

This implementation is developed in C++ with a focus on efficiency. This code is precompiled for 64 bit windows and linux and packaged in python wheels.

## Install

Installing the python bindings is as simple as
```
pip install pyslise
```

## Documentation

Full documentation can be found on 
[matslise.ugent.be](https://matslise.ugent.be/). This document contains some examples of how to use this library.

## Building from source

### Obtaining the source
To obtain the source one could simply:

```
git clone --recursive https://github.com/twist-numerical/matslise.git
```

The --recursive flag is needed to also add pybind to the project.

### Dependencies
When building from source you will always need a compiler and CMake. On Windows this means Visual Studio.

Other than those Matslise has only one required dependency: [Eigen](http://eigen.tuxfamily.org). A recent version satisfies.

If you have a package manager (apt/brew/pacman/...) Eigen can be installed via those tools. If that doesn't work installing eigen [from source](https://bitbucket.org/eigen/eigen/src/default/INSTALL) is not that difficult. Make sure you follow the steps for CMake.

### Building
```
mkdir build_dir
cd build_dir
cmake ..
```
The next thing to do is pick a target. What do you want to build?

- *matlise_test* will build and run all tests
  ```
  cmake --build . --target matslise_test
  ./test/matslise_test
  ```
- *pyslise_install* will install pyslise within the default python. To change to an other python you can configure cmake with the flag `-DPYTHON_EXECUTABLE=path/to/python`. Note that pip has to be installed for that python.
  ```
  cmake --build . --target pyslise_install
  ```
- *build_wheel* will generate a wheel that can easily be installed on other systems. This wheel will be saved in the directory wheelhouse. Besides pip the configured python also needs the wheel package.

When configuring CMake (`cmake ..`), it is possible to add some configuration options:
- `-DPYTHON_EXECUTABLE=path/to/python` enables you to build for other python installs on your system.
- `-DEigen3_DIR=path/to/eigen3/cmake` if a suitable eigen is not found automatically, one can be specified.
- `-DMATSLISE_LONG_DOUBLE=OFF` when `OFF` is changed to `ON` matslise will also be compiled for the type `long double`.
- `-DMATSLISE_QUADMATH=OFF` when `OFF` is changed to `ON` matslise will also be compiled for a quadruple precision type. For this to work Boost is needed.

For example: running the tests, also with the other floating point types:
```
cmake .. -DQUADMATH=ON -DLONG_DOUBLE=ON
cmake --build . --target matslise_test
./test/matslise_test
```

## Notes

- lapack (especially from mkl) does allocate some buffers. Because these are never freed this can cause valgrind to suspect memory leaks. See  (https://stackoverflow.com/questions/36197527/insight-as-to-why-valgrind-shows-memory-leak-for-intels-mkl-lapacke)

## Bibliography
* **[1]** Baeyens, Toon, and Marnix Van Daele. “The Fast and Accurate Computation of Eigenvalues and Eigenfunctions of Time-Independent One-Dimensional Schrödinger Equations.” Computer Physics Communications, August 26, 2020, 107568. https://doi.org/10.1016/j.cpc.2020.107568.
* **[2]** Ledoux, Veerle, and Marnix Van Daele. “MATSLISE 2.0 : A Matlab Toolbox for Sturm-Liouville Computations.” ACM TRANSACTIONS ON MATHEMATICAL SOFTWARE 42, no. 4 (2016): 18.
* **[3]** Ixaru, L. Gr. “New Numerical Method for the Eigenvalue Problem of the 2D Schrödinger Equation.” Computer Physics Communications 181 (October 1, 2010): 1738–42. https://doi.org/10.1016/j.cpc.2010.06.031.
