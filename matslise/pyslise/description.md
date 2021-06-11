# Pyslise ${PYSLISE_VERSION}

Pyslise [1] is a collection of algorithms to solve one (and two, in development) dimensional time-independent Schrödinger equations. These algorithms are based upon constant perturbation methods to efficiently solve these eigenvalue problems.

The code (and name) is based on Matslise [2]. This is a feature-rich MATLAB library for solving the one dimensional time independent Schrödinger equation.

To solve the two dimensional problem an algorithm is developed on the basis of a method proposed by Ixaru [3].

This implementation is developed in C++ with a focus on efficiency. This code is precompiled and packaged in wheels for 64 bit Linux, Windows, and Mac.


## Documentation

Full documentation can be found on 
[matslise.ugent.be](https://matslise.ugent.be/). This document contains some examples of how to use this library.

On the same page an interactive version is available.

## Examples

One dimensional problems can be tackled with:
```python
from pyslise import Pyslise
from math import pi, cos

problem = Pyslise(lambda x: 2*cos(2*x), 0, pi, tolerance=1e-6)
problem.eigenvaluesByIndex(0, 10, (0, 1), (0, 1))
```

Also two dimensional problems are possible:
```python
from pyslise import Pyslise2D

def V(x, y):
    return (1 + x**2) * (1 + y**2)

problem = Pyslise2D(V, -5.5,5.5, -5.5,5.5, tolerance=1e-6)
problem.eigenvalues(0,13)
```

## Bibliography
* **[1]** Baeyens, Toon, and Marnix Van Daele. “The Fast and Accurate Computation of Eigenvalues and Eigenfunctions of Time-Independent One-Dimensional Schrödinger Equations.” Computer Physics Communications, August 26, 2020, 107568. https://doi.org/10.1016/j.cpc.2020.107568.
* **[2]** Ledoux, Veerle, and Marnix Van Daele. “MATSLISE 2.0 : A Matlab Toolbox for Sturm-Liouville Computations.” ACM TRANSACTIONS ON MATHEMATICAL SOFTWARE 42, no. 4 (2016): 18.
* **[3]** Ixaru, L. Gr. “New Numerical Method for the Eigenvalue Problem of the 2D Schrödinger Equation.” Computer Physics Communications 181 (October 1, 2010): 1738–42. https://doi.org/10.1016/j.cpc.2010.06.031.
