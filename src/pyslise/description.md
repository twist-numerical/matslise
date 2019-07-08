# PySlise ${PYSLISE_VERSION}

PySlise is a collection of algorithms to solve the one and two dimensional time-independent Schrödinger equations. These algorithms are based upon constant perturbation methods to efficiently solve these eigenvalue problems.

The code (and name) is based on Matslise [1]. This is a feature-rich matlab library for solving the one dimensional time independent Schrödinger equation.

To solve the two dimensional problem an algorithm is developed on the basis of a method proposed by Ixaru [2].

This implementation is developed in C++ with a focus on efficiency. This code is precompiled for 64 bit windows and linux and packaged in wheels.


## Documentation

Full documentation can be found on 
[matslise.ugent.be](https://matslise.ugent.be/). This document contains some examples of how to use this library.

On the same page an interactive version is available.

## Examples

One dimensional problems can be tackled with:
```python
from pyslise import PySlise
from math import pi, cos

problem = PySlise(lambda x: 2*cos(2*x), 0, pi, tolerance=1e-5)
problem.eigenvaluesByIndex(0, 10, (0, 1), (0, 1))
```

Also two dimensional problems are possible:
```python
from pyslise import PySE2d

def V(x, y):
    return (1 + x**2) * (1 + y**2)

problem = PySE2d(V, -5.5,5.5, -5.5,5.5, tolerance=1e-5)
problem.eigenvalues(0,13)
```

## Bibliography
* **[1]** Ledoux, Veerle, and Marnix Van Daele. “MATSLISE 2.0 : A Matlab Toolbox for Sturm-Liouville Computations.” ACM TRANSACTIONS ON MATHEMATICAL SOFTWARE 42, no. 4 (2016): 18.
* **[2]** Ixaru, L. Gr. “New Numerical Method for the Eigenvalue Problem of the 2D Schrödinger Equation.” Computer Physics Communications 181 (October 1, 2010): 1738–42. https://doi.org/10.1016/j.cpc.2010.06.031.