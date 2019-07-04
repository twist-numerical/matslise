Using PySE2d
============

Ixaru's problem
***************

As a first example we'll take a look at the potential:

.. math::
  V(x, y) = (1+x^2)(1+y^2)

on the domain :math:`[-5.5; 5.5]\times[-5.5; 5.5]`.

Using ``pyslise`` this problem can be easily worked on:

.. code:: python

  from pyslise import PySE2d

  def V(x, y):
      return (1 + x**2) * (1 + y**2)

  problem = PySE2d(V, -5.5,5.5, -5.5,5.5, x_tol=1e-6, y_count=25)

To find the closest eigenvalue in the neighborhood of ``5`` one can use:

.. code:: python

  problem.findEigenvalue(5)
  # 5.52674387

``pyslise`` is able to find the closest eigenvalue to a certain guess
because the implemented algorithm is able to calculate an error-matrix
that expresses of that given guess is an eigenvalue. In turn, this matrix
can be used to improve that initial guess. Until a sufficiently accurate
estimate is found.

There is method implemented to find all eigenvalues of the two-dimensional
Schrödinger equation in a certain interval. But, as of yet, this method isn't
perfect. It is based on a few heuristics to 'guess' that all eigenvalues are
found. This heuristic is implemented in ``.findEigenvalues(Emin, Emax)``:

.. code :: python

  problem.findEigenvalues(0,13)

.. code ::

  [3.1959180850800877,
   5.526743864002774,
   5.526743877339405,
   7.55780334350954,
   8.031272354757498,
   8.444581365360518,
   9.92806092943532,
   9.928061003007299,
   11.311817072021494,
   11.31181710814099,
   12.103256481915713,
   12.201180767897501]
