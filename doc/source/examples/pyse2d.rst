Using PySE2d
#############

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

As of yet there isn't a method to automatically find eigenvalues of the two-dimensional Schrödinger equation. But ``pyslise`` is able to calculate an error-matrix that expresses of a given guess is an eigenvalue. This matrix can be used to improve that initial guess. This is what ``.findEigenvalue(E)`` does, it is able to find the closest eigenvalue to a certain guess ``E``.

To find the closest eigenvalue in the neighborhood of ``5`` one can use:

.. code:: python

  problem.findEigenvalue(5)
  # 5.52674387
