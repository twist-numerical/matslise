Klotter problem
***************

..  contents::
    :local:
    :backlinks: top

Here we will be solving the Sturm-Liouville equation

.. math::
  \frac{d}{dx}\left[p(x)\frac{dy}{dx}\right] + q(x) y  = -\lambda w(x) y
  
with

.. math::
  p(x) &= 1 \\
  q(x) &= \frac{3}{4x^2} \\
  w(x) &= \frac{64\pi^2}{9x^6}

on the interval :math:`[8/7, 8]` and boundary conditions

.. math::
    y(\frac{8}{7}) = y(8) = 0\text{.}

.. code:: python

  import pyslise
  from math import pi, cos

  def p(x):
      return 1
  
  def q(x):
      return 3/(4*x**2)
  
  def w(x):
      return 64*pi**2/(9*x**6)

  problem = pyslise.SturmLiouville(p, q, w, 8/7, 8)
  left = (0, 1) # y(8/7) = 0, y'(8/7) = 1
  right = (0, 1) # y(8) = 0, y'(8) = 1
  eigenvalues = problem.eigenvaluesByIndex(0, 10, left, right)


The variable ``eigenvalues`` contains a list of tuples. Each tuple has as
the first element the index of the eigenvalue and as second argument the eigenvalue itself.

This data can be formatted in a nice table:

.. code:: python

  print('index  eigenvalue     error')
  for index, E in eigenvalues:
      error = problem.eigenvalueError(E, left, right)
      print(f'{index:>5} {E:>11.5f} {error:>9.1e}')


.. table::
   :align: right

.. code::
  
  index      eigenvalue
      0   1.00000000000
      1   4.00000000001
      2   9.00000000001
      3  16.00000000002
      4  25.00000000003
      5  36.00000000004
      6  49.00000000005
      7  64.00000000005
      8  81.00000000008
      9 100.00000000014

As you may suspect from these results, the exact eigenvalues if this problem are squares.


