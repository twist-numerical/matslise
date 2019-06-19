Using PySlise
#############

.. code:: python

  from pyslise import PySlise
  from math import pi, cos

  problem = PySlise(lambda x: 2*cos(2*x),
                    -pi/2, pi/2, tolerance=1e-5)
  problem.computeEigenvaluesByIndex(0, 10, [0,1], [0,1])
