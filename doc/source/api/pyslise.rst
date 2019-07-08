PySlise
=======

These two classes are able to solve the one dimensional time
independent Schrödinger equation. The first one (PySlise)
solves the general case. The second class (PySliseHalf)
implements half range reduction, for this it assumes the
problem is symmetric.

.. autoclass:: pyslise.PySlise
   :members:
   :undoc-members:

   .. automethod:: __init__

Half range reduction can be applied when the problem is
symmetric. This symmetry has following requirements:

- The potential is even: V(x) = V(-x). The potential will only
  be evaluated on the positive part of the domain.
- The domain is symmetric around zero: [-max; max].
- The boundary conditions are the same left and right.

.. autoclass:: pyslise.PySliseHalf
  :members:
  :undoc-members:

   .. automethod:: __init__
