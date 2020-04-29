Pyslise
=======

These two classes are able to solve the one dimensional time
independent Schrödinger equation. The first one (Pyslise)
solves the general case. The second class (PysliseHalf)
implements half range reduction, for this it assumes the
problem is symmetric.

.. autoclass:: pyslise.Pyslise
   :members:
   :inherited-members:
   :undoc-members:
   
   .. automethod:: __init__

Half range reduction can be applied when the problem is
symmetric. This symmetry has following requirements:

- The potential is even: V(x) = V(-x). (The potential will only
  be evaluated on the positive part of the domain.)
- The domain is symmetric around zero: [-max; max].
- The boundary conditions are symmetric left and right.

.. autoclass:: pyslise.PysliseHalf
  :members:
  :inherited-members:
  :undoc-members:

   .. automethod:: __init__
