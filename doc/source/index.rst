
Welcome to PySlise's documentation!
###################################

PySlise is a Python package (written in C++) to solve Schrödinger equations. In general the time-independent Schrödinger equation is given by:

.. math::
  -\nabla^2 \varphi({\bf x}) + V({\bf x}) \varphi({\bf x}) = E \varphi({\bf x})

This package implements a constant perturbation (CP) method to efficiently solve the one-dimensional Schrödinger equation:

.. math::
  -\frac{d}{dx}\varphi(x) + V(x) \varphi(x) = E \varphi(x)

Because this CP-method is the crux of this package, the class for the one-dimensional case caries the same name: `PySlise <api/pyslise.html>`__.



.. toctree::
   :maxdepth: 2

   examples/index
   api/index



Indices and tables
------------------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
