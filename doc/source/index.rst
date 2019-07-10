
Welcome to PySlise's documentation!
###################################

PySlise is a Python package (written in C++) to solve
Schrödinger equations. In general the time-independent
Schrödinger equation is given by:

.. math::
  -\nabla^2 \varphi({\bf x}) + V({\bf x}) \varphi({\bf x}) = E \varphi({\bf x})

This package implements a constant perturbation (CP) method to
efficiently solve the one-dimensional Schrödinger equation:

.. math::
  -\frac{d^2}{dx^2}\varphi(x) + V(x) \varphi(x) = E \varphi(x)

Installation
------------

Installing PySlise is as easy as:

.. code::

  pip install pyslise

Make sure you are using a recent 64 bit version of python on a
Linux, Windows or Mac.

Documentation
-------------

.. toctree::
   :maxdepth: 2

   examples/index
   api/index

* :ref:`genindex`

Contact
-------

This package is still in it's infancy and under active development.
So naturally it will not be perfect. When you have questions about
the software, come across bugs or have suggestions for imrovements,
you can always contact us:

| Toon Baeyens
| toon.baeyens@ugent.be
| Ghent University
