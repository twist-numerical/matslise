
Welcome to the documentation of Pyslise
#######################################

Pyslise is a Python package (written in C++) to solve
Schrödinger equations. In general the time-independent
Schrödinger equation is given by:

.. math::
  -\nabla^2 \varphi({\bf x}) + V({\bf x}) \varphi({\bf x}) = E \varphi({\bf x})

This package implements a constant perturbation (CP) method to
efficiently solve the one-dimensional Schrödinger equation:

.. math::
  -\frac{d^2}{dx^2}\varphi(x) + V(x) \varphi(x) = E \varphi(x)


And also the two-dimensional Schrödinger equation (experimental):

.. math::
  -\frac{\partial^2}{\partial x^2}\varphi(x, y) -\frac{\partial^2}{\partial y^2}\varphi(x, y) + V(x, y) \varphi(x, y) = E \varphi(x, y)

Besides a python package there is also an interactive web based GUI (compiled with WebAssembly):

- `GUI for the one-dimensional problem <./ti1d>`_

(Tested in a recent version of Firefox and Chrome)


Installation
------------

Installing PySlise is as easy as:

.. code::

  pip install pyslise

Make sure you are using a recent 64 bit version of python on a
Linux, Windows or macOS (10.15).

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
the software, come across bugs or have suggestions for improvements,
you can always contact us:

| Toon Baeyens
| toon.baeyens@ugent.be
| Ghent University
