
Welcome to the documentation of Pyslise
#######################################

Pyslise is a Python package (written in C++) to solve one-dimensional time-independent Schrödinger equations and Sturm-Liouville equations.

This package implements a constant perturbation (CP) method to
efficiently solve the one-dimensional Schrödinger equation:

.. math::
  -\frac{d^2}{dx^2}\varphi(x) + V(x) \varphi(x) = E \varphi(x)



Sturm-Liouville equations

.. math::
  \frac{d}{dx}\left[p(x)\frac{dy}{dx}\right] + q(x) y  = -\lambda w(x) y

are solved by applying a Liouville transformation and then solving the Schrödinger equation.

Besides a python package there is also an interactive web based GUI (compiled with WebAssembly):

- `GUI for the one-dimensional problem <./ti1d>`_

(Tested in a recent version of Firefox and Chrome)


Documentation
-------------

This documentation consists of two main sections: a collection of examples and a reference manual.

.. toctree::
   :maxdepth: 2

   examples/index
   api/index

* :ref:`genindex`


Installation
------------

Installing PySlise is as easy as:

.. code::

  pip install pyslise

Make sure you are using a recent 64 bit version of python on a
Linux, Windows or macOS (10.15).


Contact
-------

This package is still in it's infancy and under active development.
So naturally it will not be perfect. When you have questions about
the software, come across bugs or have suggestions for improvements,
you can always contact us:

| Toon Baeyens
| toon.baeyens@ugent.be
| Ghent University
