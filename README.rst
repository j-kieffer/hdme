Description
===========

This repository contains a C library for evaluating higher-dimensional
modular equations by analytic methods. As it stands now, the following
types of modular equations in genus g=2 are supported, either over the
rationals or prime finite fields:

- Modular equations of Siegel type;

- Modular equations of Hilbert type, using parametrizations of Humbert
  surfaces for all fundamental discriminants less than 100, and
  Hilbert surfaces for all fundamental discriminants less than 20.

The choice of coordinates is dynamic, so that all corner cases (extra
endomorphisms, vanishing of I4...) can be covered.

Along the way, the library provides functionality for the evaluation
of genus 2 theta constants in quasi-linear time, and for the action of
Hecke operators on Siegel modular forms. The algorithms are described
in the associated paper `[K21]`_.

Prerequisites
=============

An installation of `Flint`_ version 3 or later.

Installation
============

::

  git clone https://github.com/j-kieffer/hdme.git && cd hdme
  ./configure                        # ./configure --help for more options
  make
  make tests                         # optional
  make check                         # optional
  make install                       # optional

Common configuration examples::

  ./configure --with-flint=/opt/homebrew --with-mpfr=/opt/homebrew --with-gmp=/opt/homebrew
  ./configure --prefix=$HOME/local --enable-debug
  ./configure --disable-openmp

Documentation
=============

A detailed documentation for each function existed in previous
versions, but is now outdated. For examples of usage, look at the test
files in the `test` subfolder of each module.

Credits
=======

- The code for fast evaluation of theta constants using Arb is
  inspired by the `cmh`_ library by Andreas Enge and Emmanuel Thomé.

- The code written by `Enea Milio`_ has been useful in the
  implementation of Mestre's algorithm and to implement the inversion
  of the Hilbert embedding.

- The design and semantics of many functions are inspired from existing
  functions in the libraries Flint and Arb, maintained by William Hart
  and Fredrik Johansson.
  
.. _[K21]: https://arxiv.org/abs/2010.10094
.. _Flint: https://flintlib.org
.. _Arb: https://arblib.org
.. _cmh: https://gitlab.inria.fr/cmh/cmh
.. _Enea Milio: https://members.loria.fr/EMilio/modular-polynomials/
