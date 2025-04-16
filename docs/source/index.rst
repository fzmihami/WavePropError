.. WavePropError documentation master file, created by
   sphinx-quickstart on Tue Apr  1 15:56:14 2025.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

WavePropError documentation
===========================

WavePropError is a hybrid Python and C++ package for simulating wave propagation by solving the the one-dimensional Nwogu's equation :cite:`nwogu1993alternative`.\
It provides different numerical options and supports diffrent order of flux reconstruction, allowing for comaparative analysis of their accuracy and performance.\
The numerical solution of the Nwogu's equation is implemented in C++ for performance, while the Python interface provides a user-friendly way to set up and run simulations.\

The package supports four numerical schemes for solving the dispersive Nwogu's equation:

- **Conservative Staggered Scheme** :cite:`stelling2003accurate, stelling2003staggered, zijlema2019role`
- **HLLC Scheme** (Harten, Lax, van Leer, and Contact) :cite:`fraccarollo1995experimental, toro2001shock`
- **HLL Scheme** (Harten-Lax-van Leer) :cite:`harten1983upstream`
- **Central Upwind Scheme**  :cite:`kurganov2001semidiscrete, kurganov2007second, chertock2015well`

.. toctree::
   :maxdepth: 1
   :caption: Contents:

   modules/installation
   modules/scheme_description
   modules/scripts
   modules/code
   modules/examples
   

References
==========

.. bibliography:: references.bib
   :style: plain
