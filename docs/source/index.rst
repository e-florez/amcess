.. AMCESS documentation master file, created by
   sphinx-quickstart on Mon Dec 13 17:36:21 2021.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. image:: _static/amcess_logo_bl.png
   :width: 400px
   :height: 200px
   :scale: 80 %
   :alt: Atomic and Molecular Cluster Energy Surface Sample (AMCESS)
   :align: center

===================================================
Atomic and Molecular Cluster Energy Surface Sampler
===================================================

AMCCES is an object-oriented Python package to sample the Potential Energy 
Surface (PES) for atomic and molecular cluster. The package provides a simple, 
lightweight, and efficient platform for exploring the PES.


It can explore the PES by sampling the potential energy surface using 

   * Simulated annealing methods

   * Bayesian methods

AMCESS calculates the energy using `PySCF <https://pyscf.org/>`_, the 
Python-based Simulations of Chemistry Framework. PySCF is an efficient 
platform for quantum chemistry calculations that can be used to simulate the
properties of molecules, crystals, and custom Hamiltonians using mean-field 
and post-mean-field methods [1]_. 



.. rubric:: Footnotes

.. [1] To ensure ease of extensibility, almost all of the features in PySCF are implemented in Python, while computationally critical parts are implemented and optimized in C. Using this combined Python/C implementation, the package is as efficient as the best existing C or Fortran based quantum chemistry programs. In addition to its core libraries, PySCF supports a rich ecosystem of Extension modules.






.. toctree::
   :maxdepth: 1
   :caption: Contents:
   :titlesonly:   

   Home <self>
   about
   support
   modules

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
