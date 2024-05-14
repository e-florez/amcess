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


.. raw:: html

   <div>
      <h1 style=font-size:45px align="center"> 
      Atomic and Molecular Cluster Energy Surface Sampler 
   </h1>
   </div>

.. toctree::
   :maxdepth: 2
   :caption: Contents:
   :titlesonly:   

   Home <self>
   About <about>
   Input <input>
   Box Bounds <boxbounds>
   Candidate Structures <candidatestructures>
   Electronic Structure <electronicstructure>
   Basic Example <quickstart>
   AMCESS API <modules>
   support

AMCESS
======

Exploration of the Potential Energy Surface (PES) of molecules or atoms clusters is
a crucial step to analyze physical–chemistry properties and processes. The Atomic and
Molecular Energy Surface Sampler (AMCESS) is an end-to-end package implemented
in Python 3.9 to generate candidate structures for the critical points sampling of the
PES. 

The amcess package uses simple input files and automates common procedures to
explore the PES using the Simulated Annealing, Simplicial Homology Global Optimiza-
tion (SHGO), and Bayesian Optimization to generate candidate structures for any kind
of critical point, such as local minima or transition states. The package also 
allows the user to perform local searches around defined regions. 

The PES is generated computing the electronic energy using standard and powerful
quantum chemistry packages such as PySCF and Psi4, also implemented in Python. 

The amcess main purpose is to be a user friendly package, easy to install, 
import, and run, available in most platforms and open-source. 
As a Python module, amcess can be integrated into any workflow. This package 
has code reviews with unit testing and continuous integration, code coverage
tools, and automatically keeps documentation up–to–date. 

AMCESS calculates the energy using `PySCF <https://pyscf.org/>`_, the 
Python-based Simulations of Chemistry Framework. PySCF is an efficient 
platform for quantum chemistry calculations that can be used to simulate the
properties of molecules and crystals using mean-field and post-mean-field 
methods [1]_. 

In addition to the above, AMCESS can also calculate the energy by using any 
input function that takes in a array of coordinates (Atom type and coordinates) 
and returns the energy of each configuration.

The basic (:doc:`examples <quickstart>`) demonstrates some of these 
features. Read more on :doc:`Installation <quickstart>`.

A previous version of AMCESS, called ASCEC :cite:`ascec` (spanish acronym 
Annealing Simulado con Energía Cuántica) was written in FORTRAN77 and 
was successfully used in the wide range of research and academic applications. 
From atomic cluster to molecular cluster, the ASCEC package has produced 
novel results (structure never seen before) published in the literature. 
Some of the most notable results are:

   1. Structural studies of the water tetramer :cite:`w4`
   2. Lithium and Bimetallic Lithium-Sodium Clusters :cite:`li-na_cluster`
   3. Structural Characterization of the (Methanol) :math:`_4` Potential Energy Surface :cite:`methanol4`
   4. Insights into the structure and stability of the carbonic acid dimer :cite:`ca_dimer`
   5. Structural studies of the water hexamer :cite:`w6`
   6. Understanding microsolvation of Li+: structural and energetical analyses :cite:`li+`
   7. Microsolvation of dimethylphosphate: a molecular model for the interaction of cell membranes with water :cite:`dmp`
   8. A combined experimental and computational study of the molecular interactions between anionic ibuprofen and water :cite:`ibuprofen`
   9. Microsolvation of methylmercury: structures, energies, bonding and NMR constants (199 Hg, 13 C and 17 O) :cite:`florez`
   10. How many water molecules does it take to dissociate HCl? :cite:`hcl`

Read more on `ASCEC publications <https://scholar.google.com/scholar?start=0&q=%22ascec%22,+annealing&hl=en&as_sdt=0,5>`_.

Availability
============

AMCESS can be easily installed with its dependencies using the pip or conda 
package managers. All source code is available under the GNU General Public 
License, version 3 from `GitLab <https://gitlab.com/ADanianZE/amcess>`_
and the `Python Package index <pypi.org/project/amcess>`_.

Participating
=============

Ask questions on the AMCESS Discord Server to talk with other users and 
developers. (In order to join our Discord server, use the invitation link 
https://discord.gg/vxQQCjpg)

Please report bugs or enhancement requests through the Issue Tracker.

AMCESS is open source and welcomes your contributions. Fork the repository on 
GitLab and submit a pull request. Participate on the developer Discord Server.



.. rubric:: Footnotes

.. [1] All of the features in PySCF are implemented in Python, while computationally critical parts are implemented and optimized in C. Using this combined Python/C implementation, the package is as efficient as the best existing C or Fortran based quantum chemistry programs.

References
==========
.. bibliography:: ./refs.bib
   :style: unsrt


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

