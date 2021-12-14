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
   Basic Example <quickstart>
   AMCESS API <modules>
   support

AMCESS
======

is an object-oriented Python package to sample the Potential Energy 
Surface (PES) for atomic and molecular cluster. The main goal of the AMCESS 
project is to provide a simple, lightweight, and efficient platform to explore
the PES. The package is designed to be used in a Python environment. 

Here we present a user-friendly, simple, and efficient Python package to
sample the Potential Energy Surface (PES) for atomic and molecular cluster.
It is a fast and affordable tool for novice users and for experts who want
to produce high-quality PES samples. 

The package only requires a a minimal input information, such as the 
number of atoms and their coordinates. Despite its click-and-go appeal, 
the package is designed to be used and customized for any specific case. 

AMCESS can explore the the potential energy surface using the following 
techniques:

   * Stochastic methods: Simulated Annealing and Parallel Tempering

      * Metropolis or Delta Energy acceptance criteria   

   * Response methods: Bayesian Optimization
   * Combinatorial topology: simplicial homology global optimization (SHGO)

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

A previous version of AMCESS, called as ASCEC :cite:`ascec` (spanish acronym 
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
https://discord.gg/AxxTsfDJwxE)

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

