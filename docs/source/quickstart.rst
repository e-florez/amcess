
Installation
************

If you have pip installed:

.. code-block:: bash

    $ pip install amcess

Or, from source code:

.. code-block:: bash

    $ git clone git@gitlab.com:ADanianZE/amcess.git
    $ cd amcess
    $ pip install -e .


Basic example
=============

In the following example the end-to-end candidate structures for a water
pentamer cluster are generated. The PES is sampled using the PySCF to 
compute the energy at Hartree-Fock level and STO basis set.

.. code-block:: python

    from amcess import Molecule, Cluster
    from amcess.search_engine import SearchConfig
    
    water_dict = {
        "atoms": [("O", 0, 0, 0), ("H", 0.587, 0.757, 0), ("H", -0.587, 0.757, 0)],
        "charge": 0,
        "multiplicity": 1,
    }

    # Create a water molecule object
    water_mol = Molecule.from_dict(water_dict)

    # Create a water pentamer cluster object
    water_pentamer = Cluster(5 * water_mol)

    # Countour Conditions: define a sphere radius (in Angstrom) for the cluster
    water_pentamer.sphere_radius = 2.0

    # define a sphere center
    water_pentamer.sphere_center = (3.0, 3.0, 3.0)

    # Simulated Annealing search   
    obj_sa = SearchConfig(water_pentamer)
    obj_sa.run(initial_temperature=1000, step=0.2, number_temperatures=10, max_cycles=20)
    
    # dual annealing search
    obj_da = SearchConfig(water_pentamer, search_methodology="dual_annealing")
    obj_da.run(maxfun=10,maxiter=10)
    
    # shgo searching
    obj_shgo = SearchConfig(water_pentamer, search_methodology="SHGO")
    obj_shgo.run(sampling_method="sobol", n=2)
    
    # Bayesian optimization search
    obj_bayesian = SearchConfig(water_pentamer, search_methodology="Bayesian")
    obj_bayesian.run(initer=10, maxiter=10)

To find out what else you can do, head over to learning AMCESS to have a look 
at the tutorials on `Binder <https://mybinder.org/v2/gl/ADanianZE%2Famcess/main?labpath=https%3A%2F%2Fgitlab.com%2FADanianZE%2Famcess%2F-%2Ftree%2Fmain%2Fworkflow>`_.


.. https://mybinder.org/v2/gl/ADanianZE%2Famcess/main?filepath=workflow%2F01_importing_atoms_and_molecules.ipynb

.. https://mybinder.org/v2/gl/ADanianZE%2Famcess/main?labpath=https%3A%2F%2Fgitlab.com%2FADanianZE%2Famcess%2F-%2Fblob%2Fmain%2Fworkflow%2F01_importing_atoms_and_molecules.ipynb

====================
Development versions
====================

To install development versions of AMCESS, you can compile it from source. 
In order to install from source, you will need a python3.9 interpreter and

    * tox==3.24.3
    * pytest==6.2.5
    * pytest-cov==2.12.1
    * coverage==5.5
    * black==21.8b0
    * flake8==3.9.2
    * flake8-black==0.2.3
    * flake8-builtins==1.5.3
    * flake8-isort==4.0.0
    * mypy==0.910 
    * sphinx==4.1.2
    * sphinx_rtd_theme==0.5.1
    * myst-parser==0.15.2
    * docutils==0.16
    * sphinxcontrib-bibtex==2.4

assuming you have already installed required dependencies, then you can compile it:

.. code-block:: bash

    $ git clone git@gitlab.com:ADanianZE/amcess.git
    $ cd amcess
    $ pip install -e .


Testing
=======

.. code-block:: bash

    $ pytest -v tests/ --cov=amcess/ --cov-fail-under 90 --cov-report term-missing

