<!-- [![Python 3.9](https://img.shields.io/badge/python-3.9-blue.svg)](https://www.python.org/downloads/release/python-390/)
![Gitlab pipeline status](https://img.shields.io/gitlab/pipeline/ADanianZE/amcess/main?style=plastic)
[![pipeline status](https://gitlab.com/ADanianZE/amcess/badges/main/pipeline.svg)](https://gitlab.com/ADanianZE/amcess/-/commits/main)
![Gitlab code coverage](https://img.shields.io/gitlab/coverage/ADanianZE/amcess/main?style=plastic)
[![Coverage Status](https://coveralls.io/repos/gitlab/ADanianZE/amcess/badge.svg?branch=main)](https://coveralls.io/gitlab/ADanianZE/amcess?branch=main)
[![coverage report](https://gitlab.com/ADanianZE/amcess/badges/main/coverage.svg)](https://gitlab.com/ADanianZE/amcess/-/commits/main)
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gl/ADanianZE%2Famcess/main?filepath=workflow%2F01_importing_atoms_and_molecules.ipynb)
![GitLab tag (latest by date)](https://img.shields.io/gitlab/v/tag/ADanianZE/amcess?style=plastic)
[![mypy: checked](http://www.mypy-lang.org/static/mypy_badge.svg)](http://mypy-lang.org/)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)
![Pod License](https://img.shields.io/badge/license-MIT-blue.svg)  -->

<!-- ![Atomic and Molecular Cluster Energy Surface Sampler](./docs/_static/amcess_logo.png) -->


<div align="center">
  <a href=http://mypy-lang.org/>
  <img src="http://www.mypy-lang.org/static/mypy_badge.svg"></a>
  <a href=https://github.com/psf/black>
  <img src="https://img.shields.io/badge/code%20style-black-000000.svg"></a>
  <br>
  <a href=https://www.python.org/downloads/release/python-390/>
  <img src=https://img.shields.io/badge/python-3.9-blue.svg></a>
  <a href=https://www.gnu.org/licenses/gpl-3.0>
  <img src="https://img.shields.io/badge/License-GPLv3-blue.svg"></a>
  <a href=https://mybinder.org/v2/gl/ADanianZE%2Famcess/main?filepath=workflow%2F01_importing_atoms_and_molecules.ipynb>
  <img src="https://mybinder.org/badge_logo.svg"></a>
  <br>  
  <a href=https://img.shields.io/gitlab/pipeline/ADanianZE/amcess/main?style=plastic>
  <img src="https://img.shields.io/gitlab/pipeline/ADanianZE/amcess/main?style=plastic"></a>
  <a href=https://gitlab.com/ADanianZE/amcess/-/commits/main>
  <img src="https://gitlab.com/ADanianZE/amcess/badges/main/pipeline.svg"></a>
  <a href=https://img.shields.io/gitlab/coverage/ADanianZE/amcess/main?style=plastic>
  <img src="https://img.shields.io/gitlab/coverage/ADanianZE/amcess/main?style=plastic"></a>
  <a href=https://coveralls.io/gitlab/ADanianZE/amcess?branch=main>
  <img src="https://coveralls.io/repos/gitlab/ADanianZE/amcess/badge.svg?branch=main"></a>
</div>

---

<div align="center">
  <a href="Atomic and Molecular Cluster Energy Surface Sampler">
  <img width="400" height="200" src="./docs/source/_static/amcess_logo.png"></a>
  <br>
</div>

---

<div align="center">
  <h1> Atomic and Molecular Cluster Energy Surface Sampler (AMCESS) </h1>
</div>

Tools to explore the Potential Energy Surface (PES) for atomic and molecular 
systems and generate candidate structures for the local and global minima.



<div align="center">
  <img width="200" height="80" src="./docs/source/_static/ibuprofen.png" VSPACE=50 HSPACE=10>
  <img width="350" height="200" src="./docs/source/_static/ibu_w6_white.gif" HSPACE=20>
  <br>
  Molecular cluster of ibuprofen and six water molecules
  [<a href="http://www.doi.org/10.1063/1.4874258">doi: 10.1063/1.4874258</a>]
  <br>
  <br>
</div>

<div align="center">
  <img width="300" height="200" src="./docs/source/_static/2d_surface_1.gif">
  <img width="300" height="200" src="./docs/source/_static/2d_surface_2.gif">
  <br>
  <img width="300" height="200" src="./docs/source/_static/2d_surface_3.gif">
  <img width="300" height="200" src="./docs/source/_static/2d_surface_4.gif">
  <br> 
  Example images taken from Deniz Yuret's blog on October 29, 2021 
  <br>
  from:
  <a href="http://www.denizyuret.com/2015/03/alec-radfords-animations-for.html">Alec Radford's animations for optimization algorithms</a>
</div>

### Description
Some of the technics implemented are:

 * Stochastic
      * Simulated Annealing (Metropolis-Monte Carlo)
      * Dual Annealing
* Bayesian Optimization

To compute the energy (cost-function), one could use some of these 
quantum chemistry packages:

- Non Relativistic: NWChem, Gaussian, DALTON,
GAMESS

- Relativistic: DIRAC, ADF

### Technical Documentation
Technical documents behind this project can be accessed [here](https://adanianze.gitlab.io/amcess).


### Requirements

First you should install the required python packages contained in the file `requirements.txt`. For developer, you should install `requirements_dev.txt`.

### Installation
AMCESS is writyten for **Python 3.9**

1. Install virtual environment:

    ```python -m venv venv```

2. Activate virtual environment:

    ```source venv/bin/activate```

3. Install the packages:
      ```
      >>>pip install -e .
      Installing collected packages: amcess
      Running setup.py develop for amcess
      Successfully installed amcess
      ```
4. Run AMCESS (check some examples below)
    
5. For developer only, install dependencies:

      ```pip install -r requirements.txt -r requirements_dev.txt```

6. Run all test:

      ``tox``

### Usage

A detail workflow is provide into `workflow` directory. It has a list of Jupyter notebook with detail examples about AMCESS tools and capabilities.

Workflow:
1. Getting starting with atoms and molecules properties. 
      - Notebook ([binder](https://mybinder.org/v2/gl/ADanianZE%2Famcess/main?filepath=workflow%2F01_importing_atoms_and_molecules.ipynb)): `01_importing_atoms_and_molecules.ipynb` 
2. Translating and rotating atoms and molecules. 
      - Notebook ([binder](https://mybinder.org/v2/gl/ADanianZE%2Famcess/main?filepath=workflow%2F02_move_rotate_molecules.ipynb)): `02_move_rotate_molecules.ipynb` 
3. Moving Molecules randomly from a Cluster. 
      - Notebook ([binder](https://mybinder.org/v2/gl/ADanianZE%2Famcess/main?filepath=workflow%2F03_move_rotate_cluster.ipynb)): `03_move_rotate_cluster.ipynb` 
4. Freezing any molecule and redefine its sphere center. 
      - Notebook ([binder](https://mybinder.org/v2/gl/ADanianZE%2Famcess/main?filepath=workflow%2F04_freeze_molecule_redefine_center.ipynb)): `04_freeze_molecule_redefine_center.ipynb` 
5. Initialize a cluster avoiding atomic overlapping
      - Notebook ([binder](https://mybinder.org/v2/gl/ADanianZE%2Famcess/main?filepath=workflow%2F05_initialize_cluster_and_move_molecule.ipynb)): `05_initialize_cluster_and_move_molecule.ipynb` 

## Support
Tell people where they can go to for help. It can be any combination of an issue tracker, a chat room, an email address, etc.

## Roadmap


## Contributing
State if you are open to contributions and what your requirements are for accepting them.

For people who want to make changes to your project, it's helpful to have some documentation on how to get started. Perhaps there is a script that they should run or some environment variables that they need to set. Make these steps explicit. These instructions could also be useful to your future self.

You can also document commands to lint the code or run tests. These steps help to ensure high code quality and reduce the likelihood that the changes inadvertently break something. Having instructions for running tests is especially helpful if it requires external setup, such as starting a Selenium server for testing in a browser.

### Licence
GNU General Public License v3 (GLPv3)

### Authors and Acknowledgment
Main authors: Alejandra Mendez, Juan Jose Aucar, Daniel Bajac, César Ibargüen, Andy Zapata, Edison Florez (_edisonffh@mail.com_)

### Project Status

Under development




---

## ASCEC (FORTRAN 77 version)

A previous version of AMCESS, called ASCEC [[1]](#1) (spanish acronym 
Annealing Simulado con Energía Cuántica) was written in FORTRAN77 and 
was successfully used in the wide range of research and academic applications. 
From atomic cluster to molecular cluster, the ASCEC package has produced 
novel results (structure never seen before) published in the literature. Read more on [ASCEC publications](https://scholar.google.com/scholar?start=0&q=%22ascec%22,+annealing&hl=en&as_sdt=0,5).


You could check the directory [ASCECV3](https://gitlab.com/ADanianZE/amcess/-/tree/main/ASCECV3)

```
ASCECV3/
|---papers/  
|---p_ascec/
|---examples/
      |---adf     
      |---dalton  
      |---g03
      |---gamess
      |---nwchem
```
### References
<div style=font-size:12px>
      <a id="1">[1]</a> 
      J Pérez and A Restrepo. Ascec v–02: annealing simulado con energía cuántica. Property, development and implementation: Grupo de Química–Física Teórica, Instituto de Química, Universidad de Antioquia: Medellín, Colombia, 2008.
<br>
</div>
