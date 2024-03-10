<!-- [![Python 3.9](https://img.shields.io/badge/python-3.9-blue.svg)](https://www.python.org/downloads/release/python-390/)
![Gitlab pipeline status](https://img.shields.io/gitlab/pipeline/ADanianZE/ascec/main?style=plastic)
[![pipeline status](https://gitlab.com/ADanianZE/ascec/badges/main/pipeline.svg)](https://gitlab.com/ADanianZE/ascec/-/commits/main)
![Gitlab code coverage](https://img.shields.io/gitlab/coverage/ADanianZE/ascec/main?style=plastic)
[![Coverage Status](https://coveralls.io/repos/gitlab/ADanianZE/ascec/badge.svg?branch=main)](https://coveralls.io/gitlab/ADanianZE/ascec?branch=main)
[![coverage report](https://gitlab.com/ADanianZE/ascec/badges/main/coverage.svg)](https://gitlab.com/ADanianZE/ascec/-/commits/main)
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gl/ADanianZE%2Fascec/main?filepath=workflow%2F01_importing_atoms_and_molecules.ipynb)
![GitLab tag (latest by date)](https://img.shields.io/gitlab/v/tag/ADanianZE/ascec?style=plastic)
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
  <a href=https://mybinder.org/v2/gl/ADanianZE%2Fascec/main?filepath=workflow%2F01_importing_atoms_and_molecules.ipynb>
  <img src="https://mybinder.org/badge_logo.svg"></a>
  <br>  
  <a href=https://img.shields.io/gitlab/pipeline/ADanianZE/ascec/main?style=plastic>
  <img src="https://img.shields.io/gitlab/pipeline/ADanianZE/ascec/main?style=plastic"></a>
  <a href=https://gitlab.com/ADanianZE/ascec/-/commits/main>
  <img src="https://gitlab.com/ADanianZE/ascec/badges/main/pipeline.svg"></a>
  <a href=https://img.shields.io/gitlab/coverage/ADanianZE/ascec/main?style=plastic>
  <img src="https://img.shields.io/gitlab/coverage/ADanianZE/ascec/main?style=plastic"></a>
  <a href=https://coveralls.io/gitlab/ADanianZE/ascec?branch=main>
  <img src="https://coveralls.io/repos/gitlab/ADanianZE/ascec/badge.svg?branch=main"></a>
</div>

---

<div align="center">
  <a href="Atomic and Molecular Cluster Energy Surface Sampler">
  <img width="400" height="200" src="https://gitlab.com/ADanianZE/amcess/-/raw/main/docs/source/_static/amcess_logo.png"></a>
  <br>
</div>

---

<div align="center">
  <h1> Atomic and Molecular Cluster Energy Surface Sampler (AMCESS) </h1>
</div>

Exploration of the Potential Energy Surface (PES) of molecular or atomic clusters is
a crucial step in analyzing physical–chemistry properties and processes. The Atomic and
Molecular Energy Surface Sampler (AMCESS) is an end-to-end package implemented
in Python 3.9 to generate candidate structures for the critical points sampling of the
PES. The amcess main purpose is to be a
user-friendly package, easy to install, import, and run, available on most platforms and
open-source. As a Python module, amcess can be integrated into any workflow. This
package has code reviews with unit testing and continuous integration, code coverage
tools, and automatically keeps documentation up–to–date.

<div align="center">
  <img width="200" height="80" src="https://gitlab.com/ADanianZE/amcess/-/raw/main/docs/source/_static/ibuprofen.png" VSPACE=50 HSPACE=10>
  <img width="350" height="200" src="https://gitlab.com/ADanianZE/amcess/-/raw/main/docs/source/_static/ibu_w6_white.gif" HSPACE=20>
  <br>
  Molecular cluster of ibuprofen and six water molecules
  [<a href="http://www.doi.org/10.1063/1.4874258">doi: 10.1063/1.4874258</a>]
  <br>
  <br>
</div>

### Description

The amcess package uses simple input files and automates common procedures to explore the PES using the Simulated Annealing, Simplicial Homology Global Optimization (SHGO), and Bayesian Optimization to generate candidate structures for any kind of critical points, such as local minima or transition states. The package also allows the user to perform local searches around defined regions. The PES is generated by computing the electronic energy using standard and powerful quantum chemistry packages such as PySCF and Psi4, also implemented in Python.

Read more about the project in the [documentation](https://github.com/e-florez/amcess/blob/main/docs/white_paper/amcess_white_paper.pdf).

### Technical Documentation

Technical documents behind this project can be accessed [here](https://adanianze.gitlab.io/amcess).

### Requirements

First, you should install the required Python packages

```
  attrs>=23.2.0
  scipy>=1.12.0
  numpy==1.23.1
  pyscf>=2.5.0
  h5py>=3.1.0
  pyberny>=0.6.3
  geomeTRIC>=0.9.7.2
  GPyOpt>=1.2.6
  pyDOE>=0.3.8
  matplotlib>=3.8.3
  matplotlib-inline>=0.1.6
```

To use Jupyter Notebook

```
  py3Dmol>=2.0.4
  notebook>=6.5.6
  notebook_shim>=0.2.4
  jupyter>=1.0.0
  ipykernel>=6.29.3
  rise>=5.7.1
```

check the file `requirements.txt`. For developers, you should install `requirements_dev.txt`.

### Installation

AMCESS is a **Python 3.9** package

1. Install virtual environment:

   `python -m venv <VirtualEnvName>`

2. Activate virtual environment:

   `source <VirtualEnvName>`/bin/activate`

3. Install the packages:

   1. after cloning this repository

      `pip install .`

   2. or directly from PyPI

      `pip install amcess`

4. Run AMCESS (check some examples below). Additionally, you can install `requirements.txt` to use Jupyter notebooks.

For developers,

1. install dependencies:

   `pip install -r requirements.txt -r requirements_dev.txt`

2. Run all test:

   `tox`

### Usage

A detailed workflow is provided in the `workflow` directory. It has a list of Jupyter notebooks with detailed examples of AMCESS tools and capabilities.

Workflow ([Binder](https://mybinder.org/v2/gl/ADanianZE%2Famcess/main?labpath=https%3A%2F%2Fgitlab.com%2FADanianZE%2Famcess%2F-%2Ftree%2Fmain%2Fworkflow)):

1. Getting started with atoms and molecules properties.
   - `01_importing_atoms_and_molecules.ipynb`
2. Translating and rotating atoms and molecules.
   - `02_move_rotate_molecules.ipynb`
3. Moving Molecules randomly from a Cluster.
   - `03_move_rotate_cluster.ipynb`
4. Freezing any molecule and redefining its sphere center.
   - `04_freeze_molecule_redefine_center.ipynb`
5. Initialize a cluster avoiding atomic overlapping
   - `05_initialize_cluster_and_move_molecule.ipynb`

## Roadmap

Some of the ideas to keep growing are:

- Integration with **RDKit** (multiple format input)
- Results: geometrical analysis (clustering, k-nearest, k-means, etc.)

## Contributing

The easiest way to get help with the project is to join the #amcess
channel on Discord.

We hang out there and you can get real-time help with your projects.

- GitHub: https://github.com/e-florez/amcess/issues

- Discord: https://discord.gg/vxQQCjpg

### License

GNU General Public License v3 (GLPv3)

### Authors and Acknowledgment

Main authors: Alejandra Mendez, Juan Jose Aucar, Daniel Bajac, César Ibargüen, Andy Zapata, Edison Florez (_edisonffh@mail.com_)

### Project Status

Under development

---

## ASCEC (FORTRAN 77 version)

A previous version of AMCESS, called ASCEC [[1]](#1) (Spanish acronym
Annealing Simulado con Energía Cuántica) was written in FORTRAN77 and
was successfully used in a wide range of research and academic applications.
From atomic cluster to molecular cluster, the ASCEC package has produced
novel results (structure never seen before) published in the literature. Read more on [ASCEC publications](https://scholar.google.com/scholar?start=0&q=%22ascec%22,+annealing&hl=en&as_sdt=0,5).

You could check the directory [ASCECV3](https://gitlab.com/ADanianZE/ascec/-/tree/main/ASCECV3)

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
