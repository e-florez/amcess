import os
import pathlib

from setuptools import setup

# =============================================================================
# CONSTANTS
# =============================================================================

PATH = pathlib.Path(os.path.abspath(os.path.dirname(__file__)))


REQUIREMENTS = [
    "attrs>=23.2.0",
    "scipy>=1.12.0",
    "numpy==1.23.1",
    "pyscf>=2.5.0",
    "h5py>=3.1.0",
    "pyberny>=0.6.3",
    "geomeTRIC>=0.9.7.2",
    "GPyOpt>=1.2.6",
    "pyDOE>=0.3.8",
    "matplotlib>=3.8.3",
    "matplotlib-inline>=0.1.6",
    "py3Dmol>=2.0.4",
    "notebook>=6.5.6",
    "notebook_shim>=0.2.4",
    "jupyter>=1.0.0",
    "ipykernel>=6.29.3",
    "rise>=5.7.1",
]

with open(PATH / "amcess" / "__init__.py") as fp:
    for line in fp.readlines():
        if line.startswith("__version__ = "):
            VERSION = line.split("=", 1)[-1].replace('"', "").strip()
            break

with open("README.md", "r") as readme:
    LONG_DESCRIPTION = readme.read()


# =============================================================================
# FUNCTIONS
# =============================================================================

setup(
    name="amcess2024",
    version="0.1.2a20",
    author="""
    Edison Florez,
    Andy Zapata,
    Daniel Bajac,
    Alejandra Mendez,
    Cesar Ibarguen,
    José Aucar
    """,
    author_email="""
    edisonffhc@gmail.com,
    danianescobarv@gmail.com
    """,
    packages=["amcess", "amcess/data"],
    install_requires=REQUIREMENTS,
    license="The GPLv3 License",
    description="Atomic and Molecular Cluster Energy Surface Sampler",
    long_description_content_type="text/markdown",
    long_description=LONG_DESCRIPTION,
    url="https://gitlab.com/ADanianZE/amcess",
    keywords=[
        "Atomic Cluster",
        "Molecular Cluster",
        "optimization",
        "Potential Energy Surface",
        "PES",
        "Monte Carlo",
        "Simulated Annealing",
        "Bayesian Optimization",
    ],
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Education",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: OS Independent",
        "Programming Language :: Python",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: Implementation :: CPython",
        "Topic :: Scientific/Engineering",
    ],
    include_package_data=True,
)
