import os
import pathlib

from setuptools import setup

# =============================================================================
# CONSTANTS
# =============================================================================

PATH = pathlib.Path(os.path.abspath(os.path.dirname(__file__)))


REQUIREMENTS = [
    "scipy==1.7.1",
    "numpy==1.21.2",
    "pyscf==1.7.6.post1",
    "h5py==2.10.0",
    "pyberny==0.6.3",
    "geomeTRIC==0.9.7.2",
    "GPyOpt==1.2.6",
    "pyDOE==0.3.8",
    "matplotlib==3.4.2",
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
    name="amcess",
    version="0.1.0",
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
    packages=["amcess"],
    install_requires=REQUIREMENTS,
    url="https://gitlab.com/ADanianZE/amcess",
    license="The GPLv3 License",
    description="Atomic and Molecular Cluster Energy Surface Sampler",
    long_description=LONG_DESCRIPTION,
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
        "License :: OSI Approved :: GPLv3 License",
        "Operating System :: OS Independent",
        "Programming Language :: Python",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: Implementation :: CPython",
        "Topic :: Scientific/Engineering",
    ],
    include_package_data=True,
)
