import os
import pathlib

from setuptools import setup

# =============================================================================
# CONSTANTS
# =============================================================================

PATH = pathlib.Path(os.path.abspath(os.path.dirname(__file__)))


REQUIREMENTS = [
    "h5py==3.1.0",
    "scipy>=1.7",
    "numpy>=1.21",
    "pyscf>=1.7.6.post1",
]

with open(PATH / "src" / "__init__.py") as fp:
    for line in fp.readlines():
        if line.startswith("__version__ = "):
            VERSION = line.split("=", 1)[-1].replace('"', "").strip()
            break

VERSION = 0.1

with open("README.md", "r") as readme:
    LONG_DESCRIPTION = readme.read()


# =============================================================================
# FUNCTIONS
# =============================================================================

setup(
    name="amcess",
    version="0.1.0",
    author="Edison Florez",
    author_email="edisonffhc@gmail.com",
    packages=["src"],
    scripts=["main.py"],
    url="http://pypi.python.org/pypi/amcess/",
    license="The MIT License",
    description="Atomic and Molecular Cluster Energy Surface Sampler",
    long_description=LONG_DESCRIPTION,
    keywords=["optimization", "PES", "Potential Energy Surface", "Monte Carlo"],
    install_requires=[
	"h5py==3.1.0",
        "scipy>=1.7",
        "numpy>=1.21",
        "pyscf>=1.7.6.post1",
    ],
)


# if __name__ == "__main__":
#     setup()
