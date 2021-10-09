from copy import deepcopy
from typing import Type

import attr
import numpy as np
from scipy.spatial.transform import Rotation

from data.atomic_data import atomic_mass


@attr.s(frozen=True)
class Atom:
    """
    Representation of an individual atomas (<element> <X> <Y> <Z>)

    Examples
    --------
    >>>Atom(element='H', x=0, y=0, z=0)
    {'element': 'H', 'x': 0, 'y': 0, 'z': 0}

    >>>Atom('F', 0, 0, 1.97)
    {'element': 'F', 'x': 0, 'y': 0, 'z': 1.97}

    Returns
    -------
    atom : object
        object like dict {'element': str, 'x': float, 'y': float, 'z': float}

    Raises
    ------
    ValueError
        format MUST be (str, float, float, float) with NOT empty filed
    """

    element: str = attr.ib()
    x: int = attr.ib()
    y: int = attr.ib()
    z: int = attr.ib()

    # ===============================================================
    # VALIDATORS
    # ===============================================================
    @element.validator
    def _check_valid_element(self, element, value):
        if not value.isalnum():
            raise ValueError(
                "\n\nMust be valid NOT empty alphanumeric character"
                f"\nyou get --> '{value}'\n"
            )

    @x.validator
    @y.validator
    @z.validator
    def _check_valid_point(self, coordinate, value):
        if not isinstance(value, (int, float)):
            raise ValueError(
                "\n\nMust be valid NOT empty float,"
                f"\nyou get --> '{value}' with type: '{type(value).__name__}'"
            )

    # ===============================================================
    # PROPERTIES
    # ===============================================================
    @property
    def atomic_mass(self) -> list:
        return atomic_mass(self.element)

    @property
    def symbol(self) -> list:
        return self.element

    # ===============================================================
    # MAGIC METHODS
    # ===============================================================
    def __str__(self):
        return str(attr.asdict(self))


# -------------------------------------------------------------
@attr.s(frozen=True)
class Molecule:
    """
    Create a Molecule that is at least ONE atom.
    The format of the INPUT coordinates must be:

    {"atoms": [(<element> <X> <Y> <Z>), (<element> <X> <Y> <Z>), ...]}

    Parameters
    ----------
    atoms : list[tuple(str, float, float, float)]
        Cartesian coordinates of each atom, by default empty list

    charge : int
        total molecular/atomic charge, by default zero (0)

    multiplicity : int
        larger than zero, by defaul one (1)
    """

    atoms = attr.ib()
    charge: int = attr.ib(default=0)
    multiplicity: int = attr.ib(default=1)

    # ===============================================================
    # VALIDATORS
    # ===============================================================
    @atoms.validator
    def _cehck_valid_atoms(self, attribute, atoms):
        for line, atom in enumerate(atoms):
            try:
                Atom(*atom)
            except (ValueError, TypeError) as err:
                raise TypeError(
                    f"\n\n{err}\ncoordinates format must be a list of tuple: "
                    "[(str, float, float, float), ...]"
                    f"\ncheck atom number {line + 1} --> {atom}\n"
                    f"from --> {atoms}\n"
                )

    @charge.validator
    def _check_valid_charge(self, attribute, charge):
        if not isinstance(charge, int):
            raise ValueError(
                "\n\ncharge must be an integer "
                f"\nyou get --> 'charge = {charge}'\n"
            )

    @multiplicity.validator
    def _check_valid_multiplicity(self, attribute, multiplicity):
        if not isinstance(multiplicity, int) or multiplicity < 1:
            raise ValueError(
                "\n\nmultiplicity must be an integer larger than 0"
                f"\nyou get --> 'multiplicity = {multiplicity}'\n"
            )

    # ===============================================================
    # CONSTRUCTORS
    # ===============================================================
    @classmethod
    def from_dict(cls, atoms_dict):
        "Dictionary type: {'atoms': [(<element> <X> <Y> <Z>), ...]}"
        if not "atoms" in atoms_dict:
            raise KeyError(
                "\n\nThe key 'atoms' is casesensitive"
                "\n{'atoms': [(<element> <X> <Y> <Z>), ...]}"
                f"\nyou get {atoms_dict}\n"
            )
        atoms = atoms_dict.get("atoms")
        charge = atoms_dict.get("charge", 0)
        multiplicity = atoms_dict.get("multiplicity", 1)
        # atoms = [Atom(*atom) for atom in atoms]
        return cls(atoms, charge, multiplicity)

    # ===============================================================
    # MAGIC METHODS
    # ===============================================================
    def __str__(self):
        return str(attr.asdict(self))


# -------------------------------------------------------------
