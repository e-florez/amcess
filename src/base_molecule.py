from copy import deepcopy

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
@attr.s(frozen=False)
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

    _atoms = attr.ib()
    charge: int = attr.ib(default=0)
    multiplicity: int = attr.ib(default=1)

    # ===============================================================
    # VALIDATORS
    # ===============================================================
    @_atoms.validator
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
        return cls(atoms, charge, multiplicity)

    # ===============================================================
    # MAGIC METHODS
    # ===============================================================
    def __str__(self):
        return str(attr.asdict(self))

    # ===============================================================
    # PROPERTIES
    # ===============================================================
    @property
    def atoms(self) -> list:
        return self._atoms

    @property
    def atomic_masses(self) -> list:
        return [Atom(*atom).atomic_mass for atom in self.atoms]

    @property
    def coordinates(self) -> list:
        return [c[1:] for c in self.atoms]

    @property
    def elements(self) -> list:
        """Show a list of unique symbols

        Returns
        -------
        list
            list of unique symbols
        """
        return list(set(self.symbols))

    @property
    def numbering_atoms(self) -> str:
        """show atom number line by line

        Returns
        -------
        str
            atom number line by line
        """
        numbered_atoms = list()
        for i in range(self.total_atoms):
            line = list(self.atoms[i])
            line[0] = f"\r  atom #{i} --> {line[0]:<6}"
            line[1] = f"{line[1]:> 15.8f}"
            line[2] = f"{line[2]:> 15.8f}"
            line[3] = f"{line[3]:> 15.8f}"
            numbered_atoms.append("".join(line))

        return "\n".join(numbered_atoms)

    @property
    def symbols(self) -> list:
        return [str(s[0]).title() for s in self.atoms]

    @property
    def total_atoms(self) -> int:
        return len(self.atoms)

    @property
    def total_mass(self) -> float:
        return sum(self.atomic_masses)

    @property
    def center_of_mass(self) -> tuple:
        """Center of mass for a N-body problem. `Jacobi coordinates`_

        Notes
        -----
            total mass for dummy atoms (not in the Periodic Table) is equal
            to ONE (1)

        Returns
        -------
        tuple : (float, float, float)
            List of N 3D tuples, where N is equal to the number of atoms

        .. _Jacobi coordinates:
            https://en.wikipedia.org/wiki/Jacobi_coordinates
        """

        total_mass = 1 if not self.total_mass else self.total_mass

        return tuple(
            np.dot(
                np.asarray(self.atomic_masses),
                np.asarray(self.coordinates),
            )
            / total_mass
        )

    @property
    def principal_axes(self) -> list:
        return [
            tuple(c)
            for c in (
                np.asarray(self.coordinates) - np.asarray(self.center_of_mass)
            )
        ]

    @property
    def xyz(self) -> str:
        """Printing Molecule coordinates using XYZ format"""
        comments = (
            f"-- charge= {self.charge:<-g} and "
            f"multiplicity= {self.multiplicity:<g} --"
        )
        write_xyz = f"""\t{self.total_atoms}\n{comments:<s}\n"""
        for atom in self.atoms:
            write_xyz += f"""{atom[0]:<6}"""
            write_xyz += f"""\t{atom[1]:> 15.8f}"""
            write_xyz += f"""\t{atom[2]:> 15.8f}"""
            write_xyz += f"""\t{atom[3]:> 15.8f}\n"""

        return write_xyz


# -------------------------------------------------------------
class Cluster(Molecule):
    """
    Create a Cluster with molecules/atoms to move and rotate
    The format of the INPUT coordinates is as follows (any):

    1. Dictionary type: {"atoms": [(<element> <X> <Y> <Z>), ...]}
    2. List type: [(<element> <X> <Y> <Z>), ...]
    3. Molecule/Cluster type (Objects)

    Parameters
    ----------
    args : List, Dict, Molecule, Cluster
        coordinates of each molecule/atom comma separates (support +,-,*)
    frozen_molecule : integer, optional
        fixing molecule to NOT move or rotate, by default NEGATIVE
        integer means all molecules can be moved freely
    sphere_radius : float, optional
        radius for the spherical boundary condition, by default None
    sphere_center : tuple, optional
        Center of the sphere, by default (0, 0, 0)

    Raises
    ------
    TypeError
        [description]
    """

    def __init__(
        self,
        *args,
        frozen_molecule: list = [],
        sphere_radius: float = None,
        sphere_center: tuple = (0, 0, 0),
    ):
        self._cluster_dict = dict()
        self._multiplicity = 1
        self._charge = 0

        # fixing molecule to NOT move or rotate
        self._frozen_molecule = frozen_molecule

        # boundary conditions
        self._sphere_radius = sphere_radius
        self._sphere_center = sphere_center


# -------------------------------------------------------------
# -------------------------------------------------------------
# -------------------------------------------------------------
def test_Molecule():

    # mol = Molecule([("a", 0, 0, 0), ("b", 10, 10, 10)], -10, 5))
    # ---
    mol = Molecule(
        [("Xe", 0, 0, 0), ("Na", 5, 5, 5), ("Na", 10, 10, 10)], -10, 5
    )

    print("\n")
    print(mol)
    print("-" * 40)
    print("atoms: :\n", mol.atoms)
    print("charge:", mol.charge)
    print("multiplicity:", mol.multiplicity)

    print("-" * 40)
    print("\nmasses: ", mol.atomic_masses)
    print("\ntotal mass: ", mol.total_mass)
    print("\nsymbols: ", mol.symbols)
    print("\ntotal atoms: ", mol.total_atoms)
    print("\nelements: ", mol.elements)

    print("-" * 40)
    print("\nnumbering atoms (index):\n", mol.numbering_atoms)
    print(mol.xyz)

    print("-" * 40)
    print("\ncenter of mass: ", mol.center_of_mass)
    print("\npp axes: ", mol.principal_axes)
