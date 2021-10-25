from copy import deepcopy
from sys import modules
from typing import overload

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
    _charge: int = attr.ib(default=0)
    _multiplicity: int = attr.ib(default=1)

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

    @_charge.validator
    def _check_valid_charge(self, attribute, charge):
        if not isinstance(charge, int):
            raise ValueError(
                "\n\ncharge must be an integer "
                f"\nyou get --> 'charge = {charge}'\n"
            )

    @_multiplicity.validator
    def _check_valid_multiplicity(self, attribute, multiplicity):
        if not isinstance(multiplicity, int) or multiplicity < 1:
            raise ValueError(
                "\n\nmultiplicity must be an integer larger than zero (0)"
                f"\nyou get --> 'multiplicity = {multiplicity}'\n"
            )

    # ===============================================================
    # CONSTRUCTORS
    # ===============================================================
    @classmethod
    def from_dict(cls, atoms_dict):
        "Dictionary type: {'atoms': [(<element> <X> <Y> <Z>), ...]}"
        if "atoms" not in atoms_dict:
            # FIXME: KeyError does not support \n
            raise TypeError(
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
    def __add__(self, other) -> object:
        return self.add_molecule(other)

    def __mul__(self, value: int):
        return value * self

    def __rmul__(self, value: int):
        """
        Replicate a molecule.
        summing or multiplying Molecule classes produce a Cluster class

        Parameters
        ----------
        value : int
            quantity to replicate Molecue

        Return
        ------
        Cluster : object
        """
        if value < 1 or not isinstance(value, int):
            raise ValueError(
                "\nMultiplier must be and integer larger than zero"
                f"\ncheck --> '{value}'"
            )

        new_cluster = deepcopy(self)
        for _ in range(value - 1):
            new_cluster = new_cluster.add_molecule(deepcopy(self))

        return new_cluster

    def __str__(self):
        return str(attr.astuple(self))

    # ===============================================================
    # PROPERTIES
    # ===============================================================
    @property
    def atoms(self) -> list:
        return self._atoms

    @atoms.setter
    def atoms(self, *args, **kwargs) -> None:
        raise AttributeError(
            "\n\nyou cannot reset 'atoms'. Consider create a new instance \n"
        )

    @property
    def atomic_masses(self) -> list:
        return [Atom(*atom).atomic_mass for atom in self.atoms]

    @property
    def charge(self) -> int:
        return self._charge

    @charge.setter
    def charge(self, new_charge) -> int:
        if not isinstance(new_charge, int):
            raise ValueError(
                "\n\ncharge must be an integer "
                f"\nyou get --> 'charge = {new_charge}'\n"
            )
        self._charge = new_charge

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
    def multiplicity(self) -> int:
        return self._multiplicity

    @multiplicity.setter
    def multiplicity(self, new_multiplicity) -> int:
        if not isinstance(new_multiplicity, int) or new_multiplicity < 1:
            raise ValueError(
                "\n\nmultiplicity must be an integer larger than zero (0)"
                f"\nyou get --> 'multiplicity = {new_multiplicity}'\n"
            )
        self._multiplicity = new_multiplicity

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
            https://en.wikipedia.org/wiki/Jacobicoordinates
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
            f"-- charge={self.charge:<-g} and "
            f"multiplicity={self.multiplicity:<g} --"
        )
        write_xyz = f"""\t{self.total_atoms}\n{comments:<s}\n"""
        for atom in self.atoms:
            write_xyz += f"""{atom[0]:<6}"""
            write_xyz += f"""\t{atom[1]:> 15.8f}"""
            write_xyz += f"""\t{atom[2]:> 15.8f}"""
            write_xyz += f"""\t{atom[3]:> 15.8f}\n"""

        return write_xyz

    def add_atoms(self, new_atoms: list) -> object:
        """adding extra atoms can NOT be MOVED or ROTATED

        Parameters
        ----------
        other : list
            cartesian coordinates; like [(<element>, <X>, <Y>, <Y>), ...]

        Returns
        -------
        Molecule : object
            a new Molecule

        Raises
        ------
        TypeError
            for anything else
        """
        if not isinstance(new_atoms, list):
            raise TypeError(
                f"\n\ncoordinates format must be a list of tuple: "
                "[(str, float, float, float), ...]"
                f"check --> \n{new_atoms}\n"
            )

        total_atoms: list = self.atoms + new_atoms
        return self.__class__(total_atoms)

    def add_molecule(self, other) -> object:
        if not isinstance(other, Molecule):
            raise TypeError(
                "\nOnly type 'Molecule', list or dict could be added"
                f"\nyou have a type: '{type(other)}', check: \n{other}"
            )
        return Cluster(self, other)

    def get_atom(self, atom: int) -> list:
        """
        Getting catesian coordinate for an atom

        Parameters
        ----------
        atom : int
            atom index

        Returns
        -------
        list
            ["element", "X", "Y", "Z"]

        Raises
        ------
        IndexError
        """
        if atom >= self.total_atoms:
            raise IndexError(
                f"\nMolecule with {self.total_atoms} total atoms "
                f"and index [0-{self.total_atoms - 1}]"
                f"\n atom index must be less than {self.total_atoms}"
                f"\nCheck! You want to get atom with index {atom}"
            )
        return self.atoms[atom]

    def remove_atom(self, atom: int) -> object:
        if atom > self.total_atoms:
            raise IndexError(
                f"\nMolecule with {self.total_atoms} total atoms "
                f"and index [0-{self.total_atoms - 1}]"
                f"\n atom index must be less than {self.total_atoms}"
                f"\nCheck! You want to remove atom with index {atom}"
            )
        new_atoms: list = self.atoms

        return self.__class__(new_atoms.pop(atom))


# -------------------------------------------------------------
class Cluster(Molecule):
    """
    Create a Cluster with molecules/atoms to move and rotate
    using spherical boundary conditions (SBC).
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
        for a wrong input argument
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

        cluster_atoms: list = list()

        # for count, mol in enumerate(args):
        for mol in args:
            size: int = len(self._cluster_dict)
            if isinstance(mol, Cluster):
                for j in mol._cluster_dict:
                    self._cluster_dict[size + j] = mol._cluster_dict[j]
                self._charge += mol.charge
                cluster_atoms += mol.atoms
                # restarting the loop
                continue
            elif isinstance(mol, Molecule):
                new_molecule = deepcopy(mol)
            elif isinstance(mol, dict):
                new_molecule = Molecule.from_dict(mol)
            elif isinstance(mol, list):
                new_molecule = Molecule(mol)
            else:
                raise TypeError(
                    "\nOnly type 'Molecule', list or dict to initialize"
                    "\n\t- Dict: {'atoms': [(<element> <X> <Y> <Z>), ...]}"
                    "\n\t- List: [(<element> <X> <Y> <Z>), ...]"
                    f"\nyou have a NOT valid '{type(mol)}', check: \n{mol}"
                )

            cluster_atoms += new_molecule.atoms
            # ! how is computed the cluster total multiplicity?
            self._charge += new_molecule.charge
            self._cluster_dict[size] = new_molecule

        # initialazing Cluster as a 'Molecule' (sum of all individual ones)
        super().__init__(
            atoms=cluster_atoms,
            charge=self.charge,
            multiplicity=self.multiplicity,
        )

    # ===============================================================
    # MAGIC METHODS
    # ===============================================================
    def __add__(self, other):
        return self.add_molecule(other)

    def __mul__(self, value: int):
        return value * self

    def __rmul__(self, value: int):
        """to replicate a molecule

        Parameters
        ----------
        value : int
            quantity to replicate Molecue

        Return
        ------
        Cluster : object
            summing or multiplying Molecule classes produce a Cluster class
        """
        if value < 1 or not isinstance(value, int):
            raise ValueError(
                "\nMultiplier must be and integer larger than zero"
                f"\ncheck --> '{value}'"
            )

        new_cluster = deepcopy(self)
        for _ in range(value - 1):
            new_cluster = new_cluster.add_molecule(deepcopy(self))

        return new_cluster

    def __str__(self):
        cluster_dict: dict = self._cluster_dict

        cluster_string: str = (
            f"Cluster of ({self.total_molecules}) molecules"
            f" and ({self.total_atoms}) total atoms\n"
        )
        for key, molecule in cluster_dict.items():
            atoms = molecule.atoms
            cluster_string += f" #{key}: molecule with {len(atoms)} atoms:\n"
            cluster_string += f"     --> atoms: {atoms}\n"
            charge = molecule.charge
            cluster_string += f"     --> charge: {charge:>+}\n"
            multiplicity = molecule.multiplicity
            cluster_string += f"     --> multiplicity: {multiplicity}\n"

        return cluster_string

    # ===============================================================
    # PROPERTIES
    # ===============================================================
    @property
    def cluster_dictionary(self) -> dict:
        return self._cluster_dict

    @property
    def frozen_molecule(self) -> int:
        return self._frozen_molecule

    @frozen_molecule.setter
    def frozen_molecule(self, values) -> None:
        if isinstance(values, list):
            self._frozen_molecule = values
        else:
            self._frozen_molecule = [values]

    @property
    def sphere_center(self) -> tuple:
        return self._sphere_center

    @sphere_center.setter
    def sphere_center(self, new_center: tuple) -> None:
        self._sphere_center = new_center

    @property
    def sphere_radius(self) -> float:
        return self._sphere_radius

    @sphere_radius.setter
    def sphere_radius(self, new_radius: float) -> None:
        self._sphere_radius = new_radius

    @property
    def total_molecules(self) -> int:
        return len(self._cluster_dict)

    def add_molecule(self, other) -> object:
        new_cluster = deepcopy(self)
        return self.__class__(
            new_cluster,
            other,
            frozen_molecule=new_cluster.frozen_molecule,
            sphere_radius=new_cluster.sphere_radius,
            sphere_center=new_cluster.sphere_center,
        )

    def get_molecule(self, molecule: int):
        if molecule not in self.cluster_dictionary:
            raise IndexError(
                f"\nMolecule with {self.total_molecules} total molecules "
                f"and index [0-{self.total_molecules - 1}]"
                f"\nmolecule index must be less than {self.total_molecules}"
                f"\nCheck! You want to get molecule with index {molecule}"
            )

        cluster_dict: dict = deepcopy(self).cluster_dictionary
        new_molecule: Molecule = cluster_dict.pop(molecule)

        return self.__class__(
            new_molecule,
            frozen_molecule=self.frozen_molecule,
            sphere_radius=self.sphere_radius,
            sphere_center=self.sphere_center,
        )

    def remove_molecule(self, molecule: int) -> object:
        if molecule not in self.cluster_dictionary:
            raise IndexError(
                f"\nMolecule with {self.total_molecules} total atoms "
                f"and index [0-{self.total_molecules - 1}]"
                f"\n molecule index must be less than {self.total_molecules}"
                f"\nCheck! You want to remove molecule with index {molecule}"
            )
        new_cluster: Cluster = deepcopy(self)
        new_cluster_dict: dict = new_cluster.cluster_dictionary
        del new_cluster_dict[molecule]

        return self.__class__(
            *new_cluster._cluster_dict.values(),
            frozen_molecule=new_cluster.frozen_molecule,
            sphere_radius=new_cluster.sphere_radius,
            sphere_center=new_cluster.sphere_center,
        )

    def rotate(
        self, molecule: int = None, x: float = 0, y: float = 0, z: float = 0
    ):
        """
        Returns a NEW Cluster Object with a ROTATED molecule (CLOCKWISE)
        around molecule internal center of mass
        """
        # avoiding to rotate a FROZEN molecule
        if molecule in self.frozen_molecule:
            return deepcopy(self)

        if not isinstance(molecule, int) and (
            0 < molecule >= self.total_molecules
        ):
            raise IndexError(
                f"\nMolecule with {self.total_molecules} total molecules "
                f"and index [0-{self.total_molecules - 1}]"
                f"\nmolecule index must be less than {self.total_molecules}"
                f"\nCheck! You want to remove molecule with index {molecule}"
            )

        molecule_to_rotate: Molecule = self._cluster_dict[molecule]
        molecule_symbols: list = molecule_to_rotate.symbols

        # avoid any rotatation attemp for a single atom system
        if len(molecule_symbols) <= 1:
            return deepcopy(self)

        molecule_center_of_mass = molecule_to_rotate.center_of_mass
        molecule_principal_axes = molecule_to_rotate.principal_axes

        # rotate around sphere center
        x, y, z = np.asarray(self.sphere_center) + np.asarray([x, y, z])

        rotation_matrix = Rotation.from_euler(
            "xyz",
            [x, y, z],
            degrees=True,
        ).as_matrix()

        rotatedcoordinates = (
            np.dot(molecule_principal_axes, rotation_matrix)
            + molecule_center_of_mass
        )

        rotated_molecule = list()
        for i, atom in enumerate(molecule_symbols):
            rotated_molecule.append(
                tuple([atom] + rotatedcoordinates[i].tolist())
            )

        new_cluster = deepcopy(self)
        new_cluster._cluster_dict[molecule] = Molecule(rotated_molecule)

        return self.__class__(
            *new_cluster._cluster_dict.values(),
            frozen_molecule=new_cluster.frozen_molecule,
            sphere_radius=new_cluster.sphere_radius,
            sphere_center=new_cluster.sphere_center,
        )

    def translate(
        self, molecule: int = None, x: float = 0, y: float = 0, z: float = 0
    ):
        """Returns a NEW Molecule Object with a TRANSLATED fragment"""
        # avoiding to rotate a FROZEN molecule
        if molecule in self.frozen_molecule:
            return deepcopy(self)

        if not isinstance(molecule, int) and (
            0 < molecule >= self.total_molecules
        ):
            raise IndexError(
                f"\nMolecule with {self.total_molecules} total molecules "
                f"and index [0-{self.total_molecules - 1}]"
                f"\nmolecule index must be less than {self.total_molecules}"
                f"\nCheck! You want to remove molecule with index {molecule}"
            )

        molecule_to_move: Molecule = self._cluster_dict[molecule]
        molecule_symbols: list = molecule_to_move.symbols

        molecule_center_of_mass = molecule_to_move.center_of_mass
        molecule_principal_axes = molecule_to_move.principal_axes

        translatedcoordinates = np.asarray(
            molecule_center_of_mass
        ) + np.asarray([x, y, z])

        # checking if the new coordinates are into the boundary conditions
        # if it is out of our sphere, we rescale it to match the sphere radius
        distance: float = np.linalg.norm(
            translatedcoordinates - np.asarray(self.sphere_center)
        )
        if self.sphere_radius and (distance > self.sphere_radius):

            max_distance: float = self.sphere_radius / np.linalg.norm(
                translatedcoordinates - np.asarray(self.sphere_center)
            )

            # rescaling to match radius
            translatedcoordinates = max_distance * translatedcoordinates + (
                1 - max_distance
            ) * np.asarray(self.sphere_center)

        translatedcoordinates = molecule_principal_axes + translatedcoordinates

        translated_molecule = list()
        for i, atom in enumerate(molecule_symbols):
            translated_molecule.append(
                tuple([atom] + translatedcoordinates[i].tolist())
            )

        new_cluster = deepcopy(self)
        new_cluster._cluster_dict[molecule] = Molecule(translated_molecule)

        return self.__class__(
            *new_cluster._cluster_dict.values(),
            frozen_molecule=new_cluster.frozen_molecule,
            sphere_radius=new_cluster.sphere_radius,
            sphere_center=new_cluster.sphere_center,
        )

    def move_molecule(
        self,
        molecule: int = 0,
        max_step: float = None,
        max_rotation: float = None,
        seed: int = None,
    ) -> object:
        """Moving (translating and rotating) randomly without overlapping
        any atom

        Parameters
        ----------
        molecule : int, optional
            molecule to move randomly, by default molecule with index zero (0)
        max_step : float, optional
            maximun value for any translation, by default None
        max_rotation : float, optional
            maximun angle fro any rotation, by default None
        seed : int, optional
            seed to initialize the random generator function, by default None

        Returns
        -------
        object : Cluster
            returns a Cluster object whit a molecule moved to a random place
            without overlapping any other

        Raises
        ------
        AttributeError : OverlappingError
            After serching for max_overlap_cycle and no place found for the
            molecule without overlapping any other
        """

        # predefine maximun closeness between any pair of atoms
        max_closeness: int = 1.0

        if not max_step or not isinstance(max_step, (int, float)):
            max_step = 1.1 * max_closeness

        if not max_rotation or not isinstance(max_rotation, (int, float)):
            max_rotation = 30

        molecule_to_move: Cluster = self.get_molecule(molecule)

        cluster_without_molecule: Cluster = self.remove_molecule(molecule)
        cluster_coordinates: Cluster = cluster_without_molecule.coordinates

        random_gen = np.random.default_rng(seed)
        max_overlap_cycle: int = 1000

        for count in range(max_overlap_cycle):

            if count % 10 == 0:
                max_step *= 1.1
                max_rotation *= 1.1

            # angle between [0, max_rotation) degrees
            rotation_x = random_gen.uniform() * max_rotation
            rotation_y = random_gen.uniform() * max_rotation
            rotation_z = random_gen.uniform() * max_rotation

            # moving between [-max_step, +max_step] Angstrom
            tranlation_x = max_step * (random_gen.uniform() - 0.5)
            tranlation_y = max_step * (random_gen.uniform() - 0.5)
            tranlation_z = max_step * (random_gen.uniform() - 0.5)

            new_molecule: Cluster = molecule_to_move.translate(
                0,
                tranlation_x,
                tranlation_y,
                tranlation_z,
            ).rotate(
                0,
                rotation_x,
                rotation_y,
                rotation_z,
            )

            molecule_coordinates: list = new_molecule.coordinates

            overlap: bool = Cluster.overlapping(
                molecule_coordinates,
                cluster_coordinates,
                max_closeness=max_closeness,
            )

            if not overlap:
                break

        cluster_dict: dict = self.cluster_dictionary
        cluster_dict[molecule] = new_molecule

        new_cluster = Cluster(
            *cluster_dict.values(),
            frozen_molecule=self.frozen_molecule,
            sphere_radius=self.sphere_radius,
            sphere_center=self.sphere_center,
        )

        if count < max_overlap_cycle - 1:
            return new_cluster
        else:
            raise AttributeError(
                "\n\n *** Overlapping Error ***"
                "\nat least one atom is overlapped with a distance"
                f" less than '{max_closeness}' Angstroms"
                "\nfor a cluster into a sphere of radius"
                f" '{self.sphere_radius}' Angstroms"
                f"\nPlease, check: \n\n{new_cluster.xyz}"
            )

    @staticmethod
    def overlapping(
        first_coordinates: list,
        second_coordinates: list,
        max_closeness: float = 0.5,
    ) -> bool:
        """pair-wise checking if any overlapping among points
        with a radius defined by `max_closeness`

        Parameters
        ----------
        first_coordinates : list
            list of tuples [(float, float, float), ...]
        second_coordinates : list
            list of tuples [(float, float, float), ...]
        max_closeness : float, optional
            maximun closeness between two pairs, by default 0.5

        Returns
        -------
        bool
            True if two point are closer than `max_closeness`
        """

        for first_atom in first_coordinates:
            for second_atom in second_coordinates:
                distance = np.linalg.norm(
                    np.asarray(first_atom) - np.asarray(second_atom)
                )

                if distance < max_closeness:
                    return True

        return False


# -------------------------------------------------------------
# -------------------------------------------------------------
# -------------------------------------------------------------
