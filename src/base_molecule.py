from copy import deepcopy

import numpy as np
from scipy.spatial.transform import Rotation

from data.atomic_data import atomic_mass


class Molecule:
    """
    Create a Molecule that is at least ONE atoms.
    The format of the INPUT coordinates is as follows (any):

    1. Dictionary type: {"atoms": [(<element> <X> <Y> <Z>), ...]}
    2. List type: [(<element> <X> <Y> <Z>), ...]
    3. Molecule type (Molecule Object)

    Parameters
    ----------
    atoms : list[tuple(str, float, float, float)]

    charge : int

    multiplicity : int
        integer larger than zero
    """

    def __init__(self, args=None):
        if not args:
            # self._atoms = [("X", 0, 0, 0)]
            self._atoms = list()
            self._charge: int = 0
            self._multiplicity: int = 1
        elif isinstance(args, Molecule):
            self._atoms = args.atoms
            self._charge: int = args.charge
            self._multiplicity: int = args.multiplicity
        elif self._atoms_format_check(args) == "is_dict":
            self._atoms = args.get("atoms")
            self._charge: int = args.get("charge", 0)
            self._multiplicity: int = args.get("multiplicity", 1)
        elif self._atoms_format_check(args) == "is_list":
            self._atoms = args
            self._charge: int = 0
            self._multiplicity: int = 1

        self._total_atoms: int = len(self._atoms)

    def __add__(self, other):
        return self.add_atoms(other)

    # ! check how to compare two Molecule objects
    # def __eq__(self, other: object) -> bool:
    #     return self.atoms == other.atoms

    def __mul__(self, value: int):
        return value * self

    def __rmul__(self, value: int):
        """to replicate a molecule

        Examples
        --------
        >>> hf_coordinates = [("H", 0, 0, 0), ("F", 1.917, 0, 0)]
        >>> hf_2 = Molecule(2 * hf_coordinates)
        >>> print(hf_3.atoms)
            4
        -- charge= 0 and multiplicity= 1 --
        H        0.00000000      0.00000000      0.00000000
        F        1.91700000      0.00000000      0.00000000
        H        0.00000000      0.00000000      0.00000000
        F        1.91700000      0.00000000      0.00000000

        Parameters
        ----------
        value : int
            quantity to replicate Molecue
        """
        if value < 1 or not isinstance(value, int):
            raise ValueError(
                "\nMultiplier must be and integer larger than zero"
                f"\ncheck --> {value}"
            )

        # ! how are multiplicity and charge computed for several atoms?
        new_molecule = deepcopy(self)
        for _ in range(value - 1):
            new_molecule += self.atoms

        return new_molecule

    def __str__(self) -> str:
        """Printing Molecule coordinates using XYZ format"""
        _comments = (
            f"-- charge= {self.charge:<-g} and "
            f"multiplicity= {self.multiplicity:<g} --"
        )
        write_xyz = f"""\t{self._total_atoms}\n{_comments:<s}\n"""
        for atom in self.atoms:
            write_xyz += f"""{atom[0]:<6}"""
            write_xyz += f"""\t{atom[1]:> .8f}"""
            write_xyz += f"""\t{atom[2]:> .8f}"""
            write_xyz += f"""\t{atom[3]:> .8f}\n"""

        return write_xyz

    def __sub__(self, other):
        return self.remove_atom(other)

    def _atoms_format_check(self, value) -> str:
        "format [(<element> <X> <Y> <Z>), ...]}"

        if isinstance(value, dict):
            try:
                atoms = value["atoms"]
            except KeyError:
                raise TypeError(
                    "\n\nMust be {'atoms': [(<element> <X> <Y> <Z>), ...]}"
                    "\nThe key 'atoms' is casesensitive"
                    f"\n--> you got {value}"
                )
            else:
                input_type = "is_dict"

        elif isinstance(value, list):
            input_type = "is_list"
            atoms = value
        else:
            raise TypeError(
                "\n Accepted either List or Dict type"
                "\n{'atoms': [(<element> <X> <Y> <Z>), ...]}"
                f"\n {type(value)}"
                f"\n {value}"
            )

        for line, atom in enumerate(atoms):
            try:
                assert len(str(atom[0]).replace(" ", ""))
                assert len([float(c) for c in atom[1:]]) == 3
            except (ValueError, AssertionError, KeyError):
                raise ValueError(
                    "\nMust be valid NOT empty element and "
                    "floats for xyz coordinates: (str, float, float, float)"
                    f"\ncheck atom number {line + 1} --> {atom}"
                )

        return input_type

    @property
    def atomic_masses(self) -> list:
        return [atomic_mass(s) for s in self.symbols]

    @property
    def atoms(self) -> list:
        return self._atoms

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
    def charge(self) -> int:
        return self._charge

    @charge.setter
    def charge(self, new_charge: int) -> None:
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
    def multiplicity(self, new_multiplicity: int) -> None:
        self._multiplicity = new_multiplicity

    @property
    def number_atoms(self) -> str:
        """show atom number line by line

        Returns
        -------
        str
            atom number line by line
        """
        numbered_atoms = list()
        for i in range(self.total_atoms):
            line = list(self.atoms[i])
            line[0] = f"\r{i:>4}: {line[0]:<6}"
            line[1] = f"\t{line[1]:> .8f}"
            line[2] = f"\t{line[2]:> .8f}"
            line[3] = f"\t{line[3]:> .8f}"
            numbered_atoms.append("".join(line))

        return "\n".join(numbered_atoms)

    @property
    def principal_axes(self) -> list:
        return [
            tuple(c)
            for c in (
                np.asarray(self.coordinates) - np.asarray(self.center_of_mass)
            )
        ]

    @property
    def symbols(self) -> list:
        return [c[0] for c in self.atoms]

    @property
    def total_atoms(self) -> int:
        return self._total_atoms

    @property
    def total_mass(self) -> float:
        """Sum atomic masses"""
        return sum(self.atomic_masses)

    @property
    def xyz(self) -> str:
        return self.__str__()

    def add_atoms(self, other):
        """adding extra atoms can NOT be MOVED or ROTATED

        Parameters
        ----------
        other : Molecue, dict, list
            with the coordinates, charge and multiplicity

        Returns
        -------
        Molecule
            a new Molecule

        Raises
        ------
        TypeError
            for anything else
        """

        # adding as a Cluster to be moved/rotated
        if isinstance(other, Cluster):
            return other.add_molecule(self)

        # adding as a Molecule add a FIXED atoms (can't be moved/rotated)
        if isinstance(other, Molecule):
            new_charge: int = self.charge + other.charge
            new_multiplicity: int = self.multiplicity
            all_atoms: list = self.atoms + other.atoms
        elif isinstance(other, dict):
            new_charge: int = self.charge + other.get("charge", 0)
            new_multiplicity: int = self.multiplicity
            all_atoms: list = self.atoms + other.get("atoms")
        elif isinstance(other, list):
            new_charge: int = 0
            new_multiplicity: int = 1
            all_atoms: list = self.atoms + other
        else:
            raise TypeError(
                "\nOnly type 'Molecule', list or dict could be added"
                f"\nyou have {type(other)}, check: "
                f"\n{other}"
            )

        return self.__class__(
            {
                "atoms": all_atoms,
                "charge": new_charge,
                "multiplicity": new_multiplicity,
            }
        )

    def get_atom(self, atom: int):
        if atom > self.total_atoms:
            raise IndexError(
                f"\nMolecule with {self.total_atoms} total atoms "
                f"and index [0-{self.total_atoms - 1}]"
                f"\n atom index must be less than {self.total_atoms}"
                f"\nCheck! You want to get atom with index {atom}"
            )

        new_atom: tuple = deepcopy(self).atoms.pop(atom)
        return self.__class__([tuple(new_atom)])

    def get_element(self, element: str):
        "Filtering all the elements by symbol"
        element_index = [
            idx for idx, e in enumerate(self.symbols) if element.title() == e
        ]

        if not element_index:
            raise ValueError(
                f"\nNot element '{element}' found in\n {self.symbols}"
            )

        new_molecule = Molecule()

        for i in range(self.total_atoms):
            if i in element_index:
                new_molecule += Molecule([self.atoms[i]])

        return self.__class__(
            {
                "atoms": new_molecule.atoms,
                "charge": self.charge,
                "multiplicity": self.multiplicity,
            }
        )

    def remove_atom(self, atom: int):
        if atom > self.total_atoms:
            raise IndexError(
                f"\nMolecule with {self.total_atoms} total atoms "
                f"and index [0-{self.total_atoms - 1}]"
                f"\n atom index must be less than {self.total_atoms}"
                f"\nCheck! You want to remove atom with index {atom}"
            )
        new_molecule = deepcopy(self)
        del new_molecule.atoms[atom]
        return self.__class__(new_molecule)

    def remove_element(self, element: str):
        "Removing all the elements by symbol"
        element_index = [
            idx for idx, e in enumerate(self.symbols) if element.title() == e
        ]

        if not element_index:
            raise ValueError(
                f"\nNot element '{element}' found in\n {self.symbols}"
            )

        new_molecule = Molecule()

        for i in range(self.total_atoms):
            if i in element_index:
                continue

            new_molecule += Molecule([self.atoms[i]])

        return self.__class__(
            {
                "atoms": new_molecule.atoms,
                "charge": self.charge,
                "multiplicity": self.multiplicity,
            }
        )


# ---------------------------------------------------------------------
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
        coordinates of each molecule/atom comma separates (support +/*)
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
        sphere_radius: float = None,
        sphere_center: tuple = (0, 0, 0),
    ):
        self._cluster_dict = dict()
        self._multiplicity = 1
        self._charge = 0

        # boundary conditions
        self._sphere_radius = sphere_radius
        self._sphere_center = sphere_center

        for i, mol in enumerate(args):
            size: int = len(self._cluster_dict)
            if isinstance(mol, Cluster):
                for j in mol.cluster_dictionary:
                    self._cluster_dict[size + j] = mol.cluster_dictionary[j]
                continue
            # else:
            try:
                molecule = Molecule(mol)
            except (TypeError, ValueError):
                raise TypeError(
                    "\nOnly type 'Molecule', list or dict to initialize"
                    "\t- Dict-> {'atoms': [(<element> <X> <Y> <Z>), ...]}"
                    "\t- List-> [(<element> <X> <Y> <Z>), ...]"
                    f"\nyou have {type(mol)}, check: {mol}"
                )
            else:
                self._charge += molecule.charge
                self._cluster_dict[size] = {
                    "atoms": molecule.atoms,
                    "charge": molecule.charge,
                    "multiplicity": molecule.multiplicity,
                }

        # initializing parent class
        self._atoms = [mol["atoms"] for mol in self._cluster_dict.values()]

        # Flatten List of Lists Using sum
        self._atoms = sum(self._atoms, [])

        super().__init__(
            {
                "atoms": self._atoms,
                "charge": self._charge,
                "multiplicity": self._multiplicity,
            }
        )

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
        """
        if value < 1 or not isinstance(value, int):
            raise ValueError(
                "\nMultiplier must be and integer larger than zero"
                f"\ncheck --> {value}"
            )

        new_cluster = deepcopy(self)
        for i in range(value - 1):
            new_cluster = new_cluster.add_molecule(deepcopy(self))

        return new_cluster

    def __str__(self) -> str:
        cluster_dict: dict = self.cluster_dictionary

        cluster_string: str = (
            f"Cluster of {self.total_molecules} molecules"
            f"and {self.total_atoms} total atoms\n"
        )
        for key, value in sorted(cluster_dict.items()):
            atoms = value["atoms"]
            cluster_string += f"molecule {key} (with {len(atoms)} atoms):\n"
            cluster_string += f"     -> atoms: {atoms}\n"
            charge = value["charge"]
            cluster_string += f"     -> charge: {charge}\n"
            multiplicity = value["multiplicity"]
            cluster_string += f"     -> multiplicity: {multiplicity}\n"

        return cluster_string

    @property
    def cluster_dictionary(self) -> dict:
        return self._cluster_dict

    @property
    def molecules(self):
        return self.cluster_dictionary

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

    @property
    def xyz(self) -> str:
        return super().__str__()

    def add_molecule(self, other):
        """adding extra molecule can be MOVED or ROTATED

        Parameters
        ----------
        other : Molecue, dict, list
            with the coordinates, charge and multiplicity

        Returns
        -------
        Cluster
            a new Molecular Cluster

        Raises
        ------
        TypeError
            for anything else
        """
        new_molecule_dict = {}

        if isinstance(other, Cluster):
            new_molecule_dict = other.cluster_dictionary
        elif isinstance(other, Molecule):
            new_molecule_dict[0] = {
                "atoms": other.atoms,
                "charge": other.charge,
                "multiplicity": other.multiplicity,
            }
        elif isinstance(other, dict):
            new_molecule_dict[0] = {
                "atoms": other.get("atoms"),
                "charge": other.get("charge", 0),
                "multiplicity": other.get("multiplicity", 1),
            }
        elif isinstance(other, list):
            new_molecule_dict[0] = {
                "atoms": other,
                "charge": 0,
                "multiplicity": 1,
            }
        else:
            raise TypeError(
                "\nOnly type 'Molecule', list or dict could be added"
                f"\nyou have {type(other)}, check: "
                f"\n{other}"
            )

        new_cluster_dict = dict(self.cluster_dictionary)

        for i in range(len(new_molecule_dict)):
            new_cluster_dict[self.total_atoms + i] = new_molecule_dict[i]

        return self.__class__(
            *new_cluster_dict.values(),
            sphere_radius=self.sphere_radius,
            sphere_center=self.sphere_center,
        )

    def get_molecule(self, molecule: int):
        if molecule > self.total_molecules:
            raise IndexError(
                f"\nMolecule with {self.total_molecules} total molecules "
                f"and index [0-{self.total_molecules - 1}]"
                f"\nmolecule index must be less than {self.total_molecules}"
                f"\nCheck! You want to get molecule with index {molecule}"
            )

        return self.__class__(
            deepcopy(self).molecules.pop(molecule),
            sphere_radius=self.sphere_radius,
            sphere_center=self.sphere_center,
        )

    def remove_molecule(self, molecule: int):
        if molecule > self.total_molecules:
            raise IndexError(
                f"\nMolecule with {self.total_molecules} total molecules "
                f"and index [0-{self.total_molecules - 1}]"
                f"\nmolecule index must be less than {self.total_molecules}"
                f"\nCheck! You want to remove molecule with index {molecule}"
            )
        new_cluster = deepcopy(self).molecules
        del new_cluster[molecule]

        return self.__class__(
            *new_cluster.values(),
            sphere_radius=self.sphere_radius,
            sphere_center=self.sphere_center,
        )

    def translate(
        self, molecule: int, x: float = 0, y: float = 0, z: float = 0
    ):
        """Returns a NEW Molecule Object with a TRANSLATED fragment"""

        if molecule > self.total_molecules:
            raise IndexError(
                f"\nMolecule with {self.total_molecules} total molecules "
                f"and index [0-{self.total_molecules - 1}]"
                f"\nmolecule index must be less than {self.total_molecules}"
                f"\nCheck! You want to remove molecule with index {molecule}"
            )

        molecule_to_move: Molecule = self.get_molecule(molecule)
        molecule_symbols: list = molecule_to_move.symbols
        molecule_coordinates: list = molecule_to_move.coordinates

        translated_coordinates = np.asarray(molecule_coordinates) + np.asarray(
            [x, y, z]
        )

        # checking if the new coordinates are into the boundary conditions
        # if it is out of our sphere, we move a on the opposite direction
        # 10% until this new move is into our sphere
        distance: float = np.linalg.norm(
            translated_coordinates - np.asarray(self.sphere_center)
        )
        if self.sphere_radius and (distance > self.sphere_radius):
            translated_coordinates = np.asarray(
                molecule_coordinates
            ) + -0.1 * np.asarray([x, y, z])

        # if self.sphere_radius and (distance > self.sphere_radius):
        #     translated_coordinates = np.asarray(
        #         molecule_coordinates
        #     ) + -0.9 * np.asarray([x, y, z])
        #     # print(self.sphere_radius, distance)
        #     while distance > self.sphere_radius:
        #         translated_coordinates = np.asarray(
        #             molecule_coordinates
        #         ) + 0.9 * np.asarray([x, y, z])

        #         distance = np.linalg.norm(
        #             translated_coordinates - np.asarray(self.sphere_center)
        #         )
        #         # print(self.sphere_radius, distance)

        translated_molecule = list()
        for i, atom in enumerate(molecule_symbols):
            translated_molecule.append(
                tuple([atom] + translated_coordinates[i].tolist())
            )

        new_molecule: dict = {
            "atoms": translated_molecule,
            "charge": molecule_to_move.charge,
            "multiplicity": molecule_to_move.multiplicity,
        }

        new_cluster = deepcopy(self).molecules
        new_cluster[molecule] = new_molecule

        return self.__class__(
            *new_cluster.values(),
            sphere_radius=self.sphere_radius,
            sphere_center=self.sphere_center,
        )

    def rotate(self, molecule: int, x: float = 0, y: float = 0, z: float = 0):
        """
        Returns a NEW Cluster Object with a ROTATED molecule (CLOCKWISE)
        around molecule internal center of mass
        """

        molecule_to_rotate: Molecule = self.get_molecule(molecule)
        molecule_symbols: list = molecule_to_rotate.symbols

        # avoid any rotatation attemp for a single atom system
        if not (len(molecule_symbols) > 1):
            return deepcopy(self)

        molecule_center_of_mass = molecule_to_rotate.center_of_mass
        molecule_principal_axes = molecule_to_rotate.principal_axes

        rotation_matrix = Rotation.from_euler(
            "xyz",
            [x, y, z],
            degrees=True,
        ).as_matrix()

        rotated_coordinates = (
            np.dot(molecule_principal_axes, rotation_matrix)
            + molecule_center_of_mass
        )

        rotated_molecule = list()
        for i, atom in enumerate(molecule_symbols):
            rotated_molecule.append(
                tuple([atom] + rotated_coordinates[i].tolist())
            )

        new_molecule: dict = {
            "atoms": rotated_molecule,
            "charge": molecule_to_rotate.charge,
            "multiplicity": molecule_to_rotate.multiplicity,
        }

        new_cluster = deepcopy(self).molecules
        new_cluster[molecule] = new_molecule

        return self.__class__(
            *new_cluster.values(),
            sphere_radius=self.sphere_radius,
            sphere_center=self.sphere_center,
        )
