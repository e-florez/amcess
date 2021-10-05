from copy import deepcopy

import attr
import numpy as np
from scipy.spatial.transform import Rotation

from data.atomic_data import atomic_mass


@attr.s
class Atom:
    """
    Representation of an individual atom. Format (<element> <X> <Y> <Z>)

    Notes
    -----
        Not to be used by itself
    """

    element = attr.ib()
    x = attr.ib(repr=True)
    y = attr.ib(repr=True)
    z = attr.ib(repr=True)

    @property
    def atom(self):
        return (self.element, self.x, self.y, self.z)


class Molecule:
    """
    Create a Molecule that is at least ONE atoms.
    The format of the INPUT coordinates is as follows (any):

    1. Dictionary type: {"atoms": [(<element> <X> <Y> <Z>), ...]}
    2. List type: [(<element> <X> <Y> <Z>), ...]
    3. Molecule type (Molecule Object)

    Parameters
    ----------
    atoms: list[tuple(str, float, float, float)]

    charge: int

    multiplicity: int
        integer larger than zero
    """

    def __init__(self, args=None):
        if not args:
            # self._atoms = [("X", 0, 0, 0)]
            self._atoms = list()
            self._charge: int = 0
            self._multiplicity: int = 1
        elif type(args).__name__ == "Molecule":
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
            )

        for line, atom in enumerate(atoms):
            try:
                assert len(str(atom[0]).replace(" ", ""))
                assert len([float(c) for c in atom[1:]]) == 3
            except (ValueError, AssertionError):
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

    def add_atoms(self, other):
        """adding extra atoms

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

        if type(other).__name__ == "Molecule":
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
        new_molecule = deepcopy(self)
        new_atom = new_molecule.atoms.pop(atom)
        return self.__class__([tuple(new_atom)])

    def get_elements(self, element: str):
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
        new_molecule = deepcopy(self)
        del new_molecule.atoms[atom]
        return self.__class__(new_molecule)

    def remove_elements(self, element: str):
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
    # def __init__(self, *args, sphere_radius, sphere_center):
    def __init__(self, *args):
        self._molecules_dict = dict()
        self._cluster_multiplicity = 1
        self._cluster_charge = 0

        for i, mol in enumerate(args):
            try:
                molecule = Molecule(mol)
            except (TypeError, ValueError):
                raise TypeError(
                    "\nOnly type 'Molecule', list or dict could be added"
                    f"\nyou have {type(mol)}, check: {mol}"
                )
            else:
                self._cluster_charge += molecule.charge
                self._molecules_dict[i] = {
                    "atoms": molecule.atoms,
                    "charge": molecule.charge,
                    "multiplicity": molecule.multiplicity,
                }

        self._cluster_atoms = [
            mol["atoms"] for mol in self._molecules_dict.values()
        ]

        # Flatten List of Lists Using sum
        self._cluster_atoms = sum(self._cluster_atoms, [])

        super().__init__(
            {
                "atoms": self._cluster_atoms,
                "charge": self._cluster_charge,
                "multiplicity": self._cluster_multiplicity,
            }
        )

    @property
    def total_molecules(self) -> int:
        return len(self._molecules_dict)


# ---------------------------------------------------------------------
# ---------------------------------------------------------------------
# ---------------------------------------------------------------------
class oldMolecule:
    """
    Create an Atomic/Molecular cluster. The formatting of the INPUT coordinates
    is as follows (any):

    1. Dictionary type: {"atoms": [(<element> <X> <Y> <Z>), ...]}
    2. List type: [(<element> <X> <Y> <Z>), ...]
    3. Molecule type (Molecule Object)

    """

    def __init__(self, *args):
        self._cluster = dict()
        self._total_fragments = 0
        self._total_atoms = 0
        self._total_fragments = len(args)

        for i, fragment in enumerate(args):
            if type(fragment) == list:
                self._cluster[i] = fragment
            elif type(fragment) == dict:
                try:
                    self._cluster[i] = fragment["atoms"]
                except KeyError as err:
                    print(
                        "\n*** ERROR *** \n"
                        + " key must be 'atoms' (casesensitive!) \n"
                        + " but you have: \n\n"
                        + "\n".join(list(fragment.keys()))
                        + "\n\n"
                    )
                    raise err
            elif type(fragment) == Molecule:
                self._cluster[i] = fragment.coordinates
            else:
                raise (
                    "\n *** ERROR ***"
                    + " only list or dict-type object must be used to "
                    + " to create a new Object. Check class help"
                )

            self._total_atoms += len(self._cluster[i])

        # checking format
        self._check_coordinates

    def __str__(self) -> str:
        """Printing Molecule coordinates using XYZ format"""
        _comments = (
            f"--system of {self.total_fragments} molecules "
            f"and {self._total_atoms} total individual atoms--"
        )
        _write_coordinates_xyz = f"""\t{self._total_atoms}\n{_comments:<s}\n"""
        for _atoms in self.coordinates:
            _write_coordinates_xyz += f"""{_atoms[0]:<s}"""
            _write_coordinates_xyz += f"""\t{_atoms[1]:> .8f}"""
            _write_coordinates_xyz += f"""\t{_atoms[2]:> .8f}"""
            _write_coordinates_xyz += f"""\t{_atoms[3]:> .8f}\n"""

        # ANSI escape codes: move cursor up one line
        # _write_coordinates_xyz += f"""\033[A"""

        return _write_coordinates_xyz

    def _check_coordinates(self) -> bool:
        """Checking XYZ format (str, float, float, float)"""

        for _line, _atoms in enumerate(self.coordinates):

            _error_message = (
                "\n *** ERROR ***\n"
                + f" line: {_line + 2} does not match format "
                + "(str, float, float, float)\n\n"
                + f" check line: {_line + 2} --> \t{_atoms}\n"
                + "...\n"
                # + str(self)
                + (
                    "\n".join(
                        [
                            ("line: " + str(i) + " | " + str(c))
                            for i, c in enumerate(str(self).split("\n"))
                        ]
                    )
                )
                + "\n...\n"
            )

            try:
                assert len(str(_atoms[0]).replace(" ", ""))
                assert len([float(c) for c in _atoms[1:]]) == 3
            except (ValueError, AssertionError) as err:
                print(_error_message)
                raise err

        return True

    def _check_fragment(self, fragment: int):
        try:
            self._fragment = int(fragment)
            self._cluster[self._fragment]
        except (KeyError, ValueError):
            print(
                f"\n* Warning! selecting a fragment: {self._fragment}"
                + "(atom/molecule) NOT found, it must be an integer in "
                + f"{list(self._cluster.keys())}"
            )
            return False
        else:
            return True

    @property
    def atomic_masses(self):
        return [atomic_mass(s) for s in self.symbols]

    @property
    def cartesian_coordinates(self):
        # ! cartesian coordinates and coordinates are different function
        # ! may cause user confusion to choose and generate logic later
        return [c[1:] for c in self.coordinates]

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

        _total_mass = 1 if not self.total_mass else self.total_mass

        self._center_of_mass = (
            np.dot(
                np.asarray(self.atomic_masses),
                np.asarray(self.cartesian_coordinates),
            )
            / _total_mass
        )

        return self._center_of_mass

    @property
    def coordinates(self):
        self._coordinates = list()
        for fragment in self._cluster:
            for atom in self._cluster[fragment]:
                self._coordinates.append(atom)

        return self._coordinates

    @property
    def principal_axes(self) -> tuple:

        self._principal_axes = np.asarray(
            self.cartesian_coordinates
        ) - np.asarray(self.center_of_mass)

        return tuple(self._principal_axes)

    @property
    def show_all(self) -> bool:
        # TODO: print everthing, symbols, mass, etc.

        print("\n")
        print("total fragments: ", self.total_fragments)
        print("total atoms: ", self.total_atoms)
        print("symbols: ", self.symbols)
        print("atomic masses: ", self.atomic_masses)
        print("total mass: ", self.total_mass)

        print("\n")
        for i, fragment in enumerate(self._cluster):
            print("-" * 50)
            print(f" fragment number: {i}")
            print("-" * 50)
            individual_fragmen = self.get_fragments(fragment)
            print(individual_fragmen)
            print("-" * 50 + "\n ***")

        return True

    @property
    def symbols(self):
        return [c[0] for c in self.coordinates]

    @property
    def total_atoms(self):
        return self._total_atoms

    @property
    def total_fragments(self):
        return self._total_fragments

    @property
    def total_mass(self) -> float:
        """Sum atomic masses"""
        return sum(self.atomic_masses)

    @property
    def xyz(self) -> str:
        return self.__str__()

    def add_fragments(self, new_coordinates):
        """Returns a NEW Molecule Object"""
        _new_cluster = dict()
        for i, fragment in enumerate(self._cluster):
            _new_cluster[i] = {"atoms": self._cluster[fragment]}

        if type(new_coordinates) == list:
            _new_cluster[len(self._cluster)] = {"atoms": new_coordinates}
        elif type(new_coordinates) == dict:
            _new_cluster[len(self._cluster)] = new_coordinates
        elif type(new_coordinates) == Molecule:
            _new_cluster[len(self._cluster)] = {
                "atoms": new_coordinates.coordinates,
            }
        else:
            print(
                "\n* Warning! fragment to add MUST be list-type "
                + " dict-type with 'coordinates' as a key. Check class help"
            )
            return deepcopy(self)

        return Molecule(*_new_cluster.values())

    def delete_fragments(self, fragment: int):
        """Returns a NEW Molecule Object"""
        self._fragment = fragment
        if not self._check_fragment(self._fragment):
            return deepcopy(self)

        _new_coordinates = deepcopy(self._cluster)
        _new_coordinates.pop(self._fragment)

        _new_cluster = dict()
        for i, fragment in enumerate(_new_coordinates):
            _new_cluster[i] = {"atoms": _new_coordinates[fragment]}

        return Molecule(*_new_cluster.values())

    def get_fragments(self, fragment):
        """Returns a NEW Object with coordinates of the selected fragment"""
        self._fragment = fragment
        if self._check_fragment(self._fragment):
            return Molecule({"atoms": self._cluster[self._fragment]})
        else:
            print(
                f"\n* Warning! selecting a fragment: {self._fragment} "
                + "(atom/molecule) must be an integer in "
                + f"{list(self._cluster.keys())}"
            )
            return deepcopy(self)

    def translate(self, fragment, x=0, y=0, z=0):
        """Returns a NEW Molecule Object with a TRANSLATED fragment"""
        self._fragment = fragment
        if not self._check_fragment(self._fragment):
            return deepcopy(self)

        _fragment_to_move: Molecule = self.get_fragments(self._fragment)
        _fragment_symbols = _fragment_to_move.symbols
        _fragment_coordinates = _fragment_to_move.cartesian_coordinates

        _translated_coordinates = np.asarray(
            _fragment_coordinates
        ) + np.asarray([x, y, z])

        _translated_fragment = list()
        for i, _atom in enumerate(_fragment_symbols):
            _translated_fragment.append(
                tuple([_atom] + _translated_coordinates[i].tolist())
            )

        return self.delete_fragments(self._fragment).add_fragments(
            _translated_fragment
        )

    def rotate(self, fragment, x=0, y=0, z=0):
        """Returns a NEW Molecule Object with a ROTATED fragment (CLOCKWISE)
        (around internal center of mass)"""
        self._fragment = fragment
        if not self._check_fragment(self._fragment):
            return deepcopy(self)

        _fragment_to_rotate: Molecule = self.get_fragments(self._fragment)
        _fragment_symbols = _fragment_to_rotate.symbols

        # avoid any rotatation attemp for a single atom system
        if not (len(_fragment_symbols) > 1):
            return deepcopy(self)

        _rotation_matrix = Rotation.from_euler(
            "xyz",
            [x, y, z],
            degrees=True,
        ).as_matrix()

        _fragment_center_of_mass = _fragment_to_rotate.center_of_mass
        _fragment_principal_axes = _fragment_to_rotate.principal_axes

        _rotated_coordinates = (
            np.dot(_fragment_principal_axes, _rotation_matrix)
            + _fragment_center_of_mass
        )

        _rotated_fragment = list()
        for i, _atom in enumerate(_fragment_symbols):
            _rotated_fragment.append(
                tuple([_atom] + _rotated_coordinates[i].tolist())
            )

        return self.delete_fragments(self._fragment).add_fragments(
            _rotated_fragment
        )


class oldCluster(oldMolecule):
    pass
