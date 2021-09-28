from copy import deepcopy

import numpy as np
from scipy.spatial.transform import Rotation

from data.atomic_data import atomic_mass


class Cluster:
    """
    Create an Atomic/Molecular cluster. The formatting of the INPUT coordinates
    is as follows (any):

    1. Dictionary type: {"coordinates": [(<element> <X> <Y> <Z>), ...]}
    2. List type: [(<element> <X> <Y> <Z>), ...]
    3. Cluster type (Cluster Object)

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
                    self._cluster[i] = fragment["coordinates"]
                except KeyError as error:
                    print(
                        f"\n*** ERROR *** \n"
                        + f" key must be 'coordinates' (casesensitive!) \n"
                        + f" but you have: \n\n"
                        + "\n".join(list(fragment.keys()))
                        + f"\n\n"
                    )
                    raise error
            elif type(fragment) == Cluster:
                self._cluster[i] = fragment.coordinates
            else:
                raise (
                    f"\n *** ERROR ***"
                    + f" only list or dict-type object must be used to "
                    + f" to create a new Object. Check class help"
                )

            self._total_atoms += len(self._cluster[i])

        # checking format
        self._check_coordinates

    @property
    def total_fragments(self):
        return self._total_fragments

    @property
    def total_atoms(self):
        return self._total_atoms

    @property
    def coordinates(self):
        self._coordinates = list()
        for fragment in self._cluster:
            for atom in self._cluster[fragment]:
                self._coordinates.append(atom)

        return self._coordinates

    @property
    def show_all(self) -> bool:
        # TODO: print everthing, symbols, mass, etc.

        print(f"\n")
        print(f"total fragments: ", self.total_fragments)
        print(f"total atoms: ", self.total_atoms)
        print(f"symbols: ", self.symbols)
        print(f"atomic masses: ", self.atomic_masses)
        print(f"total mass: ", self.total_mass)

        print(f"\n")
        for i, fragment in enumerate(self._cluster):
            print("-" * 50)
            print(f" fragment number: {i}")
            print("-" * 50)
            individual_fragmen = self.get_fragments(fragment)
            print(individual_fragmen)
            print(f"-" * 50 + "\n ***")

        return True

    @property
    def _check_coordinates(self) -> bool:
        """Checking XYZ format (str, float, float, float)"""

        for _line, _atoms in enumerate(self.coordinates):

            _error_message = (
                f"\n *** ERROR ***\n"
                + f" line: {_line + 2} does not match format (str, float, float, float)\n\n"
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
            except (ValueError, AssertionError) as error:
                print(_error_message)
                raise error

        return True

    def __str__(self) -> str:
        """Printing Cluster coordinates"""
        _write_coordinates_xyz = list()
        _write_coordinates_xyz.append(f"\t" + str(self.total_atoms))
        _write_coordinates_xyz.append(
            f" -- system of {self.total_fragments} fragments (atoms/molecules) "
            + f"and {self._total_atoms} total individual atoms --"
        )
        for _atoms in self.coordinates:
            _atoms = list(_atoms)
            _format = ""
            for i, _ in enumerate(_atoms):
                try:
                    _atoms[i] = float(_atoms[i])
                    _format += "{:> .8f}\t"
                except ValueError:
                    _atoms[i] = str(_atoms[i])
                    _format += "{:<6}\t"

            _write_coordinates_xyz.append(_format.format(*_atoms))
        return "\n".join([str(c) for c in _write_coordinates_xyz])

    @property
    def cartesian_coordinates(self):
        # ! cartesian coordinates and coordinates are different function
        # ! may cause user confusion to choose and generate logic later
        return [c[1:] for c in self.coordinates]

    @property
    def symbols(self):
        return [c[0] for c in self.coordinates]

    @property
    def atomic_masses(self):
        return [atomic_mass(s) for s in self.symbols]

    @property
    def total_mass(self) -> float:
        """Sum atomic masses"""
        return sum(self.atomic_masses)

    def _check_fragment(self, fragment: int):
        try:
            self._fragment = int(fragment)
            self._cluster[self._fragment]
        except (KeyError, ValueError):
            print(
                f"\n* Warning! selecting a fragment: {self._fragment} (atom/molecule)"
                f" NOT found, it must be an integer in {list(self._cluster.keys())}"
            )
            return False
        else:
            return True

    def get_fragments(self, fragment):
        """Returns a NEW Object with coordinates of the selected fragment"""
        self._fragment = fragment
        if self._check_fragment(self._fragment):
            return Cluster({"coordinates": self._cluster[self._fragment]})
        else:
            print(
                f"\n* Warning! selecting a fragment: {self._fragment} (atom/molecule)"
                f" must be an integer in {list(self._cluster.keys())}"
            )
            return deepcopy(self)

    def delete_fragments(self, fragment: int):
        """Returns a NEW Cluster Object"""
        self._fragment = fragment
        if not self._check_fragment(self._fragment):
            return deepcopy(self)

        _new_coordinates = deepcopy(self._cluster)
        _new_coordinates.pop(self._fragment)

        _new_cluster = dict()
        for i, fragment in enumerate(_new_coordinates):
            _new_cluster[i] = {"coordinates": _new_coordinates[fragment]}

        return Cluster(*_new_cluster.values())

    def add_fragments(self, new_coordinates):
        """Returns a NEW Cluster Object"""
        _new_cluster = dict()
        for i, fragment in enumerate(self._cluster):
            _new_cluster[i] = {"coordinates": self._cluster[fragment]}

        if type(new_coordinates) == list:
            _new_cluster[len(self._cluster)] = {"coordinates": new_coordinates}
        elif type(new_coordinates) == dict:
            _new_cluster[len(self._cluster)] = new_coordinates
        elif type(new_coordinates) == Cluster:
            _new_cluster[len(self._cluster)] = {
                "coordinates": new_coordinates.coordinates
            }
        else:
            print(
                f"\n* Warning! fragment to add MUST be list-type "
                f" dict-type with 'coordinates' as a key. Check class help"
            )
            return deepcopy(self)

        return Cluster(*_new_cluster.values())

    @property
    def center_of_mass(self) -> tuple:
        """[summary]

        Notes
        -----
            total mass for dummy atoms (not in th ePeriodic Table) is equal to ONE (1)

        Returns
        -------
        tuple : (float, float, float)
            [description]
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

######Andy centro de masa del fragmento i
    @property
    def center_of_mass_fragment(self, fragment) -> tuple:
        """[summary]

        Notes
        -----
            total mass for dummy atoms (not in th ePeriodic Table) is equal to ONE (1)

        Returns
        -------
        tuple : (float, float, float)
            [description]
        """

        self._fragment = fragment
        if not self._check_fragment(self._fragment):
            return deepcopy(self)

        _fragment_to_rotate: Cluster = self.get_fragments(self._fragment)
        _fragment_symbols = _fragment_to_rotate.symbols

        # avoid any rotatation attemp for a single atom system
        if not (len(_fragment_symbols) > 1):
            return deepcopy(self)

        _fragment_center_of_mass = _fragment_to_rotate.center_of_mass

        return _fragment_center_of_mass

    @property
    def principal_axis(self) -> tuple:

        self._principal_axis = np.asarray(self.cartesian_coordinates) - np.asarray(
            self.center_of_mass
        )

        return tuple(self._principal_axis)

    def translate(self, fragment, x=0, y=0, z=0):
        """Returns a NEW Cluster Object with a TRANSLATED fragment"""
        self._fragment = fragment
        if not self._check_fragment(self._fragment):
            return deepcopy(self)

        _fragment_to_move: Cluster = self.get_fragments(self._fragment)
        _fragment_symbols = _fragment_to_move.symbols
        _fragment_coordinates = _fragment_to_move.cartesian_coordinates

        _translated_coordinates = np.asarray(_fragment_coordinates) + np.asarray(
            [x, y, z]
        )

        _translated_fragment = list()
        for i, _atom in enumerate(_fragment_symbols):
            _translated_fragment.append(
                tuple([_atom] + _translated_coordinates[i].tolist())
            )

        return self.delete_fragments(self._fragment).add_fragments(_translated_fragment)

    def rotate(self, fragment, x=0, y=0, z=0):
        """Returns a NEW Cluster Object with a ROTATED fragment (around internal center of mass)"""
        self._fragment = fragment
        if not self._check_fragment(self._fragment):
            return deepcopy(self)

        _fragment_to_rotate: Cluster = self.get_fragments(self._fragment)
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
        _fragment_principal_axis = _fragment_to_rotate.principal_axis

        _rotated_coordinates = (
            np.dot(_fragment_principal_axis, _rotation_matrix)
            + _fragment_center_of_mass
        )

        _rotated_fragment = list()
        for i, _atom in enumerate(_fragment_symbols):
            _rotated_fragment.append(tuple([_atom] + _rotated_coordinates[i].tolist()))

        return self.delete_fragments(self._fragment).add_fragments(_rotated_fragment)
