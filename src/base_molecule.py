import context
from data.atomic_data import atomic_mass

import numpy as np
from scipy.spatial.transform import Rotation


class Molecule:
    def __init__(self, mol) -> None:
        self._molecule = mol["coordinates"]
        self._coordinates = tuple([c[1:] for c in self._molecule])
        self._total_atoms = len(self._molecule)
        self._charge = mol["charge"]
        self._multiplicity = mol["multiplicity"]

    @property
    def total_atoms(self):
        return self._total_atoms

    @property
    def charge(self):
        return self._charge

    @charge.setter
    def charge(self, charge):
        """setter to show atomic charge for a atom or molecule

        Parameters
        ----------
        charge : int
            Atomic charge

        Raises
        ------
        error : ValueError
            Atomic charge must be an integer
        """
        try:
            charge = int(charge)
        except ValueError as error:
            print(f"\n *** ERROR: Atomic/Molecular charge must be an integer \n")
            raise error

        self._charge = charge

    @property
    def multiplicity(self):
        return self._multiplicity

    @charge.setter
    def multiplicity(self, multiplicity):
        try:
            multiplicity = int(multiplicity)
            assert multiplicity > 0
        except (ValueError, AssertionError) as error:
            print(
                f"\n *** ERROR: Atomic/Molecular multiplicity must be an integer larger than zero \n"
            )
            raise error

        self._multiplicity = multiplicity

    @property
    def atoms_symbol(self):
        return tuple([atom[0] for atom in self._molecule])

    @property
    def coordinates(self):
        return self._coordinates

    @coordinates.setter
    def coordinates(self, new_coordinates):
        try:
            pass
        except:
            pass
        self._coordinates = new_coordinates

    @property
    def atomic_masses(self):
        return tuple([atomic_mass(m) for m in self.atoms_symbol])

    @property
    def total_mass(self):
        return sum(self.atomic_masses)

    @property
    def center_of_mass(self):
        """Jacobi coordinates"""
        self.cm = np.dot(np.asarray(self.atomic_masses), np.asarray(self.coordinates))
        return tuple(self.cm / self.total_mass)

    @property
    def principal_axis(self):
        return tuple(np.asarray(self.coordinates) - np.asarray(self.center_of_mass))

    def rotate(self, x=0, y=0, z=0):
        self._rotate_around_x = x
        self._rotate_around_y = y
        self._rotate_around_z = z

        self._rotation_matrix = Rotation.from_euler(
            "xyz",
            [self._rotate_around_x, self._rotate_around_y, self._rotate_around_z],
            degrees=True,
        ).as_matrix()

        self._rotate = np.dot(np.asarray(self.principal_axis), self._rotation_matrix)
        return np.asarray(self._rotate) + np.asarray(self.center_of_mass)

    def translate(self, x=0, y=0, z=0):
        self._translate_around_x = x
        self._translate_around_y = y
        self._translate_around_z = z
        self._translation_point = [
            self._translate_around_x,
            self._translate_around_y,
            self._translate_around_z,
        ]
        return np.asarray(self.coordinates) + np.asarray(self._translation_point)

    @property
    def write_molecule(self):
        self._molecule_list = list()
        for atom in range(self._total_atoms):
            self._molecule_list.append(
                tuple(
                    [self.atoms_symbol[atom]] + [str(c) for c in self.coordinates[atom]]
                )
            )
        return self._molecule_list
