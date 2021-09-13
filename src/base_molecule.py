import context
from data.atomic_data import atomic_mass


import numpy as np


class Molecule:
    def __init__(self, mol) -> None:
        self.__coordinates = mol["coordinates"]
        self.__total_atoms = len(self.__coordinates)
        self.__charge = mol["charge"]
        self.__multiplicity = mol["multiplicity"]
        # self.__atomic_symbol =
        # self.__coordinates =
        # self.__atomic_masses =
        # self.__total_mass =

    @property
    def total_atoms(self):
        return self.__total_atoms

    @property
    def charge(self):
        return self.__charge

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

        self.__charge = charge

    @property
    def multiplicity(self):
        return self.__multiplicity

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

        self.__multiplicity = multiplicity

    @property
    def atoms_symbol(self):
        return tuple([atom[0] for atom in self.__coordinates])

    @property
    def coordinates(self):
        return [c[1:] for c in self.__coordinates]

    @property
    def atomic_masses(self):
        return tuple([atomic_mass(m) for m in self.atoms_symbol])

    @property
    def total_mass(self):
        return sum(self.atomic_masses)

    @property
    def center_of_mass(self):
        """Jacobi coordinates"""
        cm = np.dot(np.asarray(self.atomic_masses), np.asarray(self.coordinates))
        return tuple(cm / self.total_mass)

    def principal_axis(self):
        pass
        # principal = list()
        # for atom in self.coordinates:


dummy = {
    "coordinates": [
        ("Xa", -1, 1, 0),
        ("Xb", 0, 0, 0),
        ("Xc", 1, 1, 0),
    ],
    "charge": 0,
    "multiplicity": 3,
}

hydrogen2 = {  # bond distance = 74 pm
    "coordinates": [
        ("H", 1, 1, 0),
        ("H", 1.74, 1, 0),
    ],
    "charge": 0,
    "multiplicity": 1,
}

water = {
    "coordinates": [
        ("H", 0.000000, 0.000000, 0.000000),
        ("H", 0.758602, 0.000000, 0.504284),
        ("H", 0.758602, 0.000000, -0.504284),
    ],
    "charge": 0,
    "multiplicity": 1,
}

lithium = {
    "coordinates": [
        ("li", 0, 0, 0),
    ],
    "charge": 1,
    "multiplicity": +1,
}

print("-" * 50)


w = Molecule(water)
li = Molecule(lithium)
du = Molecule(dummy)
h2 = Molecule(hydrogen2)

w.multiplicity = 13

w.charge

w.charge = "10"

w.atoms_symbol

# w.center_of_mass
# li.center_of_mass
# du.center_of_mass
# h2.center_of_mass


# # w_coord = w.coordinates
# w_coord = np.asarray(w.coordinates)
# # w_masses = w.atomic_masses
# w_masses = np.asarray(w.atomic_masses)
# w_total_mass = w.total_mass

# # np.array(self.atomic_masses) * np.array(self.coordinates) / self.total_mass

# cm = tuple(np.dot(w_masses, w_coord) / w_total_mass)


b = w.coordinates[0]
b = w.coordinates[2]


a = np.array(w.coordinates)

print("*" * 50)
print("END")
