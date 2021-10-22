from abc import abstractclassmethod
from attr import dataclass
from numpy import multiply


water = {
    "atoms": [
        ("O", 0, 0, 0),
        ("H", 0.58708, 0.75754, 0),
        ("H", -0.58708, 0.75754, 0),
    ],
    "charge": 0,
    "multiplicity": 1,
}


class mol:
    def __init__(self, dict):
        self.coordinates = dict.get("atoms")
        self.charge = dict.get("charge")
        self._multiplicity = dict.get("multiplicity")

    @property
    def multiplicity(self):
        return self._multiplicity

    @multiplicity.setter
    def multiplicity(self, new_multiplicity) -> None:
        if isinstance(new_multiplicity, int) and new_multiplicity > 0:
            self._multiplicity = new_multiplicity
        else:
            raise ValueError("No str")\

    @multiply.deleter

    @dataclass

    @classmethod


w = mol(water)


print("\n\n", "-" * 30)
print("charge: ", w.charge)
print("multiplicity: ", w.multiplicity)

w.charge = 10
w.multiplicity = "hola"

print("\n", "-" * 30)
print("charge: ", w.charge)
print("multiplicity: ", w.multiplicity)
