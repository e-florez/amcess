from base_molecule import Molecule


hydrogen2 = {  # bond distance = 74 pm
    "coordinates": [
        ("H", 1, 1, 0),
        ("H", 1.74, 1, 0),
    ],
    "charge": 0,
    "multiplicity": 1,
}

hydrohen_fluoride = {  # bond distance = 91.7 pm
    "coordinates": [
        ("H", 1, 1, 0),
        ("F", 1.917, 1, 0),
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

w = Molecule(water)
li = Molecule(lithium)
h2 = Molecule(hydrogen2)
hf = Molecule(hydrohen_fluoride)


print(f"\t2")
print(f" rot = 0 deg")
print("\t".join(hf.write_molecule[0]))
print("\t".join(hf.write_molecule[1]))


for i in range(8):
    ang = 45 * (i + 1)
    hf.coordinates = hf.rotate(x=0, y=0, z=45)
    hf.coordinates = hf.translate(x=0.4, y=0, z=0)

    print(f"\t2")
    print(f" rot = {ang} deg")
    print("\t".join(hf.write_molecule[0]))
    print("\t".join(hf.write_molecule[1]))


# print("*" * 50)
# print("END")
