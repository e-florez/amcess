from base_molecule import Molecule

water = {
    "coordinates": [
        ("O", 0.000000, 0.000000, 0.000000),
        ("H", 0.758602, 0.000000, 0.504284),
        ("H", 0.758602, 0.000000, -0.504284),
    ],
    "charge": 0,
    "multiplicity": 1,
}

w = Molecule(water)

print(f"\t{w.total_atoms}")
print(f" rot = 0 deg")
for i in range(w.total_atoms):
    print("\t".join(w.write_molecule[i]))


for i in range(8):
    ang = 45 * (i + 1)
    w.coordinates = w.rotate(x=0, y=0, z=45)
    w.coordinates = w.translate(x=0.4, y=0, z=0)

    print(f"\t{w.total_atoms}")
    print(f" rot = {ang} deg")
    for i in range(w.total_atoms):
        print("\t".join(w.write_molecule[i]))


# print("*" * 50)
# print("END")
