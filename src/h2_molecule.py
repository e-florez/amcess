from base_molecule import Molecule


hydrogen2 = {  # bond distance = 74 pm
    "coordinates": [
        ("H", 1, 1, 0),
        ("H", 1.74, 1, 0),
    ],
    "charge": 0,
    "multiplicity": 1,
}

h2 = Molecule(hydrogen2)


print(f"\t{h2.total_atoms}")
print(f" rot = 0 deg")
for i in range(h2.total_atoms):
    print("\t".join(h2.write_molecule[i]))


for i in range(8):
    ang = 45 * (i + 1)
    h2.coordinates = h2.rotate(x=0, y=0, z=45)
    h2.coordinates = h2.translate(x=0.4, y=0, z=0)

    print(f"\t{h2.total_atoms}")
    print(f" rot = {ang} deg")
    for i in range(h2.total_atoms):
        print("\t".join(h2.write_molecule[i]))
