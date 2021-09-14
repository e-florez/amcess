from base_molecule import Molecule


hydrogen_fluoride = {  # bond distance = 91.7 pm
    "coordinates": [
        ("H", 1, 1, 0),
        ("F", 1.917, 1, 0),
    ],
    "charge": 0,
    "multiplicity": 1,
}

hf = Molecule(hydrogen_fluoride)

print(f"\t{hf.total_atoms}")
print(f" rot = 0 deg")
for i in range(hf.total_atoms):
    print("\t".join(hf.write_molecule[i]))

for i in range(8):
    ang = 45 * (i + 1)
    hf.coordinates = hf.rotate(x=0, y=0, z=45)
    hf.coordinates = hf.translate(x=0.4, y=0, z=0)

    print(f"\t{hf.total_atoms}")
    print(f" rot = {ang} deg")
    for i in range(hf.total_atoms):
        print("\t".join(hf.write_molecule[i]))
