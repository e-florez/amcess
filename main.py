from src.base_molecule import Atom, Molecule, Cluster

# from src.base_molecule import Molecule, Cluster
from data.molecules_coordinates import (
    li,
    water,
    hydrogen2,
    hydrogen_fluoride,
    nano_boy,
    metal_complex,
)


def system_hf():
    hf = Cluster(hydrogen_fluoride)

    # print(hf)
    # print(hf.total_atoms)
    # print(hf.total_molecules)
    # print(hf.symbols)
    # print(hf.total_mass)
    # print(hf.center_of_mass)
    # print(hf.principal_axes)
    # print(hf.add_molecules(water))
    # print(hf.remove_molecule(4))

    # w = hf.add_molecule(water)
    # print(w.translate(1, x=20, y=10, z=50))
    # print(w.rotate(1, x=0, y=0, z=45))
    # print()

    for i in range(20):
        ang = 45 * (i + 1)
        print(hf)
        hf = hf.rotate(0, x=0, y=0, z=45).translate(0, x=0.3, y=0, z=0)


def system_h2():
    h2 = Cluster(hydrogen2)

    for i in range(20):
        print(h2)
        h2 = h2.rotate(0, x=0, y=0, z=45).translate(0, x=0.4, y=0, z=0)


def system_nano_boy():
    nb = Cluster(li, nano_boy)

    for i in range(2):
        angle = 90 * (-1) ** i
        print(nb)
        nb = nb.translate(1, x=0, y=0, z=0.5).rotate(1, x=angle, y=0, z=0)


def system_w_li():
    w = Cluster(li, water)

    for i in range(20):
        ang = 45 * (i + 1)
        print(w)
        w = w.translate(1, x=0.5, y=0, z=0).rotate(1, x=0, y=0, z=45)


def system_metal_complex():
    complex = Cluster(metal_complex)

    for i in range(20):
        print(complex)
        complex = complex.rotate(0, x=0, y=0, z=45).translate(0, x=2, y=0, z=0)


def system_w_li_cluster():
    # from random import random, seed, randint
    from numpy import random

    w = Cluster(
        li,
        water,
        water,
        water,
        water,
        water,
        water,
    )

    # print(w.xyz)

    total_steps = 1000

    random_gen = random.default_rng(1234)

    probability = [
        1 / (w.total_molecules - 1) if p else 0
        for p in range(w.total_molecules)
    ]

    with open("w6_li.xyz", "w") as file_xyz:
        for _ in range(total_steps):
            # molecule [0, w.total_molecules]
            mol = random_gen.choice(
                w.total_molecules,
                p=probability,
                # p=[0, 1 / 6, 1 / 6, 1 / 6, 1 / 6, 1 / 6, 1 / 6],
                # p=[0, 0, 0, 0, 0, 0, 1],
                # p=[0, 0.1, 0.3, 0.5, 0, 0, 0.1],
            )

            # angle between [0, 360)
            ax = random_gen.uniform() * 360
            ay = random_gen.uniform() * 360
            az = random_gen.uniform() * 360

            # max_t: [-2, 2]
            max_t = 2
            tx = max_t * (random_gen.uniform() - 0.5)
            ty = max_t * (random_gen.uniform() - 0.5)
            tz = max_t * (random_gen.uniform() - 0.5)

            # print(w)
            w = w.translate(mol, x=tx, y=ty, z=tz).rotate(
                mol, x=ax, y=ay, z=az
            )

            file_xyz.write(w.xyz)


def test_atom_class():

    a_dict = {
        "atoms": [
            ("a1", 0, 0, 0),
            ("a2", 0, 0, 0),
            ("a3", 0, 0, 0),
        ],
        "charge": +5,
        "multiplicity": 20,
    }

    b_list = [
        ("b1", 0, 0, 0),
        ("b1", 0, 0, 0),
        ("b3", 0, 0, 0),
    ]

    mol = Cluster(b_list)
    a = Cluster(b_list, a_dict, mol)

    print("-" * 50)
    print(a.xyz)


def test_molecule_class():
    # li2 = Molecule0(li).translate(0, 50, 0, 0)

    xe_coord = {
        "atoms": [
            ("Xe", 0, 0.000000, 0.000000),
        ],
        "charge": +1,
        "multiplicity": 2,
    }

    kf_coord = [
        ("K", 0, 0.000000, 0.000000),
        ("F", 1.96000, 0.000000, 0.000000),
    ]

    nacl_coord = {
        "atoms": [
            ("Na", 0, 0, 0),
            ("Cl", 2.386, 0, 0),
        ],
        "charge": +5,
        "multiplicity": 3,
    }

    na = Molecule(
        # li,
        # water,
        # kf_coord,
        nacl_coord,
    )

    # print(Atom("Na", 0, 0.000000, 0.000000))

    print("-" * 40)
    print(na)
    print("-" * 40)
    print("\natoms: ", na.atoms)
    print("\ncharge: ", na.charge)
    print("\nmultiplicity: ", na.multiplicity)
    print("\nsymbols: ", na.symbols)
    print("\nelements: ", na.elements)
    print("\nmasses: ", na.atomic_masses)
    print("\ntotal mass: ", na.total_mass)
    print("\ncenter of mass: ", na.center_of_mass)
    print("\npp axes: ", na.principal_axes)
    print("\nnumber atoms (index):\n", na.number_atoms)
    r = 0
    print(f"\nremoved atom {r}:\n", na.remove_atom(r))
    g = 0
    print(f"\nget atom {g}:\n", na.get_atom(g))
    print("\nadding new obj:\n", na.add_atoms(na))
    print("\nadding new dict:\n", na.add_atoms(nacl_coord))
    print("\nadding new list:\n", na.add_atoms(kf_coord))
    print("\nmagic add new obj:\n", na + na)
    print("\nmagic add new dict:\n", na + nacl_coord)
    print("\nmagic add new list:\n", na + kf_coord)
    print("\nmagic rmul obj:\n", 3 * na)
    print("\nmagic rmul new dict:\n", na + 4 * Molecule(nacl_coord))
    print("\nmagic rmul new list:\n", na + 3 * Molecule(kf_coord))
    print("\nmagic mul new list:\n", na + Molecule(kf_coord) * 2)
    print("\nremove H: \n", na.remove_element("Na"))
    print("\ngetting H: \n", na.get_element("Na"))

    print("-" * 40)
    print(na)


def test_cluster_class():
    xe_coord = {
        "atoms": [
            ("Xe", 0, 0.000000, 0.000000),
        ],
        "charge": -1,
        "multiplicity": 2,
    }

    kf_coord = [
        ("K", 0, 0.000000, 0.000000),
        ("F", 1.96000, 0.000000, 0.000000),
    ]

    nacl_coord = {
        "atoms": [
            ("Na", 0, 0, 0),
            ("Cl", 2.386, 0, 0),
        ],
        "charge": +5,
        "multiplicity": 3,
    }

    # mol = Molecule(
    #     # li,
    #     water,
    # )

    ca = Cluster(
        # li,
        # water,
        kf_coord,
        nacl_coord,
    )

    # ca = Cluster(na, xe_coord)
    # cb = Cluster(water, xe_coord)  # , kf_coord, nacl_coord)
    # cc = Cluster(
    #     Cluster(water).translate(0, x=0.09),
    #     Cluster(xe_coord),  # .translate(0, x=0.05)
    # )  # , kf_coord, nacl_coord)

    print("-" * 40)
    # print(na.cluster_dictionary)
    # print("-" * 40)
    print(ca.xyz)
    print("-" * 40)
    print("\nmolecules:\n\n", ca.molecules)
    print("\ndictionary:\n\n", ca.cluster_dictionary)
    print("\nsymbols: ", ca.symbols)
    print("\natoms: ", ca.atoms)
    print("\nmolecules: ", ca.total_molecules)
    print("\ntotal atoms: ", ca.total_atoms)
    print("\ncharge: ", ca.charge)
    print("\nmultiplicity: ", ca.multiplicity)
    print("\nelements: ", ca.elements)
    print("\nmasses: ", ca.atomic_masses)
    print("\ntotal mass: ", ca.total_mass)
    # print("\ncenter of mass: ", ca.center_of_mass)
    # print("\npp axes: ", ca.principal_axes)
    # print("\nnumber atoms (index):\n", ca.number_atoms)

    # print("-" * 40)
    # print(ca.cluster_dictionary)
    # print(na.cluster_dictionary)
    # print("\nmagic add new obj:\n", ca + na)
    # print("\nmagic add, T mol: ", (ca + na).total_molecules)
    # print("\nmagic add, T atoms: ", (ca + na).total_atoms)
    # print("\nmagic add dic:\n", (ca + na).cluster_dictionary)

    # print("cluster type: ", type(ca).__name__)
    # print("mol instance mol: ", isinstance(mol, Molecule))
    # print("mol instance cluster: ", isinstance(mol, Cluster))
    # print("ca instance mol: ", isinstance(ca, Molecule))
    # print("ca instance cluster: ", isinstance(ca, Cluster))
    # print("mol sub mol: ", issubclass(type(mol), Molecule))
    # print("mol sub cluster: ", issubclass(type(mol), Cluster))
    # print("ca sub mol: ", issubclass(type(ca), Molecule))
    # print("ca sub cluster: ", issubclass(type(ca), Cluster))

    # print("\nmagic add new dict:\n", ca + nacl_coord)
    # print("\nmagic add new list:\n", ca + kf_coord)

    # ca.frozen_molecule = 1
    # print("frozen molecule: ", ca.frozen_molecule)

    print("Translating other:\n", ca.translate(0, x=90).xyz)
    print("Translating frozen:\n", ca.translate(1, x=90).xyz)
    print("Rotating other:\n", ca.rotate(0, x=90).xyz)
    print("Rotating frozen:\n", ca.rotate(1, x=90).xyz)

    # print("-" * 40)
    # g = 2
    # print(f"\ngetting {g}: \n", ca.get_molecule(g))
    # print("\ngetting element: \n", ca.get_element("H"))
    # a = ca.get_molecule(0)
    # print(a.center_of_mass)

    # r = 2
    # print(f"\nremove {r}: \n", ca.remove_molecule(r))
    # print(f"\nremove element: \n", ca.remove_element("H"))
    # b = ca.remove_molecule(r)
    # print(b.molecules)

    # print(
    #     f"",
    #     ca.translate(m, x=10, y=20, z=30).translate(m + 3, x=10, y=20, z=30),
    # )

    # print("-" * 40)
    # print(ca)
    # print(ca.cluster_dictionary)
    # print("-" * 40)


# -------------------------------------------------------------------
def run():
    pass
    # system_h2()
    # system_hf()
    # system_nano_boy()
    # system_w_li()
    # system_w_li_cluster()
    # system_metal_complex()

    # ##test
    # test_atom_class()
    # test_molecule_class()
    test_cluster_class()


if __name__ == "__main__":
    run()
