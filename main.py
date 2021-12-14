from amcess.base_molecule import Atom, Molecule, Cluster

from data.molecules_coordinates import (
    li,
    water,
    hydrogen2,
    hydrogen_fluoride,
    nano_boy,
    metal_complex,
    ibuprofen,
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


def test_atom():
    # a = ("a", 0, 0, 0)
    # ---
    # a = ("a", "x", 0, 0)
    # a = ("`", 0, 0, 0)
    # a = ("@", 0, 0, 0)
    # a = ("a", 0, 0)
    a = ("a", 0, 0, 0)

    a = Atom(*a)

    print("\n")
    print(a)
    print("-" * 30)
    print("atoms:\n", a)


def test_molecule():

    # mol = Molecule([("a", 0, 0, 0), ("b", 10, 10, 10)], -10, 5))
    # ---
    mol = Molecule(
        [("Xe", 0, 0, 0), ("Na", 5, 5, 5), ("Na", 10, 10, 10)], -10, 5
    )

    # print("\n")
    # print(mol)
    # print("+--" * 30)
    # print("atoms: :\n", mol.atoms)
    # # mol.atoms = 1  # cannot reset atoms

    # print("charge:", mol.charge)
    # print("multiplicity:", mol.multiplicity)

    # print("+--" * 30)
    # mol.charge = 200
    # print("new charge: ", mol.charge)
    # # mol.charge = "a"  # wrong charge
    # # mol.charge = 10.2  # wrong charge
    # mol.multiplicity = 13
    # print("new multiplicity", mol.multiplicity)
    # # mol.multiplicity = -3  # wrong multiplicity
    # # mol.multiplicity = "x"  # wrong multiplicity

    # print("+--" * 30)
    # print("\nmasses: ", mol.atomic_masses)
    # print("\ntotal mass: ", mol.total_mass)
    # print("\nsymbols: ", mol.symbols)
    # print("\ntotal atoms: ", mol.total_atoms)
    # print("\nelements: ", mol.elements)

    # print("+--" * 30)
    # print("\nnumbering atoms (index):\n", mol.numbering_atoms)
    # print(mol.xyz)

    # print("+--" * 30)
    # print("\ncenter of mass: ", mol.center_of_mass)
    # print("\npp axes: ", mol.principal_axes)


def test_cluster():
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

    kf_mol = Molecule(kf_coord)
    w_mol = Molecule.from_dict(water)
    li_mol = Molecule.from_dict(li)

    # init from dict and list
    # ca = Cluster(nacl_coord, nacl_coord, nacl_coord, nacl_coord)
    # ca = Cluster(kf_coord, nacl_coord)

    # # init from Molecule
    # ca = Cluster(xe_coord)

    # # init from Molecule and dict, list
    # ca = Cluster(w_mol, kf_coord, nacl_coord)

    # print("+--" * 30)
    # print(ca)
    # print("+--" * 30)
    # print(ca.xyz)
    # print("+--" * 30)
    # print(ca.cluster_dictionary)
    # print("\nmolecules:", ca.total_molecules)
    # print("\nsymbols: ", ca.symbols)
    # print("\natoms: ", ca.atoms)
    # print("\ntotal molecules: ", ca.total_molecules)
    # print("\ntotal atoms: ", ca.total_atoms)
    # print("\ncharge: ", ca.charge)
    # print("\nmultiplicity: ", ca.multiplicity)
    # print("\nelements: ", ca.elements)
    # print("\nmasses: ", ca.atomic_masses)
    # print("\ntotal mass: ", ca.total_mass)
    # print("\ncenter of mass: ", ca.center_of_mass)
    # print("\npp axes: ", ca.principal_axes)
    # print("\nnumber atoms (index):\n", ca.numbering_atoms)

    # print("+--" * 30)
    # # print("water mol:", w_mol)
    # # print("kf mol:", kf_mol)

    # magic_add = Cluster(w_mol + kf_mol)

    # print("+--" * 30)
    # print("\nmagic add new obj:\n", magic_add)

    # print("\nmagic add new dict:\n", magic_add + nacl_coord)
    # print("\nmagic add new list:\n", magic_add + kf_coord)
    # print("\nmagic add new mol:\n", magic_add + kf_mol)

    # print("+--" * 30)
    # magic_add = Cluster(w_mol)
    # magic_add = Cluster(kf_coord)
    # magic_add = Cluster(nacl_coord)
    # magic_add = Cluster(w_mol, nacl_coord)
    # magic_add = Cluster(w_mol, kf_mol)

    # print("\nmagic rmul cluster:\n", 3 * magic_add)
    # print("\nmagic mul cluster:\n", magic_add * 3)

    # print("+--" * 30)
    # magic_add = Cluster(w_mol, kf_mol)
    # ca = magic_add + w_mol
    # print(ca.xyz)

    # print("+--" * 30)
    # g = ca.get_molecule(1)

    # m = ca.get_molecule(0)
    # print(m.coordinates)

    # r = ca.remove_molecule(0)
    # print(r.coordinates)

    # print(type(g))
    # print(g)
    # g = g.translate(0, x=100)
    # print(g)

    # ca.frozen_molecule = 1
    # # ca.frozen_molecule = [0, 1]
    # print("frozen molecule: ", ca.frozen_molecule)

    # print("No sphere radius: ", ca.sphere_radius)
    # print("frozen molecule: ", ca.frozen_molecule)
    # ca.sphere_radius = 20
    # print("choosing a sphere radius: ", ca.sphere_radius)

    # cb = ca.translate(0, x=100)
    # print("New sphere radius: ", cb.sphere_radius)
    # print("frozen molecule: ", cb.frozen_molecule)

    # print("Translating other:\n", ca.translate(0, x=90).xyz)
    # print("Translating frozen:\n", ca.translate(1, x=90).xyz)
    # print("Rotating other:\n", ca.rotate(0, z=90).xyz)
    # print("Rotating frozen:\n", ca.rotate(1, z=90).xyz)

    # print("+--" * 30)

    # print("+--" * 30)
    # g = 2
    # print(f"\ngetting {g}: \n", ca.get_molecule(g))
    # print("\ngetting element: \n", ca.get_element("H"))
    # a = ca.get_molecule(0)
    # print(a.center_of_mass)

    # m = 0
    # print(
    #     f"",
    #     ca.translate(m, x=10, y=20, z=30).translate(m + 1, x=10, y=20, z=30),
    # )

    # ca = Cluster(3 * Molecule.from_dict(water))

    # print("+--" * 30)
    # print(ca.atoms)
    # print(ca.xyz)
    # print("+--" * 30)

    # print("No sphere radius: ", ca.sphere_radius)
    # ca.sphere_radius = 1
    # print("choosing a sphere radius: ", ca.sphere_radius)

    # b = ca.move_molecule(
    #     molecule=2,
    #     max_step=None,
    #     max_rotation=None,
    #     seed=1234,
    # )
    # print(b.xyz)

    # ca = Cluster(2 * Molecule.from_dict(water))
    # print(ca.xyz)
    # # print("+--" * 30)

    # # print("sphere None: ", ca.sphere_center)
    # ca._sphere_center = (10, 0, 0)
    # # print("sphere: ", len(ca.sphere_center))

    water_molecule = [
        ("O", 0, 0, 0),
        ("H", 0.58708, 0.75754, 0),
        ("H", -0.58708, 0.75754, 0),
    ]

    w = Cluster(water_molecule)

    print(w.xyz)

    # for i in range(3):
    #     #     print(w.move_molecule(0).xyz)
    #     new_w = w.move_molecule(0)
    #     print(new_w.xyz)

    w2 = 2 * w

    # print(w2.freeze_molecule)

    w2.freeze_molecule = 0

    print(w2.freeze_molecule)
    print(w2.cluster_dictionary)

    print(w2.xyz)

    for i in range(3):
        new_w2 = w2.move_molecule(1)
        print(new_w2.xyz)

    # new_w2 = w2.initialize_cluster(max_closeness=1.0)

    # print(new_w2.xyz)

    # for i in range(5):
    #     new_w2 = new_w2.move_molecule(1, max_closeness=2)
    #     print(new_w2.xyz)

    # ca = 3 * w2

    # # let's define the spherical boundary conditions
    # ca.sphere_center = 0, 0, 0
    # w2.sphere_center = 30, 30, 30
    # ca.sphere_radius = 20

    # ini = ca.initialize_cluster(max_closeness=10)

    # print(ini.xyz)


def move_init():

    hf = {"atoms": [("H", 0, 0, 0), ("F", 1, 0, 0)]}
    cc = {"atoms": [("C", 0, 0, 0), ("C", 1.2, 0, 0)]}
    ab = Cluster(hf, cc)

    ab_xyz = ""

    ab.sphere_radius = 4

    # ab.freeze_molecule = 0

    for _ in range(100):
        mol = ab.random_generator.choice(ab.total_molecules)

        ab = ab.move_molecule(
            molecule=mol,
            max_step=3,
            max_rotation=10,
            max_closeness=0.5,
        )
        # print(ab.xyz)
        ab_xyz += ab.xyz

    file_name = "cluster_test.xyz"
    with open(file_name, "w") as f:
        f.write(ab_xyz)


def old_molecule_class():
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

    print("+--" * 30)
    print(na)
    print("+--" * 30)
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

    print("+--" * 30)
    print(na)


def old_cluster_class():
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

    print("+--" * 30)
    # print(na.cluster_dictionary)
    # print("+--" * 30)
    print(ca.xyz)
    print("+--" * 30)
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

    # print("+--" * 30)
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

    # print("+--" * 30)
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

    # print("+--" * 30)
    # print(ca)
    # print(ca.cluster_dictionary)
    # print("+--" * 30)


# # -------------------------------------------------------------------
# def run():
#     pass
#     # system_h2()
#     # system_hf()
#     # system_nano_boy()
#     # system_w_li()
#     # system_w_li_cluster()
#     # system_metal_complex()

#     # ##test
#     # test_atom()
#     # test_molecule()
#     test_cluster()


def pipeline_andy():
    from amcess.base_molecule import Cluster, Molecule
    from amcess import engine as SearchConfig

    HF = [("H", 0, 0, 0), ("F", 0.917, 0, 0)]
    HF1 = {"atoms": [("H", 0, 0, 1), ("F", 0.918, 0.0, 1)]}

    DHF = Cluster(HF, HF1)
    SC = SearchConfig(DHF)

    # print(SC.search_type)
    # SC.search_type = "Bayesian"
    # SC.search_type = "SHGO"
    # SC.search_type = "dual_annealing"
    # SC.sphere_radius = 2.0
    # SC.tolerance_contour_radius = 4.9
    # SC.run(gp_params={"initer": 3, "maxiter": 3})  # Bayesian
    # SC.run(sampling_method="sobol", n=1)
    # SC.run(maxfun=1,maxiter=1)  # !da

    SC.run(nT=1, maxCycle=3)  #


# -------------------------------------------------------------------
def pipeline_dual():
    from amcess import Cluster
    from amcess import engine

    # hf = [("H", 0, 0, 0), ("F", 0.917, 0, 0)]
    # hf_1 = {"atoms": [("H", 0, 0, 1), ("F", 0.918, 0.0, 1)]}

    # hf_dimer = Cluster(hf, hf_1)

    # water = {
    #     "atoms": [
    #         ("O", 0, 0, 0),
    #         ("H", 0.58708, 0.75754, 0),
    #         ("H", -0.58708, 0.75754, 0),
    #     ],
    #     "charge": 0,
    #     "multiplicity": 1,
    # }

    # water_dimer = Cluster(water, water)
    # water_dimer.initialize_cluster

    # print(water_dimer.xyz)

    oxygen = [("O", 0, 0, 0)]
    hydrogen = [("H", 0, 0, 0)]

    water = Cluster(oxygen, hydrogen, hydrogen)

    water.sphere_center = (0, 0, 0)
    water.sphere_radius = 1.0

    water.initialize_cluster()

    search = ["ASCEC", "dual_annealing", "SHGO", "Bayesian"]

    simulation = engine(
        water,
        # water_dimer,
        search_methodology=search[3],
    )

    # simulation.run()
    # simulation.run(T0=1000, dT=0.2, nT=5, maxCycle=10)  # ASCEC
    # simulation.run(maxfun=10, maxiter=10)  #!da
    # simulation.run(sampling_method="sobol", n=1)
    simulation.run(initer=10, maxiter=10)  # bayesian


# -------------------------------------------------------------------
def temperature_grid():
    t_ini = 10
    t_end = 1000
    n = 4

    # warmming up
    geometric = [t_ini * (t_end / t_ini) ** (i / (n - 1)) for i in range(n)]
    arithmetic = [t_ini + i * (t_end - t_ini) / (n - 1) for i in range(n)]

    # cooling down (ANNEALING)
    geometric = [t_end * (t_ini / t_end) ** (i / (n - 1)) for i in range(n)]
    arithmetic = [t_end + i * (t_ini - t_end) / (n - 1) for i in range(n)]

    print("\n")
    print(geometric)
    print("\n")
    print(arithmetic)


# -------------------------------------------------------------------
def run():
    from amcess import Cluster
    from time import time
    from amcess.search_engine import simulated_annealing

    water = {
        "atoms": [
            ("O", 0, 0, 0),
            ("H", 0.58708, 0.75754, 0),
            ("H", -0.58708, 0.75754, 0),
        ],
        "charge": 0,
        "multiplicity": 1,
    }

    water_dimer = Cluster(water, water)
    water_dimer.sphere_radius = 1.5
    water_dimer.sphere_center = (0, 0, 0)

    water_dimer = water_dimer.initialize_cluster(max_closeness=0.5)

    print("\n")
    start_time = time()

    simulation = simulated_annealing(
        cluster_settings={
            "cluster": water_dimer,
            "max_closeness": 0.5,
            "max_step": 0.6,
            "max_rotation": 30,
        },
        energy_settings={
            "energy": "pyscf",
            "hamiltonian": "b3lyp",
            "basis": "aug-cc-pvdz",
            "cpu": 4,
        },
        annealing_settings={
            "criterion": "metropolis",  # 'metropolis' or 'delta_energy'
            "file_name": "water_dimer",
            "temperature_grid": "geometric",
            "max_cycle": 1000,
        },
        temperature_settings={
            "initial_temperature": 2.0,
            "final_temperature": 700.0,
            "total_temperatures": 10,
        },
    )

    print(simulation)

    end_time = time()
    print(f"\n elapsed time: {end_time - start_time:.2f} seconds")
    print("+--" * 30)
    print("\n")


# -------------------------------------------------------------------
if __name__ == "__main__":

    # run()

    # move_init()

    # temperature_grid()

    pipeline_dual()
    # pipeline_andy()
