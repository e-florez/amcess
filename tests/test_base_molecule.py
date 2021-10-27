from _pytest.python import write_docstring
import numpy as np
import pytest

from amcess.base_molecule import Atom, Molecule, Cluster

COORDINATES = {
    "dummy": [("A", 10, 20, 30), ("X0", -0.5, 0, -10)],
    "hf": [("H", 0, 0, 0), ("F", 0.917, 0, 0)],
    "water": [
        ("O", 0, 0, 0),
        ("H", 0.58708, 0.75754, 0),
        ("H", -0.58708, 0.75754, 0),
    ],
}


# ===============================================================
# Atom class
# ===============================================================
@pytest.mark.parametrize(
    "atom, expected_result",
    [
        (("H", 0, 0, 0), {"element": "H", "x": 0, "y": 0, "z": 0}),
    ],
)
def test_atom_class(atom, expected_result):
    """Testing Atom class

    Examples
    --------
    >>>Atom('F', 0, 0, 1.97)
    {'element': 'F', 'x': 0, 'y': 0, 'z': 1.97}
    """
    input_atom = str(Atom(*atom))

    assert input_atom == str(expected_result)


@pytest.mark.parametrize(
    "wrong_atom",
    [
        (("H", 0, 0, "")),
        (("H", 0, 0, "z")),
        (("", 0, 0, 0)),
        (("'", 0, 0, 0)),
        (("@", 0, 0, 0)),
    ],
)
def test_atom_class_fails(wrong_atom):
    """Testing Atom class failures"""
    with pytest.raises(ValueError):
        Atom(*wrong_atom)


@pytest.mark.parametrize(
    "atom, expected_result",
    [
        (("O", 0, 0, 0), 16),
    ],
)
def test_atom_mass(atom, expected_result):
    """Testing Atomic mass"""
    input_atom = Atom(*atom)

    assert input_atom.atomic_mass - expected_result < 0.1


@pytest.mark.parametrize(
    "atom, expected_result",
    [
        (("Xe20", 0, 0, 0), "Xe20"),
    ],
)
def test_atom_symbol(atom, expected_result):
    """Testing Atomic mass"""
    input_atom = Atom(*atom)

    assert input_atom.symbol == expected_result


# ===============================================================
# Molecule class
# ===============================================================
def test_molecule_class_init():
    """
    Test class init
    """
    mol = Molecule([("A", 10, 20, 30), ("Z", 0, 0, 0)], -1, 10)
    assert mol.charge == -1
    assert mol.multiplicity == 10


def test_molecule_class_init_wrong_atom():
    """
    Test class init
    """
    with pytest.raises(TypeError):
        Molecule([("A", 10, 20, "")], -1, 10)
    with pytest.raises(TypeError):
        Molecule([("", 10, 20, 0)], -1, 10)
    with pytest.raises(TypeError):
        Molecule([("@", 10, 20, 0)], -1, 10)


def test_molecule_class_init_wrong_charge():
    """
    Test class init
    """
    with pytest.raises(ValueError):
        Molecule([("A", 10, 20, 0)], 0.4, 10)


def test_molecule_class_init_wrong_multiplicity():
    """
    Test class init
    """
    with pytest.raises(ValueError):
        Molecule([("A", 10, 20, 0)], -2, 0)


def test_molecule_init_from_dict():
    """
    Test class init
    """
    mol = Molecule.from_dict(
        {
            "atoms": [("A", 10, 20, 0), ("B", 0, 0, 0)],
            "charge": -2,
            "multiplicity": 10,
        }
    )

    assert mol.atoms == [("A", 10, 20, 0), ("B", 0, 0, 0)]
    assert mol.total_atoms == 2
    assert mol.charge == -2
    assert mol.multiplicity == 10


def test_molecule_init_from_dict_wrong():
    """
    Test class init
    """
    with pytest.raises(TypeError):
        Molecule.from_dict(
            {
                ".......": [("A", 10, 20, 0), ("B", 0, 0, 0)],
                "charge": -2,
                "multiplicity": 10,
            }
        )


def test_molecule_magic_add():
    """Testing add molecule method (__add__)"""
    mol1 = Molecule([("A", 0, 0, 0), ("B", 1, 1, 1)])
    mol2 = Molecule([("X", 0, 0, 0), ("Y", 1, 1, 1)])

    mol = Molecule([("A", 0, 0, 0), ("B", 1, 1, 1)]).add_molecule(mol2)
    new_mol = mol1 + mol2

    assert isinstance(new_mol, Cluster)
    assert isinstance(mol, Cluster)
    assert mol.atoms == new_mol.atoms
    assert new_mol.total_atoms == mol1.total_atoms + mol2.total_atoms
    assert new_mol.symbols == mol1.symbols + mol2.symbols
    assert new_mol.atoms == mol1.atoms + mol2.atoms


def test_molecule_magic_add_fail():
    """Testing add molecule method (__add__)"""
    mol = Molecule([("A", 0, 0, 0), ("B", 1, 1, 1)])
    with pytest.raises(TypeError):
        mol.add_molecule(Atom("H", 0, 0, 0))

    with pytest.raises(TypeError):
        mol.add_molecule("H", 0, 0, 0)


def test_molecule_magic_mul_rmul():
    """Testing magic mul"""
    mol1 = Molecule([("A", 0, 0, 0), ("B", 1, 1, 1)])
    mol2 = mol1 * 2
    mol3 = 3 * mol1
    assert mol2.atoms == 2 * mol1.atoms
    assert mol3.atoms == 3 * mol1.atoms


def test_molecule_magic_mul_rmul_fail():
    """Testing magic mul"""
    mol1 = Molecule([("A", 0, 0, 0), ("B", 1, 1, 1)])

    with pytest.raises(ValueError):
        0 * mol1
    with pytest.raises(ValueError):
        -3 * mol1


@pytest.mark.parametrize(
    "system, expected_result",
    [
        (
            [("Xe", 0, 0, 0)],
            """\t1\n-- charge=0 and multiplicity=1 --\nXe    """
            """\t     0.00000000\t     0.00000000\t     0.00000000\n""",
        )
    ],
)
def test_molecule_str(system, expected_result):
    """
    Test for list of atomic symbols
    """
    str_mol = Molecule(system)
    assert str(str_mol) == str(str_mol.xyz)
    assert str(str_mol) == expected_result


@pytest.mark.parametrize(
    "system, expected_atoms",
    [
        ("dummy", [("A", 10, 20, 30), ("X0", -0.5, 0, -10)]),
    ],
)
def test_molecule_atoms(system, expected_atoms):
    """
    Test for coordinate in XYZ format
    """
    atoms_list = Molecule(COORDINATES[system]).atoms

    assert atoms_list == expected_atoms


@pytest.mark.parametrize(
    "system, expected_atoms",
    [
        ("dummy", [("A", 10, 20, 30), ("X0", -0.5, 0, -10)]),
    ],
)
def test_molecule_setting_atoms_fail(system, expected_atoms):
    """
    Test for coordinate in XYZ format
    """
    mol = Molecule(COORDINATES[system])
    with pytest.raises(AttributeError):
        mol.atoms = expected_atoms


@pytest.mark.parametrize(
    "system, expected_masses",
    [
        ("dummy", [0, 0]),
        ("hf", [1, 19]),
        ("water", [16, 1, 1]),
    ],
)
def test_molecule_atomic_masses(system, expected_masses):
    """
    Test for atomic symbols
    """
    masses = Molecule(COORDINATES[system]).atomic_masses

    assert np.allclose(masses, expected_masses, 0.1)


@pytest.mark.parametrize(
    "system, expected_charge, expected_multiplicity",
    [
        ("dummy", 0, 1),
    ],
)
def test_molecule_default_charge_multiplicity(
    system, expected_charge, expected_multiplicity
):
    """
    Test default atomic charge and multiplicity
    """
    mol = Molecule(COORDINATES[system])

    charge = mol.charge
    multiplicity = mol.multiplicity

    assert charge == expected_charge
    assert multiplicity == expected_multiplicity


@pytest.mark.parametrize(
    "system, expected_charge",
    [
        ("dummy", -1),
        ("dummy", +5),
    ],
)
def test_molecule_setting_charge(system, expected_charge):
    """
    Test setting atomic charge
    """
    mol = Molecule(COORDINATES[system])

    mol.charge = expected_charge
    assert mol.charge == expected_charge


@pytest.mark.parametrize(
    "system, wrong_charge",
    [
        ("dummy", "x"),
        ("dummy", 1.3),
        ("dummy", ""),
    ],
)
def test_molecule_setting_charge_fails(system, wrong_charge):
    """
    Test setting atomic charge
    """
    mol = Molecule(COORDINATES[system])

    with pytest.raises(ValueError):
        mol.charge = wrong_charge


@pytest.mark.parametrize(
    "system, expected_coordinates",
    [
        ("dummy", [(10, 20, 30), (-0.5, 0, -10)]),
    ],
)
def test_molecule_coordinates(system, expected_coordinates):
    """
    Test for coordinate in XYZ format
    """
    coordinates = Molecule(COORDINATES[system]).coordinates

    assert coordinates == expected_coordinates


@pytest.mark.parametrize(
    "system, expected_elements",
    [
        ("dummy", ["A", "X0"]),
        ("hf", ["H", "F"]),
        ("water", ["O", "H"]),
    ],
)
def test_molecule_element_list(system, expected_elements):
    """
    Test for list of uniques elements
    """
    elements = Molecule(COORDINATES[system]).elements

    assert elements.sort() == expected_elements.sort()


@pytest.mark.parametrize(
    "system, expected_multiplicity",
    [
        ("dummy", 3),
        ("dummy", 15),
    ],
)
def test_molecule_setting_multiplicity(system, expected_multiplicity):
    """
    Test setting atomic multiplicity
    """
    mol = Molecule(COORDINATES[system])

    mol.multiplicity = expected_multiplicity
    assert mol.multiplicity == expected_multiplicity


@pytest.mark.parametrize(
    "system, wrong_multiplicity",
    [
        ("dummy", -1),
        ("dummy", 0),
        ("dummy", "x"),
        ("dummy", 1.3),
        ("dummy", ""),
    ],
)
def test_molecule_setting_multiplicity_fails(system, wrong_multiplicity):
    """
    Test setting atomic charge
    """
    mol = Molecule(COORDINATES[system])

    with pytest.raises(ValueError):
        mol.multiplicity = wrong_multiplicity


@pytest.mark.parametrize(
    "system, expected_result",
    [
        (
            [("Xe", 0, 0, 0)],
            """\r  atom #0 --> """
            """Xe         0.00000000     0.00000000     0.00000000""",
        )
    ],
)
def test_molecule_numbering_atoms(system, expected_result):
    """
    Test for list of atomic symbols
    """
    numbered_atoms = Molecule(system).numbering_atoms

    assert numbered_atoms == expected_result


@pytest.mark.parametrize(
    "system, expected_symbols",
    [
        ("dummy", ["A", "X0"]),
        ("hf", ["H", "F"]),
        ("water", ["O", "H", "H"]),
    ],
)
def test_molecule_atomic_symbols(system, expected_symbols):
    """
    Test for list of atomic symbols
    """
    symbols = Molecule(COORDINATES[system]).symbols

    assert symbols == expected_symbols


@pytest.mark.parametrize(
    "system, expected_atoms",
    [
        ("dummy", 2),
        ("hf", 2),
        ("water", 3),
    ],
)
def test_molecule_total_atoms(system, expected_atoms):
    """
    Test for the total number of atoms
    """
    number_of_atoms = Molecule(COORDINATES[system]).total_atoms

    assert (number_of_atoms - expected_atoms) == 0


@pytest.mark.parametrize(
    "system, expected_mass",
    [
        ("dummy", 0),
        ("hf", 20),
        ("water", 18),
    ],
)
def test_molecule_total_mass(system, expected_mass):
    """
    Test for total molecular mass
    """
    mass = Molecule(COORDINATES[system]).total_mass

    assert abs(mass - expected_mass) < 1.0e-1


@pytest.mark.parametrize(
    "system, expected_cm",
    [
        ("dummy", (0, 0, 0)),
        ("hf", (0.9, 0, 0)),
        ("water", (0, 0, 0)),
    ],
)
def test_molecule_center_of_mass(system, expected_cm):
    """
    Test center of mass
    """
    cm = Molecule(COORDINATES[system]).center_of_mass

    assert np.linalg.norm(np.asarray(cm) - expected_cm) < 0.1


@pytest.mark.parametrize(
    "system, expected_principal_axes",
    [
        ("dummy", [(10, 20, 30), (-0.5, 0, -10)]),
        ("hf", [(-0.9, 0, 0), (0.05, 0, 0)]),
        ("water", [(0, -0.09, 0), (0.6, 0.7, 0), (-0.6, 0.7, 0)]),
    ],
)
def test_molecule_principal_axes(system, expected_principal_axes):
    """
    Test principal axes using internal coordinates system.
    Number of principal axis is equal to number of atoms
    """
    pp = Molecule(COORDINATES[system]).principal_axes

    assert np.allclose(pp, expected_principal_axes, 0.1)


@pytest.mark.parametrize(
    "system, expected_result",
    [
        (
            [("Xe", 0, 0, 0)],
            """\t1\n-- charge=0 and multiplicity=1 --\nXe    """
            """\t     0.00000000\t     0.00000000\t     0.00000000\n""",
        )
    ],
)
def test_molecule_xyz_format(system, expected_result):
    """
    Test for list of atomic symbols
    """
    xyz_file = Molecule(system).xyz
    assert str(xyz_file) == expected_result


@pytest.mark.parametrize(
    "system, new_atom",
    [
        ("dummy", [("H", 0, 0, 0)]),
        ("hf", [("H", 0, 0, 0)]),
        ("water", [("H", 0, 0, 0)]),
    ],
)
def test_molecule_add_atom(system, new_atom):
    """
    Test adding a new atom
    """
    mol = Molecule(COORDINATES[system])

    # duplicating the system
    new_system = mol.add_atoms(new_atom)

    assert new_system.total_atoms == (mol.total_atoms + 1)


def test_molecule_add_atom_fail():
    """
    Test adding a new atom only as a lits
    """
    mol = Molecule([("H", 0, 0, 0)])
    with pytest.raises(TypeError):
        mol.add_atoms(("H", 0, 0, 0))
    with pytest.raises(TypeError):
        mol.add_atoms("H", 0, 0, 0)


@pytest.mark.parametrize(
    "system",
    [
        ("dummy"),
        ("hf"),
        ("water"),
    ],
)
def test_molecule_get_atom(system):
    """
    Test retriving certain atom
    """
    mol = Molecule(COORDINATES[system])

    # getting the first (index) molecule
    new_system = mol.get_atom(0)

    assert isinstance(new_system, tuple)
    assert len(new_system) == 4


def test_molecule_get_atom_fail():
    """
    Test retriving certain atom
    """
    mol = Molecule(COORDINATES["water"])

    with pytest.raises(IndexError):
        mol.get_atom(100)
    with pytest.raises(IndexError):
        mol.get_atom("@")


@pytest.mark.parametrize(
    "system",
    [
        ("dummy"),
        ("hf"),
        ("water"),
    ],
)
def test_molecule_remove_atom(system):
    """
    Test removing an atom
    """
    mol = Molecule(COORDINATES[system])

    # duplicating the system
    new_system = mol.remove_atom(0)

    assert new_system.total_atoms == mol.total_atoms - 1


def test_molecule_remove_atom_fail():
    """
    Test removing certain atom
    """
    mol = Molecule(COORDINATES["water"])

    with pytest.raises(IndexError):
        mol.remove_atom(100)
    with pytest.raises(IndexError):
        mol.remove_atom("@")


# ===============================================================
# Cluster class
# ===============================================================
def test_cluster_init_list():
    """
    Test for cluster init from list ype
    """
    w = Cluster(COORDINATES["water"])

    assert isinstance(w, Cluster)
    assert w.atoms == COORDINATES["water"]


def test_cluster_init_dict():
    """
    Test for cluster init from dictionary type
    """
    water_dict = {"atoms": COORDINATES["water"]}

    w = Cluster(water_dict)

    assert isinstance(w, Cluster)
    assert w.atoms == COORDINATES["water"]


def test_cluster_init_molecule():
    """
    Test for cluster init from Molecule class
    """
    water_dict = {"atoms": COORDINATES["water"]}

    w = Cluster(Molecule.from_dict(water_dict))

    assert isinstance(w, Cluster)
    assert w.atoms == COORDINATES["water"]


def test_cluster_init_cluster():
    """
    Test for cluster init from Cluster class
    """
    w = Cluster(COORDINATES["water"])
    w2 = Cluster(w, w)

    assert isinstance(w2, Cluster)
    assert w2.atoms == 2 * COORDINATES["water"]


def test_cluster_init_fail():
    """
    Test for cluster init raise
    """

    with pytest.raises(TypeError):
        Cluster(tuple(COORDINATES["water"]))


def test_cluster_magic_add():
    """
    Test for cluster magic method add
    """

    w = Cluster(COORDINATES["water"])
    w2 = w + w

    assert isinstance(w2, Cluster)
    assert w2.atoms == 2 * COORDINATES["water"]


def test_cluster_magic_mul_rmul():
    """
    Test for cluster magic method mul and rmul
    """

    w = Cluster(COORDINATES["water"])
    w2 = w * 2

    w6 = 3 * w2

    assert isinstance(w2, Cluster)
    assert isinstance(w6, Cluster)
    assert w2.atoms == 2 * COORDINATES["water"]
    assert w6.atoms == 6 * COORDINATES["water"]


def test_cluster_magic_rmul_fails():
    """
    Test for cluster magic method rmul fails
    """
    w = Cluster(COORDINATES["water"])
    with pytest.raises(ValueError):
        0 * w


def test_cluster_str():
    """
    Test for Cluster magic str
    """
    w = Cluster(
        {"atoms": [("H1", 0, 0, 0)]}, [("H2", 0, 0, 0), ("H3", 0, 0, 0)]
    )
    w_str = str(w)

    expected_str = (
        """Cluster of (2) molecules and (3) total atoms\n"""
        """ #0: molecule with 1 atoms:\n"""
        """     --> atoms: [('H1', 0, 0, 0)]\n"""
        """     --> charge: +0\n"""
        """     --> multiplicity: 1\n"""
        """ #1: molecule with 2 atoms:\n"""
        """     --> atoms: [('H2', 0, 0, 0), ('H3', 0, 0, 0)]\n"""
        """     --> charge: +0\n"""
        """     --> multiplicity: 1\n"""
    )

    assert w_str == expected_str


def test_cluster_dictionary():
    """
    Test for Cluster property dictionary
    """
    w = Cluster([("H", 0, 0, 0)])

    expected_dict = {0: Molecule([("H", 0, 0, 0)])}

    assert w.cluster_dictionary == expected_dict


def test_cluster_setting_frozen_molecule():
    """
    Test for Cluster setting frozen molecule from value
    """
    w = Cluster([("H1", 0, 0, 0), ("H2", 0, 0, 0)])

    w.frozen_molecule = 0

    assert w.frozen_molecule == [0]


def test_cluster_setting_frozen_molecule_list():
    """
    Test for Cluster setting frozen molecule from list
    """
    w = Cluster([("H1", 0, 0, 0), ("H2", 0, 0, 0)])

    w.frozen_molecule = [0, 1]

    assert w.frozen_molecule == [0, 1]


def test_cluster_setting_sphere_center():
    """
    Test for Cluster setting sphere center
    """
    w = Cluster([("H1", 0, 0, 0), ("H2", 0, 0, 0)])

    w.sphere_center = (-1.5, 0, 10)

    assert w.sphere_center == (-1.5, 0, 10)


def test_cluster_setting_sphere_center_fails():
    """
    Test for Cluster setting sphere center fails
    """
    w = Cluster([("H1", 0, 0, 0), ("H2", 0, 0, 0)])

    with pytest.raises(ValueError):
        w.sphere_center = (-1.5, 0)


def test_cluster_setting_sphere_radius():
    """
    Test for Cluster setting sphere radius
    """
    w = Cluster([("H1", 0, 0, 0), ("H2", 0, 0, 0)])

    w.sphere_radius = 20.2

    assert w.sphere_radius == 20.2


def test_cluster_setting_sphere_radius_fails():
    """
    Test for Cluster setting sphere radius fails
    """
    w = Cluster([("H1", 0, 0, 0), ("H2", 0, 0, 0)])

    with pytest.raises(ValueError):
        w.sphere_radius = 0.1


def test_cluster_overlap():
    """
    Testing cluster overlap check method
    """
    no_overlap = Cluster().overlapping((0, 0, 0), (1, 1, 1), 0.9)
    overlap = Cluster().overlapping((0, 0, 0), (1, 1, 1), 1.1)

    assert overlap
    assert not no_overlap  # double neg == True


def test_cluster_add_molecule():
    """
    Test cluster class add molecule
    """
    mol = Cluster(COORDINATES["water"])

    mol2 = mol.add_molecule(COORDINATES["water"])

    assert mol2.atoms == 2 * mol.atoms


@pytest.mark.parametrize(
    "system",
    [
        ("dummy"),
        ("hf"),
        ("water"),
    ],
)
def test_cluster_get_molecule(system):
    """
    Test retriving certain molecule/atom
    """
    mol = Cluster(
        COORDINATES[system],
        COORDINATES[system],
        COORDINATES[system],
    )

    # getting the first (index) molecule
    new_system = mol.get_molecule(0)

    assert new_system.atoms == COORDINATES[system]


def test_cluster_get_molecule_fails():
    """
    Test retriving certain molecule/atom
    """
    mol = Cluster(COORDINATES["water"])

    with pytest.raises(IndexError):
        mol.get_molecule(100)


def test_cluster_initialize_cluster():
    """
    Testing cluster class move (translate and rotate) molecules (if necessary)
    to avoid any atom overlapping
    """
    # a water cluster in which three molecules are at the same point
    w3 = Cluster(
        COORDINATES["water"],
        COORDINATES["water"],
        COORDINATES["water"],
        [("Xe", 10, 10, 10)],
        [("Z2", -10, -10, -10)],
    )

    w3_no_overlapping = w3.initialize_cluster

    overlap = []
    for i in w3_no_overlapping.cluster_dictionary:
        mol1 = w3_no_overlapping.cluster_dictionary[i]

        j = i + 1
        if j not in w3_no_overlapping.cluster_dictionary:
            continue

        mol2 = w3_no_overlapping.cluster_dictionary[j]

        overlap.append(
            Cluster().overlapping(mol1.coordinates, mol2.coordinates)
        )

    assert not all(overlap)


def test_cluster_move_molecule():
    """
    Testing Cluster class move_molecule method, rotating and translating
    a molecule avoiding overlaping
    """
    w5 = Cluster(
        COORDINATES["water"],
        COORDINATES["water"],
        COORDINATES["water"],
        COORDINATES["water"],
        COORDINATES["water"],
    )

    new_w5 = w5.move_molecule(1)
    new_w5 = w5.move_molecule(2)
    new_w5 = w5.move_molecule(3)
    new_w5 = w5.move_molecule(4)

    overlap = []
    for i in new_w5.cluster_dictionary:
        mol1 = new_w5.cluster_dictionary[i]

        j = i + 1
        if j not in new_w5.cluster_dictionary:
            continue

        mol2 = new_w5.cluster_dictionary[j]

        overlap.append(
            Cluster().overlapping(mol1.coordinates, mol2.coordinates)
        )

    assert not all(overlap)


def test_cluster_move_molecule_fails_too_close():
    """
    Testing Cluster class move_molecule fails
    """
    w2 = Cluster(
        COORDINATES["water"],
        COORDINATES["water"],
    )

    with pytest.raises(ValueError):
        w2.move_molecule(
            molecule=0,
            max_step=None,
            max_rotation=None,
            max_closeness=0.01,
        )


def test_cluster_move_molecule_fails_no_space():
    """
    Testing Cluster class move_molecule fails
    """
    w2 = Cluster(
        COORDINATES["water"],
        COORDINATES["water"],
    )

    w2.sphere_radius = 2

    with pytest.raises(AttributeError):
        w2.move_molecule(
            molecule=0,
            max_step=None,
            max_rotation=None,
            max_closeness=3,
        )


@pytest.mark.parametrize(
    "system",
    [
        ("dummy"),
        ("hf"),
        ("water"),
    ],
)
def test_cluster_remove_molecule(system):
    """
    Test deleting an existeing molecule/atom
    """
    mol = Cluster(
        COORDINATES[system],
        COORDINATES[system],
        COORDINATES[system],
    )

    # deleting the first (index) molecule
    new_system = mol.remove_molecule(0)

    assert (new_system.total_molecules + 1) == mol.total_molecules


def test_cluster_remove_molecule_fails():
    """
    Test removing certain molecule fails
    """
    mol = Cluster(COORDINATES["water"])

    with pytest.raises(IndexError):
        mol.remove_molecule(1)


@pytest.mark.parametrize(
    "angles, expected_rotation",
    [
        (  # single rot around x-axis "x=90, y=0, z=0",
            # where a1->a1, a2->a3, a3->-a2
            (90, 0, 0),
            [("a1", 1, 0, 0), ("a2", 0, 0, -1), ("a3", 0, 1, 0)],
        ),
        (  # single rot around y-axis "x=0, y=90, z=0",
            # where a1->a3, a2->a2, a3->-a1
            (0, 90, 0),
            [("a1", 0, 0, 1), ("a2", 0, 1, 0), ("a3", -1, 0, 0)],
        ),
        (  # single rot around z-axis "x=0, y=0, z=90",
            # where a1->-a2, a2->a1, a3->a3
            (0, 0, 90),
            [("a1", 0, -1, 0), ("a2", 1, 0, 0), ("a3", 0, 0, 1)],
        ),
        (  # MULTIPLE rot around "x=90, y=90, z=90",
            # where a1->a3, a2->a2, a3->-a1 (equivalent to y=90)
            (90, 90, 90),
            [("a1", 0, 0, 1), ("a2", 0, 1, 0), ("a3", -1, 0, 0)],
        ),
    ],
)
def test_cluster_rotate(angles, expected_rotation):
    """
    Test rotation
    """
    dummy = [("a1", 1, 0, 0), ("a2", 0, 1, 0), ("a3", 0, 0, 1)]

    ax, ay, az = angles

    mol = Cluster(dummy)
    mol_rotated = mol.rotate(0, x=ax, y=ay, z=az)

    expected_mol = Cluster(expected_rotation)

    assert np.allclose(
        mol_rotated.coordinates,
        expected_mol.coordinates,
        0.01,
    )


def test_cluster_rotate_frozen_molecule():
    """
    Test rotate a frozen molecule
    """
    dummy = [("a1", 1, 0, 0), ("a2", 0, 1, 0), ("a3", 0, 0, 1)]

    mol = Cluster(dummy, [("H", 0, 0, 0)])

    mol.frozen_molecule = 0
    mol_rotated = mol.rotate(0, x=90)

    assert mol.atoms == mol_rotated.atoms


def test_cluster_rotate_fails():
    """
    Test rotate fails molecule index does not exist
    """
    dummy = [("a1", 1, 0, 0), ("a2", 0, 1, 0), ("a3", 0, 0, 1)]

    mol = Cluster(dummy, [("H", 0, 0, 0)])

    with pytest.raises(IndexError):
        mol.rotate(-2, x=90)
    with pytest.raises(IndexError):
        mol.rotate(10, x=90)


def test_cluster_rotate_single_atom():
    """
    Test rotate a single atom
    """
    dummy = [("a1", 1, 0, 0), ("a2", 0, 1, 0), ("a3", 0, 0, 1)]

    mol = Cluster(dummy, [("H", 10, 0, 0)])

    mol_rotated = mol.rotate(1, x=90)

    assert mol.atoms == mol_rotated.atoms


@pytest.mark.parametrize(
    "steps, expected_translation",
    [
        (  # single translation on x-axis
            (10, 0, 0),
            [("a1", 11, 0, 0), ("a2", 10, 1, 0), ("a3", 10, 0, 1)],
        ),
        (  # single translation on y-axis
            (0, 5, 0),
            [("a1", 1, 5, 0), ("a2", 0, 6, 0), ("a3", 0, 5, 1)],
        ),
        (  # single translation on z-axis
            (0, 0, 2),
            [("a1", 1, 0, 2), ("a2", 0, 1, 2), ("a3", 0, 0, 3)],
        ),
        (  # MULTIPLE translation on xyz axes
            (10, 5, 2),
            [("a1", 11, 5, 2), ("a2", 10, 6, 2), ("a3", 10, 5, 3)],
        ),
    ],
)
def test_cluster_translate(steps, expected_translation):
    """
    Test translation
    """
    dummy = [("a1", 1, 0, 0), ("a2", 0, 1, 0), ("a3", 0, 0, 1)]

    tx, ty, tz = steps

    mol = Cluster(dummy)
    mol_translated = mol.translate(0, x=tx, y=ty, z=tz)

    expected_mol = Cluster(expected_translation)

    assert np.allclose(
        mol_translated.coordinates,
        expected_mol.coordinates,
        0.01,
    )


def test_cluster_translate_frozen_molecule():
    """
    Test translate a frozen molecule
    """
    dummy = [("a1", 1, 0, 0), ("a2", 0, 1, 0), ("a3", 0, 0, 1)]

    mol = Cluster(dummy, [("H", 0, 0, 0)])

    mol.frozen_molecule = 1
    mol_rotated = mol.translate(1, x=90)

    assert mol.atoms == mol_rotated.atoms


def test_cluster_translate_fails():
    """
    Test translate fails molecule index does not exist
    """
    dummy = [("a1", 1, 0, 0), ("a2", 0, 1, 0), ("a3", 0, 0, 1)]

    mol = Cluster(dummy, [("H", 0, 0, 0)])

    with pytest.raises(IndexError):
        mol.translate(-2, x=90)
    with pytest.raises(IndexError):
        mol.translate(10, x=90)
