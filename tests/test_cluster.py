import numpy as np
import pytest

from amcess.cluster import Cluster

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

    w.freeze_molecule = 0

    assert w.freeze_molecule == [0]


def test_cluster_setting_frozen_molecule_list():
    """
    Test for Cluster setting frozen molecule from list
    """
    w = Cluster([("H1", 0, 0, 0), ("H2", 0, 0, 0)])

    w.freeze_molecule = [0, 1]

    assert w.freeze_molecule == [0, 1]


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

    w3_no_overlapping = w3.initialize_cluster()

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

    mol.freeze_molecule = 0
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

    mol.freeze_molecule = 1
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


@pytest.mark.parametrize(
    "molecule1, molecule2, molecule3",
    [
        (
            [("H", 0.0, 0.0, 0.0)],
            [("H", 0.0, 0.0, 1.0), ("H", 0.78, 0.0, 1.0)],
            [("H", 0.0, 0.0, 3.0), ("H", 0.78, 0.0, 3.0)],
        ),
    ],
)
def test_the_biggest_to_initio(molecule1, molecule2, molecule3):
    """
    Test of method spherical_contour_cluster, re--sort molecule
    """
    obj_cluster = Cluster(molecule1, molecule2, molecule3)
    obj_cluster = obj_cluster.center_radius_sphere()
    assert obj_cluster.get_molecule(0).atoms == molecule2


@pytest.mark.parametrize(
    "molecule1, molecule2",
    [
        (
            [("H", 0.0, 0.0, 1.0), ("H", 0.78, 0.0, 1.0)],
            [("H", 0.0, 0.0, 3.0), ("H", 0.78, 0.0, 3.0)],
        ),
    ],
)
def test_TE_center_radius_sphere(molecule1, molecule2):
    """
    Test of TE of the tolerance argument in the method center_radius_sphere
    """
    with pytest.raises(TypeError) as e:
        obj_cluster = Cluster(molecule1, molecule2)
        obj_cluster.center_radius_sphere(add_tolerance_radius=[1])
    assert (
        str(e.value) == "\n\nThe tolerance for radius is not a float"
        f"\nplease, check: '{type([1])}'\n"
    )


@pytest.mark.parametrize(
    "molecule1, molecule2, new_seed, expected",
    [
        (
            [("H", 0.0, 0.0, 1.0), ("H", 0.78, 0.0, 1.0)],
            [("H", 0.0, 0.0, 3.0), ("H", 0.78, 0.0, 3.0)],
            2,
            0.2616121342493164,
        ),
    ],
)
def test_Cluster_seed_set(molecule1, molecule2, new_seed, expected):
    """
    Test for the setter seed
    """
    obj_cluster = Cluster(molecule1, molecule2)
    obj_cluster.seed = new_seed
    assert obj_cluster._random_gen.uniform() == expected


@pytest.mark.parametrize(
    "molecule1, molecule2, seed",
    [
        (
            [("H", 0.0, 0.0, 1.0), ("H", 0.78, 0.0, 1.0)],
            [("H", 0.0, 0.0, 3.0), ("H", 0.78, 0.0, 3.0)],
            2,
        ),
    ],
)
def test_Cluster_init_seed(molecule1, molecule2, seed):
    """
    Seed as argument in the instatiation of the Cluster class
    """
    obj_cluster = Cluster(molecule1, molecule2, seed=seed)
    assert obj_cluster._seed == seed
