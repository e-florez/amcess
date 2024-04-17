import numpy as np
import pytest

from amcess.cluster import Cluster
from amcess.cluster import Molecule

COORDINATES = {
    "dummy": [("H", 10, 20, 30), ("H", -0.5, 0, -10)],
    "hf": [("H", 0, 0, 0), ("F", 0.917, 0, 0)],
    "water": [
        ("O", 0, 0, 0),
        ("H", 0.58708, 0.75754, 0),
        ("H", -0.58708, 0.75754, 0),
    ],
}

#Missing
# amcess/cluster.py               234    11    95%
# lines missings:88, 190, 229, 440, 583-588, 690-696
# ===============================================================
# Cluster class
# ===============================================================
def test_init_():
    """
    Test for cluster init from list ype
    """
    w = Cluster(COORDINATES["water"])

    assert isinstance(w, Cluster)
    assert w.GetMolList() == COORDINATES["water"]


def test_init_dict():
    """
    Test for cluster init from dictionary type
    """
    water_dict = {"atoms": COORDINATES["water"], "charge": 0, "multiplicty": 1}

    w = Cluster(water_dict)

    assert isinstance(w, Cluster)
    assert w.GetMolList() == COORDINATES["water"]


def test_init_molecule():
    """
    Test for cluster init from Molecule class
    """
    water_dict = {"atoms": COORDINATES["water"]}

    w = Cluster(Molecule(water_dict))

    assert isinstance(w, Cluster)
    assert w.GetMolList() == COORDINATES["water"]


def test_init_cluster():
    """
    Test for cluster init from Cluster class
    """
    w = Cluster(COORDINATES["water"])
    w2 = Cluster(w, w)

    assert isinstance(w2, Cluster)
    assert w2.GetMolList() == 2 * COORDINATES["water"]


def test_init_fail():
    """
    Test for cluster init raise error
    """

    with pytest.raises(TypeError):
        Cluster(tuple(COORDINATES["water"]))


def test_magic_add():
    """
    Test for cluster magic method add
    """

    w = Cluster(COORDINATES["water"])
    w2 = w + w

    assert isinstance(w2, Cluster)
    assert w2.GetMolList() == 2 * COORDINATES["water"]


def test_magic_mul_rmul():
    """
    Test for cluster magic method mul and rmul
    """

    w = Cluster(COORDINATES["water"])
    w2 = w * 2

    w6 = 3 * w2

    assert isinstance(w2, Cluster)
    assert isinstance(w6, Cluster)
    assert w2.GetMolList() == 2 * COORDINATES["water"]
    assert w6.GetMolList() == 6 * COORDINATES["water"]


def test_magic_rmul_fails():
    """
    Test for cluster magic method rmul fails
    """
    w = Cluster(COORDINATES["water"])
    with pytest.raises(ValueError):
        0 * w


def test_str():
    """
    Test for Cluster magic str
    """
    w = Cluster(
        {"atoms": [("H", 0, 0, 0)]}, [("H", 0, 0, 0), ("H", 0, 0, 0)]
    )
    w_str = str(w)

    expected_str = (
        """Cluster of (2) molecules and (3) total atoms\n"""
        """ #0: molecule with 1 atoms:\n"""
        """     --> atoms: [('H', 0.0, 0.0, 0.0)]\n"""
        """     --> charge: +0\n"""
        """     --> multiplicity: 1\n"""
        """ #1: molecule with 2 atoms:\n"""
        """     --> atoms: [('H', 0.0, 0.0, 0.0), ('H', 0.0, 0.0, 0.0)]\n"""
        """     --> charge: +0\n"""
        """     --> multiplicity: 1\n"""
    )

    assert w_str == expected_str


def test_GetClusterDict():
    """
    Test for Cluster property dictionary
    """
    w = Cluster([("H", 0, 0, 0)])

    expected_dict = {0: Molecule([("H", 0, 0, 0)])}

    assert w.GetClusterDict() == expected_dict


def test_SetGetFreezeMol():
    """
    Test for Cluster setting frozen molecule from value
    """
    w = Cluster([("H", 0, 0, 0), ("H", 0, 0, 0)])

    w.SetFreezeMol(0)

    assert w.GetFreezeMol() == [0]


def test_SetGetFreezeMol_List():
    """
    Test for Cluster setting frozen molecule from list
    """
    w = Cluster([("H", 0, 0, 0), ("H", 0, 0, 0)])

    w.SetFreezeMol([0, 1])

    assert w.GetFreezeMol() == [0, 1]


def test_SetGetSphereCenter():
    """
    Test for Cluster setting sphere center
    """
    w = Cluster([("H", 0, 0, 0), ("H", 0, 0, 0)])

    w.SetSphereCenter((-1.5, 0, 10))

    assert w.GetSphereCenter() == (-1.5, 0, 10)


def test_SetSphereCenter_Fails():
    """
    Test for Cluster setting sphere center fails
    """
    w = Cluster([("H", 0, 0, 0), ("H", 0, 0, 0)])

    with pytest.raises(ValueError):
        w.SetSphereCenter((-1.5, 0))


def test_SetGetSphereR():
    """
    Test for Cluster setting sphere radius
    """
    w = Cluster([("H", 0, 0, 0), ("H", 0, 0, 0)])

    w.SetSphereR(20.2)

    assert w.GetSphereR() == 20.2


def test_SetSphereR_fails():
    """
    Test for Cluster setting sphere radius fails
    """
    w = Cluster([("H", 0, 0, 0), ("H", 0, 0, 0)])

    with pytest.raises(ValueError):
        w.SetSphereR(0.1)


def test_Overlap():
    """
    Testing cluster overlap check staticmethod
    """
    no_overlap = Cluster().Overlapping((0, 0, 0), (1, 1, 1), 0.9)
    overlap = Cluster().Overlapping((0, 0, 0), (1, 1, 1), 1.1)

    assert overlap
    assert not no_overlap


@pytest.mark.parametrize(
    "system",
    [
        ("dummy"),
        ("hf"),
        ("water"),
    ],
)
def test_GetMol(system):
    """
    Test retriving certain molecule/atom
    """
    mol = Cluster(
        COORDINATES[system],
        COORDINATES[system],
        COORDINATES[system],
    )

    # getting the first (index) molecule
    new_system = mol.GetMol(0)

    assert new_system.GetMolList() == COORDINATES[system]


def test_GetMol_Fails():
    """
    Test retriving certain molecule/atom, fail
    """
    mol = Cluster(COORDINATES["water"])

    with pytest.raises(KeyError):
        mol.GetMol(100)


def test_InitializeCluster():
    """
    Testing cluster class move (translate and rotate) molecules
    (if necessary) to avoid any atom overlapping
    """
    # a water cluster in which three molecules are at the same point
    w3 = Cluster(
        COORDINATES["water"],
        COORDINATES["water"],
        COORDINATES["water"],
        [("Xe", 10, 10, 10)],
        [("Zn", -10, -10, -10)],
    )

    w3_no_overlapping = w3.InitializeCluster()

    overlap = []
    for i in w3_no_overlapping.GetClusterDict():
        mol1 = w3_no_overlapping.GetClusterDict()[i]

        j = i + 1
        if j not in w3_no_overlapping.GetClusterDict():
            continue

        mol2 = w3_no_overlapping.GetClusterDict()[j]

        overlap.append(
            Cluster().Overlapping(mol1.GetAtomicCoordinates(),
                                  mol2.GetAtomicCoordinates())
        )

    assert not all(overlap)


def test_MoveMolecule():
    """
    Testing Cluster class MoveMol method, rotating
    and translating a molecule avoiding overlaping
    """
    w5 = Cluster(
        COORDINATES["water"],
        COORDINATES["water"],
        COORDINATES["water"],
        COORDINATES["water"],
        COORDINATES["water"],
    )

    new_w5 = w5.MoveMol(1)
    new_w5 = w5.MoveMol(2)
    new_w5 = w5.MoveMol(3)
    new_w5 = w5.MoveMol(4)

    overlap = []
    for i in new_w5.GetClusterDict():
        mol1 = new_w5.GetClusterDict()[i]

        j = i + 1
        if j not in new_w5.GetClusterDict():
            continue

        mol2 = new_w5.GetClusterDict()[j]

        overlap.append(
            Cluster().Overlapping(mol1.GetAtomicCoordinates(),
                                  mol2.GetAtomicCoordinates())
        )

    assert not all(overlap)


def test_MoveMolecule_Fails():
    """
    Testing Cluster class MoveMol fails
    """
    w2 = Cluster(
        COORDINATES["water"],
        COORDINATES["water"],
    )

    # ! Too Close
    with pytest.raises(ValueError):
        w2.MoveMol(
            molecule=0,
            max_step=None,
            max_rotation=None,
            max_closeness=0.01,
        )

    # ! No Space
    w2.SetSphereR(2)
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
def test_RemoveMol(system):
    """
    Test deleting an existeing molecule/atom
    """
    mol = Cluster(
        COORDINATES[system],
        COORDINATES[system],
        COORDINATES[system],
    )

    # deleting the first (index) molecule
    mol.RemoveMol(0)

    assert (mol.GetTotalMol() + 1) == 3


def test_RemoveMol_Fails():
    """
    Test removing certain molecule fails
    """
    mol = Cluster(COORDINATES["water"])

    with pytest.raises(IndexError):
        mol.RemoveMol(1)


@pytest.mark.parametrize(
    "angles, expected_rotation",
    [
        (  # single rot around x-axis "x=90, y=0, z=0",
            (90, 0, 0),
            [("H", 1, 0, 0.66666667), ("H", 0, 0, -0.33333333), ("H", 0, 1, 0.66666667)],
        ),
    ],
)
def test_RotateMol(angles, expected_rotation):
    """
    Test RotationMol
    """
    dummy = [("H", 1, 0, 0), ("H", 0, 1, 0), ("H", 0, 0, 1)]

    ax, ay, az = angles

    mol = Cluster(dummy)
    mol_rotated = mol.RotateMol(0, x=ax, y=ay, z=az)

    expected_mol = Cluster(expected_rotation)

    assert np.allclose(
        mol_rotated.GetAtomicCoordinates(),
        expected_mol.GetAtomicCoordinates(),
        0.01,
    )


def test_RotateMol_FreezeMol():
    """
    Test rotate a frozen molecule
    """
    dummy = [("H", 1, 0, 0), ("H", 0, 1, 0), ("H", 0, 0, 1)]

    mol = Cluster(dummy, [("H", 0, 0, 0)])

    mol.SetFreezeMol(0)
    mol_rotated = mol.RotateMol(0, x=90)

    print(mol.GetAtomicCoordinates()[3], mol_rotated.GetAtomicCoordinates()[3])

    assert (mol.GetAtomicCoordinates()[3].all() ==
            mol_rotated.GetAtomicCoordinates()[3].all())


def test_RotateMol_Fails():
    """
    Test RotateMol funtions, fails by molecule index does not exist
    """
    dummy = [("H", 1, 0, 0), ("H", 0, 1, 0), ("H", 0, 0, 1)]

    mol = Cluster(dummy, [("H", 0, 0, 0)])

    with pytest.raises(IndexError):
        mol.RotateMol(-2, x=90)
    with pytest.raises(IndexError):
        mol.RotateMol(10, x=90)


def test_RotateMol_single_atom():
    """
    Test RotateMol a single atom. When is selected a molecule
    composed by one atom then that doesn't rotate
    """
    dummy = [("H", 1, 0, 0), ("H", 0, 1, 0), ("H", 0, 0, 1)]

    mol = Cluster(dummy, [("H", 10, 0, 0)])

    mol_rotated = mol.RotateMol(1, x=90)

    assert (mol.GetAtomicCoordinates().all() ==
            mol_rotated.GetAtomicCoordinates().all())


@pytest.mark.parametrize(
    "steps, expected_translation",
    [
        (  # single translation on x-axis
            (10, 0, 0),
            [("H", 11, 0, 0), ("H", 10, 1, 0), ("H", 10, 0, 1)],
        ),
        (  # single translation on y-axis
            (0, 5, 0),
            [("H", 1, 5, 0), ("H", 0, 6, 0), ("H", 0, 5, 1)],
        ),
        (  # single translation on z-axis
            (0, 0, 2),
            [("H", 1, 0, 2), ("H", 0, 1, 2), ("H", 0, 0, 3)],
        ),
        (  # MULTIPLE translation on xyz axes
            (10, 5, 2),
            [("H", 11, 5, 2), ("H", 10, 6, 2), ("H", 10, 5, 3)],
        ),
    ],
)
def test_TranslateMol(steps, expected_translation):
    """
    Test TranslationMol function
    """
    dummy = [("H", 1, 0, 0), ("H", 0, 1, 0), ("H", 0, 0, 1)]

    tx, ty, tz = steps

    mol = Cluster(dummy)
    mol_translated = mol.TranslateMol(0, x=tx, y=ty, z=tz)

    expected_mol = Cluster(expected_translation)

    assert np.allclose(
        mol_translated.GetAtomicCoordinates(),
        expected_mol.GetAtomicCoordinates(),
        0.01,
    )


def test_TranslateMol_FreezeMol():
    """
    Test TranslateMol funciotn with frozen molecule
    """
    dummy = [("H", 1, 0, 0), ("H", 0, 1, 0), ("H", 0, 0, 1)]

    mol = Cluster(dummy, [("H", 0, 0, 0)])

    mol.SetFreezeMol(1)
    mol_rotated = mol.TranslateMol(1, x=90)

    assert (mol.GetAtomicCoordinates().all() ==
            mol_rotated.GetAtomicCoordinates().all())


def test_TranslateMol_Fails():
    """
    Test TranslateMol function fails, molecule index does not exist
    """
    dummy = [("H", 1, 0, 0), ("H", 0, 1, 0), ("H", 0, 0, 1)]

    mol = Cluster(dummy, [("H", 0, 0, 0)])

    with pytest.raises(IndexError):
        mol.TranslateMol(-2, x=90)
    with pytest.raises(IndexError):
        mol.TranslateMol(10, x=90)


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
def test_CalCenterRSphere(molecule1, molecule2, molecule3):
    """
    Test of CalCenterRSphere method, re--sort molecule
    """
    obj_cluster = Cluster(molecule1, molecule2, molecule3)
    obj_cluster = obj_cluster.CalCentRSphere()
    assert obj_cluster.GetMol(0).GetMolList() == molecule2


@pytest.mark.parametrize(
    "molecule1, molecule2",
    [
        (
            [("H", 0.0, 0.0, 1.0), ("H", 0.78, 0.0, 1.0)],
            [("H", 0.0, 0.0, 3.0), ("H", 0.78, 0.0, 3.0)],
        ),
    ],
)
def test_CalCenterRSphere_TE(molecule1, molecule2):
    """
    Test of TE of the tolerance argument in the CalCenterRSphere method
    """
    with pytest.raises(TypeError) as e:
        obj_cluster = Cluster(molecule1, molecule2)
        obj_cluster.CalCentRSphere(add_tolerance_radius=[1])
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
def test_SetSeed(molecule1, molecule2, new_seed, expected):
    """
    Test for the setter seed
    """
    obj_cluster = Cluster(molecule1, molecule2)
    obj_cluster.SetSeed(new_seed)
    assert obj_cluster.GetRandomGen().uniform() == expected


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
def test_GetSeed(molecule1, molecule2, seed):
    """
    Seed as argument in the instatiation of the Cluster class
    """
    obj_cluster = Cluster(molecule1, molecule2, seed=seed)
    assert obj_cluster.GetSeed() == seed
