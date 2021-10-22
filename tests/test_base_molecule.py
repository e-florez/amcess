import numpy as np
import pytest

from src.base_molecule import Cluster, Molecule

COORDINATES = {
    "dummy": [("A", 10, 20, 30), ("X0", -0.5, 0, -10)],
    "hf": [("H", 0, 0, 0), ("F", 0.917, 0, 0)],
    "water": [
        ("O", 0, 0, 0),
        ("H", 0.58708, 0.75754, 0),
        ("H", -0.58708, 0.75754, 0),
    ],
}


@pytest.mark.parametrize(
    "system, expected_atoms",
    [
        ("dummy", 2),
        ("hf", 2),
        ("water", 3),
    ],
)
def test_molecule_number_atoms(system, expected_atoms):
    """
    Test for number of atoms
    """
    number_of_atoms = Molecule(COORDINATES[system]).total_atoms

    assert (number_of_atoms - expected_atoms) == 0


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
    Test for atomic symbols
    """
    symbols = Molecule(COORDINATES[system]).symbols

    assert symbols == expected_symbols


@pytest.mark.parametrize(
    "system, expected_molecules",
    [
        ("dummy", 1),
        ("hf", 1),
        ("water", 1),
    ],
)
def test_molecule_number_molecules(system, expected_molecules):
    """
    Test for number of individual fragments (atoms/molecules)
    """
    fragments = Cluster(COORDINATES[system]).total_molecules

    assert (fragments - expected_molecules) == 0


@pytest.mark.parametrize(
    "system, expected_molecules",
    [
        (["water", "water"], 2),
        (["dummy", "hf", "water"], 3),
    ],
)
def test_molecule_number_molecules_more(system, expected_molecules):
    """
    Test for number of individual fragments (atoms/molecules)
    when more than one
    """
    new_system = [COORDINATES[s] for s in system]

    fragments = Cluster(*new_system).total_molecules

    assert (fragments - expected_molecules) == 0


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
    "system, expected_atoms",
    [
        ("dummy", [("A", 10, 20, 30), ("X0", -0.5, 0, -10)]),
    ],
)
def test_molecule_atoms(system, expected_atoms):
    """
    Test for coordinate in XYZ format
    """
    coordinates = Molecule(COORDINATES[system]).atoms

    assert coordinates == expected_atoms


@pytest.mark.parametrize(
    "system, expected_coordinates",
    [
        ("dummy", [(10, 20, 30), (-0.5, 0, -10)]),
    ],
)
def test_molecule_coordinates(system, expected_coordinates):
    """
    Test for cartesian coordinate, 3D coordinates, (x, y, z)
    """
    cartesian_coordinates = Molecule(COORDINATES[system]).coordinates

    assert cartesian_coordinates == expected_coordinates


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
    "system",
    [
        ("dummy"),
        ("hf"),
        ("water"),
    ],
)
def test_molecule_add(system):
    """
    Test adding a new molecule/atoms
    """
    mol = Cluster(COORDINATES[system])

    # duplicating the system
    new_system = mol.add_molecule(COORDINATES[system])

    assert new_system.total_atoms == (2 * mol.total_atoms)
    assert new_system.total_molecules == (2 * mol.total_molecules)


@pytest.mark.parametrize(
    "system",
    [
        ("dummy"),
        ("hf"),
        ("water"),
    ],
)
def test_molecule_delete(system):
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


@pytest.mark.parametrize(
    "system",
    [
        ("dummy"),
        ("hf"),
        ("water"),
    ],
)
def test_molecule_get(system):
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
def test_molecule_translate(steps, expected_translation):
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
def test_molecule_rotate(angles, expected_rotation):
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
