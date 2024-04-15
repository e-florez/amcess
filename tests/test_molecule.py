import numpy as np
import pytest

from amcess.molecule import Molecule

COORDINATES = {
    "dummy": [("H", 10, 20, 30), ("H", -0.5, 0, -10)],
    "hf": [("H", 0, 0, 0), ("F", 0.917, 0, 0)],
    "water": [
        ("O", 0, 0, 0),
        ("H", 0.58708, 0.75754, 0),
        ("H", -0.58708, 0.75754, 0),
    ],
}


# ===============================================================
# Molecule class
# ===============================================================
def test_molecule_class_input():
    """
    Test inputs of molecule class when is used the list format
    for molecular information
    """
    mol = Molecule([("H", 10, 20, 30), ("H", 0, 0, 0)], -1, 10)
    assert mol.GetMolCharge() == -1
    assert mol.GetMolMultiplicity() == 10


def test_molecule_class_atomic_symbols():
    """
    Test atomic symbols, this is did by rdkit
    """
    with pytest.raises(RuntimeError):
        Molecule([("A", 10, 20, "")], -1, 10)
    with pytest.raises(RuntimeError):
        Molecule([("", 10, 20, 0)], -1, 10)
    with pytest.raises(RuntimeError):
        Molecule([("@", 10, 20, 0)], -1, 10)


def test_molecule_class_charge_value():
    """
    Test to check the charge molecular
    """
    with pytest.raises(ValueError):
        Molecule([("H", 10, 20, 0)], 0.4, 10)


def test_molecule_class_multiplicity_value():
    """
    Test class init
    """
    with pytest.raises(ValueError):
        Molecule([("H", 10, 20, 0)], -2, 0)


def test_molecule_init_from_dict():
    """
    Test class init
    """
    mol = Molecule(
        {
            "atoms": [("H", 10, 20, 0), ("H", 0, 0, 0)],
            "charge": -2,
            "multiplicity": 10,
        }
    )

    assert mol.GetMolList() == [("H", 10, 20, 0), ("H", 0, 0, 0)]
    assert mol.GetNumAtoms() == 2
    assert mol.GetMolCharge() == -2
    assert mol.GetMolMultiplicity() == 10


def test_molecule_init_from_dict_wrong():
    """
    Test atomic symbols, this is did by rdkit
    """
    with pytest.raises(KeyError):
        Molecule(
            {
                ".......": [("A", 10, 20, 0), ("B", 0, 0, 0)],
                "charge": -2,
                "multiplicity": 10,
            }
        )


# #def test_molecule_magic_add():
# #    """Testing add molecule method (__add__)"""
# #    mol1 = Molecule([("A", 0, 0, 0), ("B", 1, 1, 1)])
# #    mol2 = Molecule([("X", 0, 0, 0), ("Y", 1, 1, 1)])
# #
# #    mol = Molecule([("A", 0, 0, 0), ("B", 1, 1, 1)]).add_molecule(mol2)
# #    new_mol = mol1 + mol2

# #    assert isinstance(new_mol, Cluster)
# #    assert isinstance(mol, Cluster)
# #    assert mol.atoms == new_mol.atoms
# #    assert new_mol.total_atoms == mol1.total_atoms + mol2.total_atoms
# #    assert new_mol.symbols == mol1.symbols + mol2.symbols
# #    assert new_mol.atoms == mol1.atoms + mol2.atoms


# #def test_molecule_magic_add_fail():
# #    """Testing add molecule method (__add__)"""
# #    mol = Molecule([("A", 0, 0, 0), ("B", 1, 1, 1)])
# #    with pytest.raises(TypeError):
# #        mol.add_molecule(Atom("H", 0, 0, 0))

# #    with pytest.raises(TypeError):
# #        mol.add_molecule("H", 0, 0, 0)


def test_molecule_magic_mul_rmul():
    """Testing magic mul function"""
    mol1 = Molecule([("H", 0, 0, 0), ("H", 1, 1, 1)])
    mol2 = mol1 * 2
    mol3 = 3 * mol1
    assert mol2.GetNumAtoms() == 4
    assert mol3.GetNumAtoms() == 6


def test_molecule_magic_mul_rmul_fail():
    """Testing magic mul function error"""
    mol1 = Molecule([("H", 0, 0, 0), ("H", 1, 1, 1)])

    with pytest.raises(ValueError):
        mol1 * 0
    with pytest.raises(ValueError):
        mol1 * -3


@pytest.mark.parametrize(
    "system, expected_result",
    [
        (
            [("Xe", 0, 0, 0)],
            """1\ncharge: 0 multiplicity: 1\nXe    """
            """\t     0.00000000\t     0.00000000\t     0.00000000\n""",
        )
    ],
)
def test_molecule_magic_str(system, expected_result):
    """
    Test for list of magic funtion str
    """
    str_mol = Molecule(system)
    assert str(str_mol) == str(str_mol.GetBlockXYZ())
    assert str(str_mol) == expected_result


@pytest.mark.parametrize(
    "system, expected_atoms",
    [
        ("dummy", [("H", 10, 20, 30), ("H", -0.5, 0, -10)]),
    ],
)
def test_molecule_GetMolList(system, expected_atoms):
    """
    Test GetMolList() of Molecule Class
    """
    atoms_list = Molecule(COORDINATES[system]).GetMolList()

    assert atoms_list == expected_atoms


@pytest.mark.parametrize(
    "system, expected_atoms",
    [
        ("dummy", 
         {"atoms": [("H", 10, 20, 30), ("H", -0.5, 0, -10)],
          "charge": 0, "multiplicity": 1}),
    ],
)
def test_molecule_GetMolDict(system, expected_atoms):
    """
    Test for coordinate in XYZ format
    """
    mol = Molecule(COORDINATES[system])
    assert mol.GetMolDict() == expected_atoms


@pytest.mark.parametrize(
    "system, expected_masses",
    [
        ("dummy", [1, 1]),
        ("hf", [1, 19]),
        ("water", [16, 1, 1]),
    ],
)
def test_molecule_GetAtomicMases(system, expected_masses):
    """
    Test for atomic mases
    """
    masses = Molecule(COORDINATES[system]).GetAtomicMasses()

    assert np.allclose(masses, expected_masses, 0.1)


@pytest.mark.parametrize(
    "system, expected_charge, expected_multiplicity",
    [
        ("dummy", -2, 20),
    ],
)
def test_molecule_SetMolCharge_SetMolMultiplicity(
    system, expected_charge, expected_multiplicity
):
    """
    Tests of Set atomic charge and multiplicity
    """
    mol = Molecule(COORDINATES[system])

    mol.SetMolCharge(-2)
    mol.SetMolMultiplicity(20)

    assert mol.GetMolCharge() == expected_charge
    assert mol.GetMolMultiplicity() == expected_multiplicity


@pytest.mark.parametrize(
    "system, expected_coordinates",
    [
        ("dummy", np.array([[10, 20, 30], [-0.5, 0, -10]])),
    ],
)
def test_molecule_GetAtomicCoordinates(system, expected_coordinates):
    """
    Test for coordinate in XYZ format
    """
    coordinates = Molecule(COORDINATES[system]).GetAtomicCoordinates()

    assert coordinates.all() == expected_coordinates.all()


@pytest.mark.parametrize(
    "system, expected_elements",
    [
        ("dummy", ["H", "H"]),
        ("hf", ["H", "F"]),
        ("water", ["O", "H"]),
    ],
)
def test_molecule_GetAtomicSymbols(system, expected_elements):
    """
    Test for list of uniques elements
    """
    elements = Molecule(COORDINATES[system]).GetAtomicSymbols()

    assert elements.sort() == expected_elements.sort()


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
def test_molecule_GetNumberingAtoms(system, expected_result):
    """
    Test for list of atomic symbols
    """
    numbered_atoms = Molecule(system).GetNumberingAtoms()

    assert numbered_atoms == expected_result


# @pytest.mark.parametrize(
#     "system, expected_symbols",
#     [
#         ("dummy", ["A", "X0"]),
#         ("hf", ["H", "F"]),
#         ("water", ["O", "H", "H"]),
#     ],
# )
# def test_molecule_atomic_symbols(system, expected_symbols):
#     """
#     Test for list of atomic symbols
#     """
#     symbols = Molecule(COORDINATES[system]).symbols

#     assert symbols == expected_symbols


# @pytest.mark.parametrize(
#     "system, expected_atoms",
#     [
#         ("dummy", 2),
#         ("hf", 2),
#         ("water", 3),
#     ],
# )
# def test_molecule_total_atoms(system, expected_atoms):
#     """
#     Test for the total number of atoms
#     """
#     number_of_atoms = Molecule(COORDINATES[system]).total_atoms

#     assert (number_of_atoms - expected_atoms) == 0


# @pytest.mark.parametrize(
#     "system, expected_mass",
#     [
#         ("dummy", 0),
#         ("hf", 20),
#         ("water", 18),
#     ],
# )
# def test_molecule_total_mass(system, expected_mass):
#     """
#     Test for total molecular mass
#     """
#     mass = Molecule(COORDINATES[system]).total_mass

#     assert abs(mass - expected_mass) < 1.0e-1


# @pytest.mark.parametrize(
#     "system, expected_cm",
#     [
#         ("dummy", (0, 0, 0)),
#         ("hf", (0.9, 0, 0)),
#         ("water", (0, 0, 0)),
#     ],
# )
# def test_molecule_center_of_mass(system, expected_cm):
#     """
#     Test center of mass
#     """
#     cm = Molecule(COORDINATES[system]).center_of_mass

#     assert np.linalg.norm(np.asarray(cm) - expected_cm) < 0.1


# @pytest.mark.parametrize(
#     "system, expected_principal_axes",
#     [
#         ("dummy", [(10, 20, 30), (-0.5, 0, -10)]),
#         ("hf", [(-0.9, 0, 0), (0.05, 0, 0)]),
#         ("water", [(0, -0.09, 0), (0.6, 0.7, 0), (-0.6, 0.7, 0)]),
#     ],
# )
# def test_molecule_principal_axes(system, expected_principal_axes):
#     """
#     Test principal axes using internal coordinates system.
#     Number of principal axis is equal to number of atoms
#     """
#     pp = Molecule(COORDINATES[system]).principal_axes

#     assert np.allclose(pp, expected_principal_axes, 0.1)


# @pytest.mark.parametrize(
#     "system, expected_result",
#     [
#         (
#             [("Xe", 0, 0, 0)],
#             """\t1\n-- charge=0 and multiplicity=1 --\nXe    """
#             """\t     0.00000000\t     0.00000000\t     0.00000000\n""",
#         )
#     ],
# )
# def test_molecule_xyz_format(system, expected_result):
#     """
#     Test for list of atomic symbols
#     """
#     xyz_file = Molecule(system).xyz
#     assert str(xyz_file) == expected_result


# @pytest.mark.parametrize(
#     "system, new_atom",
#     [
#         ("dummy", [("H", 0, 0, 0)]),
#         ("hf", [("H", 0, 0, 0)]),
#         ("water", [("H", 0, 0, 0)]),
#     ],
# )
# def test_molecule_add_atom(system, new_atom):
#     """
#     Test adding a new atom
#     """
#     mol = Molecule(COORDINATES[system])

#     # duplicating the system
#     new_system = mol.add_atoms(new_atom)

#     assert new_system.total_atoms == (mol.total_atoms + 1)


# def test_molecule_add_atom_fail():
#     """
#     Test adding a new atom only as a lits
#     """
#     mol = Molecule([("H", 0, 0, 0)])
#     with pytest.raises(TypeError):
#         mol.add_atoms(("H", 0, 0, 0))
#     with pytest.raises(TypeError):
#         mol.add_atoms("H", 0, 0, 0)


# @pytest.mark.parametrize(
#     "system",
#     [
#         ("dummy"),
#         ("hf"),
#         ("water"),
#     ],
# )
# def test_molecule_get_atom(system):
#     """
#     Test retriving certain atom
#     """
#     mol = Molecule(COORDINATES[system])

#     # getting the first (index) molecule
#     new_system = mol.get_atom(0)

#     assert isinstance(new_system, tuple)
#     assert len(new_system) == 4


# def test_molecule_get_atom_fail():
#     """
#     Test retriving certain atom
#     """
#     mol = Molecule(COORDINATES["water"])

#     with pytest.raises(IndexError):
#         mol.get_atom(100)
#     with pytest.raises(IndexError):
#         mol.get_atom("@")


# @pytest.mark.parametrize(
#     "system",
#     [
#         ("dummy"),
#         ("hf"),
#         ("water"),
#     ],
# )
# def test_molecule_remove_atom(system):
#     """
#     Test removing an atom
#     """
#     mol = Molecule(COORDINATES[system])

#     # duplicating the system
#     new_system = mol.remove_atom(0)

#     assert new_system.total_atoms == mol.total_atoms - 1


# def test_molecule_remove_atom_fail():
#     """
#     Test removing certain atom
#     """
#     mol = Molecule(COORDINATES["water"])

#     with pytest.raises(IndexError):
#         mol.remove_atom(100)
#     with pytest.raises(IndexError):
#         mol.remove_atom("@")
