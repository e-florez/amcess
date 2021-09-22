import pytest
from .context import src
from src.base_molecule import Cluster


@pytest.mark.parametrize(
    "input_coordinates, expected_total_mass",
    [
        ({"coordinates": [("H", 0, 0, 0), ("F", 0.917, 0, 0)]}, 20.0),
        ({"coordinates": [("H", 0, 0, 0), ("H", 0.74, 0, 0)]}, 2.0),
        # ({"coordinates": []}, 18.02),
    ],
)
def test_cluster_total_mass(input_coordinates, expected_total_mass):
    """
    Test
    """
    mol = Cluster(input_coordinates)
    molecular_mass = mol.total_mass

    # 2.016 - 2 = +0.016 > 0.001

    assert abs(molecular_mass - expected_total_mass) < 1.0e-1


@pytest.mark.parametrize(
    "input_coordinates, expected_total_atoms",
    [
        ({"coordinates": [("H", 0, 0, 0), ("F", 0.917, 0, 0)]}, 2),
        ({"coordinates": [("H", 0, 0, 0), ("H", 0.74, 0, 0)]}, 2),
        # ({"coordinates": []}, 3),
    ],
)
def test_cluster_total_atoms(input_coordinates, expected_total_atoms):
    """
    Test
    """
    mol = Cluster(input_coordinates)
    total_atoms = mol.total_atoms

    assert abs(total_atoms - expected_total_atoms) < 1.0e-10


if __name__ == "__main__":
    pytest.main(["tests"])
