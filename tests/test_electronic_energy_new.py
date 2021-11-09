import pytest
from amcess.base_molecule import Cluster
from amcess.electronic_energy_new import Electronic_energy


@pytest.mark.parametrize(
    "cluster1,cluster2, displ, angular_rotations, expected_coordinates", 
    [
        (
            [("H", 0.58708, 0.75754, 0), ("H", -0.58708, 0.75754, 0)],
            [("H", 0.58708, 0.75754, 1), ("H", -0.58708, 0.75754, 1)],
            [[0.5, 0.5, 0.5]],
            [[15, 15, 15]],
            [(1.047753097026884, 1.1487568712496972, 0.6810958364975624), 
             (-0.04775309702688413, 1.366323128750303, 0.31890416350243755), 
             (0.58708, 0.75754, 1), (-0.58708, 0.75754, 1)]
        ), 
    ]
)
def test_controled_move(cluster1, cluster2, displ, angular_rotations, 
                        expected_coordinates):

    cluster = Cluster(cluster1, cluster2)
    c = Electronic_energy(cluster)
    mole_moved = c.controled_move(displ, angular_rotations)

#    cluster.move(displacements, angular_rotations)
    assert mole_moved.coordinates == expected_coordinates


@pytest.mark.parametrize(
    "cluster1, cluster2, x, expected_energy",
    [
        (
            [("H", 0.58708, 0.75754, 0), ("H", -0.58708, 0.75754, 0)],
            [("H", 0.58708, 0.75754, 1), ("H", -0.58708, 0.75754, 1)],
            [0.5, 0.5, 0.5, 10, 10, 10],
            [-1.9068764625565566]
        )
    ]
)
def test_energy(cluster1, cluster2, expected_energy, x):

    cluster = Cluster(cluster1, cluster2)
    c = Electronic_energy(cluster)
    e = c.energy(x)
    assert e - expected_energy < 1.0e-7