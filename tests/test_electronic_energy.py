import pytest
from amcess.base_molecule import Cluster
from amcess.electronic_energy import Electronic_energy
import numpy as np


@pytest.mark.parametrize(
    "cluster1,cluster2, x, expected_coordinates",
    [
        (
            [("H", 0.58708, 0.75754, 0), ("H", -0.58708, 0.75754, 0)],
            [("H", 0.58708, 0.75754, 1), ("H", -0.58708, 0.75754, 1)],
            [0.5, 0.5, 0.5, 10, 10, 10],
            [
                (0.58708000, 0.75754000, 0.00000000),
                (-0.58708000, 0.75754000, 0.00000000),
                (1.06937737, 1.17457709, 1.61657397),
                (-0.06937737, 1.34050291, 1.38342603),
            ],
        ),
    ],
)
def test_controled_move(cluster1, cluster2, x, expected_coordinates):
    """testing controled_move function"""
    cluster = Cluster(cluster1, cluster2)

    c = Electronic_energy(cluster)
    cluster_moved = c.controled_move(x)
    assert np.allclose(cluster_moved.coordinates, expected_coordinates, 0.001)


@pytest.mark.parametrize(
    "cluster1, cluster2, x, expected_energy",
    [
        (
            [("H", 0.58708, 0.75754, 0), ("H", -0.58708, 0.75754, 0)],
            [("H", 0.58708, 0.75754, 1), ("H", -0.58708, 0.75754, 1)],
            [0.5, 0.5, 0.5, 10, 10, 10],
            [-1.9068764625565566],
        )
    ],
)
def test_energy(cluster1, cluster2, expected_energy, x):

    cluster = Cluster(cluster1, cluster2)
    c = Electronic_energy(cluster)
    e = c.energy(x)
    assert e - expected_energy < 1.0e-7
