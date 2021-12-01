import pytest
from amcess.base_molecule import Cluster
from amcess.local_minima_2 import LocalMinima, extra_functions
import numpy as np


@pytest.mark.parametrize(
    "cluster1,cluster2, expected_coordinates",
    [
        (
            [("H", 0, 0, 0), ("H", 1, 0, 0)],
            [("H", 2, 0, 0), ("F", 3, 0, 0)],
            [(2, 0, 0), (3, 0, 0), (0, 0, 0), (1, 0, 0)],
        ),
    ],
)
def test_center_sphere_mass(cluster1, cluster2, expected_coordinates):
    """testing center_sphere_mass function"""
    cluster = Cluster(cluster1, cluster2)

    cluster_moved = extra_functions(initial_cluster=cluster).center_sphere_mass
    assert np.allclose(cluster_moved.coordinates, expected_coordinates, 0.001)


@pytest.mark.parametrize(
    "cluster1,cluster2, expected_coordinates",
    [
        (
            [("H", 0, 0, 0), ("H", 1, 0, 0)],
            [("H", 2, 0, 0), ("O", 3, 0, 0), ("H", 4, 0, 0)],
            [(2, 0, 0), (3, 0, 0), (4, 0, 0), (0, 0, 0), (1, 0, 0)],
        ),
    ],
)
def test_center_sphere_atoms(cluster1, cluster2, expected_coordinates):
    """testing center_sphere_atoms function"""
    cluster = Cluster(cluster1, cluster2)

    cluster_moved = extra_functions(
        initial_cluster=cluster
    ).center_sphere_atoms
    assert np.allclose(cluster_moved.coordinates, expected_coordinates, 0.001)


@pytest.mark.parametrize(
    "cluster1,cluster2, expected_bounds",
    [
        (
            [("H", 0, 0, 0), ("H", 1, 0, 0)],
            [("H", 2, 0, 0), ("H", 3, 0, 0)],
            [
                (-4.0, 4.0),
                (-4.0, 4.0),
                (-4.0, 4.0),
                (-180, 180),
                (-180, 180),
                (-180, 180),
            ],
        ),
    ],
)
def test_bounds(cluster1, cluster2, expected_bounds):
    """testing bounds"""
    cluster = Cluster(cluster1, cluster2)

    bounds = LocalMinima(initial_cluster=cluster).bounds
    assert np.allclose(bounds, expected_bounds, 0.001)


@pytest.mark.parametrize(
    "cluster1,cluster2, expected_energy",
    [
        (
            [("H", 0, 0, 0), ("H", 1, 0, 0)],
            [("H", 2, 0, 0), ("H", 3, 0, 0)],
            [
                (-2.132227),
            ],
        ),
    ],
)
def test_shgo_optimize(cluster1, cluster2, expected_energy):
    """testing shgo optimization"""
    cluster = Cluster(cluster1, cluster2)

    minima = LocalMinima(initial_cluster=cluster).shgo_optimize()
    assert np.allclose(minima.fun, expected_energy, 0.1)


@pytest.mark.parametrize(
    "cluster1,cluster2, expected_energy",
    [
        (
            [("H", 0, 0, 0), ("H", 1, 0, 0)],
            [("H", 2, 0, 0), ("H", 3, 0, 0)],
            [
                (-2.132227),
            ],
        ),
    ],
)
def test_dual_annealing_optimize(cluster1, cluster2, expected_energy):
    """testing dual_annealing optimization"""
    cluster = Cluster(cluster1, cluster2)

    minima = LocalMinima(initial_cluster=cluster).dual_annealing_optimize()
    assert np.allclose(minima.fun, expected_energy, 0.1)


@pytest.mark.parametrize(
    "cluster1,cluster2, expected_energy",
    [
        (
            [("H", 0, 0, 0), ("H", 1, 0, 0)],
            [("H", 2, 0, 0), ("H", 3, 0, 0)],
            [
                (-2.132227),
            ],
        ),
    ],
)
def test_gaussian_optimize(cluster1, cluster2, expected_energy):
    """testing dual_annealing optimization"""
    cluster = Cluster(cluster1, cluster2)

    minima = LocalMinima(initial_cluster=cluster).gaussian_optimize()
    assert np.allclose(minima.fun, expected_energy, 0.1)
