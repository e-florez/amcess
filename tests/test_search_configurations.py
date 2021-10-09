import pytest
import sys

sys.path.append("../src")

from src.base_molecule import Cluster
from src.search_configuration import SearchConfig


@pytest.mark.parametrize(
    "cluster1, cluster2, expected_coordinates",
    [
        (
            {"atoms": [("H", 0, 0, 0), ("F", 0.917, 0, 0)]},
            {"atoms": [("H", 0, 0, 0), ("H", 0.74, 0, 0)]},
            [("H", 0, 0, 0), ("F", 0.917, 0, 0), ("H", 0, 0, 0), ("H", 0.74, 0, 0)],
        ),
    ],
)
def test_cluster_object_into_search_conf(cluster1, cluster2, expected_coordinates):
    """[summary]
    Test: Passed of object from Cluser to SearchConfig
    Args:
        cluster1 ([dict]): dictionary with symbols and coordinates of the first molecule
        cluster2 ([dict]): dictionary with symbols and coordinates of the second molecule
        expected_coordinates ([list]): expected coordinates plus symbols
    """
    search_config = SearchConfig(Cluster(cluster1, cluster2))
    assert search_config._system_object.coordinates == expected_coordinates


@pytest.mark.parametrize(
    "cluster1, cluster2, expected_basis",
    [
        (
            {"atoms": [("H", 0, 0, 0), ("F", 0.917, 0, 0)]},
            {"atoms": [("H", 0, 0, 0), ("H", 0.74, 0, 0)]},
            "sto-3g",
        ),
    ],
)
def test_cluster_object_into_search_conf(cluster1, cluster2, expected_basis):
    """[summary]
    Test: Basis set by default (sto-3g)
    Args:
        cluster1 ([dict]): dictionary with symbols and coordinates of the first molecule
        cluster2 ([dict]): dictionary with symbols and coordinates of the second molecule
        expected_basis ([string]): basis set
    """
    search_config = SearchConfig(Cluster(cluster1, cluster2))
    assert search_config._basis == expected_basis
