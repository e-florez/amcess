import pytest
from src.base_molecule import Cluster
from src.electronic_energy import hf_pyscf


@pytest.mark.parametrize(
    "cluster1, cluster2, expected_coordinates",
    [
        (
            [("H", 0, 0, 0), ("F", 0.917, 0, 0)],
            [("H", 0, 0, 0), ("H", 0.74, 0, 0)],
            [
                ("H", 0, 0, 0),
                ("F", 0.917, 0, 0),
                ("H", 0, 0, 0),
                ("H", 0.74, 0, 0),
            ],
        ),
    ],
)
def test_cluster_object_coordinates_into_search_conf(
    cluster1, cluster2, expected_coordinates
):
    """[summary]
    Test: Passed of object from Cluser to SearchConfig
    Args:
        cluster1 ([dict]): dictionary with symbols and coordinates
                            of the first molecule
        cluster2 ([dict]): dictionary with symbols and coordinates
                            of the second molecule
        expected_coordinates ([list]): expected coordinates plus symbols
    """
    search_config = hf_pyscf(Cluster(cluster1, cluster2))

    assert search_config._system_object.atoms == expected_coordinates
