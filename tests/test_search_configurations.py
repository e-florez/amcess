import pytest
from amcess.base_molecule import Cluster
from amcess.search_configuration import SearchConfig


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
    search_config = SearchConfig(Cluster(cluster1, cluster2))

    assert search_config._system_object.atoms == expected_coordinates


@pytest.mark.parametrize(
    "cluster1, cluster2, expected_coordinates",
    [
        (
            [("H", 0, 0, 0), ("F", 0.917, 0, 0)],
            [("H", 0, 0, 0), ("H", 0.74, 0, 0)],
            ["F", "H"],
        ),
    ],
)
def test_cluster_object_elements_into_search_conf(
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
    search_config = SearchConfig(Cluster(cluster1, cluster2))
    # it is necessary to sort to avoid errors because the elements attr
    # can return in different order the symbols
    assert (
        sorted(search_config._system_object.elements) == expected_coordinates
    )


@pytest.mark.parametrize(
    "cluster1, cluster2, expected_bases",
    [
        (
            [("H", 0, 0, 0), ("F", 0.917, 0, 0)],
            [("H", 0, 0, 0), ("H", 0.74, 0, 0)],
            "sto-3g",
        ),
    ],
)
def test_search_conf_basis_default(cluster1, cluster2, expected_bases):
    """[summary]
    Test: Basis set by default (sto-3g)
    Args:
        cluster1 ([dict]): dictionary with symbols and coordinates
                            of the first molecule
        cluster2 ([dict]): dictionary with symbols and coordinates
                            of the second molecule
        expected_basis ([string]): basis set
    """
    search_config = SearchConfig(Cluster(cluster1, cluster2))
    assert search_config._bases == expected_bases
