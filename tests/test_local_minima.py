import pytest
from amcess.base_molecule import Cluster
from amcess.local_minima import LocalMinima


@pytest.mark.parametrize(
    "cluster1,cluster2, expected_energy",
    [
        (
            [("O", 0, 0, 0), ("H", 0.58708, 0.75754, 0),
             ("H", -0.58708, 0.75754, 0)],
            [("O", 0, 0, 1), ("H", 0.58708, 0.75754, 1),
             ("H", -0.58708, 0.75754, 1)],
            [148],
        )
    ],
)
def test_local_minima(cluster1, cluster2, expected_energy):
    """
    Test the local minima function.
    """
    bi_water = Cluster(cluster1, cluster2)
    local_minima = LocalMinima(initial_cluster=bi_water)

    assert local_minima.shgo_optimize().fun - expected_energy[0] < 0.1
