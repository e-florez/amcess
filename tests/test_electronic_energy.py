import pytest
from amcess.base_molecule import Cluster
from amcess.electronic_energy import hf_pyscf


@pytest.mark.parametrize(
    "x0, bases, cluster1, cluster2, output, type_searching, expected_energy",
    [
        (
            [1, 2, 3, 1, 2, 3, 0.9, 0.8, 0.7, 0.09, 0.08, 0.07],
            "sto-3g",
            [("H", 0, 0, 1), ("H", 0.74, 0, 1)],
            [("H", 0, 0, 0), ("H", 0.74, 0, 0)],
            "output.xyz",
            1,
            -1.90426278138573,
        ),
    ],
)
def test_electronic_energy_with_hf_pyscf(
    x0, bases, cluster1, cluster2, output, type_searching, expected_energy
):
    """
    Test for electronic energy calculate with hf_pyscf

    Parameters
    ----------
        x0: array 1D
            Values to translate and rotate each molecule
        bases: string
            Label of the bases set
        cluster1 : dict
            Dictionary with symbols and coordinates of the
            first molecule
        cluster2 : dict
            Dictionary with symbols and coordinates of the
            second molecule
        output : string
            Name of the output xyz to save energy and coordinates
        type_searching : int
            Type of searching used to find the candidates
            structures
        expected_energy : float
            Expected electronic energy

    """
    with open(output, "w") as outxyz:
        e = hf_pyscf(
            x0, bases, Cluster(cluster1, cluster2), outxyz, type_searching
        )
    assert e - expected_energy < 1.0e-7
