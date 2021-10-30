import os

import pytest
from amcess.base_molecule import Cluster
from amcess.electronic_energy import ElectronicEnergy, hf_pyscf


@pytest.mark.parametrize(
    "x0, bases, cluster1, cluster2, output, sphere_center, sphere_radius,"
    "expected_energy",
    [
        (
            [1, 2, 3, 1, 2, 3, 0.9, 0.8, 0.7, 0.09, 0.08, 0.07],
            "sto-3g",
            [("H", 1, 1, 1), ("H", 1.74, 1, 1)],
            [("H", 0, 0, 0), ("H", 0.74, 0, 0)],
            "output.xyz",
            (0.0, 0.0, 0.0),
            1.0,
            -1.90426278138573,
        ),
    ],
)
def test_electronic_energy_with_hf_pyscf(
    x0,
    bases,
    cluster1,
    cluster2,
    output,
    sphere_center,
    sphere_radius,
    expected_energy,
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
        sphere_center : tuple
            Center mass of clusters
        sphere_radius : float
            Radius of the sphere centered in the cluster center mass
        expected_energy : float
            Expected electronic energy

    """
    with open(output, "w") as outxyz:
        e = hf_pyscf(
            x0,
            bases,
            ElectronicEnergy(
                Cluster(cluster1, cluster2),
                sphere_center,
                sphere_radius,
            ),
            outxyz,
        )
    os.remove("output.xyz")
    assert e - expected_energy < 1.0e-7


@pytest.mark.parametrize(
    "cluster1, cluster2, sphere_center, sphere_radius",
    [
        (
            [("H", 1, 1, 1), ("H", 1.74, 1, 1)],
            [("H", 0, 0, 0), ("H", 0.74, 0, 0)],
            (0.0, 0.0, 0.0),
            1.0,
        ),
    ],
)
def test_ElectronicEnergy_object_system_initial_grep(
    cluster1,
    cluster2,
    sphere_center,
    sphere_radius,
):
    """
    Test for electronic energy object save inital object system

    Parameters
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
        sphere_center : tuple
            Center mass of clusters
        sphere_radius : float
            Radius of the sphere centered in the cluster center mass
    """
    assert (
        ElectronicEnergy(
            Cluster(cluster1, cluster2),
            sphere_center,
            sphere_radius,
        )._object_system_initial
        == Cluster(cluster1, cluster2)
    )


@pytest.mark.parametrize(
    "cluster1, cluster2, cluster3, sphere_center, sphere_radius",
    [
        (
            [("H", 1, 1, 1), ("H", 1.74, 1, 1)],
            [("H", 0, 0, 0), ("H", 0.74, 0, 0)],
            [("H", 2, 2, 2), ("H", 2.74, 2, 2)],
            (0.0, 0.0, 0.0),
            1.0,
        ),
    ],
)
def test_ElectronicEnergy_object_system_initial_set(
    cluster1,
    cluster2,
    cluster3,
    sphere_center,
    sphere_radius,
):
    """
    Test for electronic energy object overwite before, initial
    and current object system when initialize new object system

    Parameters
        cluster1 : dict
            Dictionary with symbols and coordinates of the
            first molecule
        cluster2 : dict
            Dictionary with symbols and coordinates of the
            second molecule
        cluster3 : dict
            Dictionary with symbols and coordinates of the
            second molecule
        output : string
            Name of the output xyz to save energy and coordinates
        type_searching : int
            Type of searching used to find the candidates
            structures
        sphere_center : tuple
            Center mass of clusters
        sphere_radius : float
            Radius of the sphere centered in the cluster center mass

    """
    new_obj_ee = ElectronicEnergy(
        Cluster(cluster1, cluster2),
        sphere_center,
        sphere_radius,
    )
    new_obj_ee.object_system_initial = Cluster(cluster3, cluster2)
    assert (
        new_obj_ee.object_system_before,
        new_obj_ee.object_system_initial,
        new_obj_ee.object_system_current,
    ) == (
        Cluster(cluster3, cluster2),
        Cluster(cluster3, cluster2),
        Cluster(cluster3, cluster2),
    )


@pytest.mark.parametrize(
    "cluster1, cluster2, sphere_center, sphere_radius",
    [
        (
            [("H", 1, 1, 1), ("H", 1.74, 1, 1)],
            [("H", 0, 0, 0), ("H", 0.74, 0, 0)],
            (0.0, 0.0, 0.0),
            1.0,
        ),
    ],
)
def test_ElectronicEnergy_object_system_before_grep(
    cluster1,
    cluster2,
    sphere_center,
    sphere_radius,
):
    """
    Test for electronic energy object save before object system

    Parameters
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
        sphere_center : tuple
            Center mass of clusters
        sphere_radius : float
            Radius of the sphere centered in the cluster center mass
    """
    assert (
        ElectronicEnergy(
            Cluster(cluster1, cluster2),
            sphere_center,
            sphere_radius,
        )._object_system_before
        == Cluster(cluster1, cluster2)
    )


@pytest.mark.parametrize(
    "cluster1, cluster2, sphere_center, sphere_radius",
    [
        (
            [("H", 1, 1, 1), ("H", 1.74, 1, 1)],
            [("H", 0, 0, 0), ("H", 0.74, 0, 0)],
            (0.0, 0.0, 0.0),
            1.0,
        ),
    ],
)
def test_ElectronicEnergy_object_system_current_grep(
    cluster1,
    cluster2,
    sphere_center,
    sphere_radius,
):
    """
    Test for electronic energy object save current object system

    Parameters
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
        sphere_center : tuple
            Center mass of clusters
        sphere_radius : float
            Radius of the sphere centered in the cluster center mass
    """
    assert (
        ElectronicEnergy(
            Cluster(cluster1, cluster2),
            sphere_center,
            sphere_radius,
        )._object_system_current
        == Cluster(cluster1, cluster2)
    )


@pytest.mark.parametrize(
    "cluster1, cluster2, cluster3, sphere_center, sphere_radius",
    [
        (
            [("H", 1, 1, 1), ("H", 1.74, 1, 1)],
            [("H", 0, 0, 0), ("H", 0.74, 0, 0)],
            [("H", 2, 2, 2), ("H", 2.74, 2, 2)],
            (0.0, 0.0, 0.0),
            1.0,
        ),
    ],
)
def test_ElectronicEnergy_object_system_current_set(
    cluster1,
    cluster2,
    cluster3,
    sphere_center,
    sphere_radius,
):
    """
    Test for electronic energy object overwite before and
    current object system when initialize new current object

    Parameters
        cluster1 : dict
            Dictionary with symbols and coordinates of the
            first molecule
        cluster2 : dict
            Dictionary with symbols and coordinates of the
            second molecule
        cluster3 : dict
            Dictionary with symbols and coordinates of the
            second molecule
        output : string
            Name of the output xyz to save energy and coordinates
        type_searching : int
            Type of searching used to find the candidates
            structures
        sphere_center : tuple
            Center mass of clusters
        sphere_radius : float
            Radius of the sphere centered in the cluster center mass

    """
    new_obj_ee = ElectronicEnergy(
        Cluster(cluster1, cluster2),
        sphere_center,
        sphere_radius,
    )
    new_obj_ee.object_system_current = Cluster(cluster3, cluster2)
    assert (
        new_obj_ee.object_system_before,
        new_obj_ee.object_system_current,
    ) == (Cluster(cluster1, cluster2), Cluster(cluster3, cluster2))
