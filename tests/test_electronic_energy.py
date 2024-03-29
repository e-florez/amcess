from pyscf import gto
import pytest

from amcess.cluster import Cluster
from amcess.electronic_energy import ElectronicEnergy


@pytest.mark.parametrize(
    "x0, cluster1, cluster2, method_min, basis, hf, expected_energy",
    [
        (
            [1, 2, 3, 1, 2, 3, 0.9, 0.8, 0.7, 0.09, 0.08, 0.07],
            [("H", 1, 1, 1), ("H", 1.74, 1, 1)],
            [("H", 0, 0, 0), ("H", 0.74, 0, 0)],
            "dual_annealing",
            "sto-3g",
            "HF",
            -1.90426278138573,
        ),
    ],
)
def test_electronic_energy_with_hf_pyscf(
    x0,
    cluster1,
    cluster2,
    method_min,
    hf,
    basis,
    expected_energy,
):
    """
    Test for electronic energy calculate with hf_pyscf

    Parameters
    ----------
        x0: array 1D
            Values to translate and rotate each molecule
        cluster1 : dict
            Dictionary with symbols and coordinates of the
            first molecule
        cluster2 : dict
            Dictionary with symbols and coordinates of the
            second molecule
        method_min: string
            Method to minimize the energy
        basis: string
            Label of the basis set
        hf: string
            Hartree--Fock Hamiltonian
        expected_energy : float
            Expected electronic energy

    """
    e = ElectronicEnergy(
        Cluster(cluster1, cluster2),
        method_min,
        hf,
        basis,
    ).pyscf(x0)
    assert e - expected_energy < 1.0e-7


@pytest.mark.parametrize(
    "cluster1, cluster2, method_min, hf, basis",
    [
        (
            [("H", 1, 1, 1), ("H", 1.74, 1, 1)],
            [("H", 0, 0, 0), ("H", 0.74, 0, 0)],
            "dual_annealing",
            "HF",
            "sto-3g",
        ),
    ],
)
def test_ElectronicEnergy_object_system_initial_grep(
    cluster1,
    cluster2,
    method_min,
    hf,
    basis,
):
    """
    Test for electronic energy object save inital object system

    Parameters
    ----------
        cluster1 : dict
            Dictionary with symbols and coordinates of the
            first molecule
        cluster2 : dict
            Dictionary with symbols and coordinates of the
            second molecule
        method_min: string
            Method to minimize the energy
        hf: string
            Hartree--Fock Hamiltonian
        basis : string
            Label of the basis set
    """
    assert (
        ElectronicEnergy(
            Cluster(cluster1, cluster2),
            method_min,
            hf,
            basis,
        )._object_system_initial
        == Cluster(cluster1, cluster2)
    )


@pytest.mark.parametrize(
    "cluster1, cluster2, cluster3, method_min, hf, basis",
    [
        (
            [("H", 1, 1, 1), ("H", 1.74, 1, 1)],
            [("H", 0, 0, 0), ("H", 0.74, 0, 0)],
            [("H", 2, 2, 2), ("H", 2.74, 2, 2)],
            "dual_annealing",
            "HF",
            "sto-3g",
        ),
    ],
)
def test_ElectronicEnergy_object_system_initial_set(
    cluster1,
    cluster2,
    cluster3,
    method_min,
    hf,
    basis,
):
    """
    Test for electronic energy object overwite before, initial
    and current object system when initialize new object system

    Parameters
    ----------
        cluster1 : dict
            Dictionary with symbols and coordinates of the
            first molecule
        cluster2 : dict
            Dictionary with symbols and coordinates of the
            second molecule
        cluster3 : dict
            Dictionary with symbols and coordinates of the
            second molecule
        method_min: string
            Method to minimize the energy
        hf: string
            Hartree--Fock Hamiltonian
        basis : string
            Label of the basis set

    """
    new_obj_ee = ElectronicEnergy(
        Cluster(cluster1, cluster2),
        method_min,
        hf,
        basis,
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
    "cluster1, cluster2, method_min, hf, basis",
    [
        (
            [("H", 1, 1, 1), ("H", 1.74, 1, 1)],
            [("H", 0, 0, 0), ("H", 0.74, 0, 0)],
            "dual_annealing",
            "HF",
            "sto-3g",
        ),
    ],
)
def test_ElectronicEnergy_object_system_before_grep(
    cluster1,
    cluster2,
    method_min,
    hf,
    basis,
):
    """
    Test for electronic energy object save before object system

    Parameters
    ----------
        cluster1 : dict
            Dictionary with symbols and coordinates of the
            first molecule
        cluster2 : dict
            Dictionary with symbols and coordinates of the
            second molecule
        method_min: string
            Method to minimize the energy
        hf: string
            Hartree--Fock Hamiltonian
        basis   : string
            Label of the basis set
    """
    assert (
        ElectronicEnergy(
            Cluster(cluster1, cluster2),
            method_min,
            hf,
            basis,
        )._object_system_before
        == Cluster(cluster1, cluster2)
    )


@pytest.mark.parametrize(
    "cluster1, cluster2, method_min, hf, basis",
    [
        (
            [("H", 1, 1, 1), ("H", 1.74, 1, 1)],
            [("H", 0, 0, 0), ("H", 0.74, 0, 0)],
            "dual_annealing",
            "HF",
            "sto-3g",
        ),
    ],
)
def test_ElectronicEnergy_object_system_current_grep(
    cluster1,
    cluster2,
    method_min,
    hf,
    basis,
):
    """
    Test for electronic energy object save current object system

    Parameters
    ----------
        cluster1 : dict
            Dictionary with symbols and coordinates of the
            first molecule
        cluster2 : dict
            Dictionary with symbols and coordinates of the
            second molecule
        method_min: string
            Method to minimize the energy
        hf: string
            Hartree-Fock Hamiltonian
        basis : string
            Basis set used to calculate the electronic energy
    """
    assert (
        ElectronicEnergy(
            Cluster(cluster1, cluster2),
            method_min,
            hf,
            basis,
        )._object_system_current
        == Cluster(cluster1, cluster2)
    )


@pytest.mark.parametrize(
    "cluster1, cluster2, cluster3, method_min, hf, basis",
    [
        (
            [("H", 1, 1, 1), ("H", 1.74, 1, 1)],
            [("H", 0, 0, 0), ("H", 0.74, 0, 0)],
            [("H", 2, 2, 2), ("H", 2.74, 2, 2)],
            "dual_annealing",
            "HF",
            "sto-3g",
        ),
    ],
)
def test_ElectronicEnergy_object_system_current_set(
    cluster1,
    cluster2,
    cluster3,
    method_min,
    hf,
    basis,
):
    """
    Test for electronic energy object overwite before and
    current object system when initialize new current object

    Parameters
    ----------
        cluster1 : dict
            Dictionary with symbols and coordinates of the
            first molecule
        cluster2 : dict
            Dictionary with symbols and coordinates of the
            second molecule
        cluster3 : dict
            Dictionary with symbols and coordinates of the
            second molecule
        method_min: string
            Method to minimize the energy
        hf : string
            Hartree-Fock Hamiltonian
        basis : string
            Basis set used to calculate the electronic energy

    """
    new_obj_ee = ElectronicEnergy(
        Cluster(cluster1, cluster2),
        method_min,
        hf,
        basis,
    )
    new_obj_ee.object_system_current = Cluster(cluster3, cluster2)
    assert (
        new_obj_ee.object_system_before,
        new_obj_ee.object_system_current,
    ) == (Cluster(cluster1, cluster2), Cluster(cluster3, cluster2))


@pytest.mark.parametrize(
    "molecule1, method_min, hf, bases, mol_atom_input",
    [
        (
            [("H", 0.0, 0.0, 0.0), ("H", 0.00, 0.0, 0.0)],
            "dual_annealing",
            "HF",
            "sto-3g",
            "H 0 0 0; H 0 0 0",
        ),
    ],
)
def test_error_energy_hf_pyscf(
    molecule1,
    method_min,
    hf,
    bases,
    mol_atom_input,
):
    """
    Test Warning of pyscf
    """
    with pytest.warns(UserWarning) as w:
        obj_sc = ElectronicEnergy(
            Cluster(molecule1),
            method_min,
            hf,
            bases,
        )
        mol = gto.M(
            atom=mol_atom_input,
            basis=bases,
            verbose=False,
        )
        obj_sc.calculate_electronic_energy(mol)
    assert str(w) == "WarningsChecker(record=True)"


@pytest.mark.parametrize(
    "molecule1, molecule2, method_min, hf, bases",
    [
        (
            [("H", 0.0, 0.0, 0.0), ("H", 0.00, 0.0, 0.0)],
            [("H", 0.0, 0.0, 0.0), ("H", 0.00, 0.0, 0.0)],
            "dual_annealing",
            "HF",
            "sto-3g",
        ),
    ],
)
def test_metropolis_else(
    molecule1,
    molecule2,
    method_min,
    hf,
    bases,
):
    """
    Verification of else of the Metropolis method
    """
    obj_sc = ElectronicEnergy(
        Cluster(molecule1, molecule2),
        method_min,
        hf,
        bases,
    )
    obj_sc.energy_current = 1.0
    obj_sc.energy_before = 1.0
    obj_sc.metropolis()


@pytest.mark.parametrize(
    "molecule1, molecule2, method_min, dft, bases, expected",
    [
        (
            [("H", 0.0, 0.0, 0.0), ("H", 0.00, 0.0, 0.0)],
            [("H", 0.0, 0.0, 0.0), ("H", 0.00, 0.0, 0.0)],
            "dual_annealing",
            "DFT b3lyp",
            "sto-3g",
            "b3lyp",
        ),
    ],
)
def test_functional(
    molecule1,
    molecule2,
    method_min,
    dft,
    bases,
    expected,
):
    """
    Verification of functional of dft choose
    """
    obj_sc = ElectronicEnergy(
        Cluster(molecule1, molecule2),
        method_min,
        dft,
        bases,
    )
    assert obj_sc._functional == expected


@pytest.mark.parametrize(
    "molecule1, molecule2, method_min, methodology, bases",
    [
        (
            [("H", 0.0, 0.0, 0.0), ("H", 0.00, 0.0, 0.0)],
            [("H", 0.0, 0.0, 0.0), ("H", 0.00, 0.0, 0.0)],
            "dual_annealing",
            "pm",
            "sto-3g",
        ),
    ],
)
def test_method_not_implemented(
    molecule1,
    molecule2,
    method_min,
    methodology,
    bases,
):
    """
    Verification of method not implemented
    """
    with pytest.raises(ValueError) as e:
        ElectronicEnergy(
            Cluster(molecule1, molecule2),
            method_min,
            methodology,
            bases,
        )
    assert str(e.value) == "Methodology not implemented"


@pytest.mark.parametrize(
    "x0, cluster1, cluster2, method_min, basis, mp2, expected_energy",
    [
        (
            [1, 2, 3, 1, 2, 3, 0.9, 0.8, 0.7, 0.09, 0.08, 0.07],
            [("H", 1, 1, 1), ("H", 1.74, 1, 1)],
            [("H", 0, 0, 0), ("H", 0.74, 0, 0)],
            "dual_annealing",
            "sto-3g",
            "MP2",
            -2.254050875005756,
        ),
    ],
)
def test_electronic_energy_with_mp2_pyscf(
    x0,
    cluster1,
    cluster2,
    method_min,
    mp2,
    basis,
    expected_energy,
):
    """
    Test for electronic energy calculate with MP2 in the pyscf
    """
    e = ElectronicEnergy(
        Cluster(cluster1, cluster2),
        method_min,
        mp2,
        basis,
    ).pyscf(x0)
    assert e - expected_energy < 1.0e-7


@pytest.mark.parametrize(
    "x0, cluster1, cluster2, method_min, basis, ccsd, expected_energy",
    [
        (
            [1, 2, 3, 1, 2, 3, 0.9, 0.8, 0.7, 0.09, 0.08, 0.07],
            [("H", 1, 1, 1), ("H", 1.74, 1, 1)],
            [("H", 0, 0, 0), ("H", 0.74, 0, 0)],
            "dual_annealing",
            "sto-3g",
            "CCSD",
            -2.2689731245694804,
        ),
    ],
)
def test_electronic_energy_with_ccsd_pyscf(
    x0,
    cluster1,
    cluster2,
    method_min,
    ccsd,
    basis,
    expected_energy,
):
    """
    Test for electronic energy calculate with MP2 in the pyscf
    """
    e = ElectronicEnergy(
        Cluster(cluster1, cluster2),
        method_min,
        ccsd,
        basis,
    ).pyscf(x0)
    assert e - expected_energy < 1.0e-7
