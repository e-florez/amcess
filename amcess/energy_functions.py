from pyscf import lib
from pyscf import (
    gto,
    scf,
    dft,
    cc,
    mp,
)


def compute_energy(molecule, energy_function: dict):
    """Energy unit are Hartree"""

    energy_method = energy_function.get("energy", None)

    if energy_method == "pyscf":
        hamiltonian = energy_function.get("hamiltonian", "hf")
        basis_set = energy_function.get("basis_set", "3-21g")

        cpu = energy_function.get("cpu", "4")

        # cpu should be and integer but type must be
        # a STRING according to pySCF
        try:
            cpu = str(int(cpu))
        except ValueError:
            cpu = "4"

        return pyscf_energy(molecule, hamiltonian, basis_set, cpu)
    else:
        return None


# ----------------------------------------------------------------------------
def write_pyscf_input(molecule):
    input_mol = molecule.xyz.replace("\t", "").split("\n")
    input_mol = input_mol[2:]
    input_mol = ";".join(input_mol)

    return input_mol


# ----------------------------------------------------------------------------
def pyscf_energy(
    molecule,  #: Cluster,
    hamiltonian: str = "hf",
    basis_set: str = "sto-3g",
    cpu: str = "1",
) -> float:

    # define CPU number
    lib.num_threads(n=cpu)

    input_mol = write_pyscf_input(molecule)

    mol = gto.M(
        atom=input_mol,
        basis=basis_set,
        verbose=False,
    )

    if hamiltonian == "hf":
        rhf = scf.RHF(mol)
        energy = rhf.kernel()
    elif hamiltonian in ["b3lyp", "pbe0", "wb97x-d", "pbe"]:
        dft_obj = dft.RKS(mol)
        dft_obj.xc = hamiltonian
        energy = dft_obj.kernel()
    elif hamiltonian == "mp2":
        mp2 = mp.MP2(scf.RHF(mol))
        energy = mp2.kernel()
        # mf = scf.RHF(mol).run()
        # # freeze 2 core orbitals
        # pt = mp.MP2(mf).set(frozen = 2).run()
        # # freeze 2 core orbitals and 3 high lying unoccupied orbitals
        # pt.set(frozen = [0,1,16,17,18]).run()

    elif hamiltonian == "ccsd":
        cc_obj = cc.CCSD(mol)
        energy = cc_obj.kernel()
    elif hamiltonian == "ccsd(t)":
        cc_obj = cc.CCSD(mol)
        energy = cc_obj.kernel(t=True)
    elif hamiltonian == "ccsd(at)":
        cc_obj = cc.CCSD(mol)
        energy = cc_obj.kernel(at=True)
    else:
        return float("inf")

    return energy
