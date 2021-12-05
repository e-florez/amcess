def compute_energy(molecule, energy_function: dict):
    """Energy unit are Hartree"""

    # ------------------------------------------------------
    # energy setting
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


def pyscf_energy(
    molecule,  #: Cluster,
    hamiltonian: str = "hf",
    basis_set: str = "sto-3g",
    cpu: str = "1",
) -> float:
    from pyscf import lib
    from pyscf import (
        gto,
        scf,
        dft,
        cc,
        mp,
        mcscf,
        fci,
        lo,
        df,
        ao2mo,
    )

    # --------------------------------------------------

    lib.num_threads(n=cpu)

    input_mol = molecule.xyz.replace("\t", "").split("\n")
    input_mol = input_mol[2:]
    input_mol = ";".join(input_mol)

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
