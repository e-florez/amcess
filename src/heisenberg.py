import numpy as np
from pyscf import gto, scf

import src.search_configuration as SC
from src.base_molecule import Cluster


def beforeatoms(system):
    """[summary]
    Store attribute of before object made with Molecule

    Args:
        system ([type]): [description]

    Returns:
        [type]: [description]
    """
    beforeatoms.system = system or beforeatoms.system
    return beforeatoms.system


def build_input_pyscf(x, molecule_object, type_search, ncall=[0]):
    """[summary]
    Build input to pyscf

    Args:
        x ([array 1D]): possible new positions and angles.
                    dual_annealing use couchy distribution.
                    shgo #?
        system_object ([type]): [description]
        icall ([type]): [description]

    Returns:
        [type]: [description]
    """
    x_random = np.zeros((len(x[:])), dtype=float)
    if ncall[0] == 0:
        system_object = molecule_object
    else:
        system_object = beforeatoms.system
        if type_search == 1:
            # #! dual_annealing combinado con la propuesta
            # #! de juan para bounds y restando el centro
            # #!de masa al parecer evita que se alejen las
            # #! las moleculas

            mass_centers = np.zeros(
                (system_object.total_molecules, 3), dtype=float
            )

            for i in range(system_object.total_molecules):
                mass_centers[i] = system_object.get_molecule(i).center_of_mass

            # #! Es para que black no me dañe el formato 0: a 0 :
            # fmt: off
            r_random = (
                x[0: system_object.total_molecules * 3]
                - mass_centers.reshape((system_object.total_molecules * 3))
            )

            theta_random = x[
                system_object.total_molecules
                * 3: system_object.total_molecules * 6
            ]
            # fmt: on

            x_random = np.concatenate((r_random, theta_random))
        elif type_search == 2:
            # #! SHGO sobrepone las geometrías cuando resto el centor de
            # #! masas
            x_random = x

    new_geom = dict()
    for i in range(system_object.total_molecules):
        if (len(system_object.get_molecule(i).symbols)) > 1:
            # Rotate and translate
            new_geom[i] = {
                "atoms": system_object.rotate(
                    i,
                    x_random[(i + system_object.total_molecules) * 3],
                    x_random[(i + system_object.total_molecules) * 3 + 1],
                    x_random[(i + system_object.total_molecules) * 3 + 2],
                )
                .translate(
                    i,
                    x_random[i * 3],
                    x_random[i * 3 + 1],
                    x_random[i * 3 + 2],
                )
                .get_molecule(i)
                .atoms
            }
        else:
            # Translate
            new_geom[i] = {
                "atoms": system_object.get_molecule(i)
                .translate(
                    0,
                    x_random[i * 3],
                    x_random[i * 3 + 1],
                    x_random[i * 3 + 2],
                )
                .atoms
            }

    system_object = Cluster(*new_geom.values())

    # Controla que no este superpuestas las moleculas,
    # para evitar exit con pyscf y SHGO
    # En el caso que este superpuestas
    # una molecula se translada una longitud aleatoria
    # y un angulo aleatorio entre 0 a 1 en las 3 direcciones
    system_object = Cluster(*SC.overlaping(system_object).values())

    beforeatoms(system_object)
    ###
    symbols = system_object.symbols
    input_mol = "'"
    for i in range(system_object.total_atoms):
        input_mol += str(symbols[i])
        for j in range(3):
            input_mol += "  " + str(system_object.coordinates[i][j])
        if i < system_object.total_atoms - 1:
            input_mol += "; "
        else:
            input_mol += " '"
    ncall[0] += 1  # count calls
    return input_mol, system_object


def hamiltonian_pyscf(x, *args):
    """[summary]
    Calculate of electronic energy with pyscf

    Args:
        x [array 1D]: possible new positions and angles.
                    dual_annealing use couchy distribution
        shgo #?
        args [list]: basis set, Object of Molecule,
                    name output xyz, type optimizaiton
    Returns:
        e [floar]: [description]
    """

    input_pyscf, new_object = build_input_pyscf(x, args[1], args[3])

    mol = gto.M(
        atom=input_pyscf,
        basis=args[0],
    )

    e = scf.HF(mol).kernel()
    args[2].write(str(new_object.total_atoms) + "\n")
    args[2].write("Energy: " + str(e) + "\n")
    l: int = 0
    for symbols in new_object.symbols:
        args[2].write(
            str(symbols)
            + "  "
            +
            # 1 A = 1.88973 Bohr
            str(new_object.atoms[l][1] / 1.88973)
            + "  "
            + str(new_object.atoms[l][2] / 1.88973)
            + "  "
            + str(new_object.atoms[l][3] / 1.88973)
            + "\n"
        )
        l: int = l + 1

    return e
