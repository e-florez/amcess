from pyscf import gto, scf
import numpy as np

from base_molecule import *


def before_coordinates(system):
    """[summary]
    Store attribute of before object made with Molecule

    Args:
        system ([type]): [description]

    Returns:
        [type]: [description]
    """
    before_coordinates.system = system or before_coordinates.system
    return before_coordinates.system


def build_input_pyscf(x, molecule_object, ncall=[0]):
    """[summary]
    Build input to pyscf

    Args:
        x ([type]): [description]
        system_object ([type]): [description]
        icall ([type]): [description]

    Returns:
        [type]: [description]
    """
    x_random = np.zeros((len(x[:])), dtype=float)
    if ncall[0] == 0:
        system_object = molecule_object
    else:
        system_object = before_coordinates.system
        mass_centers = np.zeros((system_object._total_fragments, 3), dtype=float)
        for i in range(system_object._total_fragments):
            mass_centers[i] = system_object.get_fragments(i).center_of_mass
        x_random = x  # - mass_centers.reshape((system_object._total_fragments*3))

    new_geom = dict()
    for i in range(system_object._total_fragments):
        if (len(system_object.get_fragments(i).symbols)) > 1:
            # Rotate and translate
            #! aleja las moléculas de agua
            #! si le resto el cm a x, ahora casi no se mueven las moleculas
            # ? cuando no se alejan tiene difultad para rotar las moleculas de agua
            # ? para favorecer la disposición en que se forma el EH
            new_geom[i] = {
                "atoms": system_object.get_fragments(i)
                .rotate(0, x_random[i * 3], x_random[i * 3 + 1], x_random[i * 3 + 2])
                .translate(0, x_random[i * 3], x_random[i * 3 + 1], x_random[i * 3 + 2])
                ._coordinates
            }
        else:
            # Translate
            new_geom[i] = {
                "atoms": system_object.get_fragments(i)
                .translate(0, x_random[i * 3], x_random[i * 3 + 1], x_random[i * 3 + 2])
                ._coordinates
            }

    system_object = Molecule(*new_geom.values())
    before_coordinates(system_object)
    ###
    symbols = system_object.symbols
    input_mol = "'"
    for i in range(system_object._total_atoms):
        input_mol += str(symbols[i])
        for j in range(3):
            input_mol += "  " + str(system_object.cartesian_coordinates[i][j])
        if i < system_object._total_atoms - 1:
            input_mol += "; "
        else:
            input_mol += " '"
    ncall[0] += 1  # count calls
    return input_mol, system_object


def hamiltonian_pyscf(x, *args):
    """[summary]
    Calculate of electronic energy with pyscf

    Args:
        x ([type]): [description]

    Returns:
        [type]: [description]
    """

    input_pyscf, new_object = build_input_pyscf(x, args[1])

    mol = gto.M(
        atom=input_pyscf,
        basis=args[0],
    )

    print(new_object._total_atoms)
    e = scf.HF(mol).kernel()
    l = 0
    for symbols in new_object.symbols:
        print(
            symbols,
            "  ",
            new_object._coordinates[l][1] / 1.88973,  # 1 A = 1.88973 Bohr
            "  ",
            new_object._coordinates[l][2] / 1.88973,
            "  ",
            new_object._coordinates[l][3] / 1.88973,
        )
        l += 1
    l = 0
    # print("distance ",np.sqrt(DX+DY+DZ), e)
    # exit()
    return e
