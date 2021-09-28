from pyscf import gto, scf
import numpy as np

from base_molecule import *
#from search_configuration import SearchConfig


def build_input_pyscf(x, system_object, icall):
#Evitar lo siguiente cuando se hace el primer llamado, para cuando x0 = coordinates

    fragment_center_of_mass = np.zeros((2,3),dtype=float)
    fragment: Cluster = system_object.get_fragments(0)
    fragment_center_of_mass[0] = fragment.center_of_mass
    fragment: Cluster = system_object.get_fragments(1)
    fragment_center_of_mass[1] = fragment.center_of_mass

    #if x.all() != np.asarray(system_object.center_of_mass).reshape(-1).all():
    x_random = np.zeros((len(x[:])),dtype=float)
    x_random = x - np.asarray(fragment_center_of_mass).reshape(-1)
    print("x: ",x_random)
    #x_random = x - np.asarray(system_object.cartesian_coordinates).reshape(-1)
    #else:
    #    x_random = x
    ###
    for i in range(system_object._total_fragments):
    #Translate
        print("i ",i)
        system_object._coordinates = system_object.translate(i,x_random[i*3],x_random[i*3 + 1],x_random[i*3 + 2])
        if (len(system_object.get_fragments(i).symbols)) > 1:
        #Rotate
            system_object._coordinates = system_object.rotate(i,x_random[i*3],x_random[i*3 + 1],x_random[i*3 + 2])
    print(system_object.cartesian_coordinates[0][0])
    exit()
    #
    ###
    symbols = system_object.symbols
    input_mol = "'"
    for i in range(system_object._total_atoms):
        input_mol += str(symbols[i])
        for j in range(3):
            input_mol += "  " + str(x[i*3+j])
        if i < system_object._total_atoms - 1:
            input_mol +=  "; "
        else:
            input_mol +=  " '"
    print(input_mol)
    exit()
    return input_mol

def hamiltonian_pyscf(x, *args):

    #input_pyscf = build_input_pyscf(x, args[1], args[2])
    #mol = gto.M(atom = input_pyscf, basis = args[0],)
    mol = gto.M(atom = [['H', x[0], x[1], x[2]],
                        ['H', x[3], x[4], x[5]],
                       ], basis = args[0],)
    #mol = gto.M(atom = [['H', x[0], x[1], x[2]],
    #                    ['O', x[3], x[4], x[5]],
    #                   ['H', x[6], x[7], x[8]]], basis = 'sto-3g',)
    DX = x[0]-x[3]
    DX *= DX
    DY = x[1]-x[4]
    DY *= DY
    DZ = x[2]-x[5]
    DZ *= DZ

    e = scf.HF(mol).kernel()
    #print("distance ", e)
    print("distance ",np.sqrt(DX+DY+DZ), e)
    #exit()
    return e
    #return 0.0