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
    x_random = x #- np.asarray(fragment_center_of_mass).reshape(-1)
    new_geom = []
    for i in range(system_object._total_fragments):
        if (len(system_object.get_fragments(i).symbols)) > 1:
        #Rotate and translate
            new_geom.append(system_object.translate(i,x_random[i*3],x_random[i*3 + 1],x_random[i*3 + 2])\
            .rotate(i,x_random[i*3],x_random[i*3 + 1],x_random[i*3 + 2]))
        else:
        #Translate
            new_geom.append(system_object.translate(i,x_random[i*3],x_random[i*3 + 1],x_random[i*3 + 2]))

    w1 = {"coordinates": new_geom[0]}
    w2 = {"coordinates": new_geom[1]}
    system_object = Cluster(w1, w2)

    ###
    symbols = system_object.symbols
    input_mol = "'"
    for i in range(system_object._total_atoms):
        input_mol += str(symbols[i])
        for j in range(3):
            input_mol += "  " + str(system_object.cartesian_coordinates[i][j])
        if i < system_object._total_atoms - 1:
            input_mol +=  "; "
        else:
            input_mol +=  " '"
    #print(input_mol)
    #exit()
    return input_mol

def hamiltonian_pyscf(x, *args):

    #input_pyscf = build_input_pyscf(x, args[1], args[2])

    x_random = x #- np.asarray(fragment_center_of_mass).reshape(-1)
    system_object = args[1]
    new_geom = []
    for i in range(system_object._total_fragments):
        if (len(system_object.get_fragments(i).symbols)) > 1:
        #Rotate and translate
            new_geom.append(system_object.translate(i,x_random[i*3],x_random[i*3 + 1],x_random[i*3 + 2])\
            .rotate(i,x_random[i*3],x_random[i*3 + 1],x_random[i*3 + 2]))
        else:
        #Translate
            new_geom.append(system_object.translate(i,x_random[i*3],x_random[i*3 + 1],x_random[i*3 + 2]))
    #objectaa = print(*new_geom, sep = ", ")
    w1 = {}
    w2 = {}
    w1['coordinates'] = new_geom[0]
    w2['coordinates'] = new_geom[1]
    system_object = Cluster(w1, w2)

    ###
    symbols = system_object.symbols
    input_mol = "'"
    for i in range(system_object._total_atoms):
        input_mol += str(symbols[i])
        for j in range(3):
            input_mol += "  " + str(system_object.cartesian_coordinates[i][j])
        if i < system_object._total_atoms - 1:
            input_mol +=  "; "
        else:
            input_mol +=  " '"

    mol = gto.M(atom = input_mol, basis = args[0],)
    #mol = gto.M(atom = [['H', x[0], x[1], x[2]],
    #                    ['H', x[3], x[4], x[5]],
    #                   ], basis = args[0],)

    #DX = x[0]-x[3]
    #DX *= DX
    #DY = x[1]-x[4]
    #DY *= DY
    #DZ = x[2]-x[5]
    #DZ *= DZ
    print("   ",system_object._total_atoms)
    e = scf.HF(mol).kernel()
    l = 0
    for symbols in system_object.symbols:
        print(symbols,"  ",system_object._coordinates[l][1]/1.88973, #1 A = 1.88973 Bohr
        "  ",system_object._coordinates[l][2]/1.88973,
        "  ",system_object._coordinates[l][3]/1.88973)
        l += 1
    l = 0
    # args[1] = Cluster(w1, w2) no se deja actualizar el objeto 2h2o
    #print("distance ", e)
    #print("distance ",np.sqrt(DX+DY+DZ), e)
    #exit()
    return e
    #return 0.0