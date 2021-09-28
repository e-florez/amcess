from pyscf import gto, scf
import numpy as np
#from search_configuration import SearchConfig


def build_input_pyscf(mass_center, atomic_symb_principal_axes):
#    M = 1.008*2+15.999
#    rH1a = atomic_symb_principal_axes[6][0] + mass_center[0:3]/M
#    rO1 = atomic_symb_principal_axes[6][1] + mass_center[0:3]/M
#    rH2a = atomic_symb_principal_axes[6][2] + mass_center[0:3]/M
#    rH1b = atomic_symb_principal_axes[6][3] + mass_center[3:6]/M
#    rO2 = atomic_symb_principal_axes[6][4] +  mass_center[3:6]/M
#    rH2b = atomic_symb_principal_axes[6][5] + mass_center[3:6]/M

#    print(6,"\n")
#    print('H ',rH1a[0], rH1a[1], rH1a[2])
#    print('O ',rO1[0], rO1[1], rO1[2])
#    print('H ',rH2a[0], rH2a[1], rH2a[2])
#    print('H ',rH1b[0], rH1b[1], rH1b[2])
#    print('O ',rO2[0], rO2[1], rO2[2])
#    print('H ',rH2b[0], rH2b[1], rH2b[2])

#    mol = str(atomic_symb_principal_axes[0]) + "  " + str(rH1a[0]) + "  " +\
#    str(rH1a[1]) + "  " + str(rH1a[2]) + ";" + "  "
#    mol = mol + str(atomic_symb_principal_axes[1]) + "  " + str(rO1[0]) + "  " +\
#    str(rO1[1]) + "  " + str(rO1[2]) + ";" + "  "
#    mol = mol + str(atomic_symb_principal_axes[2]) + "  " + str(rH2a[0]) + "  " +\
#    str(rH2a[1]) + "  " + str(rH2a[2]) + ";" + "  "
#    mol = mol + str(atomic_symb_principal_axes[3]) + "  " + str(rH1b[0]) + "  " +\
#    str(rH1b[1]) + "  " + str(rH1b[2]) + ";" + "  "
#    mol = mol + str(atomic_symb_principal_axes[4]) + "  " + str(rO2[0]) + "  " +\
#    str(rO2[1]) + "  " + str(rO2[2]) + ";" + "  "
#    mol = mol + str(atomic_symb_principal_axes[5]) + "  " + str(rH2b[0]) + "  " +\
#    str(rH2b[1]) + "  " + str(rH2b[2])
#    mol = "'" + mol + " '"

    print(6,"\n")
    print('H ', mass_center[0], mass_center[1], mass_center[2])
    print('O ', mass_center[3], mass_center[4], mass_center[5])
    print('H ', mass_center[6], mass_center[7], mass_center[8])
    print('H ', mass_center[9], mass_center[10], mass_center[11])
    print('O ', mass_center[12], mass_center[13], mass_center[14])
    print('H ', mass_center[15], mass_center[16], mass_center[17])

    #mol = str(atomic_symb_principal_axes[0]) + "  " + str(mass_center[0]) + "  " +\
    #str(mass_center[1]) + "  " + str(mass_center[2]) + ";" + "  "
    #mol = mol + str(atomic_symb_principal_axes[1]) + "  " + str(mass_center[3]) + "  " +\
    #str(mass_center[4]) + "  " + str(mass_center[5]) + ";" + "  "
    #mol = mol + str(atomic_symb_principal_axes[2]) + "  " + str(mass_center[6]) + "  " +\
    #str(mass_center[7]) + "  " + str(mass_center[8]) + ";" + "  "
    #mol = mol + str(atomic_symb_principal_axes[3]) + "  " + str(mass_center[9]) + "  " +\
    #str(mass_center[10]) + "  " + str(mass_center[11]) + ";" + "  "
    #mol = mol + str(atomic_symb_principal_axes[4]) + "  " + str(mass_center[12]) + "  " +\
    #str(mass_center[13]) + "  " + str(mass_center[14]) + ";" + "  "
    #mol = mol + str(atomic_symb_principal_axes[5]) + "  " + str(mass_center[15]) + "  " +\
    #str(mass_center[16]) + "  " + str(mass_center[17])
    #mol = "'" + mol + " '"

    #print(mol)
    return mol

def hamiltonian_pyscf(x, *args):
    #input_pyscf = build_input_pyscf(mases_centers, args)
    #mol = gto.M(atom = input_pyscf, basis = 'sto-3g',)
    mol = gto.M(atom = [['H', x[0], x[1], x[2]],
                        ['O', x[3], x[4], x[5]],
                       ['H', x[6], x[7], x[8]]], basis = 'sto-3g',)
    DX = x[0]-x[3]
    DX *= DX
    DY = x[1]-x[4]
    DY *= DY
    DZ = x[2]-x[5]
    DZ *= DZ

    D1X = x[6]-x[3]
    D1X *= D1X
    D1Y = x[7]-x[4]
    D1Y *= D1Y
    D1Z = x[8]-x[5]
    D1Z *= D1Z

    e = scf.HF(mol).kernel()
    print("distance ",np.sqrt(DX+DY+DZ), e, np.sqrt(D1X+D1Y+D1Z))
    #exit()
    return e
    #return 0.0