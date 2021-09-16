from pyscf import gto, scf

def heisenberg(x, T=None, *args):
    mol = gto.M(atom = 'H 0.0 x[0] x[1]; H 0 0 0', basis = 'sto-3g',)
    e = scf(mol)
    return e.kernel()