from pyscf import gto, scf, dft
import attr
from copy import deepcopy


from amcess.base_molecule import Cluster


@attr.s
class Electronic_energy:
    """Functionn to calculate the energy of a cluster rotated and translated
    in a controled way.

    You can introduce the basis set you want to use and the method, which is
    HF for defect, or can introduce any functional for the DFT calculation.

    Returns
    -------
    [Electronic_energy object]]
    """
    initial_cluster = attr.ib(default=None, type=Cluster)
    basis = attr.ib(default='sto-3g')
    method = attr.ib(default='hf')

    def controled_move(self, displacements, angular_rotations):
        """Function to move molecules of a cluster in a controled way

        Parameters
        ----------
        displacements : [list]
            [Displacements in x,y,z directions]
        angular_rotations : [list]
            [description]

        Returns
        -------
        [Cluster]
            [The cluster with the new positions]
        """
        new_cluster = deepcopy(self.initial_cluster)
        for i in range(new_cluster.total_molecules-1):
            new_cluster = new_cluster.translate(
                i, displacements[i][0],
                displacements[i][1],
                displacements[i][2]).rotate(
                i, angular_rotations[i][0],
                angular_rotations[i][0],
                angular_rotations[i][0])

        return new_cluster

    def energy(self, x):
        """Function to calculate the energy of a cluster

        Parameters
        ----------
        x : [numpy array]
            [The array with the displacements and rotations]

        Returns
        -------
        [type]
            [The energy of the cluster]
        """
        cluster = deepcopy(self.initial_cluster)
        displacements = []
        angular_rotations = []
        for i in range(cluster.total_molecules-1):
            displacements.append(x[i*3:i*3+3])
            angular_rotations.append(x[(
                cluster.total_molecules-1+i)*3:
                (cluster.total_molecules+i)*3])

        candidate_to_test = self.controled_move(
                                displacements, angular_rotations)
        mol = gto.Mole()
        mol.fromstring(str(candidate_to_test.xyz))
        mol.build(basis=self.basis)
        if self.method == 'hf':
            mf = scf.RHF(mol)
            mf.conv_tol = 1e-10
            mf.conv_tol_grad = 1e-10
            mf.max_cycle = 200
        elif self.method == 'b3lyp':
            mf = dft.RKS(mol)
            mf.xc = 'b3lyp'
            mf.conv_tol = 1e-10
            mf.conv_tol_grad = 1e-10
            mf.max_cycle = 200
        mf.kernel()
        return mf.energy_tot()
