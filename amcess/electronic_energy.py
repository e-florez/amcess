import numpy as np
from pyscf import gto, scf, dft, mp, cc

from amcess.base_molecule import Cluster


class ElectronicEnergy:
    def __init__(
        self,
        object_system: object,
        search_type: str,
        methodology: str,
        basis_set: str,
        seed: int = None,
    ) -> None:
        """
        Class to calculate electronic energy

        Attributes
        ----------
        object_system : object
            Object initialized with Molecule or Cluster class
        """

        self._object_system_initial = object_system
        self._object_system_before = object_system
        self._object_system_current = object_system

        self._search_type = search_type

        self._method = methodology.split()[0]
        if len(methodology.split()) > 1:
            self._functional = methodology.split()[1]
        self._basis_set = basis_set

        self._move_seed = seed

        self.store_structures = []

        if self._search_type != "ASCEC ":
            mol = gto.M(
                atom=self.input_atom_mol_pyscf(),
                basis=self._basis_set,
                verbose=False,
            )
            self._e0 = self.calculate_electronic_energy(mol)
            self.energy_before = self._e0

    # ===============================================================
    # Decorators
    # ===============================================================
    def build_input_gto_pyscf(func_energy):
        def new_input_gto(self, x):
            """
            Build input to pyscf

            Parameters
            ----------
            x : array 1D
                possible new positions and angles.
            system_object : object Cluster
                Object initialized with Molecule or Cluster class

            Returns
            -------
            input_gto_pyscf: list
                Atom's symbols and coordinates
            system_object: Cluster
                Cluster objects
            """
            system_object = self._object_system_current

            # ------------------------------------------------------------
            # Rotate and translate each molecule into object Cluster
            new_geom = dict()
            new_geom[0] = {"atoms": system_object.get_molecule(0).atoms}
            for i in range(system_object.total_molecules - 1):
                new_geom[i + 1] = {
                    "atoms": system_object.translate(
                        i + 1,
                        x=x[i * 3],
                        y=x[i * 3 + 1],
                        z=x[i * 3 + 2],
                    )
                    .rotate(
                        i + 1,
                        x=x[(i + system_object.total_molecules - 1) * 3],
                        y=x[(i + system_object.total_molecules - 1) * 3 + 1],
                        z=x[(i + system_object.total_molecules - 1) * 3 + 2],
                    )
                    .get_molecule(i + 1)
                    .atoms
                }

            # ------------------------------------------------------------
            # New object Cluster with new geometries
            self.object_system_current = Cluster(
                *new_geom.values(),
                sphere_radius=system_object._sphere_radius,
                sphere_center=system_object._sphere_center
            )
            # ------------------------------------------------------------
            # Build input to pyscf
            self.input_atom_mol_pyscf()

            return func_energy(self, x)

        return new_input_gto

    # ===============================================================
    # PROPERTIES
    # ===============================================================
    @property
    def object_system_initial(self):
        return self._object_system_initial

    @object_system_initial.setter
    def object_system_initial(self, new_object_system):
        (
            self._object_system_initial,
            self._object_system_before,
            self._object_system_current,
        ) = (
            new_object_system,
            new_object_system,
            new_object_system,
        )

    @property
    def object_system_before(self):
        return self._object_system_before

    @property
    def object_system_current(self):
        return self._object_system_current

    @object_system_current.setter
    def object_system_current(self, new_object_system):
        self._object_system_before = self._object_system_current
        self._object_system_current = new_object_system

    # ===============================================================
    # Methods
    # ===============================================================
    def input_atom_mol_pyscf(self):
        """
        Build a portion of the input for the gto object of pyscf
            'X 0.0 0.0 0.0; X 0.0 0.0 1.0'

        Returns
        -------
            input_gto_pyscf: list
                Atom's symbols and coordinates
        """

        symbols = self._object_system_current.symbols
        self.input_gto_pyscf = "'"
        for i in range(self._object_system_current.total_atoms):
            self.input_gto_pyscf += str(symbols[i])
            for j in range(3):
                self.input_gto_pyscf += "  " + str(
                    self._object_system_current.coordinates[i][j]
                )
            if i < self._object_system_current.total_atoms - 1:
                self.input_gto_pyscf += "; "
            else:
                self.input_gto_pyscf += " '"
        return self.input_gto_pyscf

    def write_to_file(self, filename):
        """
        Write all accepted structures to a file

        Parameters
        ----------
        filename: str
            File name where is save structure and energy
        """

        n_atoms = len(self.store_structures[0]) - 1
        with open(filename, "w") as f:
            for system in self.store_structures:
                f.write(str(n_atoms) + "\n")
                f.write("Energy: " + str(system[0]) + "\n")
                for terms in system[1:]:
                    f.write(" ".join([str(x) for x in terms]) + "\n")

    def store_structure(self):
        """
        Store the accept systems in a list of lists of the energy more
        a tuples with the coordinates

        [[Energy, ('X', 0., 0., 0.), ('Y', 1., 0., 0.)], [ ... ], ...]
        """
        self.store_structures.append(
            [self.energy_current] + self._object_system_current.atoms
        )

    def calculate_electronic_energy(self, mol):
        """
        Calculate electronic energy with pyscf

        Parameters
        ----------
        mol: object
            gto pyscf object

        Returns
        -------
        Electronic energy
        """
        try:
            if self._method == "HF":
                return scf.RHF(mol).kernel()
            elif self._method == "DFT":
                dft_call = dft.RKS(mol)
                dft_call.xc = self._functional
                return dft_call.kernel()
            elif self._method == "MP2":
                mf = scf.RHF(mol).run()
                energy = mp.MP2(mf).run()
                return energy.e_tot
            elif self._method == "CCSD":
                mf = scf.RHF(mol).run()
                energy = cc.CCSD(mf).run()
                return energy.e_tot
            else:
                raise ValueError("Methodology not implemented")
        except (UserWarning, np.linalg.LinAlgError):
            print("*** Exception in SCF Calculation \n")
            return float("inf")

    def metropolis(self):
        """
        Metroplois algorithm (symmetric proposal distribution).
        If structure is accepted, it will add into the list of
        lists self.store_structures

        """
        if self.energy_current < self.energy_before:
            print("New configuration ACCEPTED: lower energy")
            self.energy_before = self.energy_current
            self.store_structure()
        else:
            RE = self.energy_current / self.energy_before
            if np.random.random(1)[0] <= RE:
                print("New configuration ACCEPTED: Metropolis criterion")
                self.energy_before = self.energy_current
                self.store_structure()

    @build_input_gto_pyscf
    def pyscf(self, x):
        """
        Calculate of electronic energy with pyscf

        Parameters
        ----------
        x : array 1D
            Possible new positions and angles

        Returns
        -------
        Electronic energy

        """
        # ------------------------------------------------------
        # Build input of gto pyscf
        mol = gto.M(
            atom=self.input_gto_pyscf,
            basis=self._basis_set,
            verbose=False,
        )

        # ------------------------------------------------------
        # Calculate electronic energy with pyscf
        self.energy_current = self.calculate_electronic_energy(mol)

        if self._search_type != "ASCEC":
            # --------------------------------------------------
            # Metroplis
            self.metropolis()

        return self.energy_current
