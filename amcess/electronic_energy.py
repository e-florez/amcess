import numpy as np
from pyscf import cc, dft, gto, mp, scf

from amcess.cluster import Cluster


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

        self._initial_system = object_system
        self._before_system = object_system
        self._current_system = object_system

        self._search_type = search_type

        self._method = methodology.split()[0]
        if len(methodology.split()) > 1:
            self._functional = methodology.split()[1]
            
        self._basis_set = basis_set

        self._move_seed = seed

        self.store_structures = []

        if self._search_type != "ASCEC ":
            mol = gto.M(
                atom=self.InputMolPyscf(),
                basis=self._basis_set,
                verbose=False,
            )
            self._e0 = self.SCF(mol)  # calculate_electronic_energy(mol)
            self.energy_before = self._e0

    # ===============================================================
    # Decorators
    # ===============================================================
    def BuildInputPyscf(func_energy):
        def NewInputPyscf(self, x):
            """
            Build input to pyscf

            .. rubric:: Parameters

            x : array 1D
                possible new positions and angles.
            system_object : object Cluster
                Object initialized with Molecule or Cluster class

            .. rubric:: Returns

            input_gto_pyscf: list
                Atom's symbols and coordinates
            system_object: Cluster
                Cluster objects
            """
            system_object = self.GetCurrentSystem()

            # ------------------------------------------------------------
            # Rotate and translate each molecule into object Cluster
            new_geom = dict()
            new_geom[0] = {"atoms": system_object.GetMol(0).GetAtoms()}
            for i in range(system_object.GetNumMols() - 1):
                new_geom[i + 1] = {
                    "atoms": system_object.TranslateMol(
                        i + 1,
                        x=x[i * 3],
                        y=x[i * 3 + 1],
                        z=x[i * 3 + 2],
                    )
                    .RotateMol(
                        i + 1,
                        x=x[(i + system_object.GetNumMols() - 1) * 3],
                        y=x[(i + system_object.GetNumMols() - 1) * 3 + 1],
                        z=x[(i + system_object.GetNumMols() - 1) * 3 + 2],
                    )
                    .GetMol(i + 1)
                    .GetAtoms()
                }

            # ------------------------------------------------------------
            # New object Cluster with new geometries
            self.current_system = Cluster(
                *new_geom.values(),
                sphere_radius=system_object.GetSphereR(),
                sphere_center=system_object.GetSphereCenter()
            )
            # ------------------------------------------------------------
            # Build input to pyscf
            self.InputMolPyscf()

            return func_energy(self, x)

        return NewInputPyscf

    # ===============================================================
    # PROPERTIES
    # ===============================================================
    def GetInitialSystem(self):
        return self._initial_system

    def GetCurrentSystem(self):
        return self._current_system

    def GetBeforeSystem(self):
        return self._before_system

    def SetInitialSystem(self, new_initial_system):
        (
            self._initial_system,
            self._before_system,
            self._current_system,
        ) = (
            new_initial_system,
            new_initial_system,
            new_initial_system,
        )

    def SetCurrentSystem(self, new_current_system):
        self.SetBeforeSystem(self._current_system)
        self._current_system = new_current_system
    
    def SetBeforeSystem(self, new_before_system):
        self._before_system = new_before_system

    # ===============================================================
    # Methods
    # ===============================================================
    def InputMolPyscf(self):
        """
        Build a portion of the input for the gto object of pyscf
            'X 0.0 0.0 0.0; X 0.0 0.0 1.0'

        .. rubric:: Returns

        input_gto_pyscf: list
            Atom's symbols and coordinates
        """

        symbols = self.GetCurrentSystem().GetAtomicSymbols()
        self.input_gto_pyscf = "'"
        for i in range(self.GetCurrentSystem().GetTotalAtoms()):
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

    def WriteOutPut(self, filename):
        """
        Write all accepted structures to a file

        .. rubric:: Parameters

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

        .. rubric:: Parameters

        mol: object
            gto pyscf object

        .. rubric:: Returns

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

    @BuildInputPyscf
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
