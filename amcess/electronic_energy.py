import numpy as np
from pyscf import gto, scf

from amcess.base_molecule import Cluster


class ElectronicEnergy:
    def __init__(
        self,
        object_system: object,
        search_type: str,
        sphere_center: tuple,
        sphere_radius: float,
        basis_set: str,
        max_closeness: float = 1.0,
        seed: int = None,
    ) -> None:
        """
        Class to calculate electronic energy

        Attributes
        ----------
            molecule_object : object
                Object initialized with Molecule or Cluster class
            sphere_center : list
                Center of the sphere where evolve the system
            sphere_radius : float
                Radius of the sphere where should evolve the system
        """

        self._object_system_initial = object_system
        self._object_system_before = object_system
        self._object_system_current = object_system

        self._search_type = search_type
        self._sphere_center = sphere_center
        self._sphere_radius = sphere_radius

        self._basis_set = basis_set

        self._max_closeness = max_closeness
        self._move_seed = seed

        self.store_structures = []

        if self._search_type != "ASCEC ":
            mol = gto.M(
                atom=self.input_atom_mol_pyscf(),
                basis=self._basis_set,
                verbose=False,
            )
            self._e0 = self.calculate_electronic_e(mol)
            self.energy_before = self._e0

    # ===============================================================
    # Decorators
    # ===============================================================
    def build_input_pyscf(func_energy):
        def new_input(self, x):
            """
            Build input to pyscf

            Parameters
            ----------
                x : array 1D
                    possible new positions and angles.
                system_object : object
                    Object initialized with Molecule or Cluster class
                icall : integer
                    number of call

            Returns
            -------
                imput_mol: list
                    list of atoms and coordinates
                system_object: Cluster
                    Cluster objects
            """
            system_object = self._object_system_current
            # Rotate and translate
            new_geom = dict()
            # ! Martin
            # ! tener cuidado como indexo el diccionario
            new_geom[0] = {"atoms": system_object.get_molecule(0).atoms}
            for i in range(system_object.total_molecules - 1):
                new_geom[i + 1] = {
                    "atoms": system_object.move_molecules(
                        i + 1,
                        (
                            x[i * 3],
                            x[i * 3 + 1],
                            x[i * 3 + 2],
                        ),
                        (
                            x[(i + system_object.total_molecules - 1) * 3],
                            x[(i + system_object.total_molecules - 1) * 3 + 1],
                            x[(i + system_object.total_molecules - 1) * 3 + 2],
                        ),
                        self._max_closeness,
                        self._move_seed,
                    )
                    .get_molecule(i + 1)
                    .atoms
                }

            self.object_system_current = Cluster(
                *new_geom.values(),
                sphere_radius=self._sphere_radius,
                sphere_center=self._sphere_center
            )

            self.input_atom_mol_pyscf()

            return func_energy(self, x)

        return new_input

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
        Write the current system to a file
        """
        # Build input to pyscf
        symbols = self._object_system_current.symbols
        self.input_mol = "'"
        for i in range(self._object_system_current.total_atoms):
            self.input_mol += str(symbols[i])
            for j in range(3):
                self.input_mol += "  " + str(
                    self._object_system_current.coordinates[i][j]
                )
            if i < self._object_system_current.total_atoms - 1:
                self.input_mol += "; "
            else:
                self.input_mol += " '"
        return self.input_mol

    def write_to_file(self):
        """
        Write the current system to a file

        Parameters
        ----------
            filename: str
                File name where is save structure and energy
        """
        n_atoms = len(self.store_structures[0]) - 1
        with open("configurations.txt", "w") as f:
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
        pass

    def calculate_electronic_e(self, mol):
        """
        Calculate electronic energy

        Parameters
        ----------
            mol: object
                pyscf object

        Returns
        -------
            e: float
                Electronic energy
        """
        try:
            return scf.HF(mol).kernel()
        except (UserWarning, np.linalg.LinAlgError):
            print("Error in pyscf")  # usar un warning
            return float("inf")

    def metropolis(self):
        """[summary]
        Metroplois algorithm (symmetric proposal distribution).
        If structure is accepted, it will add into file xyz

        Parameters
        ----------
            filename: str
                File name where is save structure and energy
        """
        if self.energy_current < self.energy_before:
            print("Accept 1")
            self.energy_before = self.energy_current
            self.store_structure
        else:
            RE = self.energy_current / self.energy_before
            if np.random.random(1)[0] <= RE:
                print("Accept 2")
                self.energy_before = self.energy_current
                self.store_structure

    @build_input_pyscf
    def hf_pyscf(self, x):
        """
        Calculate of electronic energy with pyscf

        Parameters
        ----------
            x : array 1D
                Possible new positions and angles
            args : list
                basis set, Cluster Object, name output xyz

        Returns
        -------
            e : float
                electronic energy

        """
        # Build input to pyscf
        mol = gto.M(
            atom=self.input_mol,
            basis=self._basis_set,
            verbose=False,
        )

        # Calculate electronic energy
        self.energy_current = self.calculate_electronic_e(mol)

        if self._search_type != "ASCEC":
            # Metroplis
            self.metropolis()

        return self.energy_current


# ! Martin
# ! hacer la visualización desde el archivo
# ! evaluar lo métodos para almacenar los archivos
