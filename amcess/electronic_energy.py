import numpy
from pyscf import gto, scf

from amcess.base_molecule import Cluster


class ElectronicEnergy:
    def __init__(
        self,
        object_system: object,
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

        self._sphere_center = sphere_center
        self._sphere_radius = sphere_radius

        self._basis_set = basis_set

        self._max_closeness = max_closeness
        self._move_seed = seed

    # ===============================================================
    # Decorators
    # ===============================================================
    def build_input_pyscf(func_energy):
        def new_input(self, x_random, *args, ncall=[0]):
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
            new_geom[0] = {"atoms": system_object.get_molecule(0).atoms}
            for i in range(system_object.total_molecules - 1):
                new_geom[i + 1] = {
                    "atoms": system_object.move_molecules(
                        i + 1,
                        (
                            x_random[i * 3],
                            x_random[i * 3 + 1],
                            x_random[i * 3 + 2],
                        ),
                        (
                            x_random[
                                (i + system_object.total_molecules - 1) * 3
                            ],
                            x_random[
                                (i + system_object.total_molecules - 1) * 3 + 1
                            ],
                            x_random[
                                (i + system_object.total_molecules - 1) * 3 + 2
                            ],
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

            # Build input to pyscf
            symbols = self._object_system_current.symbols
            input_mol = "'"
            for i in range(self._object_system_current.total_atoms):
                input_mol += str(symbols[i])
                for j in range(3):
                    input_mol += "  " + str(
                        self._object_system_current.coordinates[i][j]
                    )
                if i < self._object_system_current.total_atoms - 1:
                    input_mol += "; "
                else:
                    input_mol += " '"

            args = (args[0], args[1], input_mol, system_object)
            return func_energy(self, x_random, *args, ncall=[0])

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
    @build_input_pyscf
    def energy_hf_pyscf(self, x, *args, ncall=[0]):
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
            atom=args[2],
            basis=self._basis_set,
        )

        try:
            e = scf.HF(mol).kernel()
        except (UserWarning, numpy.linalg.LinAlgError):
            print("Error in pyscf")
            e = float("inf")

        new_object = args[3]
        args[1].write(str(new_object.total_atoms) + "\n")
        args[1].write("Energy: " + str(e) + "\n")
        l: int = 0
        for symbols in new_object.symbols:
            args[1].write(
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

        ncall[0] += 1  # count calls
        return e


def hf_pyscf(x, *args):
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
    return args[0].energy_hf_pyscf(x, *args, ncall=[0])
