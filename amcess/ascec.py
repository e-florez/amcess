import numpy as np
from scipy import constants
from amcess.electronic_energy import ElectronicEnergy


class Ascec(ElectronicEnergy):
    """
    ASCEC algorithm
    """

    def __init__(
        self,
        object_system: object,
        search_type: str,
        sphere_center: tuple,
        sphere_radius: float,
        basis_set: str,
        call_function: int,
        bounds: list,
        max_closeness: float = 1.0,
        seed: int = None,
        T0: float = 1000.0,
        nT: int = 100,
        dT: float = 0.1,
        maxCycle: int = 3000,
    ):
        """
        Initialize the Ascec class

        Parameters
        ----------
            call_function : callable
                The function to calculate electronic energy
            bounds : array, float
                The bounds of the search space
            T0 : float
                Initial temperature
            nT : int
                Number of temperature steps
            dT : float
                Temperature step
            maxCylce : int
                Maximum number of cycles
        """
        super().__init__(
            object_system,
            search_type,
            sphere_center,
            sphere_radius,
            basis_set,
            max_closeness,
            seed,
        )
        self._object_system_current = object_system
        self._bounds = bounds
        self._call_function = call_function

        self._T0 = T0
        self._nT = nT
        self._dT = dT
        self._maxCylce = maxCycle

        self.store_structures = []

        # initial energy
        self.electronic_e(np.zeros(len(bounds)))
        self._e0 = self.energy_current
        self.e_before = self._e0

    # ===============================================================
    # Methods
    # ===============================================================
    def electronic_e(self, x):
        """
        Evaluate the electronic energy

        Parameters
        ----------
            function : callable
                The function to calculate electronic energy
            x : array, float
                Value to move the molecules, in the 1D array

        Returns
        -------
            Electronic energy of the new configuration in
            the attribute self.electronic_e
        """
        if self._call_function == 1:
            self.energy_current = self.hf_pyscf(x)

    def random_mov(self, n):
        """
        Randomly move the molecules

        Parameters
        ----------
            n : int
                dimension of the 1D array

        Returns
        -------
            x : array, float
                Random value to move the molecules, in the 1D array
        """
        # ! Editar para que el rotar y la translacion esten asociados con los bounds
        return np.random.rand(n)

    def ascec_criterion(self, T):
        """[summary]
        ASCEC criterion for acceptance, based in Markov Chain Monte Carlo

        Parameters
        ----------
        x : array, float
            Value of each coordinate in the 1D array
        e : float
            Value of the cost function
        T : float
            Annealing temperature
        """
        KB: float = 3.166811563e-6  # Eh/K
        if self.energy_current < self.e_before:
            return True
        else:
            DE = (
                np.abs(self.energy_current - self.e_before)
                / self.energy_current
            )
            # TKb = T * constants.k  # Boltzmann constant [J/K]
            TKb = T * KB  # Boltzmann constant [Eh/K]
            exp = np.exp(-DE / TKb)
            if DE < exp:
                print("DE < Boltzmann Poblation ", DE)
                return True

    def ascec_run(self):
        """
        Run ASCEC algoritm
        """
        iT = 0
        T = self._T0
        while iT <= self._nT:
            count = 0
            while count <= self._maxCylce:
                # 3 values to translate and another 3 to rotate
                x = self.random_mov(len(self._bounds))
                self.electronic_e(x)
                self.store_structure()
                accept = self.ascec_criterion(T)
                if accept:
                    self.e_before = self.energy_current
                    count = self._maxCylce + 1
                    self.store_structure()
                    print("Accept in the Temperature ", T)
                count += 1
            T = T - T * self._dT
            iT += 1
