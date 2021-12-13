import numpy as np
from amcess.electronic_energy import ElectronicEnergy

from copy import deepcopy


class Ascec(ElectronicEnergy):
    """
    ASCEC algorithm
    """

    def __init__(
        self,
        object_system: object,
        search_type: str,
        methodology: str,
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
            methodology,
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
        energy :
            Electronic energy of the new configuration
        """
        if self._call_function == 1:
            self.energy_current = self.pyscf(x)

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
        translate = np.random.uniform(
            low=-self._object_system_current._sphere_radius,
            high=self._object_system_current._sphere_radius,
            size=(int(n / 2),),
        )
        rotate = np.random.uniform(low=-180.0, high=180.0, size=(int(n / 2),))
        return np.concatenate((translate, rotate))

    def ascec_criterion(self, T):
        """
        ASCEC criterion for acceptance, based in Markov Chain Monte Carlo

        Parameters
        ----------
        e : float
            Value of the cost function

        T : float
            Annealing temperature
        """
        KB: float = 3.166811563e-6  # Eh/K
        accepted = False
        lower_energy = False
        if self.energy_current < self.e_before:
            accepted = True
            lower_energy = True
        else:
            DE = self.energy_current - self.e_before
            TKb = T * KB  # Boltzmann constant [Eh/K]
            boltzmann = np.exp(-DE / TKb)
            delta_energy = np.abs(DE / self.energy_current)

            # random = np.random.uniform(0, 1)
            # if boltzmann > random:
            if boltzmann > delta_energy:
                # print(f"Boltzmann Accepted {boltzmann:.3e}")
                accepted = True

        return accepted, lower_energy

    def ascec_run(self):
        """
        Run ASCEC algoritm
        """

        iT = 0
        T = self._T0
        configurations_accepted = 0
        while iT <= self._nT:
            count = 0
            while count <= self._maxCylce:

                print(
                    f"\r Current temperature {T:7.2f} K, progress:"
                    f" {100*iT/self._nT:.2f}%, with "
                    f" {configurations_accepted:3d} configurations accepted"
                    f" (cycle {count:>4d}/{self._maxCylce:<4d})",
                    end="",
                )

                # 3 values to translate and another 3 to rotate
                x = self.random_mov(len(self._bounds))
                self.electronic_e(x)
                # self.store_structure()
                accepted, lower_energy = self.ascec_criterion(T)
                if accepted:
                    configurations_accepted += 1
                    self.e_before = self.energy_current
                    self.store_structure()
                    # print("Configuration Temperature ", T)

                    if lower_energy:
                        count = float("inf")

                count += 1
            T = T - T * self._dT
            iT += 1
