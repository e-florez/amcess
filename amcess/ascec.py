import numpy as np

from amcess.electronic_energy import ElectronicEnergy


class Ascec(ElectronicEnergy):
    """
    reference:
    ASCEC algorithm
    J. F. Pérez, C. Z. Hadad, A. Restrepo.
    Int. J. Quantum Chem. 2008, 108, 1653--1659
    """

    def __init__(
        self,
        object_system: object,
        search_type: str,
        methodology: str,
        basis_set: str,
        program: str,
        bounds: list,
        seed: int = 1234567890,
        initial_temperature: float = 1000.0,
        number_temperatures: int = 100,
        step_temperature: float = 0.1,
        maxCycles: int = 3000,
    ) -> None:
        """
        Initialize the Ascec class

        Parameters
        ----------
        call_function : callable
            The function to calculate electronic energy
        bounds : array, float
            The bounds of the search space
        initial_temperature : float
            Initial temperature
        number_temperatures : int
            Number of temperature steps
        step_temperature : float
            Temperature step
        maxCylces : int
            Maximum number of cycles
        """
        super().__init__(
            object_system,
            search_type,
            methodology,
            basis_set,
            seed,
        )
        self._bounds = bounds
        self._call_program = program

        self._initial_temperature = initial_temperature
        self._number_temperature = number_temperatures
        self._step_temperature = step_temperature
        self._maxCylces = maxCycles

        # initial energy
        self.electronic_e(np.zeros(len(bounds)))
        self._e0 = self.energy_current
        self.e_before = self._e0

    # ===============================================================
    # Methods
    # ===============================================================
    def electronic_e(self, x):
        """Evaluate the electronic energy

        .. rubric:: Parameters

        x : array, float
            Value to move the molecules, in the 1D array

        .. rubric:: Returns

        energy :
            Electronic energy of the new configuration
        """
        if self._call_program == "pyscf":
            self.energy_current = self.pyscf(x)

    def random_mov(self, n):
        """
        Randomly move the molecules

        .. rubric:: Parameters

        n : int
            dimension of the 1D array

        .. rubric:: Returns

        x : array, float
            Random value to move the molecules, in the 1D array
        """
        translate = np.random.uniform(
            low=self._bounds[0][0],
            high=-self._bounds[0][0],
            size=(int(n / 2),),
        )
        rotate = np.random.uniform(low=-180.0, high=180.0, size=(int(n / 2),))
        return np.concatenate((translate, rotate))

    def ascec_criterion(self, temperature):
        """
        ASCEC criterion for acceptance, based in Markov Chain Monte Carlo

        .. rubric:: Parameters

        temperature : float
            Annealing temperature

        .. rubric:: Returns

        accepted : boolean
            True if the configuration is accepted
        lower_energy : boolean
            True if the new energy is lower than the previous one
        """

        KB: float = 3.166811563e-6  # Eh/K
        accepted = False
        lower_energy = False
        # ------------------------------------------------------------
        # 1) Accepted: if the new energy is lower than the previous one
        if self.energy_current < self.e_before:
            accepted = True
            lower_energy = True
        else:
            DE = self.energy_current - self.e_before
            TKb = temperature * KB
            boltzmann = np.exp(-DE / TKb)
            # 2) Accepted: if DE is lower than the Boltzmann distribution
            if boltzmann > np.abs(DE):
                accepted = True

        return accepted, lower_energy

    def RunASCEC(self):
        """
        Run ASCEC algoritm.
        """
        iT = 0
        temperature = self._initial_temperature
        configurations_accepted = 0
        while iT <= self._number_temperature:
            count = 0
            while count <= self._maxCylces:
                # --------------------------------------------------------------
                # Information about the before cycle
                print(
                    f"\r Current temperature {temperature:7.2f} K, progress:"
                    f" {100*iT/self._number_temperature:.2f}%, with "
                    f" {configurations_accepted:3d} configurations accepted"
                    f" (cycle {count:>4d}/{self._maxCylces:<4d})",
                    end="",
                )
                # ------------------------------------------------------------
                # Generate 3 random values to translate and other 3 to rotate
                x = self.random_mov(len(self._bounds))
                # ------------------------------------------------------------
                # Electronic energy calculation
                self.electronic_e(x)
                # ------------------------------------------------------------
                # ASCEC criterion
                accepted, lower_energy = self.ascec_criterion(temperature)
                if accepted:
                    # --------------------------------------------------------
                    # -- Counter of accepted configurations
                    configurations_accepted += 1
                    # -- Store the new energy
                    self.e_before = self.energy_current
                    # -- Store the new configuration
                    self.store_structure()
                    # --------------------------------------------------------
                    # Skip to the next temperature
                    if lower_energy:
                        count = float("inf")
                # ------------------------------------------------------------
                # Counter of cycles
                count += 1
            # ------------------------------------------------------------
            # -- Update the temperature
            temperature = temperature - temperature * self._step_temperature
            # -- Counter of temperature steps
            iT += 1
