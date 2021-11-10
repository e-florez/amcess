import numpy as np
from scipy import constants


class Ascec:
    """
    ASCEC algorithm
    """

    def __init__(
        self,
        call_function: callable,
        bounds,
        T0: float = 1000.0,
        nT: int = 100,
        dT: float = 0.1,
        maxCylce: int = 3000,
        args: tuple = (),
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
            args : tuple
                Any additional fixed parameters needed to completely specify
        """

        self._bounds = bounds
        self._call_function = call_function

        self._T0 = T0
        self._nT = nT
        self._dT = dT
        self._maxCylce = maxCylce

        self.args = args

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
            Electronic energy of the new configuration into
            the attribute self.electronic_e
        """

        self.energy_current = self._call_function(x, *self.args)

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

        return np.random.rand(n)

    def ascec_criterion(self, T):
        """[summary]
        Evaluate ASCEC criterion for acceptance

        Parameters
        ----------
        x : array, float
            Value of each coordinate in the 1D array
        e : float
            Value of the cost function
        T : float
            Annealing temperature
        args :
            Any additional fixed parameters needed to completely specify
            the objective function.
        """

        if self.energy_current < self.e_before:
            return True
        else:
            DE = self.energy_current - self.e_before
            TKb = T * constants.k  # Boltzmann constant [J/K]
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
                x = self.random_mov(
                    (self.args[0]._object_system_current.total_molecules - 1)
                    * 6
                )
                self.electronic_e(x)
                accept = self.ascec_criterion(T)
                if accept:
                    self.e_before = self.energy_current
                    count = self._maxCylce + 1
                    print("Accept in the Temperature ", T)
                count += 1
            T = T - T * self._dT
            iT += 1
