import numpy as np
from scipy import constants


class ascec(object):
    """
    stores or starts variables about ascec
    """
    def __init__(self, ascec, *args):
        """
        Store the boolean B. When B is True is used the ASCEC criterio
        in the acceptance

        Parameters
        ----------
            B : Boolean
                True activate ASCEC criterio

        Returns
        -------
            B : Boolean
        """
        self.ascec = ascec

    def ascec_criterion(self, e_before, e_now, t):
        """[summary]
        Evaluate ASCEC criterion for acceptance

        Parameters
        ----------
        x : array, float
            Value of each coordinate in the 1D array
        e : float
            Value of the cost function
        t : float
            Annealing temperature
        args :
            Any additional fixed parameters needed to completely specify
            the objective function.
        """

        DE = e_before - e_now
        if DE < 0.0000000001:
            print("DE < 0", DE)
        else:
            TKb = t*constants.k  # Boltzmann constant [J/K]
            exp = np.exp(-DE*constants.h/TKb)
            if DE < exp:
                print("DE < Boltzmann Poblation ", DE)