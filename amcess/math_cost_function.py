import numpy as np


def Rastrigin(x, T=None, *args):
    """
    Rastrigin function

    f(x) = An + Sum_i^n [x_i*x_i - Acos(2pix_i)]
    A = 10 y x e [-5.12,5.12]

    bounds
    -5.12,5.12

    Parameters
    ----------
        x : array, float
            Value of each coordinate in the 1D array
        T : float
            Annealing temperature
        args :
            Any additional fixed parameters needed to completely specify
            the objective function.

    Returns
    -------
        Value of cost function in x

    Note
    ----
    Reference:
        [1] Rastrigin, L. A. "Systems of extremal control." Mir, Moscow (1974).
    """
    return np.sum(x * x - 10 * np.cos(2 * np.pi * x)) + 10 * np.size(x)


def Himmelblau(x, T=None, *args):
    """
    Himmelblau function

    bounds
    -5.0, 5.0

    Roots
    Himmelbau(3.0,2.0)=0.0,
    Himmelbau( − 2.805118 , 3.131312 )   = 0.0
    Himmelbau( − 3.779310 , − 3.283186 ) = 0.0
    Himmelbau( 3.584428 , − 1.848126 )   = 0.0

    Parameters
    ----------
        x : array, float
            Value of each coordinate in the 1D array
        T : float
            Annealing temperature
        args :
            Any additional fixed parameters needed to completely specify
            the objective function.

    Returns
    -------
        Value of cost function in x

    Note
    ----
    Reference:
        [1] Himmelblau, D. (1972). Applied Nonlinear Programming. McGraw-Hill.
        ISBN 0-07-028921-2.
    """
    # !function is not good to analyze ascec criterio,
    # !because it's very uniform
    return (x[0] ** 2 + x[1] - 11) ** 2 + (x[0] + x[1] ** 2 - 7) ** 2
