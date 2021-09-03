import types
from scipy.optimize import dual_annealing


def solve_dual_annealing(function, bounds):
    """Find the global minimum of a function using Dual Annealing [1].

        Example:
            lambda x: x ** 2
            lambda x: (x[0] ** 2 + x[1] - 11) ** 2

    Parameters
    ----------
    function : callable
        The objective function to be minimized. Must be in the form f(x, *args),
        where x is the argument in the form of a 1-D array and args is a tuple
        of any additional fixed parameters needed to completely specify the function

    bounds : sequence, shape (n, 2)
        Bounds for variables. (min, max) pairs for each element in x, defining bounds
        for the objective function parameter.

    Returns
    -------
    Optimize Result
        The optimization result represented as a OptimizeResult object. Important
        attributes are: x the solution array, fun the value of the function at the
        solution, and message which describes the cause of the termination.
        See OptimizeResult for a description of other attributes.

    Exception
    ---------


    Note
    ----
        References [1]: Tsallis C. Possible generalization of Boltzmann-Gibbs
        statistics. Journal of Statistical Physics, 52, 479-487 (1998).

        Check `scipy.optimize.dual_annealing` official documentation
    """

    if not bounds:
        raise Exception(f"\n\n *** ERROR: No 'bounds' founds\n")

    if not isinstance(function, types.FunctionType):
        raise Exception(f"\n\n *** ERROR: No 'function' found\n")

    result = dual_annealing(
        function,
        bounds=bounds,
        maxiter=1000,
        initial_temp=2000,
    )

    return result.x
