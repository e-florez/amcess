import types

from scipy.optimize import dual_annealing


def solve_dual_annealing(function, bounds):
    """
    Find the global minimum of a function using Dual Annealing [1].

    Example
    -------
    >>> print(solve_dual_annealing(
                function=lambda x: (x[0] ** 2 + x[1] - 11) ** 2,
                bounds=[(-5, 5), (-5, 5)],
                )
            )
    [0.0 -0.2]

    Parameters
    ----------
    function : callable
        The objective function to be minimized.

    bounds : sequence, shape (n, 2) -> n: target function dimension
        A (min, max) pairs for each element in x, defining bounds
        for the objective function parameter.

    Returns
    -------
    Optimize Result
        The optimization result represented as a OptimizeResult object.
        Important attributes are: x the solution array, fun the value of the
        function at the solution, and message which describes the cause of
        the termination.
        See OptimizeResult for a description of other attributes.

    Raises
    ------
    Exception
        if not ``function`` are set for passed in as parameters




    References
    ----------
    .. [1]: Tsallis C. Possible generalization of Boltzmann-Gibbs
    statistics. Journal of Statistical Physics, 52, 479-487 (1998).

    Check ``scipy.optimize.dual_annealing`` official documentation
    """

    try:
        assert isinstance(
            function, types.FunctionType
        ), f"\n\n *** ERROR: function: '{function}' has not FunctionType\n"
    except AssertionError as error:
        print(error)
        return

    result = dual_annealing(
        function,
        bounds=bounds,
        maxiter=1000,
        initial_temp=2000,
    )

    return result.x


# fun = lambda x: (x[0] ** 2 + x[1] - 11) ** 2 + (x[0] + x[1] ** 2 - 7) ** 2
# bounds = [(-5, 5), (-5, 5)]


# a = solve_dual_annealing(fun, bounds)
