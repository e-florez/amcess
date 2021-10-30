import pytest
import numpy as np
from amcess.minimization import solve_gaussian_processes


@pytest.mark.parametrize(
    "cost_function, bounds, initer, maxiter, expected_minima",
    [
        (
            lambda x: (x + 1) ** 2,
            [(-10, 10)],
            3,
            10,
            [0]
        ),
        (
            lambda x: x ** 2,
            [(-10, 10)],
            3,
            10,
            [0],
        ),
        # (
        #     lambda x: (x[0,0] ** 2 + x[0,1] - 11) ** 2
        #     + (x[0,0] + x[0,1] ** 2 - 7) ** 2,
        #     [(-5, 5), (-5, 5)],
        #     30,
        #     50,
        #     [0],
        # ),
    ],
)
def test_solve_gaussian_processes(cost_function,
                                  bounds,
                                  initer,
                                  maxiter,
                                  expected_minima):
    """Test for a function to find root using Gaussian Processes
    """
    seed = 666
    gp_params = {'initer': initer, 'maxiter': maxiter}
    opt = solve_gaussian_processes(cost_function,
                                   bounds,
                                   seed=seed,
                                   gp_params=gp_params)
    assert (opt.fun - expected_minima) < 1e-6
