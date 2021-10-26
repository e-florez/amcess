import pytest
import sys
sys.path.append("../src")

import numpy as np
import src.m_GPyOpt
import src.math_cost_function

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

def test_solve_gaussian_processes(cost_function, bounds, initer, maxiter, expected_minima):
    """Test for a function to find root using Gaussian Processes
    """
    seed = np.random.seed(666)
    gp_params = {'initer': initer}
    opt = src.m_GPyOpt.solve_gaussian_processes(cost_function, bounds, maxiter, seed, gp_params)
    assert (opt.fun - expected_minima) < 1e-6

