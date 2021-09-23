import pytest

from .context import src
from src.dual_annealing import solve_dual_annealing


@pytest.mark.parametrize(
    "cost_function, bounds, expected_minima",
    [
        (
            lambda x: (x + 1) ** 2,
            [(-10, 10)],
            [0],
        ),
        (
            lambda x: x ** 2,
            [(-10, 10)],
            [0],
        ),
        (
            lambda x: (x[0] ** 2 + x[1] - 11) ** 2 + (x[0] + x[1] ** 2 - 7) ** 2,
            [(-5, 5), (-5, 5)],
            [0],
        ),
    ],
)
def test_solve_dual_annealing(cost_function, bounds, expected_minima):
    """
    Test for a function to find root using Dual Annealing procedure
    """

    points = solve_dual_annealing(cost_function, bounds)

    global_minima = cost_function(points)

    assert ((global_minima - expected_minima) < 1.0e-6).all()
