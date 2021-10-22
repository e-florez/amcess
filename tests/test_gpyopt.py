import pytest
import sys
sys.path.append("../src")

import src.m_GPyOpt
import src.math_cost_function

def test_gpyopt():
    bounds = [(2,4), (1,3)]
    initer = 10
    maxiter = 40
    seed = 666
    gp_params = {'initer': initer}
    opt = src.m_GPyOpt.solve_gaussian_processes(src.math_cost_function.Himmelblau, bounds, maxiter, seed, gp_params)

    assert opt.fun < 1e-1