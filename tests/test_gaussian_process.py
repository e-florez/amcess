import math

import pytest

from amcess.cluster import Cluster
from amcess.search_engine import SearchConfig
from amcess.gaussian_process import GPyOpt_formatted_bounds
from amcess.gaussian_process import define_optimization_args
from amcess.gaussian_process import define_run_optimization_args
from amcess.gaussian_process import define_run_parallel_optimization_args
from amcess.gaussian_process import Objective


# Bayesian functions tests
@pytest.mark.parametrize(
    "bounds",
    [
        [
            (-2.3260049557867735, 2.3260049557867735),
            (-2.3260049557867735, 2.3260049557867735),
            (-2.3260049557867735, 2.3260049557867735),
            (0, math.pi),
            (0, math.pi),
            (0, math.pi),
        ],
    ],
)
def test_Bayesian_bounds_formatter(bounds):
    """
    Test Bayesian bounds formatter
    """
    keys = ("name", "type", "domain", "dimensionality")
    xbounds = GPyOpt_formatted_bounds(bounds)
    assert all(k in xbounds[0] for k in keys)


@pytest.mark.parametrize(
    "kargs",
    [
        {},
    ],
)
def test_Bayesian_optimization_args(kargs):
    """
    Test Bayesian minimal optimization parameters
    """
    keys = ("initial_design", "optimize_restarts", "xi", "MCMC")
    opt_keys = define_optimization_args(**kargs)
    assert all(k in opt_keys for k in keys)


@pytest.mark.parametrize(
    "kargs",
    [
        {"initial_design": "random"},
    ],
)
def test_Bayesian_optimization_args_with_input(kargs):
    """
    Test Bayesian minimal optimization parameters with input
    parameters
    """
    keys = ("initial_design", "optimize_restarts", "xi", "MCMC")
    opt_keys = define_optimization_args(**kargs)
    assert all(k in opt_keys for k in keys)


@pytest.mark.parametrize(
    "kargs",
    [
        {},
    ],
)
def test_Bayesian_run_optimization_args(kargs):
    """
    Test Bayesian minimal run optimization parameters
    """
    keys = ("save_models_parameters", "evaluations_file", "models_file")
    opt_keys = define_run_optimization_args(**kargs)
    assert all(k in opt_keys for k in keys)


@pytest.mark.parametrize(
    "kargs",
    [
        {"save_models_parameters": False},
    ],
)
def test_Bayesian_run_optimization_args_with_input(kargs):
    """
    Test Bayesian minimal run optimization parameters with input
    parameters
    """
    keys = ("save_models_parameters", "evaluations_file", "models_file")
    opt_keys = define_run_optimization_args(**kargs)
    assert all(k in opt_keys for k in keys)


@pytest.mark.parametrize(
    "num_cores",
    [
        (2),
    ],
)
def test_Bayesian_run_parallel_optimization_args(num_cores):
    """
    Test Bayesian parallel optimization parameters
    """
    key = "num_cores"
    opt_keys = define_run_parallel_optimization_args(num_cores=num_cores)
    assert key in opt_keys


@pytest.mark.parametrize(
    "molecule1, molecule2, search_meth, initer",
    [
        (
            [("H", 0.0, 0.0, 0.0), ("H", 0.78, 0.0, 0.0)],
            [("H", 0.0, 0.0, 1.0), ("H", 0.78, 0.0, 1.0)],
            "Bayesian",
            3,
        ),
    ],
)
def test_Bayesian_input_maxeval(molecule1, molecule2, search_meth, initer):
    """
    Test ValueError for maxeval
    """
    with pytest.raises(ValueError) as e:
        SearchConfig(Cluster(molecule1, molecule2), search_methodology=search_meth).run(
            initer=3
        )
    assert str(e.value) == "maxiter not defined in Bayesian Optimization"


@pytest.mark.parametrize(
    "molecule1, molecule2, search_meth, maxiter",
    [
        (
            [("H", 0.0, 0.0, 0.0), ("H", 0.78, 0.0, 0.0)],
            [("H", 0.0, 0.0, 1.0), ("H", 0.78, 0.0, 1.0)],
            "Bayesian",
            3,
        ),
    ],
)
def test_Bayesian_input_initer(molecule1, molecule2, search_meth, maxiter):
    """
    Test ValueError for maxeval
    """
    with pytest.raises(ValueError) as e:
        SearchConfig(Cluster(molecule1, molecule2), search_methodology=search_meth).run(
            maxiter=3
        )
    assert str(e.value) == "initer not defined in Bayesian Optimization"


@pytest.mark.parametrize(
    "x",
    [
        (0),
    ],
)
def test_Bayesian_Objective(x):
    """
    Test NotImplementedError in Objective class
    """
    with pytest.raises(NotImplementedError) as e:
        Objective().evaluate(x=x)
    assert str(e.value) == ""
