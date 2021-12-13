import time

import GPyOpt
import numpy as np
from scipy.optimize import OptimizeResult


def GPyOpt_formatted_bounds(bounds):
    """
    Create dictionary with GPyOpt format for bounds
    """
    keys = ["name", "type", "domain", "dimensionality"]
    xbounds = list()
    for i, bound in enumerate(bounds):
        values = ["x" + str(i), "continuous", bound, 1]
        xbounds.append(dict(zip(keys, values)))
    return xbounds


def define_optimization_args(**kargs):
    """
    Define arguments for optimization. If no values are given by user,
    it returns default values.

    Parameters
    ----------
    kargs : dictionary
        parameters given by user

    Returns
    -------
    dictionary
        output with model parameters
    """
    default_opt = {
        "initial_design": "latin",
        "optimize_restarts": 5,
        "xi": 0.001,
        "MCMC": None,
    }
    opt_args = dict()
    for key, value in default_opt.items():
        opt_args[key] = value
        if key in kargs.keys():
            opt_args[key] = kargs[key]
    return opt_args


def define_run_optimization_args(**kargs):
    """
    Define arguments for run_optimization method. If no values are
    given by user, it returns default values.

    Parameters
    ----------
    kargs : dictionary
        parameters given by user

    Returns
    -------
    dictionary
        output with model parameters
    """
    default_runopt = {
        "save_models_parameters": False,
        "evaluations_file": None,
        "models_file": None,
    }
    runopt_args = {}
    for key, value in default_runopt.items():
        runopt_args[key] = value
        if key in kargs.keys():
            runopt_args[key] = kargs[key]
    return runopt_args


def define_run_parallel_optimization_args(**kargs):
    """
    Define arguments for parallel computation. If no values are given,
    it returns default values.

    Parameters
    ----------
    kargs : dictionary
        parameters given by user

    Returns
    -------
    dictionary
        output with parallel computation parameters
    """
    default_paropt = {"num_cores": 1}
    par_kargs = {}
    for key, value in default_paropt.items():
        par_kargs[key] = value
        if key in kargs.keys():
            par_kargs[key] = kargs[key]
    return par_kargs


# #! Classes copied from GPyOpt/core/objective.py to can pass the
# #! electronic energy object in the energy evaluation
class Objective(object):
    """
    General class to handle the objective function internally.
    """

    def evaluate(self, x):
        raise NotImplementedError("")


class SingleObjective_Edited(Objective):
    """
    Class to handle problems with one single objective function.

    Parameters
    ----------
    param func:
        objective function.
    param batch_size:
        size of the batches (default, 1)
    param num_cores:
        number of cores to use in the process of evaluating
        the objective (default, 1).
    param objective_name:
        name of the objective function.
    param batch_type:
        Type of batch used. Only 'synchronous' evaluations are possible
        at the moment.
    param space:
        Not in use.

    Notes
    -----
    the objective function should take 2-dimensional numpy arrays
    as input and outputs. Each row should contain a location (in the case of
    the inputs) or a function evaluation (in the case of the outputs).

    """

    def __init__(
        self,
        func,
        args,
        num_cores=1,
        objective_name="no_name",
        batch_type="synchronous",
        space=None,
    ):
        self.func = func
        self.args = args
        self.n_procs = num_cores
        self.num_evaluations = 0
        self.space = space
        self.objective_name = objective_name
        self.parallel_error = False

    def evaluate(self, x):
        """
        Performs the evaluation of the objective at x.
        """

        if self.n_procs == 1:
            f_evals, cost_evals = self._eval_func(x)
        else:
            if not self.parallel_error:
                print(
                    "Parallel computation not implemented.\n"
                    "Fall back to single process!"
                )
                self.parallel_error = True
            f_evals, cost_evals = self._eval_func(x)

        return f_evals, cost_evals

    def _eval_func(self, x):
        """
        Performs sequential evaluations of the function at x (single
        location or batch). The computing time of each evaluation is
        also provided.
        """
        cost_evals = []
        f_evals = np.empty(shape=[0, 1])

        for i in range(x.shape[0]):
            st_time = time.time()
            # rlt = self.func(np.atleast_2d(x[i]), self.args)
            rlt = self.func(np.ravel(x[i]), *self.args)
            f_evals = np.vstack([f_evals, rlt])
            cost_evals += [time.time() - st_time]
        return f_evals, cost_evals


def solve_gaussian_processes(
    func,
    bounds,
    object_system: object = None,
    basis_set: str = None,
    args: tuple = (),
    seed=None,
    initer: int = None,
    maxiter: int = None,
    **kargs
):
    """
    Find the global minimum of a function using Bayesian Optimization
    with Gaussian Processes [1].

    Parameters
    ----------
    func : callable
        The objective function to be minimized. Must be in the form
        `f(x, *args)`, where x is the argument in the form of a 1-D array and
        args is a tuple of any additional fixed parameters needed to
        completely specify the function
    bounds : sequence, shape (n, 2)
        Bounds for variables. (min, max) pairs for each element in x,
        defining bounds for the objective function parameter.
    args : tuple
        Basis set, Cluster object, and name output file.
    iseed : None, int
        If seed is None ...
    initer :
        Number of initial evaluations (prior)
    maxiter :
        Maximum number of iterations
    kargs : dict
        Dictionary with Gaussian process parameters.

    Notes
    -----
    For more specific parameters, see `GPyOpt`_ official documentation

    [1] Gaussian Processes for Machine Learning. C. E. Rasmussen and
        C. K. I. Williams. MIT Press, 2006.

    .. _GPyOpt:
        https://sheffieldml.github.io/GPyOpt/
    """

    # Check input parameters
    if not initer:
        raise ValueError("initer not defined in Bayesian Optimization")
    if not maxiter:
        raise ValueError("maxiter not defined in Bayesian Optimization")
    if seed is None:
        seed = np.random.seed(1)

    # Define optimization parameters
    opt_kargs = define_optimization_args(**kargs)
    run_kargs = define_run_optimization_args(**kargs)
    par_kargs = define_run_parallel_optimization_args(**kargs)
    kargs.update(opt_kargs)

    # Define search space
    xbounds = GPyOpt_formatted_bounds(bounds)
    space = GPyOpt.Design_space(space=xbounds)

    # Define initial evaluations (random seed must be fixed before)
    np.random.seed(seed)
    initial_design = GPyOpt.experiment_design.initial_design(
        kargs["initial_design"], space, initer
    )

    # define model and acquisition function
    acquisition_opt = GPyOpt.optimization.AcquisitionOptimizer(space)
    if kargs["MCMC"]:
        model = GPyOpt.models.GPModel_MCMC(exact_feval=True, verbose=False)
        acquisition = GPyOpt.acquisitions.AcquisitionEI_MCMC(
            model, space, acquisition_opt
        )
    else:
        model = GPyOpt.models.GPModel(
            exact_feval=True,
            optimize_restarts=kargs["optimize_restarts"],
            verbose=False,
            ARD=True,
        )
        acquisition = GPyOpt.acquisitions.AcquisitionEI(
            model, space, acquisition_opt, jitter=kargs["xi"]
        )

    # Define evaluation method and objective function
    objective = SingleObjective_Edited(func, args, **par_kargs)
    evaluator = GPyOpt.core.evaluators.Sequential(acquisition)

    opt = GPyOpt.methods.ModularBayesianOptimization(
        model, space, objective, acquisition, evaluator, initial_design
    )
    opt.run_optimization(max_iter=maxiter, **run_kargs)

    # Setting the OptimizeResult values
    optimize_res = OptimizeResult()
    optimize_res.setdefault("X", None)
    optimize_res.setdefault("Y", None)
    optimize_res.setdefault("plot_acquisition", None)
    optimize_res.setdefault("plot_convergence]", None)
    optimize_res.success = True
    optimize_res.status = 0
    optimize_res.x = opt.x_opt
    optimize_res.fun = opt.fx_opt
    optimize_res.nfev = maxiter
    optimize_res.update(
        {
            "X": opt.X,
            "Y": opt.Y,
            "plot_acquisition": opt.plot_acquisition,
            "plot_convergence": opt.plot_convergence,
        }
    )

    return optimize_res
