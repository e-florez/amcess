import time

import GPyOpt
import numpy as np
from GPyOpt.util.general import spawn
from scipy.optimize import OptimizeResult


def GPyOpt_formatted_bounds(bounds):
    """
    Create dictionary with GPyOpt format for bounds
    """
    return [
        dict(
            zip(
                ["name", "type", "domain", "dimensionality"],
                ["x" + str(i), "continuous", bound, 1],
            )
        )
        for i, bound in enumerate(bounds)
    ]


def define_run_optimization_args(gp_params):
    """
    Define arguments for run_optimization method. If no values are
    given by user, it returns default values.

    Parameters
    ----------
    gp_params : dictionary
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
        if key in gp_params.keys():
            runopt_args[key] = gp_params[key]
        else:
            runopt_args[key] = value

    return runopt_args


# #! Classes copied from GPyOpt/core/objective.py to can pass the
# #! electronic energy object in the energy evaluation
class Objective(object):
    """
    General class to handle the objective function internally.
    """

    def evaluate(self, x):
        raise NotImplementedError()


class SingleObjective_Edited(Objective):
    """
    Class to handle problems with one single objective function.

    param func: objective function.
    param batch_size: size of the batches (default, 1)
    param num_cores: number of cores to use in the process of evaluating
    the objective (default, 1).
    param objective_name: name of the objective function.
    param batch_type: Type of batch used. Only 'synchronous' evaluations
    are possible at the moment.
    param space: Not in use.

    .. Note:: the objective function should take 2-dimensional numpy arrays
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

    def evaluate(self, x):
        """
        Performs the evaluation of the objective at x.
        """

        if self.n_procs == 1:
            f_evals, cost_evals = self._eval_func(x)
        else:
            try:
                f_evals, cost_evals = self._syncronous_batch_evaluation(x)
            # #! it is add E722 to avoid the error by :: to tox.init for flake8
            except:
                if not hasattr(self, "parallel_error"):
                    print(
                        "Error in parallel computation.\n"
                        "Fall back to single process!"
                    )
                else:
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

    def _syncronous_batch_evaluation(self, x):
        """
        Evaluates the function a x, where x can be a single location
        or a batch. The evaluation is performed in parallel according
        to the number of accessible cores.
        """

        from multiprocessing import Pipe, Process

        # --- parallel evaluation of the function
        # #! it is add E203 to avoid the error by :: to tox.init for flake8
        divided_samples = [x[i :: self.n_procs] for i in range(self.n_procs)]
        pipe = [Pipe() for i in range(self.n_procs)]
        proc = [
            Process(target=spawn(self._eval_func), args=(c, k))
            for k, (p, c) in zip(divided_samples, pipe)
        ]
        [p.start() for p in proc]
        [p.join() for p in proc]

        # --- time of evaluation is set to constant (=1).
        # This is one of the hypothesis of synchronous batch methods.
        f_evals = np.zeros((x.shape[0], 1))
        cost_evals = np.ones((x.shape[0], 1))
        i = 0
        for (p, c) in pipe:
            f_evals[i :: self.n_procs] = p.recv()[0]  # throw away costs
            i += 1
        return f_evals, cost_evals


def solve_gaussian_processes(
    func,
    bounds,
    args: tuple,
    seed=None,
    gp_params={},
):
    """
    Find the global minimum of a function using Bayesian Optimization
    with Gaussian Processes [1].

            Example:
                lambda x: x ** 2
                lambda x: (x[0] ** 2 + x[1] - 11) ** 2

        Parameters
        ----------
        func : callable
            The objective function to be minimized. Must be in the form
            f(x, *args), where x is the argument in the form of a 1-D array and
            args is a tuple of any additional fixed parameters needed to
            completely specify the function

        bounds : sequence, shape (n, 2)
            Bounds for variables. (min, max) pairs for each element in x,
            defining bounds for the objective function parameter.

        args : tuple
            Basis set, Cluster object, and name output file.

        iseed : None, int

            If seed is None ...

        gp_params : dictionary

            Dictionary with Gaussian process parameters
            - initer : Number of initial evaluations (prior)
            - maxiter : Maximum number of iterations
            For more specific parameters, see [2]

        Returns
        -------
    [1] "Gaussian Processes for Machine Learning" C. E. Rasmussen and
    C. K. I. Williams. MIT Press, 2006.

    [2] GPyOpt official documentation: https://sheffieldml.github.io/GPyOpt/
    """

    if not bounds:
        raise Exception("\n\n *** ERROR: No 'bounds' founds\n")
    if seed is None:
        seed = np.random.seed(1)

    lu = list(zip(*bounds))
    lower, upper = np.array(lu[0]), np.array(lu[1])

    # Checking bounds are valid
    if (
        np.any(np.isinf(lower))
        or np.any(np.isinf(upper))
        or np.any(np.isnan(lower))
        or np.any(np.isnan(upper))
    ):
        raise ValueError("Some bounds values are inf values or nan values")
    # Checking that bounds are consistent
    if not np.all(lower < upper):
        raise ValueError("Bounds are not consistent min < max")
    # Checking that bounds are the same length
    if not len(lower) == len(upper):
        raise ValueError("Bounds do not have the same dimensions")

    # Check values in optimization_parameters:
    gpyopt_keys = gp_params.keys()
    if "initer" not in gpyopt_keys:
        gp_params["initer"] = 3 * len(bounds)
    if "maxiter" not in gpyopt_keys:
        gp_params["maxiter"] = 10 * len(bounds)
    if "initial_design" not in gpyopt_keys:
        gp_params["initial_design"] = "latin"
    if "optimize_restarts" not in gpyopt_keys:
        gp_params["optimize_restarts"] = 5
    if "xi" not in gpyopt_keys:
        gp_params["xi"] = 0.001
    if "MCMC" not in gpyopt_keys:
        gp_params["MCMC"] = None

    # Define run_optimization parameters
    run_optimization_args = define_run_optimization_args(gp_params)

    # Define search space
    xbounds = GPyOpt_formatted_bounds(bounds)
    space = GPyOpt.Design_space(space=xbounds)

    # Define function to be minimized
    # objective = GPyOpt.core.task.SingleObjective(func)

    objective = SingleObjective_Edited(func, args)

    # Define initial evaluations (random seed must be fixed before)
    np.random.seed(seed)
    initial_design = GPyOpt.experiment_design.initial_design(
        gp_params["initial_design"], space, gp_params["initer"]
    )

    # define model and acquisition function
    acquisition_opt = GPyOpt.optimization.AcquisitionOptimizer(space)
    if gp_params["MCMC"]:
        model = GPyOpt.models.GPModel_MCMC(exact_feval=True, verbose=False)
        acquisition = GPyOpt.acquisitions.AcquisitionEI_MCMC(
            model, space, acquisition_opt
        )
    else:
        model = GPyOpt.models.GPModel(
            exact_feval=True,
            optimize_restarts=gp_params["optimize_restarts"],
            verbose=False,
            ARD=True,
        )
        acquisition = GPyOpt.acquisitions.AcquisitionEI(
            model, space, acquisition_opt, jitter=gp_params["xi"]
        )
    evaluator = GPyOpt.core.evaluators.Sequential(acquisition)

    opt = GPyOpt.methods.ModularBayesianOptimization(
        model, space, objective, acquisition, evaluator, initial_design
    )
    opt.run_optimization(
        max_iter=gp_params["maxiter"], **run_optimization_args
    )

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
    optimize_res.nfev = gp_params["maxiter"]
    optimize_res.update(
        {
            "X": opt.X,
            "Y": opt.Y,
            "plot_acquisition": opt.plot_acquisition,
            "plot_convergence": opt.plot_convergence,
        }
    )

    return optimize_res
