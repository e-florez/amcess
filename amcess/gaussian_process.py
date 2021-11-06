import GPyOpt
import numpy as np
from scipy.optimize import OptimizeResult


def GPyOpt_formatted_bounds(bounds):
    """
    Create dictionary with GPyOpt format for bounds
    """
    return [dict(zip(['name', 'type', 'domain', 'dimensionality'],
                     ['x' + str(i), 'continuous', bound, 1]))
            for i, bound in enumerate(bounds)]


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
    default_runopt = {'save_models_parameters': False,
                      'evaluations_file': None,
                      'models_file': None}

    runopt_args = {}
    for key, value in default_runopt.items():
        if key in gp_params.keys():
            runopt_args[key] = gp_params[key]
        else:
            runopt_args[key] = value

    return runopt_args


def solve_gaussian_processes(func,
                             bounds,
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
    if 'initer' not in gpyopt_keys:
        gp_params['initer'] = 3*len(bounds)
    if 'maxiter' not in gpyopt_keys:
        gp_params['maxiter'] = 10*len(bounds)
    if 'initial_design' not in gpyopt_keys:
        gp_params['initial_design'] = 'latin'
    if 'optimize_restarts' not in gpyopt_keys:
        gp_params['optimize_restarts'] = 5
    if 'xi' not in gpyopt_keys:
        gp_params['xi'] = 0.001
    if 'MCMC' not in gpyopt_keys:
        gp_params['MCMC'] = None

    # Define run_optimization parameters
    run_optimization_args = define_run_optimization_args(gp_params)

    # Define search space
    xbounds = GPyOpt_formatted_bounds(bounds)
    space = GPyOpt.Design_space(
                        space=xbounds
                        )

    # Define function to be minimized
    objective = GPyOpt.core.task.SingleObjective(func)

    # Define initial evaluations (random seed must be fixed before)
    np.random.seed(seed)
    initial_design = GPyOpt.experiment_design.initial_design(
                        gp_params['initial_design'],
                        space,
                        gp_params['initer']
                        )

    # define model and acquisition function
    acquisition_opt = GPyOpt.optimization.AcquisitionOptimizer(
                        space
                        )
    if gp_params['MCMC']:
        model = GPyOpt.models.GPModel_MCMC(
                        exact_feval=True,
                        verbose=False
                        )
        acquisition = GPyOpt.acquisitions.AcquisitionEI_MCMC(
                        model,
                        space,
                        acquisition_opt
                        )
    else:
        model = GPyOpt.models.GPModel(
                        exact_feval=True,
                        optimize_restarts=gp_params['optimize_restarts'],
                        verbose=False,
                        ARD=True
                        )
        acquisition = GPyOpt.acquisitions.AcquisitionEI(
                        model,
                        space,
                        acquisition_opt,
                        jitter=gp_params['xi']
                        )
    evaluator = GPyOpt.core.evaluators.Sequential(
                        acquisition
                        )

    opt = GPyOpt.methods.ModularBayesianOptimization(
                        model,
                        space,
                        objective,
                        acquisition,
                        evaluator,
                        initial_design
                        )
    opt.run_optimization(
                        max_iter=gp_params['maxiter'],
                        **run_optimization_args
                        )

    # Setting the OptimizeResult values
    optimize_res = OptimizeResult()
    optimize_res.setdefault('X', None)
    optimize_res.setdefault('Y', None)
    optimize_res.setdefault('plot_acquisition', None)
    optimize_res.setdefault('plot_convergence]', None)
    optimize_res.success = True
    optimize_res.status = 0
    optimize_res.x = opt.x_opt
    optimize_res.fun = opt.fx_opt
    optimize_res.nfev = gp_params['maxiter']
    optimize_res.update({'X': opt.X,
                         'Y': opt.Y,
                         'plot_acquisition': opt.plot_acquisition,
                         'plot_convergence': opt.plot_convergence})

    return optimize_res
