import GPyOpt
import numpy as np
from scipy._lib._util import check_random_state
from scipy.optimize import OptimizeResult

def GPyOpt_formatted_bounds(bounds):
    """ Create dictionary with GPyOpt format for bounds """
    return [dict(zip(['name', 'type', 'domain', 'dimensionality'], 
                     ['x' + str(i), 'continuous', bound, 1])) for i, bound in enumerate(bounds)]

def solve_gaussian_processes(
    func, 
    bounds,
    NT,
    seed=None,
    gp_parameters={},
    ):
    """Find the global minimum of a function using Bayesian Optimization 
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

        seed : {None, int, numpy.random.Generator,

            numpy.random.RandomState}, optional
            If seed is None (or np.random), the numpy.random.RandomState
            singleton is used. If seed is an int, a new RandomState instance
            is used, seeded with seed. If seed is already a Generator or
            RandomState instance then that instance is used. Specify seed
            for repeatable minimizations. The random numbers generated
            with this seed only affect the visiting distribution function
            and new coordinates generation.
        
        NT : int, optional
            The maximum number of global search iterations. Default value
            is 1000.

        Returns
        -------
    [1] "Gaussian Processes for Machine Learning" C. E. Rasmussen and 
    C. K. I. Williams. MIT Press, 2006.

    Check https://sheffieldml.github.io/GPyOpt/ official documentation
"""

    if not bounds:
        raise Exception("\n\n *** ERROR: No 'bounds' founds\n")

    if not NT:
        NT = 10*len(bounds)

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
    bo_keys = gp_parameters.keys()
    if 'initer' not in bo_keys: gp_parameters['initer'] = 3*len(bounds)
    if 'initial_design' not in bo_keys: gp_parameters['initial_design'] = 'latin'
    if 'optimize_restarts' not in bo_keys: gp_parameters['optimize_restarts'] = 5
    if 'xi' not in bo_keys: gp_parameters['xi'] = 0.001
    if 'save_models_parameters' not in bo_keys: gp_parameters['save_models_parameters'] = False
    if 'evaluations_file' not in bo_keys: gp_parameters['evaluations_file'] = None
    if 'models_file' not in bo_keys: gp_parameters['models_file'] = None
    if 'MCMC' not in bo_keys: gp_parameters['MCMC'] = None
    # Initialization of random Generator for reproducible runs if seed provided
    rand_state = check_random_state(seed)

    # Define search space
    xbounds = GPyOpt_formatted_bounds(bounds)
    space = GPyOpt.Design_space(space = xbounds)
    acquisition_optimizer = GPyOpt.optimization.AcquisitionOptimizer(space)

    objective = GPyOpt.core.task.SingleObjective(func)
    initial_design = GPyOpt.experiment_design.initial_design(gp_parameters['initial_design'], space, gp_parameters['initer'])
    # define model and acquisition function
    if gp_parameters['MCMC']:
        model = GPyOpt.models.GPModel_MCMC(exact_feval=True, verbose=False)
        acquisition=GPyOpt.acquisitions.AcquisitionEI_MCMC(model, space, acquisition_optimizer)
    else:
        model = GPyOpt.models.GPModel(exact_feval=True,optimize_restarts=gp_parameters['optimize_restarts'], verbose=False, ARD=True)
        acquisition=GPyOpt.acquisitions.AcquisitionEI(model, space, acquisition_optimizer, jitter=gp_parameters['xi'])
    evaluator = GPyOpt.core.evaluators.Sequential(acquisition)

    # OptimizeResult object to be returned
    optimize_res = OptimizeResult()

    opt = GPyOpt.methods.ModularBayesianOptimization(model, space, objective, acquisition, evaluator, initial_design)
    opt.run_optimization(max_iter = NT, 
                         save_models_parameters=gp_parameters['save_models_parameters'], 
                         evaluations_file=gp_parameters['evaluations_file'], 
                         models_file=gp_parameters['models_file'])

    # Setting the OptimizeResult values
    optimize_res.success = True
    optimize_res.status = 0
    optimize_res.x = opt.x_opt
    optimize_res.fun = opt.fx_opt
    optimize_res.nfev = NT
    
    return optimize_res