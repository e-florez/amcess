from src.base_molecule import Cluster
import types
# ipdb para debugger python
# import ipdb; ipdb.set_trace()
# Object into _dual_annealing.py
from scipy.optimize._dual_annealing import ObjectiveFunWrapper
from scipy.optimize._dual_annealing import LocalSearchWrapper
from scipy.optimize._dual_annealing import EnergyState
from scipy.optimize._dual_annealing import VisitingDistribution
from scipy.optimize._dual_annealing import StrategyChain

from scipy.optimize import OptimizeResult
from scipy.optimize import minimize
from scipy.special import gammaln
from scipy._lib._util import check_random_state

from scipy import constants

import numpy as np  # 1.21.2
from math_cost_function import *
from ascec_criterion import *
from heisenberg import *
#from base_molecule import Cluster
from move_fragments import *
# pytest 6.2.5

def solve_dual_annealing(func, bounds, ascec_activation=True,
                         there_is_molecule = None,
                         system_object = None,
                         seed=None,
                         NT=1000, T0=5230.0, dT=2e-5, mxcycle=10000000.0,
                         local_search_options={},
                         no_local_search=False, visit_regions=2.62,
                         accept=-5.0,
                         x0=None, args=(), callback=None):
    """Find the global minimum of a function using Dual Annealing [1].

        Example:
            lambda x: x ** 2
            lambda x: (x[0] ** 2 + x[1] - 11) ** 2


    This funtions is a copy tthe function dual_annealing to have more control
    about variables #!(Acá debería ir algo sobre la licencia, porqué se copío
    #! las lineas de la función dual_annealing)

    Parameters
    ----------
    func : callable
        The objective function to be minimized. Must be in the form
        f(x, *args), where x is the argument in the form of a 1-D array and
        args is a tuple of any additional fixed parameters needed to
        completely specify the function

    bounds : sequence, shape (n, 2)
        Bounds for variables. (min, max) pairs for each element in x, defining
        bounds for the objective function parameter.

    seed : {None, int, numpy.random.Generator,

        numpy.random.RandomState}, optional
        If seed is None (or np.random), the numpy.random.RandomState singleton
        is used. If seed is an int, a new RandomState instance is used, seeded
        with seed. If seed is already a Generator or RandomState instance then
        that instance is used. Specify seed for repeatable minimizations. The
        random numbers generated with this seed only affect the visiting
        distribution function and new coordinates generation.

    NT : int, optional
        The maximum number of global search iterations. Default value is 1000.

    T0 : float, optional
        The initial temperature, use higher values to facilitates a wider
        search of the energy landscape, allowing dual_annealing to escape
        local minima that it is trapped in. Default value is 5230. Range
        is (0.01, 5.e4].

    dT : float, optional
        During the annealing process, temperature is decreasing, when it
        reaches initial_temp * restart_temp_ratio, the reannealing process is
        triggered. Default value of the ratio is 2e-5. Range is (0, 1).

    mxcyle : int, optional
        Soft limit for the number of objective function calls. If the
        algorithm is in the middle of a local search, this number will
        be exceeded, the algorithm will stop just after the local search
        is done. Default value is 1e7.

    no_local_search : bool, optional
        If no_local_search is set to True, a traditional Generalized Simulated
        Annealing will be performed with no local search strategy applied.

    local_search_options : dict, optional
        Extra keyword arguments to be passed to the local minimizer (minimize).
        Some important options could be: method for the minimizer method to
        use and args for objective function additional arguments.

#* Funcionan los siguientes local_search_options
#   BFGS o L-BFGS-B o CG o TNC o SLSQP (no usa Hess)
#   COBYLA, Nelder-Mead (no usa jac y hess), Powell (no usa jac y hess)
#todos BFGS, L-BFGS-B o SLSQP se usan por defecto
#!  trust_constr/ncg Whenever the gradient is estimated via finite-differences,
#!  we require the Hessian to be estimated using one of the quasi-Newton
#!  strategies.

    visit_regions : float, optional
        Parameter for visiting distribution. Default value is 2.62. Higher
        values give the visiting distribution a heavier tail, this makes
        the algorithm jump to a more distant region. The value range is
        (1, 3].

    accept : float, optional
        Parameter for acceptance distribution. It is used to control the
        probability of acceptance. The lower the acceptance parameter,
        the smaller the probability of acceptance. Default value is -5.0
        with a range (-1e4, -5].

    callback : callable, optional
        A callback function with signature ``callback(x, f, context)``,
        which will be called for all minima found.
        ``x`` and ``f`` are the coordinates and function value of the
        latest minimum found, and ``context`` has value in [0, 1, 2], with the
        following meaning:

            - 0: minimum detected in the annealing process.
            - 1: detection occurred in the local search process.
            - 2: detection done in the dual annealing process.

        If the callback implementation returns True, the algorithm will stop.

    x0 : ndarray, shape(n,), optional
        Coordinates of a single N-D starting point.

    Returns
    -------
    Optimize Result
        The optimization result represented as a OptimizeResult object.
        Important attributes are: x the solution array, fun the value of the
        function at the solution, and message which describes the cause of
        the termination. See OptimizeResult for a description of other
        attributes.

    Exception
    ---------


    Note
    ----
        References [1]: Tsallis C. Possible generalization of Boltzmann-Gibbs
        statistics. Journal of Statistical Physics, 52, 479-487 (1998).

        Check `scipy.optimize.dual_annealing` official documentation
    """

    if not bounds:
        raise Exception("\n\n *** ERROR: No 'bounds' founds\n")

    if not isinstance(func, types.FunctionType):
        raise Exception("\n\n *** ERROR: No 'function' found \n", func)

    if x0 is not None and not len(x0) == len(bounds):
        raise ValueError('Bounds size does not match x0')

    lu = list(zip(*bounds))
    lower = np.array(lu[0])
    upper = np.array(lu[1])

    # Check that restart temperature ratio is correct
    if dT <= 0. or dT >= 1.:
        raise ValueError('Restart temperature ratio has to be in range (0, 1)')
    # Checking bounds are valid
    if (np.any(np.isinf(lower)) or np.any(np.isinf(upper)) or np.any(
            np.isnan(lower)) or np.any(np.isnan(upper))):
        raise ValueError('Some bounds values are inf values or nan values')
    # Checking that bounds are consistent
    if not np.all(lower < upper):
        raise ValueError('Bounds are not consistent min < max')
    # Checking that bounds are the same length
    if not len(lower) == len(upper):
        raise ValueError('Bounds do not have the same dimensions')
    # Wrapper for the objective function
    func_wrapper = ObjectiveFunWrapper(func, mxcycle, *args)
    # Wrapper ascec
    ascec_wrapper = Ascec(ascec_activation)
    # Wrapper fot the minimizer
    minimizer_wrapper = LocalSearchWrapper(
        bounds, func_wrapper, **local_search_options)
    # Initialization of random Generator for reproducible runs if seed provided
    rand_state = check_random_state(seed)
    # Initialization of the energy state
    energy_state = EnergyState(lower, upper, callback)
    # Save attibutes as: energy_state.ebest, energy_state.current_energy, ...
    energy_state.reset(func_wrapper, rand_state, x0)
    # Minimum value of annealing temperature reached to perform
    # re-annealing
    temperature_restart = T0 * dT

    if there_is_molecule == 0:
        # VisitingDistribution instance (Class used to generate new coordinates)
        #! VisitingDistribution is based Cauchy-Lorentz distribution, se podría
        #! usar para cluster de átomo, ya que con moléculas las puede destruir,
        #! porqué no diferencia entre átomos y molécula.
        #! Si usara el cm como sln, tampoco se podría debido a que también
        #! necesitaría los ejes principales, cada vez que se edite el cm
        #? Una posible sln sería usar la clase Molecule
        visit_dist = VisitingDistribution(lower, upper, visit_regions, rand_state)
    else:
        #visit_dist = Cluster(args)
        visit_dist = MoveFragments(system_object, lower, upper, visit_regions, rand_state)

    # Strategy chain instance (Markov Chain, call: VisitingDistribution)
    strategy_chain = StrategyChain(accept, visit_dist, func_wrapper,
                                   minimizer_wrapper, rand_state,
                                   energy_state)
    need_to_stop = False
    iteration = 0
    message = []
    # OptimizeResult object to be returned
    optimize_res = OptimizeResult()
    optimize_res.success = True
    optimize_res.status = 0

    t1 = np.exp((visit_regions - 1) * np.log(2.0)) - 1.0
    # Run the search loop
    while(not need_to_stop):
        for i in range(NT):
            # Compute temperature for this step
            s = float(i) + 2.0
            t2 = np.exp((visit_regions - 1) * np.log(s)) - 1.0
            temperature = T0 * t1 / t2
            if iteration >= NT:
                message.append("Maximum number of iteration reached")
                need_to_stop = True
                break
            # Need a re-annealing process?
            if temperature < temperature_restart:
                energy_state.reset(func_wrapper, rand_state)
                break
            # starting strategy chain
            val = strategy_chain.run(i, temperature) #acá se llama a visit_dist.visiting
            #print("val, Temperature : ", energy_state.ebest,
            #      energy_state.current_energy, temperature)

            # ASCEC
            if ascec_wrapper.ascec is True:
                e_before = energy_state.ebest
                e_now = energy_state.current_energy
                ascec_wrapper.ascec_criterion(e_before, e_now, temperature)

            if val is not None:
                message.append(val)
                need_to_stop = True
                optimize_res.success = False
                break
            # Possible local search at the end of the strategy chain
            if not no_local_search:
                val = strategy_chain.local_search()
                if val is not None:
                    message.append(val)
                    need_to_stop = True
                    optimize_res.success = False
                    break
            iteration += 1

    # Setting the OptimizeResult values
    optimize_res.x = energy_state.xbest
    optimize_res.fun = energy_state.ebest
    optimize_res.nit = iteration
    optimize_res.nfev = func_wrapper.nfev
    optimize_res.njev = func_wrapper.ngev
    optimize_res.nhev = func_wrapper.nhev
    optimize_res.message = message
    print("NT ",NT)
    return optimize_res


#bounds = [(-5.0, 5.0), (-5.0, 5.0)]
# bounds = [(-5.12, 5.12)]
#seed = 333
#ascec_activation = True
#if ascec_activation is True:
#    ascec_wrapper = ascec(ascec_activation)

#search_config = solve_dual_annealing(Himmelblau, bounds, NT=3,
#                                     no_local_search=False, visit_regions=2.9)

#print(search_config)
