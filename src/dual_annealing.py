import types
from scipy.optimize import dual_annealing
from scipy import constants
import numpy as np # 1.19.1

def store_T(T):
    """
    Store annealing temperature because when is used the local_search is
    not return temperature by dual_annealing

    Args:
        T : float
            Annealing temperature

    Returns:
        T : float
            Store temperature
    """
    store_T.T = T or store_T.T
    return store_T.T
class ascec(object):
    """
    stores or starts varaibles about ascec
    """
    def __init__(self,ascec,*args):
        """
        Store the boolean B. When B is True is used the ASCEC criterio
        in the acceptance

        Parameters
        ----------
            B : Boolean
                True activate ASCEC criterio

        Returns
        -------
            B : Boolean
        """
        self.ascec = ascec

    def ascec_energy(self,E=None):
        """
        Store the Energy of structure accept

        Parameters
        ----------
            E : Float
                Energy of structure accept

        Returns
        -------
            E : Float
                Energy of structure accept
        """
        if E is not None:
            self.E = E or self.E
            return self.E

    def ascec_criterion(self,x,e,T,it):
        """[summary]
        Evaluate ASCEC criterion for acceptance

        Parameters
        ----------
        x : array, float
            Value of each coordinate in the 1D array
        e : float
            Value of the cost function
        T : float
            Annealing temperature
        args :
            Any additional fixed parameters needed to completely specify
            the objective function.

        """

        if T is not None:
            store_T(T)
        elif T is None:
            T = store_T.T

        if it > 1:
            DE = self.E - e
            if DE < 0.0000000001:
                print("DE < 0",DE)
                ascec_wrapper.ascec_energy(e)
            else:
                TKb = T*constants.k #Boltzmann constant [J/K]
                exp = np.exp(-DE*constants.h/TKb)
                if DE < exp:
                    print("DE < Boltzmann Poblation ",DE)
                    ascec_wrapper.ascec_energy(e)
        elif it == 1:
            ascec_wrapper.ascec_energy(e)

def Rastrigin(x,T=None,*args):
    """
    Rastrigin function

    f(x) = An + Sum_i^n [x_i*x_i - Acos(2pix_i)]
    A = 10 y x e [-5.12,5.12]

    bounds
    -5.12,5.12

    Parameters
    ----------
        x : array, float
            Value of each coordinate in the 1D array
        T : float
            Annealing temperature
        args :
            Any additional fixed parameters needed to completely specify
            the objective function.

    Returns
    -------
        Value of cost function in x

    Note
    ----
    Reference:
        [1] Rastrigin, L. A. "Systems of extremal control." Mir, Moscow (1974).
    """
    if ascec_wrapper.ascec is True:
        e = np.sum(x*x - 10*np.cos(2*np.pi*x)) + 10*np.size(x)
        ascec_wrapper.ascec_criterion(x,e,T,args[0])
        return e
    return np.sum(x*x - 10*np.cos(2*np.pi*x)) + 10*np.size(x)

def Himmelblau(x, T=None, *args):
    """
    Himmelblau function

    bounds
    -5.0, 5.0

    Roots
    Himmelbau(3.0,2.0)=0.0,\quad
    Himmelbau( − 2.805118 , 3.131312 )   = 0.0
    Himmelbau( − 3.779310 , − 3.283186 ) = 0.0
    Himmelbau( 3.584428 , − 1.848126 )   = 0.0

    Parameters
    ----------
        x : array, float
            Value of each coordinate in the 1D array
        T : float
            Annealing temperature
        args :
            Any additional fixed parameters needed to completely specify
            the objective function.

    Returns
    -------
        Value of cost function in x

    Note
    ----
    Reference:
        [1] Himmelblau, D. (1972). Applied Nonlinear Programming. McGraw-Hill.
        ISBN 0-07-028921-2.
    """
    if T is not None:
        store_T(T)
    elif T is None:
        T = store_T.T

    if ascec_wrapper.ascec is True:
        e=(x[0]**2 + x[1] - 11)**2 + (x[0] + x[1]**2 -7)**2
        if args[0] > 1:
            DE = ascec_wrapper.E - e
            if DE < 0.2 and DE > -0.0000000001:
                print("DE < 0",DE,x[:])
                ascec_wrapper.ascec_energy(e)
        elif args[0] == 1:
            ascec_wrapper.ascec_energy(e)
        return e
#!function is not good to analyze ascec criterio, because it's very uniform
    return (x[0]**2 + x[1] - 11)**2 + (x[0] + x[1]**2 -7)**2

def solve_dual_annealing(function, bounds, seed=None, NT=1000, T0=5230.0,
                        dT=2e-5, mxcycle=10000000.0, local_search_options={},
                        no_local_search=False, visit_regions=2.62, accept=-5.0):
    """Find the global minimum of a function using Dual Annealing [1].

        Example:
            lambda x: x ** 2
            lambda x: (x[0] ** 2 + x[1] - 11) ** 2

    Parameters
    ----------
    function : callable
        The objective function to be minimized. Must be in the form f(x, *args),
        where x is the argument in the form of a 1-D array and args is a tuple
        of any additional fixed parameters needed to completely specify the function

    bounds : sequence, shape (n, 2)
        Bounds for variables. (min, max) pairs for each element in x, defining bounds
        for the objective function parameter.

    seed : {None, int, numpy.random.Generator,

        numpy.random.RandomState}, optional
    If seed is None (or np.random), the numpy.random.RandomState singleton is used.
    If seed is an int, a new RandomState instance is used, seeded with seed. If seed
    is already a Generator or RandomState instance then that instance is used.
    Specify seed for repeatable minimizations. The random numbers generated with this
    seed only affect the visiting distribution function and new coordinates generation.

    NT : int, optional
    The maximum number of global search iterations. Default value is 1000.

    T0 : float, optional
    The initial temperature, use higher values to facilitates a wider search of the
    energy landscape, allowing dual_annealing to escape local minima that it is
    trapped in. Default value is 5230. Range is (0.01, 5.e4].

    dT : float, optional
    During the annealing process, temperature is decreasing, when it reaches
    initial_temp * restart_temp_ratio, the reannealing process is triggered. Default
    value of the ratio is 2e-5. Range is (0, 1).

    mxcyle : int, optional
    Soft limit for the number of objective function calls. If the algorithm is in the
    middle of a local search, this number will be exceeded, the algorithm will stop just
    after the local search is done. Default value is 1e7.

    no_local_search : bool, optional
    If no_local_search is set to True, a traditional Generalized Simulated Annealing
    will be performed with no local search strategy applied.

    local_search_options : dict, optional
    Extra keyword arguments to be passed to the local minimizer (minimize). Some
    important options could be: method for the minimizer method to use and args for
    objective function additional arguments.

#* Funcionan los siguientes local_search_options
#   BFGS o L-BFGS-B o CG o TNC o SLSQP (no usa Hess)
#   COBYLA, Nelder-Mead (no usa jac y hess), Powell (no usa jac y hess)
#todos BFGS, L-BFGS-B o SLSQP se usan por defecto
#!  trust_constr/ncg Whenever the gradient is estimated via finite-differences,
#!  we require the Hessian to be estimated using one of the quasi-Newton strategies.

    visit_regions : float, optional
    Parameter for visiting distribution. Default value is 2.62. Higher values give
    the visiting distribution a heavier tail, this makes the algorithm jump to a more
    distant region. The value range is (1, 3].

    accept : float, optional
    Parameter for acceptance distribution. It is used to control the probability of
    acceptance. The lower the acceptance parameter, the smaller the probability of
    acceptance. Default value is -5.0 with a range (-1e4, -5].

    Returns
    -------
    Optimize Result
        The optimization result represented as a OptimizeResult object. Important
        attributes are: x the solution array, fun the value of the function at the
        solution, and message which describes the cause of the termination.
        See OptimizeResult for a description of other attributes.

    Exception
    ---------


    Note
    ----
        References [1]: Tsallis C. Possible generalization of Boltzmann-Gibbs
        statistics. Journal of Statistical Physics, 52, 479-487 (1998).

        Check `scipy.optimize.dual_annealing` official documentation
    """

    if not bounds:
        raise Exception(f"\n\n *** ERROR: No 'bounds' founds\n")

    if not isinstance(function, types.FunctionType):
        raise Exception(f"\n\n *** ERROR: No 'function' found\n")

    #* Store initial temperature
    store_T(T0)

    result = dual_annealing(
        function,
        bounds=bounds,
        maxiter=NT,
        initial_temp=T0,
        restart_temp_ratio=dT,
        seed=seed,
        maxfun=mxcycle,
        local_search_options=local_search_options,
        no_local_search=no_local_search,
        visit=visit_regions,
        accept=accept
    )

    return result
bounds = [(-5.0,5.0),(-5.0,5.0)]
#bounds = [(-5.12,5.12)]
seed = 333
ascec_activation = True
if ascec_activation is True:
    ascec_wrapper = ascec(ascec_activation)

search_config = solve_dual_annealing(Himmelblau, bounds,
                no_local_search=True,visit_regions=2.9)

print(search_config)