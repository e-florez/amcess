from scipy.optimize import shgo
import attr
import numpy as np
from amcess.base_molecule import Cluster
from amcess.electronic_energy_new import Electronic_energy


class LocalMinima(Electronic_energy):
    """[summary]

    Parameters
    ----------
    Electronic_energy : [type]
        [description]

    Returns
    -------
    [type]
        [description]
    """
    initial_cluster = attr.ib(default=None, type=Cluster)
    basis = attr.ib(default='sto-3g')
    method = attr.ib(default='hf')

    def shgo_optimize(self, n_sampling=None, numb_iters=None,
                      sampling_method_=None):
        """Function to optimize the electronic energy of a cluster.

        Parameters
        ----------
        n_sampling : [integer], optional
            [Number of sampling points used in the construction of the
             simplicial complex.], by default 32
        numb_iters : [int], optional
            [Number of iterations used in the construction of the
            simplicial complex], by default 3
        sampling_method_ : [str], optional
            [Current built in sampling method options are halton, sobol and
            simplicial. The default simplicial provides the theoretical
            guarantee of convergence to the global minimum in finite time.
            halton and sobol method are faster in terms of sampling point
            generation at the cost of the loss of guaranteed convergence.],
            by default sobol

        Returns
        -------
        [OptimizeResult object]
            [The optimization result represented as a OptimizeResult object.
            Important attributes are: x the solution array corresponding to the
            global minimum, fun the function output at the global solution,
            xl an ordered list of local minima solutions, funl the function
            output at the corresponding local solutions]
        """
        if n_sampling is None:
            n_sampling = 16
        if numb_iters is None:
            numb_iters = 1
        if sampling_method_ is None:
            sampling_method_ = 'sobol'
        energie = self.energy
        sphere_radius = self.initial_cluster.total_atoms*0.5*0.5
        print("sphere radiuos", sphere_radius)
        discretization = sphere_radius / 1
        bounds = np.array([(-discretization, discretization),
                           (-discretization, discretization),
                           (-discretization, discretization),
                           (0, 360), (0, 360), (0, 360)])
        shgo_tolerance_dict = {"ftol": 1e-6}  # Default is 1e-12
        shgo_options_dict = {"options": shgo_tolerance_dict}
        bounds_full = np.repeat(bounds, self.initial_cluster.total_molecules-1,
                                axis=0)
        print("bounds full ", len(bounds_full),
              self.initial_cluster.total_molecules, bounds_full)
        minimos = shgo(energie, bounds=bounds_full,
                       minimizer_kwargs=shgo_options_dict,
                       sampling_method=sampling_method_, n=n_sampling,
                       iters=numb_iters)

        return minimos
