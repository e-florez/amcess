from scipy.optimize import shgo
from scipy.optimize import dual_annealing
import attr
import numpy as np
from amcess.base_molecule import Cluster
from amcess.electronic_energy_2 import Electronic_energy
from amcess.gaussian_process import solve_gaussian_processes


@attr.s
class extra_functions(Electronic_energy):
    """Extra functions used for ordering the cluster introduced by the user"""

    initial_cluster = attr.ib(
        default=None,
        type=Cluster,
        validator=attr.validators.instance_of(Cluster),
    )
    minima_energy = []
    minima_energy_x = []

    @property
    def center_sphere_mass(self):
        """Put the most massive molecule as the first molecule, in order to
        move the others molecules in the cluster around there.
        """

        max_mass = 0

        for i in range(self.initial_cluster.total_molecules):
            if self.initial_cluster.get_molecule(i).total_mass > max_mass:
                max_mass_order = i

        new_cluster_center = self.initial_cluster.get_molecule(max_mass_order)
        self.initial_cluster = self.initial_cluster.remove_molecule(
            max_mass_order
        )
        self.initial_cluster = new_cluster_center + self.initial_cluster
        return self.initial_cluster

    @property
    def center_sphere_atoms(self):
        """Put the molecule with most atoms as the first molecule, in order to
        move the others molecules in the cluster around there.
        """

        max_atom = 0

        for i in range(self.initial_cluster.total_molecules):
            if len(self.initial_cluster.get_molecule(i).atoms) > max_atom:
                max_atom_order = i

        new_cluster_center = self.initial_cluster.get_molecule(max_atom_order)
        self.initial_cluster = self.initial_cluster.remove_molecule(
            max_atom_order
        )
        self.initial_cluster = new_cluster_center + self.initial_cluster
        return self.initial_cluster

    def callback_dual_annealing(self, x, func, context):
        """This function is created in order to get local minima of the function
        that is being optimized.


        Parameters
        ----------
        x : [type]
            [values that are being optimized]
        func : [type]
            [values of the function that is being optimized]
        context : [type]
            [0 or 1, 1 is that the x is a local minimum of the function]

        Raises
        ------
        ValueError
            [description]
        """
        # print(context)
        if context == 1:
            func = np.around(func, decimals=8)
            if self.minima_energy.count(func) == 0:
                self.minima_energy.append(func)
                self.minima_energy_x.append(x)


@attr.s
class LocalMinima(extra_functions):
    """Class that contains methods to find the local minima of a molecular
    system
    """

    initial_cluster = attr.ib(
        default=None,
        type=Cluster,
        validator=attr.validators.instance_of(Cluster),
    )
    tolerance = attr.ib(default=2, type=float)
    basis = attr.ib(
        default="sto-3g", type=str, validator=attr.validators.instance_of(str)
    )
    method = attr.ib(default="hf", type=str)
    ordening_cluster = attr.ib(
        default="mass",
        type=str,
        validator=attr.validators.in_(["mass", "atoms"]),
    )

    def __attrs_post_init__(self):
        """Function to initialize the class"""
        if self.ordening_cluster == "mass":
            self.change_center_sphere = self.center_sphere_mass
        elif self.ordening_cluster == "atoms":
            self.change_center_sphere = self.center_sphere_atoms

    @property
    def bounds(self):
        """Define the bonds of wich the optimization method will search the
        minumus.
        The bound_translate, is defined as the distance between the center of
        mass of the target molecule and each molecule.
        """

        cluster = self.change_center_sphere
        radius_center = cluster.get_molecule(0).center_of_mass
        bound_translate = []
        for i in range(1, cluster.total_molecules, 1):
            distance_between_molecules = np.linalg.norm(
                np.asarray(radius_center)
                - np.asarray(cluster.get_molecule(i).center_of_mass)
            )
            self.sphere_radius = distance_between_molecules + self.tolerance
            mol_translate = [
                (-self.sphere_radius, self.sphere_radius),
                (-self.sphere_radius, self.sphere_radius),
                (-self.sphere_radius, self.sphere_radius),
            ]
            bound_translate = bound_translate + mol_translate

        bound_rotate = [
            (-180, 180),
            (-180, 180),
            (-180, 180),
        ] * (cluster.total_molecules - 1)

        self._bounds = bound_translate + bound_rotate
        # print("los bounds de los métodos son", self._bounds)
        return self._bounds

    def shgo_optimize(self, n=None, iters=None, sampling_method=None):

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
        if n is None:
            n = 16
        if iters is None:
            iters = 1
        if sampling_method is None:
            sampling_method = "sobol"

        energie = self.energy

        minimos = shgo(
            energie,
            bounds=self.bounds,
            iters=iters,
            n=n,
            sampling_method=sampling_method,
        )
        with open("configurations.xyz", "a") as f:
            for i in range(len(minimos.funl)):
                cluster_in_minima = self.controled_move(minimos.xl[i])
                f.write(str(len(self.initial_cluster.atoms)) + "\n")
                f.write("Energy: " + str(minimos.funl[i]) + "\n")
                for terms in cluster_in_minima.atoms:
                    f.write(" ".join([str(x) for x in terms]) + "\n")

        return minimos

    def dual_annealing_optimize(
        self, num_runs: int = None, maxiter: int = None
    ):
        """Function to find the mimina of a molecular system using the the
        dual annealing optimization procedure

        Parameters
        ----------
        num_runs : [int] number of times that the dual_annealing procces will
        be repeated. Default is 5

        maxiter : [int] maximum number of iterations that the dual_annealing
        procces will be repeated. Default is 5

        Returns
        -------
        [OptimizeResult object]
        """
        if num_runs is None:
            num_runs = 2
        if maxiter is None:
            maxiter = 100

        if not isinstance(num_runs, (int)):
            raise ValueError("\n\n num_runs must be an integer")

        for i in range(num_runs):
            energie = self.energy
            minimo = dual_annealing(
                energie,
                bounds=self.bounds,
                maxiter=maxiter,
                callback=self.callback_dual_annealing,
            )

            with open("configurations.xyz", "a") as f:
                for i in range(len(self.minima_energy)):
                    f.write(str(len(self.initial_cluster.atoms)) + "\n")
                    f.write("Energy: " + str(self.minima_energy[i]) + "\n")
                    cluster_in_minima = self.controled_move(
                        self.minima_energy_x[i]
                    )
                    for terms in cluster_in_minima.atoms:
                        f.write(" ".join([str(x) for x in terms]) + "\n")

        return minimo

    def gaussian_optimize(self):
        """Function to find the mimina of a molecular system using the the
        gaussian procces optimization procedure
        """
        energie = self.energy
        minimo = solve_gaussian_processes(energie, bounds=self.bounds)
        return minimo


#    def ascec(self):
#        """Function to find the mimina of a molecular system using the the
#        Andy's optimization procedure
#        """
#        bounds = self.bounds

#        minimo = Ascec(
#            object_system=self.initial_cluster,
#            search_type="ASCEC",
#            sphere_center=self.initial_cluster.get_molecule(0).center_of_mass,
#            sphere_radius=self.sphere_radius,
#            basis_set=self.basis,
#            call_function=1,
#            bounds=self.bounds,
#        )
#        return minimo
