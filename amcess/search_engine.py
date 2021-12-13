from copy import deepcopy

import numpy as np
from scipy.optimize import shgo
from scipy.optimize import dual_annealing


from amcess.ascec import Ascec
from amcess.base_molecule import Cluster
from amcess.electronic_energy import ElectronicEnergy
from amcess.gaussian_process import solve_gaussian_processes

METHODS = {
    "ASCEC": Ascec,
    "dual_annealing": dual_annealing,
    "SHGO": shgo,
    "Bayesian": solve_gaussian_processes,
}


def lennard_jones(r, epsilon=1.0, sigma=1.0) -> float:
    """Lennard-Jones potential

    Parameters
    ----------
    r : float
        Distance between two atoms.
    epsilon : float, optional
        Depth of potential well, by default 1.0
    sigma : float, optional
        Width of potential well, by default 1.0

    Returns
    -------
    float
        Lennard-Jones potential
    """

    return 4 * epsilon * ((sigma / r) ** 12 - (sigma / r) ** 6)


COST_FUNCTIONS = {
    "pyscf": ElectronicEnergy.pyscf,
    "Lennard_Jones": lennard_jones,
}


class SearchConfig:
    """
    Interface to articulate the cluster object with type of search
    and the calculation of energy

    Parameters
    ----------
    system_object : object
        Object made with the Cluster class
    search_methodology : int
        Integer associated with type of searching
    basis : string
        Label of basis set
    program_electronic_structure : int
        Integer associated with the program to make the
        electronic structure calculations
    outxyz : string
        Name of the output xyz with coordinates of the
        configurations accepts

    Returns
    -------
    Output xyz with coordinates and electronic structure

    Raises
    ------
    TypeError
        System_object isn't difinite
        AttributeError system_object isn't difinite as
        an object Cluster
    """

    def __init__(
        self,
        system_object: Cluster = None,
        search_methodology: str = "ASCEC",
        methodology: str = "HF",
        basis: str = "sto-3g",
        tolerance_contour_radius: float = 1.0,
        outxyz: str = "configurations.xyz",
        cost_function="pyscf",
    ) -> None:

        # Verfication and assigment of variables (type, value)
        self.system_object = system_object
        self.search_type = search_methodology
        self.methodology = methodology
        self.basis_set = basis
        self.output_name = outxyz
        self.tolerance_contour_radius = tolerance_contour_radius
        self.func_cost = cost_function

        # Check Overlaping
        self._system_object.initialize_cluster()

        # Build bounds, format for scipy functions
        if system_object._sphere_radius is None:
            self.spherical_contour_cluster()
        else:
            # ! arreglar esto, al fin que va a pasar con sphere_center y radius
            # ! y eso articular con bounds
            self._sphere_center = system_object.sphere_center
            self._sphere_radius = system_object.sphere_radius
            new_radius_t = (
                self._sphere_radius
            )  # self._tolerance_contour_radius + new_radius

            bound_translate = [
                (-new_radius_t, new_radius_t),
                (-new_radius_t, new_radius_t),
                (-new_radius_t, new_radius_t),
            ]
            bound_rotate = [(-180, 180), (-180, 180), (-180, 180)]

            bound_translate = bound_translate * (
                self._system_object.total_molecules - 1
            )
            bound_rotate = bound_rotate * (
                self._system_object.total_molecules - 1
            )

            self._bounds = bound_translate + bound_rotate

    # ===============================================================
    # Decorators
    # ===============================================================

    def bounds_sphere_change(function_change_radius):
        def new_bounds(self, new_radius):
            """
            Define the bounds for the optimization algorithm

            Returns
            -------
            bounds : list
                Bounds for the optimization algorithm
            """
            new_radius_t = self._tolerance_contour_radius + new_radius

            bound_translate = [
                (-new_radius_t, new_radius_t),
                (-new_radius_t, new_radius_t),
                (-new_radius_t, new_radius_t),
            ]
            bound_rotate = [(-180, 180), (-180, 180), (-180, 180)]

            bound_translate = bound_translate * (
                self._system_object.total_molecules - 1
            )
            bound_rotate = bound_rotate * (
                self._system_object.total_molecules - 1
            )

            self._bounds = bound_translate + bound_rotate

            return function_change_radius(self, new_radius)

        return new_bounds

    # ===============================================================
    # PROPERTIES
    # ===============================================================
    @property
    def system_object(self):
        return self._system_object

    @system_object.setter
    def system_object(self, new_object):
        if new_object is None:
            raise TypeError("System_object isn't difinite\n" "It's NoneType")
        if not isinstance(new_object, Cluster):
            raise TypeError(
                "System_object isn't difinite as an object Cluster\n"
                f"please, check:\n'{new_object}'"
            )
        self._system_object = new_object

    @property
    def bounds(self):
        return self._bounds

    @bounds.setter
    def bounds(self, new_bounds):
        if len(new_bounds) != len(self._bounds):
            raise ValueError(
                "\n\nArray dimensions insufficient: "
                f"\ndimensions of old bounds: '{len(self._bounds)}'\n"
                f"\ndimensions of new bounds: '{len(new_bounds)}'\n"
            )

        self._bounds = new_bounds

    @property
    def output_name(self):
        return self._output_name

    @output_name.setter
    def output_name(self, new_name_output):
        if not isinstance(new_name_output, str):
            raise TypeError(
                "\n\nThe new name to output is not a string"
                f"\nplease, check: '{type(new_name_output)}'\n"
            )

        self._output_name = new_name_output

    @property
    def search_type(self):
        return self._search_methodology

    @search_type.setter
    def search_type(self, change_search_methodology):
        if not isinstance(change_search_methodology, str):
            raise TypeError(
                "\n\nThe new search methodology is not a string"
                f"\nplease, check: '{type(change_search_methodology)}'\n"
            )
        if change_search_methodology not in METHODS and not callable(
            change_search_methodology
        ):
            available = list(METHODS.keys())
            raise ValueError(f"Invalid value. options are: {available}")

        self._search_methodology = change_search_methodology

    @property
    def methodology(self):
        return self._methodology

    @methodology.setter
    def methodology(self, new_methodology):
        if not isinstance(new_methodology, str):
            raise TypeError(
                "\n\nThe new name to basis set is not a string"
                f"\nplease, check: '{type(new_methodology)}'\n"
            )

        self._methodology = new_methodology

    @property
    def basis_set(self):
        return self._basis_set

    @basis_set.setter
    def basis_set(self, new_basis_set):
        if not isinstance(new_basis_set, str):
            raise TypeError(
                "\n\nThe new name to basis set is not a string"
                f"\nplease, check: '{type(new_basis_set)}'\n"
            )

        self._basis_set = new_basis_set

    @property
    def tolerance_contour_radius(self):
        return self._tolerance_contour_radius

    @tolerance_contour_radius.setter
    def tolerance_contour_radius(self, new_tol_radius: float):
        if not isinstance(new_tol_radius, float):
            raise TypeError(
                "\n\nThe new tolerance radius is not a float"
                f"\nplease, check: '{type(new_tol_radius)}'\n"
            )
        self._tolerance_contour_radius = new_tol_radius

    @property
    def sphere_center(self) -> tuple:
        return self._sphere_center

    @sphere_center.setter
    def sphere_center(self, new_center: tuple) -> None:
        if not isinstance(new_center, tuple):
            raise TypeError(
                "\n\nThe Sphere center must be a tuple with three elements: "
                "(float, float, float)"
                f"\nplease, check: '{type(new_center)}'\n"
            )
        if len(new_center) != 3:
            raise ValueError(
                "\n\nThe Sphere center must be a tuple with three elements: "
                "(float, float, float)"
                f"\nplease, check: '{new_center}'\n"
            )

        self._sphere_center = new_center

    @property
    def sphere_radius(self) -> float:
        return self._sphere_radius

    @sphere_radius.setter
    @bounds_sphere_change
    def sphere_radius(self, new_radius: float) -> None:
        if not isinstance(new_radius, float):
            raise TypeError(
                "\n\nThe Sphere  Radius must be a float"
                f"\nplease, check: '{type(new_radius)}'\n"
            )
        self._sphere_radius = new_radius + self._tolerance_contour_radius
        if self._sphere_radius <= self._tolerance_contour_radius:
            raise ValueError(
                "\n\nThe Sphere Radius more tolerance must be larger than 1 A"
                f"\nplease, check: '{new_radius}'\n"
            )

    @property
    def func_cost(self):
        return self._func_cost

    @func_cost.setter
    def func_cost(self, new_func_cost):
        self._func_cost = new_func_cost

    # ===============================================================
    # Methods
    # ===============================================================

    def spherical_contour_cluster(self, new_tol: float = None):
        """
        Define a spherical outline that contains our cluster

        Parameters
        ----------
        tolerance : float
            Tolerance with the radius between the mass center to the
            furthest atom

        Returns
        -------
        sphere_center : tuple
            Mass center of the biggest molecule
        sphere_radius : float
            Radius between the sphere center to the furthest atom

        """
        if new_tol is not None:
            self._tolerance_contour_radius = new_tol

        max_distance_cm = 0.0
        molecule = 0
        max_atoms = 0

        # The biggest molecule
        for i in range(self._system_object.total_molecules):
            if self._system_object.get_molecule(i).total_atoms > max_atoms:
                max_atoms = self._system_object.get_molecule(i).total_atoms
                molecule = i

        self._sphere_center = self._system_object.get_molecule(
            molecule
        ).center_of_mass

        # Move the biggest molecule to initio in the cluster object,
        # if is necessary
        if molecule != 0:
            new_geom = dict()
            for i in range(self._system_object.total_molecules):
                if i == 0:
                    new_geom[i] = self._system_object.get_molecule(molecule)
                elif i == molecule:
                    new_geom[i] = self._system_object.get_molecule(0)
                else:
                    new_geom[i] = self._system_object.get_molecule(i)

            self.system_object = Cluster(
                *new_geom.values(), sphere_center=self._sphere_center
            )

        # Radius between the sphere center to the furthest atom
        for xyz in self._system_object.coordinates:
            temp_r = np.linalg.norm(
                np.asarray(self._sphere_center) - np.asarray(xyz)
            )
            if temp_r > max_distance_cm:
                max_distance_cm = temp_r

        self.sphere_radius = max_distance_cm

    def run(self, **kwargs):
        """
        Alternative to execute the searching methodologies in METHODS

        Parameters
        ----------
        kwargs : dict
            Dictionary with the parameters to be used in the
            search methodologies
        """
        func = (
            self._search_methodology
            if callable(self._search_methodology)
            else METHODS[self._search_methodology]
        )

        if self._search_methodology == "ASCEC":
            print("*** Minimization: ASCEC ***")
            self._search = func(
                object_system=self._system_object,
                search_type=self._search_methodology,
                sphere_center=self._sphere_center,
                sphere_radius=self._sphere_radius,
                methodology=self._methodology,
                basis_set=self._basis_set,
                call_function=1,
                bounds=self._bounds,
                **kwargs,
            )
            self._search.ascec_run()
            self._search.write_to_file(self.output_name)
        else:
            if self._search_methodology == "dual_annealing":
                print("*** Minimization: Dual Annealing ***")
            if self._search_methodology == "SHGO":
                print("*** Minimization: SHGO from Scipy ***")
            if self._search_methodology == "Bayesian":
                print("*** Minimization: Bayesian ***")

            if self._search_methodology != "ASCEC":
                obj_ee = ElectronicEnergy(
                    self._system_object,
                    self._search_methodology,
                    self._sphere_center,
                    self._sphere_radius,
                    self._methodology,
                    self._basis_set,
                )

                if self._func_cost == "pyscf":
                    cost_func = obj_ee.pyscf
                else:
                    cost_func = COST_FUNCTIONS[self._func_cost]

            self._search = func(
                cost_func,
                bounds=self._bounds,
                **kwargs,
            )
            obj_ee.write_to_file(self.output_name)
