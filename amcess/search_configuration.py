import numpy as np
import scipy
from scipy.optimize import shgo

from amcess.base_molecule import Cluster
from amcess.electronic_energy import ElectronicEnergy, hf_pyscf
from amcess.gaussian_process import solve_gaussian_processes
from amcess.m_dual_annealing import solve_dual_annealing


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
        search_methodology: int = 1,
        basis: str = "sto-3g",
        program_electronic_structure: int = 1,
        tolerance_contour_radius: float = 1.0,
        outxyz: str = "configurations.xyz",
    ) -> None:

        # Verfication and assigment of variables (type, value)
        self.system_object = system_object
        self.search_type = search_methodology
        self.basis_set = basis
        self.cost_function_number = program_electronic_structure
        self.output_name = outxyz
        self.tolerance_contour_radius = tolerance_contour_radius
        self.cost_function_number = program_electronic_structure

        # Check Overlaping
        self._system_object.initialize_cluster()

        # Build bounds, format for scipy functions
        if system_object._sphere_radius is None:
            self.spherical_contour_cluster()

        self._func = self.program_cost_function(
            self._program_calculate_cost_function
        )

        self._obj_ee = ElectronicEnergy(
            self._system_object,
            self._sphere_center,
            self._sphere_radius,
            self._basis_set,
        )

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
            bound_rotate = [(0, scipy.pi), (0, scipy.pi), (0, scipy.pi)]

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
        if not isinstance(change_search_methodology, int):
            raise TypeError(
                "\n\nThe new search methodology is not an integer"
                f"\nplease, check: '{type(change_search_methodology)}'\n"
            )
        if change_search_methodology > 3:
            raise ValueError(
                "\n\nThe search methodology is associated with a integer \n"
                "1 -> Dual Annealing \n"
                "2 -> SHGO \n"
                "3 -> Bayessiana \n"
                f"\nplease, check: '{type(change_search_methodology)}'\n"
            )

        self._search_methodology = change_search_methodology

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
        if not isinstance(new_radius, (int, float)):
            raise TypeError(
                "\n\nThe Sphere  Radius must be a float or int"
                f"\nplease, check: '{type(new_radius)}'\n"
            )
        self._sphere_radius = new_radius + self._tolerance_contour_radius
        if self._sphere_radius <= self._tolerance_contour_radius:
            raise ValueError(
                "\n\nThe Sphere Radius more tolerance must be larger than 1 A"
                f"\nplease, check: '{new_radius}'\n"
            )

    @property
    def cost_function_ee(self):
        return self._func

    @property
    def cost_function_number(self):
        return self._program_calculate_cost_function

    @cost_function_number.setter
    def cost_function_number(self, new_func):
        if not isinstance(new_func, int):
            raise TypeError(
                "\n\nThe new cost function is not a integer"
                f"\nplease, check: '{type(new_func)}'\n"
            )
        elif new_func > 1:
            raise ValueError(
                "\n\nThe new cost function is not implemeted "
                "\n 1 -> Hartree Fock into pyscf"
                f"\nplease, check: '{new_func}'\n"
            )

        self._program_calculate_cost_function = new_func
        self._func = self.program_cost_function(new_func)

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

    def program_cost_function(self, _program_calculate_cost_function):
        """
        Assign the name of the cost function, which is associated with
        determined program to calculate the energy of the electronic
        structure

        Parameters
        ----------
            _program_calculate_cost_function : int
                Integer associated with the program to calculate the cost
                and methodology (Hamiltonian, Functional, etc)

        Returns
        -------
            called
            name of the function cost which associated with a specify
            program and methodology (Hamiltonian, Functional, etc)

        """
        if _program_calculate_cost_function == 1:
            print(
                "\n\n"
                "*** Cost function is Hartree--Fock implemented into pyscf ***"
                "\n\n"
            )
            return hf_pyscf

    def run(self, **kwargs):
        """
        Alternative to execute the searching methodologies

        Parameters
        ----------
            **kwargs : dict
                Dictionary with the parameters to be used in the
                search methodologies
        """
        if self._search_methodology == 1:
            self.da(**kwargs)
        if self._search_methodology == 2:
            self.shgo(**kwargs)
        if self._search_methodology == 3:
            self.bayesian(**kwargs)

    def da(self, **kwargs):
        """
        Execute solve dual annealing to search candidate structure
        and open output file

        Parameters
        ----------
            **kwargs : dict
                Dictionary with the parameters to be used in the
                dual annealing methodology
        """
        print("*** Minimization: Dual Annealing ***")
        with open(self._output_name, "w") as outxyz:
            self._search = solve_dual_annealing(
                self._func,
                self._bounds,
                self._system_object,
                args=(self._obj_ee, outxyz),
                **kwargs,
            )

    def shgo(self, **kwargs):
        """
        Execute solve shgo to search candidate structure
        and open output file

        Parameters
        ----------
            **kwargs : dict
                Dictionary with the parameters to be used in the
                shgo methodology
        """
        print("*** Minimization: SHGO from Scipy ***")
        with open(self._output_name, "w") as outxyz:
            self._search = shgo(
                self._func,
                bounds=self._bounds,
                args=(self._obj_ee, outxyz),
                **kwargs,
            )

    def bayesian(self, **kwargs):
        """
        Execute solve Bayesian to search candidate structure
        and open output file

        Parameters
        ----------
            **kwargs : dict
                Dictionary with the parameters to be used in the
                Bayesian methodology
        """
        print("*** Minimization: Bayesian ***")
        with open(self._output_name, "w") as outxyz:
            self._search = solve_gaussian_processes(
                self._func,
                bounds=self._bounds,
                args=(self._obj_ee, outxyz),
                **kwargs,
            )
