from scipy.optimize import dual_annealing, shgo

from amcess2024.ascec import Ascec
from amcess2024.base_molecule import Cluster
from amcess2024.electronic_energy import ElectronicEnergy
from amcess2024.gaussian_process import solve_gaussian_processes

METHODS = {
    "ASCEC": Ascec,
    "dual_annealing": dual_annealing,
    "SHGO": shgo,
    "Bayesian": solve_gaussian_processes,
}


class SearchConfig:
    """
    Interface to articulate the cluster object with type of search
    and the calculation of energy

    .. rubric:: Parameters

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

    .. rubric:: Returns

    Output xyz with coordinates and electronic structure

    .. rubric:: Raises

    TypeError
        System_object isn't define. AttributeError system_object isn't
        define as an object Cluster
    """

    def __init__(
        self,
        system_object: Cluster = None,
        search_methodology: str = "ASCEC",
        methodology: str = "HF",
        basis: str = "sto-3g",
        outxyz: str = "configurations.xyz",
        cost_function="pyscf",
    ) -> None:
        # ---------------------------------------------------------------
        # Verfication and instantiation (type, value)
        # -- Cluster Object
        #    Calculate center and radius sphere when are null
        if system_object._sphere_radius is None:
            self.system_object = system_object.center_radius_sphere()
        else:
            self.system_object = system_object
        # -- Search Methodology: ASCEC, SHGO, dual_annealing, Bayesian
        self.search_type = search_methodology
        # -- Methodology: HF, DFT, MP2, etc.
        self.methodology = methodology
        # -- Basis Set: sto-3g, 6-31g, 6-31g**, etc.
        self.basis_set = basis
        # -- Output name: xyz
        self.output_name = outxyz
        # -- Cost function: pyscf, Lennard_Jones
        self.func_cost = cost_function
        # ---------------------------------------------------------------
        # Build bounds, format for scipy functions
        sphere_radius = self._system_object._sphere_radius
        # -- translate bounds
        bound_translate = [
            (-sphere_radius, sphere_radius),
            (-sphere_radius, sphere_radius),
            (-sphere_radius, sphere_radius),
        ]
        # -- rotate bounds
        bound_rotate = [(-180, 180), (-180, 180), (-180, 180)]
        # -- Multiply bounds by the amount of molecules
        bound_translate = (  # noqa
            self._system_object.total_molecules - 1
        ) * bound_translate

        bound_rotate = bound_rotate * (self._system_object.total_molecules - 1)
        # -- concatenate bounds
        self._bounds = bound_translate + bound_rotate

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
                "\n\nThe new name to methodology is not a string"
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
    def func_cost(self):
        return self._func_cost

    @func_cost.setter
    def func_cost(self, new_func_cost):
        if not isinstance(new_func_cost, str):
            raise TypeError(
                "\n\nThe new cost function is not a string"
                f"\nplease, check: '{type(new_func_cost)}'\n"
            )

        self._func_cost = new_func_cost

    # ===============================================================
    # Methods
    # ===============================================================

    def run(self, **kwargs):
        """
        Alternative to execute the searching methodologies in METHODS

        .. rubric:: Parameters

        kwargs : dict
            Dictionary with the parameters to be used in the search
            methodologies
        """
        # ---------------------------------------------------------------
        # Choose the search methodologies
        func = (
            self._search_methodology
            if callable(self._search_methodology)
            else METHODS[self._search_methodology]
        )
        # ---------------------------------------------------------------
        # Execute the search methodologies
        if self._search_methodology == "ASCEC":
            print("*** Minimization: ASCEC ***")
            self._search = func(
                object_system=self._system_object,
                search_type=self._search_methodology,
                methodology=self._methodology,
                basis_set=self._basis_set,
                program=self._func_cost,
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
                    self._methodology,
                    self._basis_set,
                )

            cost_func = obj_ee.pyscf

            self._search = func(
                cost_func,
                bounds=self._bounds,
                **kwargs,
            )
            obj_ee.write_to_file(self.output_name)
