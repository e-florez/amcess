from scipy.optimize import dual_annealing, shgo  # type: ignore

from amcess.ascec import Ascec
from amcess.cluster import Cluster
from amcess.electronic_energy import ElectronicEnergy
# from amcess.gaussian_process import solve_gaussian_processes

METHODS = {
    "ASCEC": Ascec,
    "dual_annealing": dual_annealing,
    # "SHGO": shgo,
    # "Bayesian": solve_gaussian_processes,
}


class SearchConfig:
    """
    Interface to articulate the cluster object with type of search
    and the calculation of energy

    .. rubric:: Parameters

    system_object : object
        Object made with the Cluster class
    search_type : int
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
        search_type: str = "ASCEC",
        methodology: str = "HF",
        basis: str = "sto-3g",
        outxyz: str = "configurations.xyz",
        cost_fun="pyscf",
        bounds=None,
    ) -> None:
        # ---------------------------------------------------------------
        # Verfication and instantiation (type, value)
        # -- Cluster Object
        #    Calculate center and radius sphere when are null
        self._system_object = system_object
        if self._system_object.GetSphereR() is None and bounds is None:
            self._system_object = system_object.CalCentRSphere()
        # -- Search Methodology: ASCEC, SHGO, dual_annealing, Bayesian
        self._search_type = search_type
        # -- Methodology: HF, DFT, MP2, etc.
        self._methodology = methodology
        # -- Basis Set: sto-3g, 6-31g, 6-31g**, etc.
        self._basis_set = basis
        # -- Output name: xyz
        self._output_name = outxyz
        # -- Cost function: pyscf, Lennard_Jones
        self._cost_fun = cost_fun
        # ---------------------------------------------------------------
        if bounds is None:
            # Build bounds, format for scipy functions
            sphere_radius = self._system_object.GetSphereR()
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
                self._system_object.GetNumMols() - 1
            ) * bound_translate

            bound_rotate = bound_rotate * (self._system_object.GetNumMols() - 1)
            # -- concatenate bounds
            self._bounds = bound_translate + bound_rotate
        else:
            self._bounds = bounds

    # ===============================================================
    # PROPERTIES
    # ===============================================================
    # ! Getter
    def GetBasisSet(self):
        """Basis set"""
        return self._basis_set

    def GetBounds(self):
        """System limit"""
        return self._bounds

    def GetCostFun(self):
        """Cost function"""
        return self._cost_fun

    def GetMethodology(self):
        """Hamiltonian or energy methodology"""
        return self._methodology

    def GetOutputName(self):
        """Output file name"""
        return self._output_name

    def GetSearchType(self):
        """Search/optimization type"""
        return self._search_type

    def GetSystemObject(self):
        """System"""
        return self._system_object

    # ! Setter
    def SetBasisSet(self, new_basis_set):
        """Basis set"""
        if not isinstance(new_basis_set, str):
            raise TypeError(
                "\n\nThe new name to basis set is not a string"
                f"\nplease, check: '{type(new_basis_set)}'\n"
            )

        self._basis_set = new_basis_set

    def SetCostFun(self, new_cost_fun):
        """Cost function"""
        if not isinstance(new_cost_fun, str):
            raise TypeError(
                "\n\nThe new cost function is not a string"
                f"\nplease, check: '{type(new_cost_fun)}'\n"
            )

        self._cost_func = new_cost_fun

    def SetBounds(self, new_bounds):
        """System limit"""
        if len(new_bounds) != len(self._bounds):
            raise ValueError(
                "\n\nArray dimensions insufficient: "
                f"\ndimensions of old bounds: '{len(self._bounds)}'\n"
                f"\ndimensions of new bounds: '{len(new_bounds)}'\n"
            )

        self._bounds = new_bounds

    def SetMethodology(self, new_methodology):
        """Hamiltonian or energy methodology"""
        if not isinstance(new_methodology, str):
            raise TypeError(
                "\n\nThe new name to methodology is not a string"
                f"\nplease, check: '{type(new_methodology)}'\n"
            )

        self._methodology = new_methodology

    def SetOutputName(self, new_name_output):
        """Output file name"""
        if not isinstance(new_name_output, str):
            raise TypeError(
                "\n\nThe new name to output is not a string"
                f"\nplease, check: '{type(new_name_output)}'\n"
            )

        self._output_name = new_name_output

    def SetSearchType(self, change_search_type):
        """Search/optimization type"""
        if not isinstance(change_search_type, str):
            raise TypeError(
                "\n\nThe new search methodology is not a string"
                f"\nplease, check: '{type(change_search_type)}'\n"
            )
        if change_search_type not in METHODS and not callable(change_search_type):
            available = list(METHODS.keys())
            raise ValueError(f"Invalid value. options are: {available}")

        self._search_type = change_search_type

    def SetSystemObject(self, new_object):
        """System"""
        if new_object is None:
            raise TypeError("System_object isn't difinite\n" "It's NoneType")
        if not isinstance(new_object, Cluster):
            raise TypeError(
                "System_object isn't difinite as an object Cluster\n"
                f"please, check:\n'{new_object}'"
            )
        self._system_object = new_object

    # ===============================================================
    # Methods
    # ===============================================================

    def RunSearch(self, **kwargs):
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
            self.GetSearchType()
            if callable(self.GetSearchType())
            else METHODS[self.GetSearchType()]
        )
        # ---------------------------------------------------------------
        # Execute the search methodologies
        if self.GetSearchType() == "ASCEC":
            print("*** Minimization: ASCEC ***")
            self._search = func(
                object_system=self._system_object,
                search_type=self.GetSearchType(),
                methodology=self.GetMethodology(),
                basis_set=self.GetBasisSet(),
                program=self.GetCostFun(),
                bounds=self.GetBounds(),
                **kwargs,
            )
            self._search.RunASCEC()
            self._search.WriteOutPut(self.GetOutputName())
        else:
            if self.GetSearchType() == "dual_annealing":
                print("*** Minimization: Dual Annealing ***")
            # if self.GetSearchType() == "SHGO":
            #    print("*** Minimization: SHGO from Scipy ***")
            # if self.GetSearchType() == "Bayesian":
            #    print("*** Minimization: Bayesian ***")

            if self.GetSearchType() != "ASCEC":
                obj_ee = ElectronicEnergy(
                    self._system_object,
                    self.GetSearchType(),
                    self._methodology,
                    self._basis_set,
                )

            cost_func = obj_ee.Pyscf # ! with parenthesis return error

            self._search = func(
                cost_func,
                bounds=self.GetBounds(),
                **kwargs,
            )
            obj_ee.WriteOutPut(self.GetOutputName())
