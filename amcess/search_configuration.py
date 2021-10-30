import random
import sys

import scipy
from scipy.optimize import shgo

from amcess.electronic_energy import ElectronicEnergy, hf_pyscf
from amcess.base_molecule import Cluster
from amcess.m_dual_annealing import solve_dual_annealing


def overlaping(_system_object):
    """
    Confirm that the molecules aren't overlaping

    Parameters
    ----------
        _system_object : object
            Cluster object

    Returns
    -------
        fragments : dictionary
            Input for Cluster class

    Print warning when there is overlaping
    """
    import warnings  # !Falta definir el tipo de warning

    fragments = dict()
    for i in range(_system_object.total_molecules - 1):
        for j in range(i + 1, _system_object.total_molecules):
            fragments[i] = _system_object.get_molecule(i)
            equal_arrays = (
                _system_object.get_molecule(i).center_of_mass
                == _system_object.get_molecule(j).center_of_mass
            )
            if equal_arrays:
                r = random.random()
                fragments[j] = (
                    Cluster(_system_object.get_molecule(j))
                    .translate(0, r, r, r)
                    .rotate(0, r, r, r)
                ).get_molecule(0)
                message = (
                    "Center of Mass of the fragments "
                    + str(i)
                    + " and "
                    + str(j)
                    + " are overlapping. Then "
                    + str(j)
                    + " is move "
                    + str(r)
                )
                warnings.warn(message)
            else:
                fragments[j] = _system_object.get_molecule(j)
    return fragments


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
        bases : string
            Label of bases set
        program_electronic_structure : int
            Integer associated with the program to make the
            electronic structure calculations
        outxyz : string
            Name of the output xyz with coordinates of the
            configurations accepts

    Returns
    -------
        Output xyz with coordinates and electronic structure
    """

    def __init__(
        self,
        system_object=None,
        search_methodology=1,
        bases="sto-3g",
        program_electronic_structure=1,
        tolerance_contour_radius=1,
        outxyz="configurations.xyz",
    ) -> None:

        if system_object is None:
            raise ValueError(
                "AttributeError system_object isn't difinite\n" "It's None"
            )
        elif not isinstance(system_object, object):
            raise ValueError(
                "AttributeError system_object isn't difinite\n" "It's None"
            )
        elif system_object.total_molecules == 1:
            raise ValueError(
                "System of study most have AT LEAST TWO FRAGMENTS"
            )

        self._system_object = system_object

        self._search_methodology = search_methodology

        self._bases = bases
        self._program_calculate_cost_function = program_electronic_structure

        self._tolerance_contour_radius = tolerance_contour_radius

        self._outxyz = outxyz

        # Check Overlaping
        self._system_object.initialize_cluster()

        # Build bounds, format for scipy functions
        if system_object._sphere_radius is None:
            system_object.spherical_contour_cluster(tolerance_contour_radius)

        bound_translate = [
            (-system_object.sphere_radius, system_object.sphere_radius),
            (-system_object.sphere_radius, system_object.sphere_radius),
            (-system_object.sphere_radius, system_object.sphere_radius),
        ]
        bound_rotate = [(0, scipy.pi), (0, scipy.pi), (0, scipy.pi)]
        bound_translate = bound_translate * self._system_object.total_molecules
        bound_rotate = bound_rotate * self._system_object.total_molecules
        self._bounds = bound_translate + bound_rotate

        self._func = self.program_cost_function(
            self._program_calculate_cost_function
        )

        self._obj_ee = ElectronicEnergy(self._system_object)

    # ===============================================================
    # PROPERTIES
    # ===============================================================
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

        self._outxyz = new_name_output

    @property
    def search_name(self):
        return self._search_methodology

    @search_name.setter
    def search_name(self, change_search_methodology):
        if not isinstance(change_search_methodology, int):
            raise TypeError(
                "\n\nThe search methodology is associated with a integer \n"
                "1 -> Dual Annealing \n"
                "2 -> SHGO \n"
                "3 -> Bayessiana \n"
                f"\nplease, check: '{type(change_search_methodology)}'\n"
            )

        self._search_methodology = change_search_methodology

    @property
    def bases_set(self):
        return self._bases_set

    @bases_set.setter
    def bases_set(self, new_bases_set):
        if not isinstance(new_bases_set, str):
            raise TypeError(
                "\n\nThe new name to output is not a string"
                f"\nplease, check: '{type(new_bases_set)}'\n"
            )

        self._bases_set = new_bases_set

    @property
    def radius_contour(self):
        return self._tolerance_contour_radius

    @radius_contour.setter
    def radius_contour(self, new_tol_radius):
        if not isinstance(new_tol_radius, float):
            raise TypeError(
                "\n\nThe new tolerance for contour radius is not a float"
                f"\nplease, check: '{type(new_tol_radius)}'\n"
            )

        self._tolerance_contour_radius = new_tol_radius

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
                f"\nplease, check: '{type(new_func)}'\n"
            )

        self._func = self.program_cost_function(new_func)

    # ===============================================================
    # Methods
    # ===============================================================

    def program_cost_function(self, _program_calculate_cost_function):
        """
        Assign the name of the cost function, which is associated with
        determined program to calculate the energy of the electronic
        structure

        Parameters
        ----------
            _program_calculate_cost_function ([type]): [description]

        Returns
        -------
            called
            name of the function cost which associated with a specify
            program

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
        with open(self._outxyz, "w") as outxyz:
            self._search = solve_dual_annealing(
                self._func,
                self._bounds,
                self._system_object,
                args=(
                    self._bases,
                    self._obj_ee,
                    outxyz,
                    self._search_methodology,
                ),
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
        with open(self._outxyz, "w") as outxyz:
            self._search_methodology = 2

            self._search = shgo(
                self._func,
                bounds=self._bounds,
                sampling_method="sobol",
                args=(
                    self._bases,
                    self._obj_ee,
                    outxyz,
                    self._search_methodology,
                ),
                **kwargs,
            )
