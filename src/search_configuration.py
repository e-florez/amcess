import sys
from scipy.optimize import shgo
import random

from m_dual_annealing import *


def overlaping(_system_object):
    """[summary]
    Confirm that the molecules aren't overlaping

    Args:
        _system_object ([Object]): Molecule object

    Returns:
        fragments [dictionary]: Input for Molecule Class

    Print warning when there is overlaping
    """
    import warnings  #! Falta definir el tipo de warning

    fragments = dict()
    for i in range(_system_object._total_fragments - 1):
        for j in range(i + 1, _system_object._total_fragments):
            fragments[i] = {"atoms": _system_object.get_fragments(i)._coordinates}
            comparison = (
                _system_object.get_fragments(i).center_of_mass
                == _system_object.get_fragments(j).center_of_mass
            )
            equal_arrays = comparison.all()
            if equal_arrays:
                r = random.random()
                fragments[j] = {
                    "atoms": _system_object.get_fragments(j)
                    .translate(0, r, r, r)
                    .rotate(0, r, r, r)
                    ._coordinates
                }
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
                fragments[j] = {"atoms": _system_object.get_fragments(j)._coordinates}
    return fragments


class SearchConfig:
    def __init__(
        self,
        system_object=None,
        search_methodology=1,
        basis="sto-3g",
        program_electronic_structure=1,
        outxyz="configurations.xyz",
    ) -> None:

        if system_object is None:
            sys.exit(
                "AttributeError system_object isn't a object of Molecule. It's None"
            )

        self._system_object = system_object

        #
        self._search_methodology = search_methodology
        self._search_name = self.search_name(self._search_methodology)

        # cost function
        self._basis = basis
        self._program_calculate_cost_function = program_electronic_structure
        self._func = self.program_cost_function(self._program_calculate_cost_function)

        # archivo de salida xyz con todas las configuraciones
        self._outxyz = outxyz

        # Siguiendo propuesta de Juan, con ediciones, para definir el bounds
        # al parcer con este bounds no se alejan las moleculas en shgo pero
        # con dual_annealing no se evita todavía que se alejen
        sphere_radius = self._system_object._total_atoms * 1.5 * 0.5
        discretization = sphere_radius / 1.6
        bound_translate = [
            (-discretization, discretization),
            (-discretization, discretization),
            (-discretization, discretization),
        ]
        bound_rotate = [(0, 360), (0, 360), (0, 360)]
        bound_translate = bound_translate * self._system_object._total_fragments
        bound_rotate = bound_rotate * self._system_object._total_fragments
        self._bounds = bound_translate + bound_rotate

        # verificar superposición de las moleculas
        self._system_object = Molecule(*overlaping(self._system_object).values())

    def bounds(self):
        return self._bounds

    def search_name(self, search_name):
        if self._search_methodology == 1:
            return "dual_annealing from Scipy"
        if self._search_methodology == 2:
            return "shgo from Scipy"
        if self._search_methodology == 3:
            return "Bayesiana"

    def program_cost_function(self, _program_calculate_cost_function):
        if _program_calculate_cost_function == 1:
            return hamiltonian_pyscf

    def run(self, **kwargs):
        if self._search_methodology == 1:
            self.DA(**kwargs)
        if self._search_methodology == 2:
            print("*** Calculo Realizado con SHGO from Scipy ***")
            self.SHGO(**kwargs)

    def DA(self, **kwargs):
        with open(self._outxyz, "w") as outxyz:
            self._search = solve_dual_annealing(
                self._func,
                self._bounds,
                self._system_object,
                args=(
                    self._basis,
                    self._system_object,
                    outxyz,
                    self._search_methodology,
                ),
                **kwargs
            )

    def SHGO(self, **kwargs):
        with open(self._outxyz, "w") as outxyz:

            self._search_methodology = 2
            shgo_tolerance_dict = {"ftol": 1e-6}
            shgo_options_dict = {"options": shgo_tolerance_dict}

            self._search = shgo(
                self._func,
                bounds=self._bounds,
                minimizer_kwargs=shgo_options_dict,
                sampling_method="sobol",
                args=(
                    self._basis,
                    self._system_object,
                    outxyz,
                    self._search_methodology,
                ),
                **kwargs
            )
