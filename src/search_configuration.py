from m_dual_annealing import *
from scipy.optimize import shgo


class SearchConfig:
    def __init__(self, system_object, search_setting) -> None:
        self._system_object = system_object
        self._bounds = search_setting["bounds"]
        self._search_methodology = search_setting["search_methodology"]
        self._basis = search_setting["basis"]
        self._there_is_molecule = search_setting[
            "type of fragments"
        ]  # 0: todos átomos, 1: hay moléculas
        # cost function
        self._program_calculate_cost_function = search_setting["program_cost_function"]

        self._func = self.program_cost_function(self._program_calculate_cost_function)
        self._search_name = self.search_name(self._search_methodology)

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
            self._search = solve_dual_annealing(
                self._func,
                self._bounds,
                self._system_object,
                self._there_is_molecule,
                args=(self._basis, self._system_object),  # 0 asociar al llamado
                **kwargs
            )
        if self._search_methodology == 2:
            print("*** Calculo Realizado con SHGO from Scipy ***")
            self._search = shgo(
                self._func,
                self._bounds,
                args=(self._basis, self._system_object),  # 0 asociar al llamado
                **kwargs
            )
