from m_dual_annealing import *

class SearchConfig:
    def __init__(self, system_object, search_setting) -> None:
        self._system_object = system_object
        self._bounds = search_setting["bounds"]
        self._search_methodology = search_setting["search_methodology"]
        self._basis = search_setting["basis"]
        self._there_is_molecules = search_setting["type of fragments"] #0: todos átomos, 1: hay moléculas
        #cost function
        self._program_calculate_cost_function = search_setting["program_cost_function"]

        self._func = self.program_cost_function(self._program_calculate_cost_function)
        self._search_name = self.search_name(self._search_methodology)

        #por trabajar
        self._ascec_activation = False
        self._seed = None
        self._NT = None
        self._T0 = None
        self._dT = None
        self._mxcycle = None
        self._local_search_options = None
        self._no_local_search = None
        self._visit_regions = None
        self._accept = None
        self._x0 = None
        self._args = None
        self._callback = None

    def bounds(self):
        return self._bounds

    def search_name(self, search_name):
        if self._search_methodology == 0:
            return "dual_annealing from Scipy"
        if self._search_methodology == 1:
            return "Bayesiana"

    def program_cost_function(self, _program_calculate_cost_function):
        if _program_calculate_cost_function == 0 or _program_calculate_cost_function == None:
            return hamiltonian_pyscf
#        if _program_calculate_cost_function == 1:
#            return hamiltonian_horton

    def run(self, **kwargs):
        if self._search_methodology == 0:
            self._search = solve_dual_annealing(self._func, self._bounds,
                                                self._ascec_activation,
                                                self._there_is_molecules,
                                                self._system_object,
                                                **kwargs)
#        if self._se... :
#            ... = abc(x,y,z, ...)
#        if self._.....:
#           ... = bayesiana(x,y,z, ....)