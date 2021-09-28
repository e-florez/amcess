from m_dual_annealing import *

class SearchConfig:
    def __init__(self, system_object, search_setting) -> None:
        self._system_object = system_object
        self._bounds = search_setting["bounds"]
        self._search_methodology = search_setting["search_methodology"]
        self._basis = search_setting["basis"]
        self._there_is_molecule = search_setting["type of fragments"] #0: todos átomos, 1: hay moléculas
        #cost function
        self._program_calculate_cost_function = search_setting["program_cost_function"]
        self._func = self.program_cost_function(self._program_calculate_cost_function)
        self._search_name = self.search_name(self._search_methodology)

    def bounds(self):
        return self._bounds

    def search_name(self, search_name):
        if self._search_methodology == 1:
            #args dual_annealing
            self._NT = 1000
            self._T0 = 5230.0
            self._dT = 2e-5
            self._mxcycle = 10000000.0
            self._visit_regions = 2.62
            self._accept = -5.0
            self._no_local_search = False
            self._args = (self._basis, self._system_object, 0) #0 asociar al llamado
            self._local_search_options = {}
            self._x0 = [0.00000000e+00,  0.00000000e+00, -9.46852623e-06, 1., 1., 0.99999053]
            self._callback = None
            self._ascec_activation = False
            self._seed = None
            return "dual_annealing from Scipy"
        if self._search_methodology == 2:
            return "Bayesiana"

    def program_cost_function(self, _program_calculate_cost_function):
        if _program_calculate_cost_function == 1:
            return hamiltonian_pyscf
#        if _program_calculate_cost_function == 1:
#            return hamiltonian_horton

    def run(self, **kwargs):
        if self._search_methodology == 1:
            #self._search = solve_dual_annealing(self._func, self._bounds,
            #                                    self._ascec_activation,
            #                                    self._there_is_molecule,
            #                                    self._system_object,
            #                                    **kwargs)
            self._search = solve_dual_annealing(self._func, self._bounds,
                                                self._ascec_activation,
                                                self._there_is_molecule,
                                                self._system_object,
                                                self._seed,
                                                self._NT, self._T0,
                                                self._dT, self._mxcycle,
                                                self._local_search_options,
                                                self._no_local_search,
                                                self._visit_regions,
                                                self._accept,
                                                self._x0, self._args,
                                                self._callback)


#        if self._se... :
#            ... = abc(x,y,z, ...)
#        if self._.....:
#           ... = bayesiana(x,y,z, ....)