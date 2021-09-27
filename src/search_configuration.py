from m_dual_annealing import *

class SearchConfig:
    def __init__(self, system_object, search_setting) -> None:
        self._system_object = system_object
        self._bounds = search_setting["bounds"]
        self._search_methodology = search_setting["search_methodology"]
        self._basis = search_setting["basis"]
        self._there_is_molecules = ["type of fragments"] #0: todos átomos, 1: hay moléculas
        self._search_name = ["search_methodology_name"]
        #cost function
        self._func = heisenberg
        #
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

    def search_name(self):
        if self._search_methodology == 0:
            self._search_name = "dual_annealing from Scipy"
        if self._search_methodology == 1:
            self._search_name = "Bayesiana"

    def run(self, **kwargs):
        if self._search_methodology == 0:
            self._search = solve_dual_annealing(self._func, self._bounds,
                                                self._ascec_activation,
                                                self._there_is_molecules,
                                                self._system_object,
                                                **kwargs)