import context
from m_dual_annealing import *

class SearchConfig:
    def __init__(self, search_setting) -> None:
        self._bounds = search_setting["bounds"]
        self._search_methodology = search_setting["search_methodology"]
        self._search_name = ["search_methodology_name"]
        self._search = None

    def bounds(self):
        return self._bounds

    def search_name(self):
        if self._search_methodology == 0:
            self._search_name = "dual_annealing from Scipy"

    def run(self):
        self._search = solve_dual_annealing(heisenberg, self._bounds, NT=3,
                                     no_local_search=False, visit_regions=2.9)