from m_dual_annealing import *
from scipy.optimize import shgo


class SearchConfig:
    def __init__(self, system_object, search_setting) -> None:
        self._system_object = system_object
        #self._bounds = search_setting["bounds"]
        self._search_methodology = search_setting["search_methodology"]
        self._basis = search_setting["basis"]
        self._there_is_molecule = search_setting[
            "type of fragments"
        ]  # 0: todos átomos, 1: hay moléculas
        # cost function
        self._program_calculate_cost_function = search_setting["program_cost_function"]
        self._func = self.program_cost_function(
            self._program_calculate_cost_function)
        self._search_name = self.search_name(self._search_methodology)

        # archivo de salida xyz con todas las configuraciones
        self._outxyz = 'configurations.xyz'

        # Siguiendo propuesta de Juan, con ediciones, para definir el bounds
        # al parcer con este bounds no se alejan las moleculas en shgo pero
        # con dual_annealing no se evita todavía que se alejen
        sphere_radius = self._system_object._total_atoms*1.5*0.5
        discretization = sphere_radius / 1.6
        bound_translate = [(-discretization, discretization),
                           (-discretization, discretization),
                           (-discretization, discretization),
                           ]
        bound_rotate = [(0, 360), (0, 360), (0, 360)]
        bound_translate = bound_translate * self._system_object._total_fragments
        bound_rotate = bound_rotate * self._system_object._total_fragments
        self._bounds = bound_translate + bound_rotate

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
        with open(self._outxyz, "w") as outxyz:
            if self._search_methodology == 1:
                self._search = solve_dual_annealing(
                    self._func,
                    self._bounds,
                    self._system_object,
                    self._there_is_molecule,
                    args=(self._basis, self._system_object,
                          outxyz, self._search_methodology),
                    **kwargs
                )
            if self._search_methodology == 2:
                print("*** Calculo Realizado con SHGO from Scipy ***")
                shgo_tolerance_dict = {"ftol": 1e-6}
                shgo_options_dict = {"options": shgo_tolerance_dict}

                self._search = shgo(
                    self._func,
                    bounds=self._bounds,
                    minimizer_kwargs=shgo_options_dict,
                    sampling_method='sobol',
                    args=(self._basis, self._system_object,
                          outxyz, self._search_methodology),
                    **kwargs
                )
