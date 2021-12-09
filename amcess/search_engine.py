from copy import deepcopy

import numpy as np
import scipy
from scipy.optimize import shgo
from scipy.optimize import dual_annealing

from .energy_functions import compute_energy

from amcess.ascec import Ascec
from amcess.base_molecule import Cluster
from amcess.electronic_energy import ElectronicEnergy
from amcess.gaussian_process import solve_gaussian_processes

METHODS = {
    "ASCEC": Ascec,
    "dual_annealing": dual_annealing,
    "SHGO": shgo,
    "Bayesian": solve_gaussian_processes,
}


def lennard_jones(r, epsilon=1.0, sigma=1.0):
    """
    Lennard-Jones potential.

    Parameters
    ----------
    r : float
        Distance between two atoms.
    epsilon : float
        Depth of potential well.
    sigma : float
        Width of potential well.

    Returns
    -------
    float
        Lennard-Jones potential.
    """
    return 4 * epsilon * ((sigma / r) ** 12 - (sigma / r) ** 6)


COST_FUNCTIONS = {
    "pyscf": ElectronicEnergy.pyscf,
    "Lennard_Jones": lennard_jones,
}


def simulated_annealing(
    cluster_settings: dict,
    energy_settings: dict,
    annealing_settings: dict,
    temperature_settings: dict,
):

    # Boltzmann constant
    # KB: float = 1.38064852e-23  # J/K
    KB: float = 3.166811563e-6  # Eh/K

    # ------------------------------------------------------
    # cluster setting
    if not cluster_settings:
        raise ValueError("cluster_settings is empty")

    cluster = cluster_settings.get("cluster", None)

    if not cluster:
        raise ValueError("molecular cluster is empty")

    max_closeness = cluster_settings.get("max_closeness", 0.8)
    max_step = cluster_settings.get("max_step", 0.5)
    max_rotation = cluster_settings.get("max_rotation", 10)

    # ------------------------------------------------------
    # Annealing setting
    annealing_criterion = annealing_settings.get("criterion", "metropolis")
    file_name = annealing_settings.get("file_name", "annealing")
    temperature_grid = annealing_settings.get("temperature_grid", "arithmetic")
    max_cycle = annealing_settings.get("max_cycle", 1000)

    # ------------------------------------------------------
    # Temperature setting
    initial_temperature = temperature_settings.get("initial_temperature", 100)
    final_temperature = temperature_settings.get("final_temperature", 500)
    total_temperatures = temperature_settings.get("total_temperatures", 100)

    # temperatures progession (cooling down for Annealing, going tf -> to)
    if temperature_grid == "arithmetic":
        temperature_list = [
            final_temperature
            + (initial_temperature - final_temperature)
            * (i / (total_temperatures - 1))
            for i in range(total_temperatures)
        ]
    elif temperature_grid == "geometric":
        temperature_list = [
            final_temperature
            * (initial_temperature / final_temperature)
            ** (i / (total_temperatures - 1))
            for i in range(total_temperatures)
        ]

    # ------------------------------------------------------
    # Creating output files
    configuration_file: str = file_name + "_configuration.xyz"
    with open(configuration_file, "w") as f:
        f.write("")

    energy_data_file: str = file_name + "_energy.dat"
    with open(energy_data_file, "w") as f:
        # --------0123456789012345678901234567890123456789012345678901
        f.write("#Accepted conf #    Temperature [K]     Energy [Eh]\n")

    lowest_energy_file: str = file_name + "_lowest_energy.xyz"
    with open(lowest_energy_file, "w") as f:
        f.write("")

    # ------------------------------------------------------
    # Initial configuration
    current_configuration = deepcopy(cluster)
    random_gen = current_configuration.random_generator

    energy_evaluations = 1
    current_energy = compute_energy(
        molecule=current_configuration,
        energy_function=energy_settings,
    )

    configuration_accepted = False
    lowest_energy = current_energy
    total_accepted = 0
    total_rejected = 0
    configuration_number = 0

    # ------------------------------------------------------
    # Printing setting information
    print("#" + "-" * 80)
    print(
        "Simulated Annealing:"
        f"\n\tCooling progression = {temperature_grid}"
        f"\n\tInitial temperature= {initial_temperature} K"
        f"\n\tFinal temperature= {final_temperature} K"
        f"\n\tTotal number of temperatures = {total_temperatures} points"
        f"\n\tAcceptance criterion = {annealing_criterion}"
        f"\n\tMaximum search per temperature = {max_cycle}"
    )
    print("#" + "-" * 80)
    print(
        f"Molecular Cluster:"
        f"\n\tOutput file name (prefix) = {file_name}"
        f"\n\tLowest energy configuration file = {lowest_energy_file}"
        f"\n\tAccepted configuration file = {configuration_file}"
        f"\n\tEnergy data file = {energy_data_file}"
        f"\n\tSphere radius = {current_configuration.sphere_radius}"
        f"\n\tSphere center = {current_configuration.sphere_center}"
        f"\n\tMaximum translation = {max_step}"
        f"\n\tMaximum rotation = {max_rotation}"
        f"\n\tMaximum closeness between molecules = {max_closeness}"
        f"\n\tNumber of molecules = {current_configuration.total_molecules}"
        f"\n\tNumber of atoms = {current_configuration.total_atoms}"
    )

    print("#" + "-" * 80)
    print("Energy function:")
    for k, v in energy_settings.items():
        print(f"\t{k} = {v}")

    print("#" + "-" * 80)
    print("Running simulated annealing...")

    # ------------------------------------------------------
    # Annealing
    for count, Temperature in enumerate(temperature_list):
        print(
            f"\rTemperature: {Temperature:.3g} K\t"
            f"progress {count:d}/{total_temperatures:d}"
            f" ({100.0*(count/total_temperatures):.3f}%)",
            end="",
        )

        for _ in range(max_cycle):
            configuration_number += 1

            random_choice = random_gen.choice(
                current_configuration.total_molecules
            )

            next_configuration = current_configuration.move_molecule(
                molecule=random_choice,
                max_closeness=max_closeness,
                max_step=max_step,
                max_rotation=max_rotation,
            )

            next_energy = compute_energy(
                molecule=next_configuration,
                energy_function=energy_settings,
            )
            energy_evaluations += 1

            if next_energy < current_energy:
                configuration_accepted = True
                accepted = "lower energy configuration"
            else:
                boltzmann_factor = np.exp(
                    -(next_energy - current_energy) / Temperature / KB
                )

                if annealing_criterion == "metropolis":
                    # Metropolis criterion
                    criterion = random_gen.uniform(0, 1)
                elif annealing_criterion == "delta_energy":
                    # ASCEC original criterion
                    criterion = (
                        np.abs(next_energy - current_energy) / next_energy
                    )

                if boltzmann_factor > criterion:
                    configuration_accepted = True
                    accepted = (
                        f"Boltzmann after {configuration_number} evaluations"
                    )

            if configuration_accepted:
                total_accepted += 1
                current_energy = next_energy
                current_configuration = next_configuration

                with open(configuration_file, "a") as f:
                    f.write(f"\t{next_configuration.total_atoms:d}\n")
                    f.write(f"energy: {next_energy:+.20e} Eh\n")
                    f.write(next_configuration.write_atoms)

                with open(energy_data_file, "a") as f:
                    f.write(
                        f"{configuration_number:<20d}"
                        f"{Temperature:<18.3f} "
                        f"{next_energy:<+.20e}"
                        f"\t{accepted}\n"
                    )

                if next_energy < lowest_energy:
                    lowest_energy = next_energy

                    with open(lowest_energy_file, "w") as f:
                        f.write(f"\t{next_configuration.total_atoms:d}\n")
                        f.write(f"energy: {next_energy:+.20e} Eh\n")
                        f.write(next_configuration.write_atoms)

                break
        else:
            total_rejected += 1
            msg1 = "Rejected"
            with open(energy_data_file, "a") as f:
                f.write(
                    f"{configuration_number:<20d}"
                    f"{Temperature:<18.3f} "
                    f"{msg1:<20s}\n"
                )
    result_message: str = (
        f"\n\n----------------------------------------------------------------"
        f"\nTotal configuration accepted: {total_accepted}"
        f"\nwe evaluate the electronic energy: {energy_evaluations} times"
        f"\nLowest energy: {lowest_energy} Eh"
        f"\nAccepted ratio: {total_accepted/energy_evaluations:.2f}"
        "\n\nSIMULATED ANNEALING DONE \n\n"
    )

    return result_message


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
        basis : string
            Label of basis set
        program_electronic_structure : int
            Integer associated with the program to make the
            electronic structure calculations
        outxyz : string
            Name of the output xyz with coordinates of the
            configurations accepts

    Returns
    -------
        Output xyz with coordinates and electronic structure

    Raises
    ------
        TypeError
            System_object isn't difinite
            AttributeError system_object isn't difinite as
            an object Cluster
    """

    def __init__(
        self,
        system_object: Cluster = None,
        search_methodology: str = "ASCEC",
        methodology: str = "HF",
        basis: str = "sto-3g",
        # program_electronic_structure: int = 1,
        tolerance_contour_radius: float = 1.0,
        outxyz: str = "configurations.xyz",
        # costo_function="pyscf_HF",
        program_cost_function="pyscf",
    ) -> None:

        # Verfication and assigment of variables (type, value)
        self.system_object = system_object
        self.search_type = search_methodology
        self.methodology = methodology
        self.basis_set = basis
        # self.cost_function_number = program_electronic_structure
        self.output_name = outxyz
        self.tolerance_contour_radius = tolerance_contour_radius
        # self.func_costo = costo_function
        self.cost_function = program_cost_function
        # Check Overlaping
        self._system_object.initialize_cluster()

        # Build bounds, format for scipy functions
        if system_object._sphere_radius is None:
            self.spherical_contour_cluster()
        else:
            # ! arreglar esto, al fin que va a pasar con sphere_center y radius
            # ! y eso articular con bounds
            self._sphere_center = system_object.sphere_center
            self._sphere_radius = system_object.sphere_radius
            new_radius_t = (
                self._sphere_radius
            )  # self._tolerance_contour_radius + new_radius

            bound_translate = [
                (-new_radius_t, new_radius_t),
                (-new_radius_t, new_radius_t),
                (-new_radius_t, new_radius_t),
            ]
            bound_rotate = [(-180, 180), (-180, 180), (-180, 180)]

            bound_translate = bound_translate * (
                self._system_object.total_molecules - 1
            )
            bound_rotate = bound_rotate * (
                self._system_object.total_molecules - 1
            )

            self._bounds = bound_translate + bound_rotate

    # ===============================================================
    # Decorators
    # ===============================================================

    def bounds_sphere_change(function_change_radius):
        def new_bounds(self, new_radius):
            """
            Define the bounds for the optimization algorithm

            Returns
            -------
                bounds : list
                    Bounds for the optimization algorithm
            """
            new_radius_t = self._tolerance_contour_radius + new_radius

            bound_translate = [
                (-new_radius_t, new_radius_t),
                (-new_radius_t, new_radius_t),
                (-new_radius_t, new_radius_t),
            ]
            bound_rotate = [(-180, 180), (-180, 180), (-180, 180)]

            bound_translate = bound_translate * (
                self._system_object.total_molecules - 1
            )
            bound_rotate = bound_rotate * (
                self._system_object.total_molecules - 1
            )

            self._bounds = bound_translate + bound_rotate

            return function_change_radius(self, new_radius)

        return new_bounds

    def init_electronic_energy(function_minimization):
        def wrapper(self, **kwargs):
            if self._search_methodology != "ASCEC":
                self._obj_ee = ElectronicEnergy(
                    self._system_object,
                    self._search_methodology,
                    self._sphere_center,
                    self._sphere_radius,
                    self._methodology,
                    self._basis_set,
                )
            return function_minimization(self, **kwargs)

        return wrapper

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
                "\n\nThe new name to basis set is not a string"
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
    def tolerance_contour_radius(self):
        return self._tolerance_contour_radius

    @tolerance_contour_radius.setter
    def tolerance_contour_radius(self, new_tol_radius: float):
        if not isinstance(new_tol_radius, float):
            raise TypeError(
                "\n\nThe new tolerance radius is not a float"
                f"\nplease, check: '{type(new_tol_radius)}'\n"
            )
        self._tolerance_contour_radius = new_tol_radius

    @property
    def sphere_center(self) -> tuple:
        return self._sphere_center

    @sphere_center.setter
    def sphere_center(self, new_center: tuple) -> None:
        if not isinstance(new_center, tuple):
            raise TypeError(
                "\n\nThe Sphere center must be a tuple with three elements: "
                "(float, float, float)"
                f"\nplease, check: '{type(new_center)}'\n"
            )
        if len(new_center) != 3:
            raise ValueError(
                "\n\nThe Sphere center must be a tuple with three elements: "
                "(float, float, float)"
                f"\nplease, check: '{new_center}'\n"
            )

        self._sphere_center = new_center

    @property
    def sphere_radius(self) -> float:
        return self._sphere_radius

    @sphere_radius.setter
    @bounds_sphere_change
    def sphere_radius(self, new_radius: float) -> None:
        if not isinstance(new_radius, float):
            raise TypeError(
                "\n\nThe Sphere  Radius must be a float"
                f"\nplease, check: '{type(new_radius)}'\n"
            )
        self._sphere_radius = new_radius + self._tolerance_contour_radius
        if self._sphere_radius <= self._tolerance_contour_radius:
            raise ValueError(
                "\n\nThe Sphere Radius more tolerance must be larger than 1 A"
                f"\nplease, check: '{new_radius}'\n"
            )

    # @property
    # def cost_function_number(self):
    #    return self._program_calculate_cost_function

    @property
    def cost_function(self):
        return self._cost_function

    @cost_function.setter
    def cost_function(self, new_cost_function):
        self._cost_function = (
            new_cost_function
            if callable(new_cost_function)
            else COST_FUNCTIONS[new_cost_function]
        )

    # ===============================================================
    # Methods
    # ===============================================================

    def spherical_contour_cluster(self, new_tol: float = None):
        """
        Define a spherical outline that contains our cluster

        Parameters
        ----------
            tolerance : float
                Tolerance with the radius between the mass center to the
                furthest atom

        Returns
        -------
            sphere_center : tuple
                Mass center of the biggest molecule
            sphere_radius : float
                Radius between the sphere center to the furthest atom

        """
        if new_tol is not None:
            self._tolerance_contour_radius = new_tol

        max_distance_cm = 0.0
        molecule = 0
        max_atoms = 0

        # The biggest molecule
        for i in range(self._system_object.total_molecules):
            if self._system_object.get_molecule(i).total_atoms > max_atoms:
                max_atoms = self._system_object.get_molecule(i).total_atoms
                molecule = i

        self._sphere_center = self._system_object.get_molecule(
            molecule
        ).center_of_mass

        # Move the biggest molecule to initio in the cluster object,
        # if is necessary
        if molecule != 0:
            new_geom = dict()
            for i in range(self._system_object.total_molecules):
                if i == 0:
                    new_geom[i] = self._system_object.get_molecule(molecule)
                elif i == molecule:
                    new_geom[i] = self._system_object.get_molecule(0)
                else:
                    new_geom[i] = self._system_object.get_molecule(i)

            self.system_object = Cluster(
                *new_geom.values(), sphere_center=self._sphere_center
            )

        # Radius between the sphere center to the furthest atom
        for xyz in self._system_object.coordinates:
            temp_r = np.linalg.norm(
                np.asarray(self._sphere_center) - np.asarray(xyz)
            )
            if temp_r > max_distance_cm:
                max_distance_cm = temp_r

        self.sphere_radius = max_distance_cm

    # def program_cost_function(self, _program_calculate_cost_function):
    #     """
    #     Assign the name of the cost function, which is associated with
    #     determined program to calculate the energy of the electronic
    #     structure

    #     Parameters
    #     ----------
    #         _program_calculate_cost_function : int
    #             Integer associated with the program to calculate the cost
    #             and methodology (Hamiltonian, Functional, etc)

    #     """
    #     if _program_calculate_cost_function == 1:
    #         print(
    #             "\n\n"
    #             "*** Cost function is Hartree--Fock implemented into pyscf ***"
    #             "\n\n"
    #         )

    @init_electronic_energy
    def run(self, **kwargs):
        """
        Alternative to execute the searching methodologies in METHODS

        Parameters
        ----------
            **kwargs : dict
                Dictionary with the parameters to be used in the
                search methodologies
        """
        func = (
            self._search_methodology
            if callable(self._search_methodology)
            else METHODS[self._search_methodology]
        )

        if self._search_methodology == "ASCEC":
            print("*** Minimization: ASCEC ***")
            self._search = func(
                object_system=self._system_object,
                search_type=self._search_methodology,
                sphere_center=self._sphere_center,
                sphere_radius=self._sphere_radius,
                methodology=self._methodology,
                basis_set=self._basis_set,
                call_function=1,
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

            self._search = func(
                # self._func,
                self._cost_function,
                bounds=self._bounds,
                **kwargs,
            )
            self._obj_ee.write_to_file(self.output_name)
