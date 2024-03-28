from copy import deepcopy

import numpy as np
from scipy.spatial.transform import Rotation
from rdkit.Chem.rdchem import Mol

from amcess.molecule import Molecule


class Cluster(Molecule):
    """
    Create a Cluster with molecules/atoms to move and rotate
    using spherical boundary conditions (SBC).
    The format of the INPUT coordinates is as follows (any):

    1. Dictionary type: {"atoms": [(<element> <X> <Y> <Z>), ...]}
    2. List type: [(<element> <X> <Y> <Z>), ...]
    3. Molecule/Cluster type (Objects)

    .. rubric:: Parameters

    args : List, Dict, Molecule, Cluster
        coordinates of each molecule/atom comma separates (support +,-,*)
    freeze_molecule : integer, optional
        fixing molecule to NOT move or rotate, by default NEGATIVE
        integer means all molecules can be moved freely
    sphere_radius : float, optional
        radius for the spherical boundary condition, by default None
    sphere_center : tuple, optional
        Center of the sphere, by default (0, 0, 0)
    seed : int, optional
        seed to initialize the random generator function, by default None

    .. rubric:: Raises

    TypeError
        for a wrong input argument
    """

    def __init__(
        self,
        *args,
        freeze_molecule: list = None,
        sphere_radius: float = None,
        sphere_center: tuple = (0, 0, 0),
        seed: int = None,
    ):
        self._cluster_dict = dict()
        self._multiplicity = 1
        self._charge = 0

        # fixing molecule to NOT move or rotate
        # initialize with an empty list
        self._freeze_molecule = (  # noqa
            [] if freeze_molecule is None else freeze_molecule
        )

        self._sphere_radius = sphere_radius
        self._sphere_center = sphere_center

        # initialize random generator
        if not seed:
            self._seed = np.random.randint(0, 999999)
        else:
            self._seed: int = seed

        # ----------------------------------------------------
        # attrs post-initialization

        cluster_atoms: list = list()

        # for count, mol in enumerate(args):
        for mol in args:
            size: int = len(self._cluster_dict)
            if isinstance(mol, Cluster):
                for j in mol._cluster_dict:
                    self._cluster_dict[size + j] = mol._cluster_dict[j]
                self._charge += mol.GetMolCharge
                cluster_atoms += mol.GetMolList
                # restarting the loop
                continue
            elif isinstance(mol, Molecule):
                new_molecule = mol
            elif isinstance(mol, Mol):
                new_molecule = Molecule(mol)
            elif isinstance(mol, dict):
                new_molecule = Molecule(mol)
            elif isinstance(mol, list):
                new_molecule = Molecule(mol)
            else:
                raise TypeError(
                    "\nOnly type 'Molecule', list or dict to initialize"
                    "\n\t- Dict: {'atoms': [(<element> <X> <Y> <Z>), ...]}"
                    "\n\t- List: [(<element> <X> <Y> <Z>), ...]"
                    f"\nyou have a NOT valid '{type(mol)}', check: \n{mol}"
                )

            cluster_atoms += new_molecule.GetMolList #new_molecule.atoms
            # ! how is computed the cluster total multiplicity?
            self._charge += new_molecule.GetMolCharge   #charge
            self._cluster_dict[size] = new_molecule

        # initializing Cluster as a 'Molecule' (sum of all individual ones)
        super().__init__(
            atoms=cluster_atoms,
            charge=self._charge,
            multiplicity=self._multiplicity,
        )

    # ===============================================================
    # MAGIC METHODS
    # ===============================================================
    def __add__(self, other):
        """Adding two molecules/clusters, return a new Cluster object"""
        # ! Martin
        # ! ver en que puede se diferenciar de Molecule
        # * la idea es que debe ser similar para que ambos
        # * creen una nueva instancia nueva de Cluster
        new_cluster: dict = {}
        for count, m in self.GetClusterDict.items():
            new_cluster[count] = m.GetMolDict
        new_cluster[self.GetTotalMol] = Molecule(other).GetMolDict
        return Cluster(*new_cluster.values())

    def __mul__(self, value: int):
        """multiply the cluster by a number"""
        return value * self

    def __rmul__(self, value: int):
        """to replicate a molecule

        Parameters
        ----------
        value : int
            quantity to replicate Molecue

        Return
        ------
        Cluster : object
            summing or multiplying Molecule classes produce a Cluster class
        """
        if value < 1 or not isinstance(value, int):
            raise ValueError(
                "\nMultiplier must be and integer larger than zero"
                f"\ncheck --> '{value}'"
            )

        new_cluster = self
        for _ in range(value - 1):
            new_cluster += new_cluster

        return new_cluster

    def __str__(self):
        """print the cluster"""
        cluster_dict: dict = self._cluster_dict

        cluster_string: str = (
            f"Cluster of ({self.GetTotalMol}) molecules"
            f" and ({self.GetNumAtoms()}) total atoms\n"
        )
        for key, molecule in cluster_dict.items():
            atoms = molecule.GetMolList
            cluster_string += f" #{key}: molecule with {len(atoms)} atoms:\n"
            cluster_string += f"     --> atoms: {atoms}\n"
            charge = molecule.GetMolCharge
            cluster_string += f"     --> charge: {charge:>+}\n"
            multiplicity = molecule.GetMolMultiplicity
            cluster_string += f"     --> multiplicity: {multiplicity}\n"

        return cluster_string

    # ===============================================================
    # PROPERTIES
    # ===============================================================
    @property
    def GetClusterDict(self) -> dict:
        """return the cluster dictionary"""
        return self._cluster_dict
    
    @property
    def GetClusterList(self) -> list:
        """return the cluster list"""
        return [mol for mol in self._cluster_dict.values()]

    @property
    def GetFreezeMol(self) -> int:
        """return a list with freezed molecules"""
        return self._freeze_molecule

    @GetFreezeMol.setter
    def SetFreezeMol(self, values) -> None:
        """set the freeze molecules"""
        if isinstance(values, list):
            self._freeze_molecule = values
        elif isinstance(values, int):
            self._freeze_molecule = [values]
        else:
            raise TypeError(f"values {type(values)} can only be a list or int")

    @property
    def GetRandomGen(self) -> np.random.Generator:
        """return the random generator"""
        # self._random_gen: np.random.Generator = np.random.default_rng(seed)
        return np.random.default_rng(self.GetSeed)

    @property
    def GetSeed(self) -> int:
        """return the seed for the random generator"""
        return self._seed

    @GetSeed.setter
    def SetSeed(self, new_seed: int) -> None:
        """set the seed for the random generator"""
        self._seed = new_seed
        self._random_gen = np.random.default_rng(new_seed)

    @property
    def GetSphereCenter(self) -> tuple:
        """return the sphere center for the Cluster boundary conditions"""
        return self._sphere_center

    @GetSphereCenter.setter
    def SetSphereCenter(self, new_center: tuple) -> None:
        """set the sphere center for the Cluster boundary conditions"""
        if len(new_center) != 3:
            raise ValueError(
                "\n\nThe Sphere center must be a tuple with three elements: "
                "(float, float, float)"
                f"\nplease, check: '{new_center}'\n"
            )

        self._sphere_center = new_center

    @property
    def GetSphereR(self) -> float:
        """return the sphere radius for the Cluster boundary conditions"""
        return self._sphere_radius

    @GetSphereR.setter
    def SetSphereR(self, new_radius: float) -> None:
        """set the sphere radius for the Cluster boundary conditions"""
        if not isinstance(new_radius, (int, float)) or new_radius < 0.9:
            raise ValueError(
                "\n\nThe Sphere  Radius must be larger than 1 Angstrom"
                f"\nplease, check: '{new_radius}'\n"
            )

        self._sphere_radius = new_radius

    @property
    def GetTotalMol(self) -> int:
        """return the total number of molecules in the cluster"""
        return len(self._cluster_dict)

    # ===============================================================
    # METHODS
    # ===============================================================
    def GetMol(self, molecule) -> Mol:
        """"""
        return self.GetClusterDict[molecule]

    @staticmethod
    def Overlapping(
        first_coordinates: list,
        second_coordinates: list,
        max_closeness: float = 1.0,
    ) -> bool:
        """pair-wise checking if any overlapping among points
        with a radius defined by `max_closeness`

        .. rubric:: Parameters

        first_coordinates : list
            list of tuples [(float, float, float), ...]
        second_coordinates : list
            list of tuples [(float, float, float), ...]
        max_closeness : float, optional
            maximun closeness between two pairs, by default 1.0

        .. rubric:: Returns

        bool
            True if two point are closer than `max_closeness`
        """
        # ! Martin
        # ! itertools para optimizar los for
        for first_atom in first_coordinates:
            for second_atom in second_coordinates:
                distance = np.linalg.norm(
                    np.asarray(first_atom) - np.asarray(second_atom)
                )

                if distance < max_closeness:
                    return True

        return False

    def InitializeCluster(self, max_closeness: float = 1.0) -> object:
        """Create a new cluster object which any atom is overlapped

        .. rubric:: Parameters

        max_closeness : float, optional
            maximun closeness between two pairs, by default 1.0

        .. rubric:: Returns

        Cluster : object
            returns a new Cluster object
        """
        # center of mass coordinates
        sc_x = self.GetSphereCenter[0]
        sc_y = self.GetSphereCenter[1]
        sc_z = self.GetSphereCenter[2]

        # initializing a new cluster moving the first molecule
        # to the center of the cluster sphere

        molecule = Cluster(self.GetClusterDict[0])
        new_cluster = molecule.TranslateMol(0, sc_x, sc_y, sc_z)

        for i in range(1, self.GetTotalMol):
            # moving the next single molecule into the cluster sphere
            molecule = Cluster(self.GetClusterDict[i]).TranslateMol(0, sc_x, sc_y, sc_z)
            if Cluster.Overlapping(  # noqa
                molecule.GetAtomicCoordinates, new_cluster.GetAtomicCoordinates
                #molecule.coordinates, new_cluster.coordinates
            ):
                new_cluster += molecule
                new_cluster = new_cluster.MoveMol(
                    i,
                    max_step=None,
                    max_rotation=None,
                    max_closeness=max_closeness,
                )
            else:
                new_cluster += molecule

        return self.__class__(
            *new_cluster.GetClusterDict.values(),
            freeze_molecule=self.GetFreezeMol,
            sphere_radius=self.GetSphereR,
            sphere_center=self.GetSphereCenter,
        )

    def MoveMol(
        self,
        molecule: int = 0,
        max_step: float = None,
        max_rotation: float = None,
        max_closeness: int = 1.0,
    ) -> object:
        """Moving (translating and rotating) randomly without overlapping
        any atom

        .. rubric:: Parameters

        molecule : int, optional
            molecule to move randomly, by default molecule with index zero (0)
        max_step : float, optional
            maximun value for any translation, by default None
        max_rotation : float, optional
            maximun angle fro any rotation, by default None
        max_closeness : float
            maximun closeness between any pair of atoms, by default 1.0 A

        .. rubric:: Returns

        object : Cluster
            returns a Cluster object whit a molecule moved to a random place
            without overlapping any other

        .. rubric:: Raises

        AttributeError : OverlappingError
            After serching for max_overlap_cycle and no place found for the
            molecule without overlapping any other
        """
        if not isinstance(max_closeness, (int, float)) or max_closeness < 0.1:
            raise ValueError(
                "\n\n Maximun closeness between any pair of atom must be"
                f" larger than '0.1' Angstrom\nPlease, check '{max_closeness}'"
            )

        if not max_step or not isinstance(max_step, (int, float)):
            max_step = 1.1 * max_closeness

        if not max_rotation or not isinstance(max_rotation, (int, float)):
            max_rotation = 30

        molecule_to_move: Cluster = Cluster(self.GetClusterDict[molecule])

        cluster_without_molecule: Cluster = self.RemoveMol(molecule)
        cluster_coordinates = cluster_without_molecule.GetAtomicCoordinates

        random_gen: np.random.Generator = self.GetRandomGen

        max_overlap_cycle: int = 10000

        for count in range(max_overlap_cycle):

            if count % 10 == 0:
                max_step *= 1.1
                max_rotation *= 1.1

            # angle between [0, max_rotation) degrees
            rotation_x = random_gen.uniform(-1, 1) * max_rotation
            rotation_y = random_gen.uniform(-1, 1) * max_rotation
            rotation_z = random_gen.uniform(-1, 1) * max_rotation

            # moving between [-max_step, +max_step] Angstrom
            tranlation_x = max_step * random_gen.uniform(-1, 1)
            tranlation_y = max_step * random_gen.uniform(-1, 1)
            tranlation_z = max_step * random_gen.uniform(-1, 1)

            new_molecule: Cluster = molecule_to_move.TranslateMol(
                0,
                tranlation_x,
                tranlation_y,
                tranlation_z,
            ).RotateMol(
                0,
                rotation_x,
                rotation_y,
                rotation_z,
            )

            molecule_coordinates: list = new_molecule.GetAtomicCoordinates

            overlap: bool = Cluster.Overlapping(
                molecule_coordinates,
                cluster_coordinates,
                max_closeness=max_closeness,
            )

            if not overlap:
                break
        # if overlapping and max_overlap_cycle reached
        else:
            raise AttributeError(
                "\n\n *** Overlapping Error ***"
                "\nat least one atom is overlapped with a distance"
                f" less than '{max_closeness}' Angstroms"
                "\nfor a cluster into a sphere of radius"
                f" '{self.GetSphereR}' Angstroms"
                # f"\nPlease, check: \n\n{self.xyz}"
            )

        cluster_dict: dict = deepcopy(self.GetClusterDict)
        cluster_dict[molecule] = new_molecule

        return self.__class__(
            *cluster_dict.values(),
            freeze_molecule=self.GetFreezeMol,
            sphere_radius=self.GetSphereR,
            sphere_center=self.GetSphereCenter,
        )

    def RemoveMol(self, molecule: int) -> object:
        """Removing molecule from cluster"""
        if molecule not in self.GetClusterDict:
            raise IndexError(
                f"\nMolecule with {self.GetTotalMol} total atoms "
                f"and index [0-{self.GetTotalMol - 1}]"
                f"\n molecule index must be less than {self.GetTotalMol}"
                f"\nCheck! You want to remove molecule with index {molecule}"
            )
        new_cluster: Cluster = self
        new_cluster_dict: dict = new_cluster.GetClusterDict
        del new_cluster_dict[molecule]
        # ! Martin
        # ! tener cuidado con el self porque se puede compartir Molecule
        # ! o Cluster
        # ! Preferible usar Cluster(...)
        return self.__class__(
            *new_cluster._cluster_dict.values(),
            freeze_molecule=new_cluster.GetFreezeMol,
            sphere_radius=new_cluster.GetSphereR,
            sphere_center=new_cluster.GetSphereCenter,
        )

    def RotateMol(  # noqa
        self, molecule: int = None, x: float = 0, y: float = 0, z: float = 0
    ):
        """
        Returns a NEW Cluster Object with a ROTATED molecule (CLOCKWISE)
        around molecule internal center of mass
        """
        # avoiding to rotate a FROZEN molecule
        if molecule in self.GetFreezeMol:
            return self

        if (
            not isinstance(molecule, int)
            or molecule >= self.GetTotalMol
            or molecule < 0
        ):
            raise IndexError(
                f"\nMolecule with {self.GetTotalMol} total molecules "
                f"and index [0-{self.GetTotalMol - 1}]"
                f"\nmolecule index must be less than {self.GetTotalMol}"
                f"\nCheck! You want to remove molecule with index {molecule}"
            )

        molecule_to_rotate: Molecule = self.GetClusterDict[molecule]
        molecule_symbols: list = molecule_to_rotate.GetAtomicSymbols

        # avoid any rotatation attemp for a single atom system
        if len(molecule_symbols) <= 1:
            return self

        molecule_center_of_mass = molecule_to_rotate.GetMolCM
        molecule_principal_axes = molecule_to_rotate.GetMolPrincipalAxes

        # rotate around sphere center
        x, y, z = np.asarray(self.GetSphereCenter) + np.asarray([x, y, z])

        rotation_matrix = Rotation.from_euler(
            "xyz",
            [x, y, z],
            degrees=True,
        ).as_matrix()

        rotatedcoordinates = (
            np.dot(molecule_principal_axes, rotation_matrix)  # noqa
            + molecule_center_of_mass
        )

        rotated_molecule = list()
        for i, atom in enumerate(molecule_symbols):
            rotated_molecule.append(  # noqa
                tuple([atom] + rotatedcoordinates[i].tolist())
            )

        new_cluster = self
        new_cluster.GetClusterDict[molecule] = Molecule(rotated_molecule)

        return self.__class__(
            *new_cluster._cluster_dict.values(),
            freeze_molecule=new_cluster.GetFreezeMol,
            sphere_radius=new_cluster.GetSphereR,
            sphere_center=new_cluster.GetSphereCenter,
        )

    def TranslateMol(  # noqa
        self, molecule: int = None, x: float = 0, y: float = 0, z: float = 0
    ):
        """Returns a NEW Molecule Object with a TRANSLATED fragment"""
        # avoiding to rotate a FROZEN molecule
        if molecule in self.GetFreezeMol:
            return deepcopy(self)

        if (
            not isinstance(molecule, int)
            or molecule >= self.GetTotalMol
            or molecule < 0
        ):
            raise IndexError(
                f"\nMolecule with {self.GetTotalMol} total molecules "
                f"and index [0-{self.GetTotalMol - 1}]"
                f"\nmolecule index must be less than {self.GetTotalMol}"
                f"\nCheck! You want to remove molecule with index {molecule}"
            )

        molecule_to_move: Molecule = self._cluster_dict[molecule]
        molecule_symbols: list = molecule_to_move.GetAtomicSymbols

        molecule_center_of_mass = molecule_to_move.GetMolCM
        molecule_principal_axes = molecule_to_move.GetMolPrincipalAxes

        translatedcoordinates = np.asarray(  # noqa
            molecule_center_of_mass
        ) + np.asarray([x, y, z])

        # checking if the new coordinates are into the boundary conditions
        # if it is out of our sphere, we rescale it to match the sphere radius
        distance: float = np.linalg.norm(
            translatedcoordinates - np.asarray(self.GetSphereCenter)
        )
        if self.GetSphereR and (distance > self.GetSphereR):

            max_distance: float = self.GetSphereR / np.linalg.norm(
                translatedcoordinates - np.asarray(self.GetSphereCenter)
            )

            # rescaling to match radius
            translatedcoordinates = max_distance * translatedcoordinates + (
                1 - max_distance
            ) * np.asarray(self.GetSphereCenter)

        translatedcoordinates = molecule_principal_axes + translatedcoordinates

        translated_molecule = list()
        for i, atom in enumerate(molecule_symbols):
            translated_molecule.append(
                tuple([atom] + translatedcoordinates[i].tolist())
            )

        new_cluster = self
        new_cluster._cluster_dict[molecule] = Molecule(translated_molecule)

        return self.__class__(
            *new_cluster._cluster_dict.values(),
            freeze_molecule=new_cluster.GetFreezeMol,
            sphere_radius=new_cluster.GetSphereR,
            sphere_center=new_cluster.GetSphereCenter
        )

    def CalCentRSphere(self, add_tolerance_radius: float = 1.0):
        """
        Define a spherical outline that contains our cluster

        .. rubric:: Parameters

        add_tolerance_radius : float
            Tolerance with the radius between the mass center to the
            furthest atom

        .. rubric:: Returns

        sphere_center : tuple
            Mass center of the biggest molecule
        sphere_radius : float
            Radius between the sphere center to the furthest atom

        """
        # ----------------------------------------------------------------
        # Verfication
        if not isinstance(add_tolerance_radius, float):
            raise TypeError(
                "\n\nThe tolerance for radius is not a float"
                f"\nplease, check: '{type(add_tolerance_radius)}'\n"
            )
        # ----------------------------------------------------------------
        # Initialize cluster to avoid overlaping, then can calculate of
        # radius
        new_cluster: Cluster = self.InitializeCluster()
        # ---------------------------------------------------------------
        maximum_r_cm = 0.0
        molecule = 0
        max_atoms = 0
        # ---------------------------------------------------------------
        # The biggest molecule
        molecules_number: int = new_cluster.GetTotalMol
        for i in range(molecules_number):
            if new_cluster.GetClusterDict[i].GetNumAtoms() > max_atoms: 
                max_atoms = new_cluster.GetClusterDict[i].GetNumAtoms()
                molecule = i
        # ---------------------------------------------------------------
        # Define sphere center above the cm of the biggest molecule
        center = new_cluster.GetClusterDict[molecule].GetMolCM 
        # ---------------------------------------------------------------
        # Radius between the sphere center to the furthest atom
        for xyz in new_cluster.GetAtomicCoordinates:
            temporal_r = np.linalg.norm(
                np.asarray(new_cluster.GetSphereCenter) - np.asarray(xyz)
            )
            if temporal_r > maximum_r_cm:
                maximum_r_cm = temporal_r
        # ---------------------------------------------------------------
        # Move the biggest molecule to the first position in the cluster
        # object, if is necessary
        if molecule != 0:
            new_geom = dict()
            for i in range(molecules_number):
                if i == 0:
                    new_geom[i] = new_cluster.GetClusterDict[molecule]
                elif i == molecule:
                    new_geom[i] = new_cluster.GetClusterDict[0]
                else:
                    new_geom[i] = new_cluster.GetClusterDict[i]
            # ---------------------------------------------------------------
            # Instantiation of Cluster object with radius and center sphere
            return new_cluster.__class__(
                *new_geom.values(),
                sphere_center=center,
                sphere_radius=maximum_r_cm,
            )
        else:
            return new_cluster.__class__(
                *new_cluster.GetClusterDict.values(),
                sphere_center=center,
                sphere_radius=maximum_r_cm,
            )