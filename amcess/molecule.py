from copy import deepcopy

import attr
from pathlib import Path

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdDetermineBonds
from rdkit.Chem.rdchem import Mol


from amcess.atom import Atom


# This dictionary is a menu to call  RDKit's function according 
# input's format that contains the molecular information
EXT_FILE: dict[str] = {'.mol': Chem.rdmolfiles.MolFromMolFile,
                       '.mol2': Chem.rdmolfiles.MolFromMol2File, 
                       '.xyz': Chem.rdmolfiles.MolFromXYZFile, 
                       '.tpl': Chem.rdmolfiles.MolFromTPLFile, 
                       '.png': Chem.rdmolfiles.MolFromPNGFile, 
                       '.pdb': Chem.rdmolfiles.MolFromPDBFile, 
                       '.mrv': Chem.rdmolfiles.MolFromMrvFile
                       }

@attr.s(frozen=False)
class Molecule(Mol):
    """
    This class inherits attributes of the Mol class from RDKit, 
    for more information:
        *) https://www.rdkit.org/docs/source/rdkit.Chem.rdchem.html#rdkit.Chem.rdchem.Mol
    
    !Class description:
    Create a Molecule that is at least ONE atom.

    The format of the INPUT coordinates must be:

    {"atoms": [(<element> <X> <Y> <Z>), (<element> <X> <Y> <Z>), ...]}

    .. rubric:: Parameters

    atoms : list[tuple(str, float, float, float)]
        Cartesian coordinates of each atom, by default empty list

    charge : int
        total molecular/atomic charge, by default zero (0)

    multiplicity : int
        larger than zero, by defaul one (1)
    """

    #! attributes initial of Molecule class, some they have default value
    _atoms = attr.ib()
    _charge: int = attr.ib(default=0)
    _multiplicity: int = attr.ib(default=1)
    _file: bool = attr.ib(default=False)
    _addHs: bool = attr.ib(default=False)
    _removeHs: bool = attr.ib(default=False)

    # ===============================================================
    # VALIDATORS
    # ===============================================================
    def _init_mol_rdkit(self, mol):
        """Build Molecule object using Mol class from RDKit"""
        super().__init__(mol)
    
    def _check_atom(self, line, atom, atoms):
        """
        Check atomic information
        """
        try:
            #! Check atom informations
            Atom(*atom)
        except (ValueError, TypeError) as err:
            raise TypeError(
                f"\n\n{err}\ncoordinates format must be a dict: "
                "{'atoms':[(str, float, float, float), ...], ...}"
                f"\ncheck atom number {line + 1} --> {atom}\n"
                f"from --> {atoms}\n"
            )

    def _check_atoms_dict(self, attribute, atoms):
        """
        Check that information in dict is Ok
        Example:
        {"atoms": [(<element> <X> <Y> <Z>), ...], "charge": 0, "multiplicty": 1}
        """
        for line, atom in enumerate(atoms["atoms"]):
            self._check_atom(line, atom, atoms)

        total_atoms: int = len(atoms["atoms"]) 
        block_xyz: str = f"""{total_atoms}\n\n"""
        for atom in atoms["atoms"]:
            block_xyz += f"""{atom[0]:<6}"""
            block_xyz += f"""\t{atom[1]:> 15.8f}"""
            block_xyz += f"""\t{atom[2]:> 15.8f}"""
            block_xyz += f"""\t{atom[3]:> 15.8f}\n"""

        mol: Mol = Chem.rdmolfiles.MolFromXYZBlock(block_xyz)
        rdDetermineBonds.DetermineConnectivity(mol)
        self._init_mol_rdkit(mol)

    def _check_atoms_list(self, attribute, atoms):
        """
        Check that information in list is Ok
        Example:
        [(<element> <X> <Y> <Z>), ...]
        """
        for line, atom in enumerate(atoms):
            self._check_atom(line, atom, atoms)

        total_atoms: int = len(atoms) 
        block_xyz: str = f"""{total_atoms}\n\n"""
        for atom in atoms:
            block_xyz += f"""{atom[0]:<6}"""
            block_xyz += f"""\t{atom[1]:> 15.8f}"""
            block_xyz += f"""\t{atom[2]:> 15.8f}"""
            block_xyz += f"""\t{atom[3]:> 15.8f}\n"""

        mol: Mol = Chem.rdmolfiles.MolFromXYZBlock(block_xyz)
        rdDetermineBonds.DetermineConnectivity(mol)
        self._init_mol_rdkit(mol)

    def _mol_charge_multiplicity(self, mol):
        """Save charge and multiplicty in variables of object"""
        self._charge = Chem.rdmolops.GetFormalCharge(mol)
        self._multiplicity = Descriptors.NumRadicalElectrons(mol) + 1 

    def _check_atoms_smiles(self, attribute, atoms):
        """
        Check that the smiles is Ok
        Example:
        'CC'
        """
        try:
            Chem.MolFromSmiles(atoms, sanitize=False)
        except (ValueError, TypeError) as err:
            raise TypeError(f"\n\n{err}\n atoms must be a smiles: "
                            " 'CCO' ")
        mol: Mol = Chem.MolFromSmiles(atoms)
        mol = Chem.AddHs(mol, explicitOnly=self._addHs)
        # NOTE: Explanation of EmbedMolecule process
        #       https://www.rdkit.org/docs/GettingStartedInPython.html#working-with-3d-molecules
        AllChem.EmbedMolecule(mol)
        self._mol_charge_multiplicity(mol)
        self._init_mol_rdkit(mol)

    def _check_atoms_file(self, attribute, file):
        """
        Check that molecular information is Ok and if the file exists
        Inputs formats allowed: .xyz, .mol, .mol2, .tpl, png, .pdb, .mrv
        """
        if not Path(file).exists():
            raise ValueError(
                f"The file {file} doesn't exist"
            )
        if Path(file).suffix.lower() not in EXT_FILE.keys():
            raise TypeError(
                f"File with extension {Path(file).suffix} can't be reading"
            )
        if Path(file).suffix.lower() in ['.xyz', '.png', '.tpl']:
            mol: Mol = EXT_FILE[Path(file).suffix.lower()](file)
            if Path(file).suffix.lower() in ['.xyz']:
                rdDetermineBonds.DetermineConnectivity(mol)    
        else:
            mol = EXT_FILE[Path(file).suffix.lower()](file, removeHs=self._removeHs)
        
        if self._addHs:
            if Path(file).suffix.lower() == '.mol2':
                raise TypeError('RDKit have problem to add H to mol from'
                                'mol2 file')                
            mol = Chem.AddHs(mol, addCoords=True)
        self._mol_charge_multiplicity(mol)
        self._init_mol_rdkit(mol)

    @_atoms.validator
    def _cehck_valid_atoms(self, attribute, atoms):
        """check if the atoms are valid"""
        if isinstance(atoms, (Molecule, Chem.rdchem.Mol)):
            self._init_mol_rdkit(atoms)
        elif self._file:
            self._check_atoms_file(attribute, atoms)
        elif isinstance(atoms, list):
            self._check_atoms_list(attribute, atoms)
        elif isinstance(atoms, dict):
            self._check_atoms_dict(attribute, atoms)
        elif isinstance(atoms, str):
            self._check_atoms_smiles(attribute, atoms)
        else:
            raise TypeError(
                    "\ncoordinates format must be a list of tuple"
                    " or str (Smiles):\n[(str, float, float, float), ...]"
                    "\n'CCO' "
                )

    @_addHs.validator
    def _cehck_valid_addHs(self, attribute, addHs):
        if not isinstance(addHs, bool):
            raise ValueError(
                "\n\naddHs must be an bool "  # noqa
                f"\nyou get --> 'addHs = {addHs}'\n"
            )

    @_charge.validator
    def _check_valid_charge(self, attribute, charge):
        """check if the charge is valid"""
        if not isinstance(charge, int):
            raise ValueError(
                "\n\ncharge must be an integer "  # noqa
                f"\nyou get --> 'charge = {charge}'\n"
            )

    @_multiplicity.validator
    def _check_valid_multiplicity(self, attribute, multiplicity):
        """check if the multiplicity is valid"""
        if not isinstance(multiplicity, int) or multiplicity < 1:
            raise ValueError(
                "\n\nmultiplicity must be an integer larger than zero (0)"
                f"\nyou get --> 'multiplicity = {multiplicity}'\n"
            )

    # ===============================================================
    # MAGIC METHODS
    # ===============================================================
    def __add__(self, other) -> object:
        """Magic method '__add__' to add two molecules, return a new one"""
        return Molecule(self.GetMolList() + other.GetMolList(),
                        self.GetMolCharge() + other.GetMolCharge(),
                        (self.GetMolMultiplicity() + other.GetMolMultiplicity())
                        - 1)

    def __mul__(self, value: int):
        """Magic method '__mul__' to multiply a molecule by a number"""
        return value * self

    def __rmul__(self, value: int):
        """
        Replicate a molecule.
        summing or multiplying Molecule classes produce a Cluster class

        Parameters
        ----------
        value : int
            quantity to replicate Molecue

        Return
        ------
        Cluster : object
        """
        if not isinstance(value, int) or value < 1:
            raise ValueError(
                "\nMultiplier must be and integer larger than zero"
                f"\ncheck --> '{value}'"
            )

        tcoordinates: list = [at for i in range(value)
                              for at in self.GetMolList()]
        tcharge: int = 0
        tmultiplicity: int = 0
        for i in range(value):
            tcharge += self.GetMolCharge()
            tmultiplicity += self.GetMolMultiplicity()
        tmultiplicity -= value - 1

        return Molecule(tcoordinates, tcharge, tmultiplicity)
        
    def __str__(self):
        """Magic method '__str__' to print the Molecule in XYZ format"""
        return self.GetBlockXYZ()

    # ===============================================================
    # PROPERTIES
    # ===============================================================
    #! Getter
    def GetAtomicSymbols(self) -> list[str]:
        """Return the list of atoms"""
        return [a.GetSymbol() for a in self.GetAtoms()]

    def GetAtomicNumbers(self) -> list[int]:
        """Return the list of atoms"""
        return [a.GetAtomicNum() for a in self.GetAtoms()]

    def GetBlockXYZ(self) -> str:
        """Printing Molecule coordinates using XYZ format"""
        write_coordinates: str = ""
        comment: str = f"charge: {self.GetMolCharge()} multiplicity:{self.GetMolMultiplicity()}"
        write_coordinates += f"{len(self.GetAtomicSymbols())}\n{comment}\n"
        for atom, xyz in zip(self.GetAtoms(), self.GetConformer().GetPositions()):
            write_coordinates += f"""{atom.GetSymbol():<6}"""
            write_coordinates += f"""\t{xyz[0]:> 15.8f}"""
            write_coordinates += f"""\t{xyz[1]:> 15.8f}"""
            write_coordinates += f"""\t{xyz[2]:> 15.8f}\n"""

        return write_coordinates

    def GetAtomicMasses(self) -> list[float]:
        """Atomic mass of the molecule"""
        return [atom.GetMass() for atom in self.GetAtoms()]

    def GetMolCharge(self) -> int:
        """Total molecular/atomic charge"""
        return self._charge

    def GetAtomicCoordinates(self) -> list[float]:
        """Return the list of coordinates"""
        return self.GetConformer().GetPositions()

    def GetMolMultiplicity(self) -> int:
        """Return the multiplicity"""
        return self._multiplicity

    def GetMolList(self) -> list[str, float]:
        """Return the dict atoms"""
        list_atoms: list(tuple(str, float)) = [tuple([a.GetSymbol()] + list(xyz))
                for a, xyz in zip(self.GetAtoms(), self.GetConformer().GetPositions())]
        return list_atoms

    def GetMolDict(self) -> dict:
        """Return the dict atoms"""
        return {"atoms": self.GetMolList(),
                "charge": self.GetMolCharge(),
                "multiplicity": self.GetMolMultiplicity()}

    def GetMolMass(self) -> float:
        """Return the total mass of the molecule"""
        return sum(self.GetAtomicMasses())

    def GetMolCM(self) -> tuple[float]:
        """Center of mass for a N-body problem. `Jacobi coordinates`_

        .. rubric:: Notes

        total mass for dummy atoms (not in the Periodic Table) is equal
        to ONE (1)

        .. rubric:: Returns

        tuple : (float, float, float)
            List of N 3D tuples, where N is equal to the number of atoms

        .. _Jacobi coordinates:
            https://en.wikipedia.org/wiki/Jacobicoordinates
        """
        # total_mass = 1 if not self.total_mass else self.total_mass
        return tuple(
            np.dot(
                np.asarray(self.GetAtomicMasses()),
                np.asarray(self.GetAtomicCoordinates()),
            )
            / self.GetMolMass()
        )

    def GetMolPrincipalAxes(self) -> list[float]:
        """Principal axes for according to Jacobi coordinates"""
        return [
            tuple(c)
            for c in (  # noqa
                np.asarray(self.GetAtomicCoordinates())
                - np.asarray(self.GetMolCM())
            )
        ]

    #! Setter
    def SetMolCharge(self, new_charge) -> int:
        """Set the total molecular/atomic charge"""
        if not isinstance(new_charge, int):
            raise ValueError(
                "\n\ncharge must be an integer "
                f"\nyou get --> 'charge = {new_charge}'\n"
            )
        self._charge = new_charge

    def SetMolMultiplicity(self, new_multiplicity) -> int:
        """Set the multiplicity"""
        if not isinstance(new_multiplicity, int) or new_multiplicity < 1:
            raise ValueError(
                "\n\nmultiplicity must be an integer larger than zero (0)"
                f"\nyou get --> 'multiplicity = {new_multiplicity}'\n"
            )
        self._multiplicity = new_multiplicity


    #! Para que?
    # @property
    # def GetNumberingAtoms(self) -> str:
    #     """show atom number line by line

    #     .. rubric:: Returns

    #     str
    #         atom number line by line
    #     """
    #     numbered_atoms = list()
    #     for a, xyz in zip(self.GetAtoms(), self.GetCoordinates):
    #         line = list(self.atoms[i])
    #         line[0] = f"\r  atom #{i} --> {line[0]:<6}"
    #         line[1] = f"{line[1]:> 15.8f}"
    #         line[2] = f"{line[2]:> 15.8f}"
    #         line[3] = f"{line[3]:> 15.8f}"
    #         numbered_atoms.append("".join(line))
    #     return "\n".join(numbered_atoms)

    def AddAtoms(self, new_atoms: list, attribute: None = None) -> object:
        """adding extra atoms can NOT be MOVED or ROTATED

        .. rubric:: Parameters

        other : list
            cartesian coordinates; like [(<element>, <X>, <Y>, <Y>), ...]

        .. rubric:: Returns

        Molecule : object
            a new Molecule

        .. rubric:: Raises

        TypeError
            for anything else
        """
        if not isinstance(new_atoms, list):
            raise TypeError(
                f"\n\ncoordinates format must be a list of tuple: "
                "[(str, float, float, float), ...]"
                f"check --> \n{new_atoms}\n"
            )

        atoms: list = self.GetMolList() + new_atoms

        return self._cehck_valid_atoms(attribute, atoms)

    # def add_molecule(self, other) -> object:
    #     """adding molecule return a new Cluster object"""
    #     if not isinstance(other, dict):    #Molecule):
    #         raise TypeError(
    #             "\nOnly type 'Molecule', list or dict could be added"
    #             f"\nyou have a type: '{type(other)}', check: \n{other}"
    #         )
    #     return Cluster(self, other)

    def GetAtomWithIndix(self, atom: int) -> list:
        """
        Getting catesian coordinate for an atom

        .. rubric:: Parameters

        atom : int
            atom index

        .. rubric:: Returns

        tuple
            ("element", "X", "Y", "Z")

        .. rubric:: Raises

        IndexError
        """
        if not isinstance(atom, int) or atom >= self.GetNumAtoms():
            raise IndexError(
                f"\nMolecule with {self.GetNumAtoms()} total atoms "
                f"and index [0-{self.GetNumAtoms() - 1}]"
                f"\n atom index must be less than {self.GetNumAtoms()}"
                f"\nCheck! You want to get atom with index {atom}"
            )
        return self.GetMolList()[atom]

    def RemoveAtom(self, atom: int, attribute: None = None) -> object:
        """remove one atom"""
        if not isinstance(atom, int) or atom >= self.GetNumAtoms():
            raise IndexError(
                f"\nMolecule with {self.GetNumAtoms()} total atoms "
                f"and index [0-{self.GetNumAtoms() - 1}]"
                f"\n atom index must be less than {self.GetNumAtoms()}"
                f"\nCheck! You want to remove atom with index '{atom}'"
            )

        atoms: list = self.GetMolList()

        del atoms[atom]

        return self._cehck_valid_atoms(attribute, atoms)
