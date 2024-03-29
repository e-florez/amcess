# from copy import deepcopy
# from rdkit.Chem import AllChem  # type: ignore
import attr
from pathlib import Path

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors


from amcess.atom import Atom
from rdkit.Chem import Descriptors  # type: ignore


# NOTE: Format allow in rdkit
EXT_FILE: dict[str] = {'.mol': Chem.rdmolfiles.MolFromMolFile,
                       '.mol2': Chem.rdmolfiles.MolFromMol2File, 
                       '.xyz': Chem.rdmolfiles.MolFromXYZFile, 
                       '.tpl': Chem.rdmolfiles.MolFromTPLFile, 
                       '.png': Chem.rdmolfiles.MolFromPNGFile, 
                       '.pdb': Chem.rdmolfiles.MolFromPDBFile, 
                       '.mrv': Chem.rdmolfiles.MolFromMrvFile
                       }

@attr.s(frozen=False)
class Molecule:
    """
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

    _atoms = attr.ib()
    _charge: int = attr.ib(default=0)
    _multiplicity: int = attr.ib(default=1)
    _file: bool = attr.ib(default=False)
    _addHs: bool = attr.ib(default=False)
    _removeHs: bool = attr.ib(default=False)

    # ===============================================================
    # VALIDATORS
    # ===============================================================
    def _check_atoms_list(self, attribute, atoms):
        for line, atom in enumerate(atoms):
            try:
                Atom(*atom)
            except (ValueError, TypeError) as err:
                raise TypeError(
                    f"\n\n{err}\ncoordinates format must be a list of tuple: "
                    "[(str, float, float, float), ...]"
                    f"\ncheck atom number {line + 1} --> {atom}\n"
                    f"from --> {atoms}\n"
                )
    def _molecule_list_building(self, mol):
        self._atoms = [tuple([a.GetSymbol()] + list(xyz))
                for a, xyz in zip(mol.GetAtoms(), mol.GetConformer().GetPositions())]
        self._charge = Chem.rdmolops.GetFormalCharge(mol)
        self._multiplicity = Descriptors.NumRadicalElectrons(mol) + 1 

    def _check_atoms_smiles(self, attribute, atoms):
        try:
            Chem.MolFromSmiles(atoms, sanitize=False)
        except (ValueError, TypeError) as err:
            raise TypeError(
                f"\n\n{err}\n atoms must be a smiles: "
                " 'CCO' "
            )
        mol = Chem.MolFromSmiles(atoms)
        mol = Chem.AddHs(mol, explicitOnly=self._addHs)
        # NOTE: Explanation of EmbedMolecule process
        #       https://www.rdkit.org/docs/GettingStartedInPython.html#working-with-3d-molecules
        AllChem.EmbedMolecule(mol)
        self._molecule_list_building(mol)

    def _check_atoms_file(self, attribute, file):
        if not Path(file).exists():
            raise ValueError(
                f"The file {file} doesn't exist"
            )
        if Path(file).suffix.lower() not in EXT_FILE.keys():
            raise TypeError(
                f"File with extension {Path(file).suffix} can't be reading"
            )
        if Path(file).suffix.lower() in ['.xyz', '.png', '.tpl']:
            mol = EXT_FILE[Path(file).suffix.lower()](file)    
        else:
            mol = EXT_FILE[Path(file).suffix.lower()](file, removeHs=self._removeHs)
        if self._addHs:
            if Path(file).suffix.lower() == '.mol2':
                raise TypeError('RDKit have problem to add H to mol from'
                                'mol2 file')                
            mol = Chem.AddHs(mol, addCoords=True)
        self._molecule_list_building(mol)

    @_atoms.validator
    def _cehck_valid_atoms(self, attribute, atoms):
        """check if the atoms are valid"""
        if self._file:
            self._check_atoms_file(attribute, atoms)
        elif isinstance(atoms, list):
            self._check_atoms_list(attribute, atoms)
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
    # CONSTRUCTORS
    # ===============================================================
    @classmethod
    def from_dict(cls, atoms_dict):
        "Dictionary type: {'atoms': [(<element> <X> <Y> <Z>), ...]}"
        if "atoms" not in atoms_dict:
            # FIXME: KeyError does not support \n
            raise TypeError(
                "\n\nThe key 'atoms' is casesensitive"
                "\n{'atoms': [(<element> <X> <Y> <Z>), ...]}"
                f"\nyou get {atoms_dict}\n"
            )
        atoms = atoms_dict.get("atoms")
        charge = atoms_dict.get("charge", 0)
        multiplicity = atoms_dict.get("multiplicity", 1)
        return cls(atoms, charge, multiplicity)

    # ===============================================================
    # MAGIC METHODS
    # ===============================================================
    def __add__(self, other) -> object:
        """Magic method '__add__' to add two molecules, return a new one"""
        return Molecule(
            self.atoms + other.atoms,
            self.charge + self.charge,
            (self.multiplicity + other.multiplicity) - 1,
        )

    # def __mul__(self, value: int):
    #    """Magic method '__mul__' to multiply a molecule by a number"""
    #    return value * self

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

        tcoordinates: list = [at for i in range(value) for at in self.atoms]
        tcharge: int = 0
        tmultiplicity: int = 0
        for i in range(value):
            tcharge += self.charge
            tmultiplicity += self.multiplicity
        tmultiplicity -= value - 1

        return Molecule(tcoordinates, tcharge, tmultiplicity)

    def __str__(self):
        """Magic method '__str__' to print the Molecule in XYZ format"""
        return self.xyz

    # ===============================================================
    # PROPERTIES
    # ===============================================================
    @property
    def atoms(self) -> list:
        """Return the list of atoms"""
        return self._atoms

    @atoms.setter
    def atoms(self, *args, **kwargs) -> None:
        """Set the list of atoms"""
        raise AttributeError(
            "\n\nyou cannot reset 'atoms'. Consider create a new instance \n"
        )

    @property
    def write_atoms(self) -> str:
        """Printing Molecule coordinates using XYZ format"""
        write_coordinates = ""
        for atom in self.atoms:
            write_coordinates += f"""{atom[0]:<6}"""
            write_coordinates += f"""\t{atom[1]:> 15.8f}"""
            write_coordinates += f"""\t{atom[2]:> 15.8f}"""
            write_coordinates += f"""\t{atom[3]:> 15.8f}\n"""

        return write_coordinates

    @property
    def atomic_masses(self) -> list:
        """Atomic mass of the molecule"""
        return [Atom(*atom).atomic_mass for atom in self.atoms]

    @property
    def charge(self) -> int:
        """Total molecular/atomic charge"""
        return self._charge

    @charge.setter
    def charge(self, new_charge) -> int:
        """Set the total molecular/atomic charge"""
        if not isinstance(new_charge, int):
            raise ValueError(
                "\n\ncharge must be an integer "
                f"\nyou get --> 'charge = {new_charge}'\n"
            )
        self._charge = new_charge

    @property
    def coordinates(self) -> list:
        """Return the list of coordinates"""
        return [c[1:] for c in self.atoms]

    @property
    def elements(self) -> list:
        """Show a list of unique symbols

        .. rubric:: Returns

        list
            list of unique symbols
        """
        return list(set(self.symbols))

    @property
    def multiplicity(self) -> int:
        """Return the multiplicity"""
        return self._multiplicity

    @multiplicity.setter
    def multiplicity(self, new_multiplicity) -> int:
        """Set the multiplicity"""
        if not isinstance(new_multiplicity, int) or new_multiplicity < 1:
            raise ValueError(
                "\n\nmultiplicity must be an integer larger than zero (0)"
                f"\nyou get --> 'multiplicity = {new_multiplicity}'\n"
            )
        self._multiplicity = new_multiplicity

    @property
    def molecule(self) -> dict:
        """Return the dict atoms"""
        return {
            "atoms": self.atoms,
            "charge": self.charge,
            "multiplicity": self.multiplicity,
        }

    @property
    def numbering_atoms(self) -> str:
        """show atom number line by line

        .. rubric:: Returns

        str
            atom number line by line
        """
        numbered_atoms = list()
        for i in range(self.total_atoms):
            line = list(self.atoms[i])
            line[0] = f"\r  atom #{i} --> {line[0]:<6}"
            line[1] = f"{line[1]:> 15.8f}"
            line[2] = f"{line[2]:> 15.8f}"
            line[3] = f"{line[3]:> 15.8f}"
            numbered_atoms.append("".join(line))

        return "\n".join(numbered_atoms)

    @property
    def symbols(self) -> list:
        """Return the list of symbols"""
        return [str(s[0]).title() for s in self.atoms]

    @property
    def total_atoms(self) -> int:
        """Return the total number of atoms"""
        return len(self.atoms)

    @property
    def total_mass(self) -> float:
        """Return the total mass of the molecule"""
        return sum(self.atomic_masses)

    @property
    def center_of_mass(self) -> tuple:
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

        total_mass = 1 if not self.total_mass else self.total_mass

        return tuple(
            np.dot(
                np.asarray(self.atomic_masses),
                np.asarray(self.coordinates),
            )
            / total_mass
        )

    @property
    def principal_axes(self) -> list:
        """Principal axes for according to Jacobi coordinates"""
        return [
            tuple(c)
            for c in (  # noqa
                np.asarray(self.coordinates) - np.asarray(self.center_of_mass)
            )
        ]

    @property
    def xyz(self) -> str:
        """Printing Molecule coordinates using XYZ format"""
        comments = (
            f"-- charge={self.charge:<-g} and "
            f"multiplicity={self.multiplicity:<g} --"
        )
        write_xyz = f"""\t{self.total_atoms}\n{comments:<s}\n"""
        for atom in self.atoms:
            write_xyz += f"""{atom[0]:<6}"""
            write_xyz += f"""\t{atom[1]:> 15.8f}"""
            write_xyz += f"""\t{atom[2]:> 15.8f}"""
            write_xyz += f"""\t{atom[3]:> 15.8f}\n"""

        return write_xyz

    def add_atoms(self, new_atoms: list) -> object:
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

        total_atoms: list = self.atoms + new_atoms
        return self.__class__(total_atoms)

    # def add_molecule(self, other) -> object:
    #     """adding molecule return a new Cluster object"""
    #     if not isinstance(other, dict):    #Molecule):
    #         raise TypeError(
    #             "\nOnly type 'Molecule', list or dict could be added"
    #             f"\nyou have a type: '{type(other)}', check: \n{other}"
    #         )
    #     return Cluster(self, other)

    def get_atom(self, atom: int) -> list:
        """
        Getting catesian coordinate for an atom

        .. rubric:: Parameters

        atom : int
            atom index

        .. rubric:: Returns

        list
            ["element", "X", "Y", "Z"]

        .. rubric:: Raises

        IndexError
        """
        if not isinstance(atom, int) or atom >= self.total_atoms:
            raise IndexError(
                f"\nMolecule with {self.total_atoms} total atoms "
                f"and index [0-{self.total_atoms - 1}]"
                f"\n atom index must be less than {self.total_atoms}"
                f"\nCheck! You want to get atom with index {atom}"
            )
        return self.atoms[atom]

    def remove_atom(self, atom: int) -> object:
        """remove one atom"""
        if not isinstance(atom, int) or atom >= self.total_atoms:
            raise IndexError(
                f"\nMolecule with {self.total_atoms} total atoms "
                f"and index [0-{self.total_atoms - 1}]"
                f"\n atom index must be less than {self.total_atoms}"
                f"\nCheck! You want to remove atom with index '{atom}'"
            )

        new_atoms: list = list(self.atoms)

        del new_atoms[atom]

        return self.__class__(new_atoms)
