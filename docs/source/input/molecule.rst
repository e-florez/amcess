Molecule class
--------------

The Molecule class inherits everything from RDKit's Mol class because it has the properties that are needed within the AMCESS program. 
Additionally, the Molecule class inherits everything from the Atom class, when the objects are created.

RDKit's Mol class documentation is located in the following website:
.. centered:: `<https://www.rdkit.org/docs/source/rdkit.Chem.rdchem.html\#rdkit.Chem.rdchem.Mol>`_


RDKit's ol class is inside of the Chem package (`rdkit.Chem <https://www.rdkit.org/docs/source/rdkit.html>`_), 
because it belongs to the rdchem module (`rdkit.Chem.rdchem <https://www.rdkit.org/docs/source/rdkit.Chem.html>`_).

The next RDKit's classes can be utility:

* `class rdkit.Chem.rdchem.RWMol((object)self, (Mol)m) <https://www.rdkit.org/docs/source/rdkit.Chem.rdchem.html\#rdkit.Chem.rdchem.RWMol>`_: This class is a more-performant version of the EditableMolecule class in that it is a live molecule and shares the interface from the Mol class. All changes are performed without the need to create a copy of the molecule using GetMol() (this is still available, however). n.b. Eventually this class may become a direct replacement for EditableMol.
* `class rdkit.Chem.rdchem.MolBundle((object)self) <https://www.rdkit.org/docs/source/rdkit.Chem.rdchem.html\#rdkit.Chem.rdchem.MolBundle>`_: A class for storing groups of related molecules.
* `class rdkit.Chem.rdchem.PeriodicTable <https://www.rdkit.org/docs/source/rdkit.Chem.rdchem.html\#rdkit.Chem.rdchem.PeriodicTable>`_: A class which stores information from the Periodic Table (GetAtomicWeight, GetAtomicNumber, GetElementSymbol, GetElementName, GetRvdw (van der Waals radius), GetRCovalent (covalent radius), GetDefaultValence, GetValenceList, GetNOuterElecs (number of valence electrons), GetMostCommonIsotope, GetMostCommonIsotopeMass, GetRb0, GetAbundanceForIsotope, GetMassForIsotope).

Constructor
^^^^^^^^^^^

Due to the diversity of ways to obtain the structural information of a molecule, the attr package is used to define the Molecule class constructor. The :ref:`python-m1` shows the arguments of the Molecule class, accompanied by their respective type and default values.

.. _python-m1:
.. code-block:: python
    :caption: Listing 6 : Molecule class' Arguments.
    :linenos:

    _atoms = attr.ib()
    _charge: int = attr.ib(default=0)
    _multiplicity: int = attr.ib(default=1)
    _file: bool = attr.ib(default=False)
    _addHs: bool = attr.ib(default=False)
    _removeHs: bool = attr.ib(default=False)

The :ref:`python-m2` shows the validator functions of the Molecule class. Keep in mind that the line number is not necessarily equal to the respective one in the .py.

.. _python-m2:
.. code-block:: python
    :caption: Listing 7 : Validators functions of Molecule class
    :linenos:
    :emphasize-lines: 6

    @_atoms.validator
    def _cehck_valid_atoms(self, attribute, atoms):
        """check if the atoms are valid"""

        if isinstance(atoms, (Molecule, Chem.rdchem.Mol)):
            self._init_mol_rdkit(atoms)
            self._check_atom(atoms.GetMolList())
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
        """check if _addHs is a bool"""
        if not isinstance(addHs, bool):
            raise ValueError(
                "\n\naddHs must be an bool "  # noqa
                f"\nyou get --> 'addHs = {addHs}'\n"
            )

    @_file.validator
    def _cehck_valid_file(self, attribute, file):
        """check if _file is a bool"""
        if not isinstance(file, bool):
            raise ValueError(
                "\n\nfile must be an bool "  # noqa
                f"\nyou get --> 'file = {file}'\n"
            )

    @_charge.validator
    def _check_valid_charge(self, attribute, charge):
        """check if the charge is valid"""
        if not isinstance(charge, int):
            raise ValueError(
                "\n\ncharge must be an integer "  # noqa
                f"\nyou get --> 'charge = {charge}'\n"
            )
        if charge != 0:
            self._charge = charge

    @_multiplicity.validator
    def _check_valid_multiplicity(self, attribute, multiplicity):
        """check if the multiplicity is valid"""
        if not isinstance(multiplicity, int) or multiplicity < 1:
            raise ValueError(
                "\n\nmultiplicity must be an integer larger than zero (0)"
                f"\nyou get --> 'multiplicity = {multiplicity}'\n"
            )
        if multiplicity != 1:
            self._multiplicity = multiplicity

The **_init_mol_rdkit** function is responsible of the inheritance from the RDKit's Mol class. The mol variable must be any of the file formats shown in lines 6-8 (`mol <https://www.rdkit.org/docs/source/rdkit.Chem.rdmolfiles.html\#rdkit.Chem.rdmolfiles.MolFromMolFile>`_, `mol2 <https://www.rdkit.org/docs/source/rdkit.Chem.rdmolfiles.html\#rdkit.Chem.rdmolfiles.MolFromMol2File>`_, `xyz <https://www.rdkit.org/docs/source/rdkit.Chem.rdmolfiles.html\#rdkit.Chem.rdmolfiles.MolFromXYZFile>`_, `tpl <https://www.rdkit.org/docs/source/rdkit.Chem.rdmolfiles.html\#rdkit.Chem.rdmolfiles.MolFromTPLFile>`_, `png <https://www.rdkit.org/docs/source/rdkit.Chem.rdmolfiles.html\#rdkit.Chem.rdmolfiles.MolFromPNGFile>`_, `pdb <https://www.rdkit.org/docs/source/rdkit.Chem.rdmolfiles.html\#rdkit.Chem.rdmolfiles.MolFromPDBFile>`_, `mrv <https://www.rdkit.org/docs/source/rdkit.Chem.rdmolfiles.html\#rdkit.Chem.rdmolfiles.MolFromMrvFile>`_) or `SMILE notation <https://www.rdkit.org/docs/source/rdkit.Chem.rdmolfiles.html\#rdkit.Chem.rdmolfiles.MolFromSmiles>`_. When the argument received by the Molecule class is a Molecule or Mol object, the **_init_mol_rdkit** is called in line 6 of :ref:`python-m2`.

.. _python-m3:
.. code-block:: python
    :caption: Listing 8 : Inheritance from RDKit's mol class to AMCESS's Molecule class. Additionally, the files formats allow is shown.
    :linenos:
    :emphasize-lines: 2-8

    EXT_FILE: dict[str] = {
    ".mol": Chem.rdmolfiles.MolFromMolFile,
    ".mol2": Chem.rdmolfiles.MolFromMol2File,
    ".xyz": Chem.rdmolfiles.MolFromXYZFile,
    ".tpl": Chem.rdmolfiles.MolFromTPLFile,
    ".png": Chem.rdmolfiles.MolFromPNGFile,
    ".pdb": Chem.rdmolfiles.MolFromPDBFile,
    ".mrv": Chem.rdmolfiles.MolFromMrvFile,
    }
    
    def _init_mol_rdkit(self, mol):
        """Build Molecule object using Mol class from RDKit"""
        super().__init__(mol)

The following functions check the format of molecular information input. The main function is **_check_atom**, because it checks the input format of the AMCESS's Atom class (:ref:`python-a1`). This is necessary because the Molecule object's fundamental property is a list of Atom objects.

The **_check_atom** function is called from the **_check_atoms_dict** (line=27), **_check_atoms_list** (line=53), **_check_atoms_smiles** (line=89), and **_check_atoms_file** (line=124) functions. These functions check the dictionary, list, SMILES, and files formats (:ref:`python-m3`), respectively.

The **_check_atoms_smiles** (line=82, 86) and **_check_atoms_file** (líneas=127) functions admit that a molecule can be defined with only the heavy atoms, i. e., without the Hs explicits, because these functions add Hs by mean of **ff**.

The atomic coordinates in the **_check_atoms_smiles** (lines=84-87) and **_check_atoms_file** (lines=119-122) functions are obtained with **Mol.GetComformer().GetPositions()**. However, this procedure is not used in the Getters or Setters methods, for design reasons, since as mentioned, Molecule class objects are lists of AMCESS's Atom objects.

.. _python-m4:
.. code-block:: python
    :caption: Listing 9 : Functions that check the formats of input molecular information for Molecule class.
    :linenos:
    :emphasize-lines: 27, 53, 84-87, 89, 119-122, 124 

    def _check_atom(self, atoms) -> None:
        """
        Check atomic information
        """

        self._molecule: list = []
        for line, atom in enumerate(atoms):
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

            self._molecule.append(Atom(*atom))

    def _check_atoms_dict(self, attribute, atoms):
        """
        Check that information in dict is Ok
        Example:
        {"atoms": [(<element> <X> <Y> <Z>), ...], "charge": 0, "multiplicty": 1}
        """
        self._check_atom(atoms["atoms"])

        if "charge" in atoms.keys():
            self._check_valid_charge(attribute, atoms["charge"])

        if "multiplicity" in atoms.keys():
            self._check_valid_multiplicity(attribute, atoms["multiplicity"])

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
        self._check_atom(atoms)

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

    def _check_atoms_smiles(self, attribute, atoms):
        """
        Check that the smiles is Ok
        Example:
        'CC'
        """
        try:
            Chem.MolFromSmiles(atoms, sanitize=False)
        except (ValueError, TypeError) as err:
            raise TypeError(f"\n\n{err}\n atoms must be a smiles: " " 'CCO' ")

        mol: Mol = Chem.MolFromSmiles(atoms)
        mol = Chem.AddHs(mol, explicitOnly=self._addHs)
        # NOTE: Explanation of EmbedMolecule process
        #       https://www.rdkit.org/docs/GettingStartedInPython.html#working-with-3d-molecules
        AllChem.EmbedMolecule(mol)

        atoms = [
            tuple(a.GetSymbol()) + tuple(r)
            for a, r in zip(mol.GetAtoms(), mol.GetConformer().GetPositions())
        ]

        self._check_atom(atoms)

        self._mol_charge_multiplicity(mol)
        self._init_mol_rdkit(mol)

    def _check_atoms_file(self, attribute, file):
        """
        Check that molecular information is Ok and if the file exists
        Inputs formats allowed: .xyz, .mol, .mol2, .tpl, png, .pdb, .mrv
        """
        if not Path(file).exists():
            raise ValueError(f"The file {file} doesn't exist")
        if Path(file).suffix.lower() not in EXT_FILE.keys():
            raise TypeError(
                f"File with extension {Path(file).suffix}" + "can't be reading"
            )
        if Path(file).suffix.lower() in [".xyz", ".png", ".tpl"]:
            mol: Mol = EXT_FILE[Path(file).suffix.lower()](file)
            if Path(file).suffix.lower() in [".xyz"]:
                rdDetermineBonds.DetermineConnectivity(mol)
        else:
            mol = EXT_FILE[Path(file).suffix.lower()](file, removeHs=self._removeHs)

        if self._addHs:
            if Path(file).suffix.lower() == ".mol2":
                raise TypeError(
                    "RDKit have problem to add H to mol from " + "mol2 file"
                )
            mol = Chem.AddHs(mol, addCoords=True)

        atoms = [
            tuple(a.GetSymbol()) + tuple(r)
            for a, r in zip(mol.GetAtoms(), mol.GetConformer().GetPositions())
        ]

        self._check_atom(atoms)

        self._mol_charge_multiplicity(mol)
        self._init_mol_rdkit(mol)

Finally, the **_mol_charge_multiplicity** function checks the _charge and _multiplicity properties. This function is used when the molecule input is a file or SMILES.

.. _python-m5:
.. code-block:: python
    :caption: Listing 10 : Function that reviews the _charge and _multiplicity properties.
    :linenos:

    def _mol_charge_multiplicity(self, mol):
        """Save charge and multiplicty in variables of object"""
        self._charge = Chem.rdmolops.GetFormalCharge(mol)
        self._multiplicity = Descriptors.NumRadicalElectrons(mol) + 1

Properties
^^^^^^^^^^

RDKit's Mol class is composed for he next properties:

.. tabularcolumns:: p{0.12\linewidth}p{0.196\linewidth}p{0.30\linewidth}p{0.30\linewidth}
.. table:: Table 3 : Properties from RDKit's Mol class.
   :name: tab_m1
   :widths: 30, 40, 20, 10
   :class: longtable
   :align: center
   :width: 66%

   +------------------------+----------------------------+------------------+-------------------------+
   |**Propety**             |  **Argument**              | **Return**       |  **Description**        |
   +========================+============================+==================+=========================+
   |.. centered::                             **Getters**                                             |
   +------------------------+----------------------------+------------------+-------------------------+
   | GetAromaticAtoms       |                            | _ROQAtomSeq      | Aromatic atoms          |
   +------------------------+----------------------------+------------------+-------------------------+
   | GetAtoms               |                            | Iterador         | Atoms objects           |
   +------------------------+----------------------------+------------------+-------------------------+
   | GetAtomsMatchingQuery  | qa:QueryAtom               | _ROQAtomSeq      | Atoms belong to query   |
   +------------------------+----------------------------+------------------+-------------------------+
   | GetAtomWithIdx         | idx: int                   | Atom             | Atom object             |
   +------------------------+----------------------------+------------------+-------------------------+
   | GetBondBetweenAtoms    |  idx1:int, idx2:int        | Bond             | Bonds between atoms     |
   +------------------------+----------------------------+------------------+-------------------------+
   | GetBondWithIdx         | idx: int                   |                  | Bond selection          |
   +------------------------+----------------------------+------------------+-------------------------+
   | GetBonds               |                            | Iterador         | Bonds                   |
   +------------------------+----------------------------+------------------+-------------------------+
   | GetBoolProp            | key:str                    | bool             | Property with bool value|
   +------------------------+----------------------------+------------------+-------------------------+
   | GetConformer           |                            | Conformer        | Conformer selection     |
   +------------------------+----------------------------+------------------+-------------------------+
   | GetConformers          |                            | _ROConformerSeq  | Conformers              |
   +------------------------+----------------------------+------------------+-------------------------+
   | GetDoubleProp          | key:str                    | double           |Property with float value|
   +------------------------+----------------------------+------------------+-------------------------+
   | GetIntProp             | key:str                    | int              |Property with int value  |
   +------------------------+----------------------------+------------------+-------------------------+
   | GetNumAtoms            | onlyExplicit:bool(T)       | int              | Atoms number            |
   +------------------------+----------------------------+------------------+-------------------------+
   | GetNumBonds            | onlyHeavy:bool(T)          | int              | Bonds number            |
   +------------------------+----------------------------+------------------+-------------------------+
   | GetNumConformers       |                            | int              | Conformers number       |
   +------------------------+----------------------------+------------------+-------------------------+
   | GetNumHeavyAtoms       |                            | int              | Heavy atoms number      |
   +------------------------+----------------------------+------------------+-------------------------+
   | GetProp                | key:str,                   | object           | Property value          |
   +------------------------+----------------------------+------------------+-------------------------+
   |                        | autoConvert:bool(F)        |                  |                         |
   +------------------------+----------------------------+------------------+-------------------------+
   | GetPropNames           | includePrivate:bool(F),    | tuple            |Properties names         |
   +------------------------+----------------------------+------------------+-------------------------+
   |                        | includeComputed:bool(F)    |                  |                         |
   +------------------------+----------------------------+------------------+-------------------------+
   |  GetPropsAsDict        | includePrivate:bool(F),    | dict             |Properties in a dict     |
   +------------------------+----------------------------+------------------+-------------------------+
   |                        | includeComputed:bool(F),   |                  |                         |
   +------------------------+----------------------------+------------------+-------------------------+
   |                        |autoConvertStrings:bool(T)  |                  |                         |
   +------------------------+----------------------------+------------------+-------------------------+
   | GetRingInfo            |                            | RingInfo         | Ring's information      |
   +------------------------+----------------------------+------------------+-------------------------+
   | GetStereoGroups        |                            | StereoGroup_vect  | Lis: stereochemistry   |
   +------------------------+----------------------------+------------------+-------------------------+
   | GetSubstructMatch      | query:Mol,                 | tuple            | ID: substrucure's atoms |
   +------------------------+----------------------------+------------------+-------------------------+
   |                        | useChirality:bool(F),      |                  |                         |
   +------------------------+----------------------------+------------------+-------------------------+
   |                        |useQueryQueryMatches:bool(F)|                  |                         |
   +------------------------+----------------------------+------------------+-------------------------+
   | GetSubstructMatches    |query:Mol, uniquify:bool(T),| tuple            |ID: substructures' atoms |
   +------------------------+----------------------------+------------------+-------------------------+
   |                        |useChirality:bool(F),       |                  |                         |
   +------------------------+----------------------------+------------------+-------------------------+
   |                        |useQueryQueryMatches:bool(F)|                  |                         |
   +------------------------+----------------------------+------------------+-------------------------+
   |                        |maxMatches:int(1000)        |                  |                         |
   +------------------------+----------------------------+------------------+-------------------------+
   | GetUnsignedProp        | key:str                    | int              |Property with value > 0  |
   +------------------------+----------------------------+------------------+-------------------------+
   |.. centered:: **Setters**                                                                         |
   +------------------------+----------------------------+------------------+-------------------------+
   | SetBoolProp            | key:str, value:bool,       |                  | Property with bool value|  
   +------------------------+----------------------------+------------------+-------------------------+
   |                        |computed:bool(F)            |                  |                         |
   +------------------------+----------------------------+------------------+-------------------------+
   | SetDoubleProp          | key:str, value:double,     |                  |Property with float value|
   +------------------------+----------------------------+------------------+-------------------------+
   |                        |     computed:bool(F)       |                  |                         |
   +------------------------+----------------------------+------------------+-------------------------+
   | SetIntProp             | key:str, value:int,        |                  |Property with int value  |
   +------------------------+----------------------------+------------------+-------------------------+
   |                        |computed:bool(F)            |                  |                         |
   +------------------------+----------------------------+------------------+-------------------------+
   | SetProp                | key:str, value:str,        |                  |Property with str value  |
   +------------------------+----------------------------+------------------+-------------------------+
   |                        |computed:bool(F)            |                  |                         |
   +------------------------+----------------------------+------------------+-------------------------+
   | SetUnsignedProp        | key:str, value:str,        |                  |Property with value > 0  | 
   +------------------------+----------------------------+------------------+-------------------------+
   |                        | computed:bool(F)           |                  |                         |
   +------------------------+----------------------------+------------------+-------------------------+
   |.. centered :: **Generals**                                                                       |
   +------------------------+----------------------------+------------------+-------------------------+
   | AddConformer           | conf:Conformer,            | int              | Add Conformer           |
   +------------------------+----------------------------+------------------+-------------------------+
   |                        | assignId:bool(F)           |                  |                         |
   +------------------------+----------------------------+------------------+-------------------------+
   | ClearComputedProps     | includeRings:bool(T)       |                  | Remove properties       |
   +------------------------+----------------------------+------------------+-------------------------+
   | ClearProp              | key:str                    |                  | Clean properties        |
   +------------------------+----------------------------+------------------+-------------------------+
   | Compute2DCoords        | canonOrient:bool(T),       | int              | Compute atomic coord.   |
   +------------------------+----------------------------+------------------+-------------------------+
   |                        | clearConfs:bool(T),        |                  |                         |
   +------------------------+----------------------------+------------------+-------------------------+
   |                        | coordMap:dict,             |                  |                         |
   +------------------------+----------------------------+------------------+-------------------------+
   |                        | nFlipsPerSample:int(0),    |                  |                         |
   +------------------------+----------------------------+------------------+-------------------------+
   |                        | nSample:int(0),            |                  |                         |
   +------------------------+----------------------------+------------------+-------------------------+
   |                        | sampleSeed:int(0),         |                  |                         |
   +------------------------+----------------------------+------------------+-------------------------+
   |                        | permuteDeg4Nodes:bool(F),  |                  |                         |
   +------------------------+----------------------------+------------------+-------------------------+
   |                        | bondLength:float(-1.),     |                  |                         |
   +------------------------+----------------------------+------------------+-------------------------+
   |                        | forceRDKit:bool(F),        |                  |                         |
   +------------------------+----------------------------+------------------+-------------------------+
   |                        | useRingTemplates:bool(F)   |                  |                         |
   +------------------------+----------------------------+------------------+-------------------------+
   |ComputeGasteigerCharges | nIter:int(12),             |                  |Compute Gasteiger charge*|
   +------------------------+----------------------------+------------------+-------------------------+
   |                        |throwOnParamFailure:bool(F) |                  |                         |
   +------------------------+----------------------------+------------------+-------------------------+
   | Debug                  |  useStdout:bool(T)         |                  | Debugging information   |
   +------------------------+----------------------------+------------------+-------------------------+
   | HasProp                | key:str                    | int              | Query property          |
   +------------------------+----------------------------+------------------+-------------------------+
   | HasQuery               |                            | bool             | Query exist             |
   +------------------------+----------------------------+------------------+-------------------------+
   | HasSubstructMatch      | query:Mol,                 | bool             | Query substructure      |
   +------------------------+----------------------------+------------------+-------------------------+
   |                        | recursionPossible:bool(T), |                  |                         |
   +------------------------+----------------------------+------------------+-------------------------+
   |                        | useChirality:bool(F),      |                  |                         |
   +------------------------+----------------------------+------------------+-------------------------+
   |                        |useQueryQueryMatches:bool(F)|                  |                         |
   +------------------------+----------------------------+------------------+-------------------------+
   |NeedsUpdatePropertyCache|                            | bool             | Check valence           |
   +------------------------+----------------------------+------------------+-------------------------+
   | RemoveAllConformers    |                            |                  | Remove conformers       |
   +------------------------+----------------------------+------------------+-------------------------+
   | RemoveConformer        | id:int                     |                  | Remove conformer        |
   +------------------------+----------------------------+------------------+-------------------------+
   | ToBinary               |                            | object           | Binary representation   |
   +------------------------+----------------------------+------------------+-------------------------+
   | UpdatePropertyCache    | strict:bool(T)             |                  | Updated property        |
   +------------------------+----------------------------+------------------+-------------------------+

The AMCESS's Molecule class also has the next properties. Some properties represent shortcuts that combine methods from the Mol class with methods from RDKit's Atom class.

.. tabularcolumns:: p{0.12\linewidth}p{0.196\linewidth}p{0.30\linewidth}p{0.30\linewidth}
.. table:: Table 4 : Properties only of AMCESS's Molecule class.
   :name: tab_m2
   :widths: 30, 40, 20, 10
   :class: longtable
   :align: center
   :width: 66%

   +------------------------+----------------------------+------------------+-------------------------+
   |**Propety**             |  **Argument**              | **Return**       |  **Description**        |
   +========================+============================+==================+=========================+
   |.. centered::                             **attrs**                                               |
   +------------------------+----------------------------+------------------+-------------------------+
   | __attrs_init__         |  atoms:list                |                  | Re--init Molecule object|
   +------------------------+----------------------------+------------------+-------------------------+
   |.. centered::                             **Getters**                                             |
   +------------------------+----------------------------+------------------+-------------------------+
   | GetAtom                | atom:int                   | str              | Atom's properties       |
   +------------------------+----------------------------+------------------+-------------------------+
   | GetAtoms               |                            | list             | Atom objects            |
   +------------------------+----------------------------+------------------+-------------------------+
   | GetAtomicMasses        |                            | list             | Atomic masses           |
   +------------------------+----------------------------+------------------+-------------------------+
   | GetAtomicNumbers       |                            | list             | Atomuic numbers         |
   +------------------------+----------------------------+------------------+-------------------------+
   | GetAtomicSymbols       |                            | list             | Atomic Symbols          |
   +------------------------+----------------------------+------------------+-------------------------+
   | GetBlockXYZ            |                            | str              | Molecule object to  XYZ |
   +------------------------+----------------------------+------------------+-------------------------+
   | GetMolCoord            |                            | list             | Atomic Coordinates      |
   +------------------------+----------------------------+------------------+-------------------------+
   | GetNumberingAtoms      |                            | str              | Chain of coordinates    |
   +------------------------+----------------------------+------------------+-------------------------+
   | GetMolCharge           |                            | int              | Charge                  |
   +------------------------+----------------------------+------------------+-------------------------+
   | GetMolCM               |                            | tuple            | Mass center             |
   +------------------------+----------------------------+------------------+-------------------------+
   | GetMolDict             |                            | dict             | Molecule object to dict |
   +------------------------+----------------------------+------------------+-------------------------+
   | GetMolList             |                            | list             | Molecule object to list |
   +------------------------+----------------------------+------------------+-------------------------+
   | GetMolMass             |                            | float            | Molecular mass          |
   +------------------------+----------------------------+------------------+-------------------------+
   | GetMolMultiplicity     |                            | int              | Multiplicity            |
   +------------------------+----------------------------+------------------+-------------------------+
   | GetMolPrincipalAxes    |                            | list             | Principal axes          |
   +------------------------+----------------------------+------------------+-------------------------+
   |.. centered::                             **Setters**                                             |
   +------------------------+----------------------------+------------------+-------------------------+
   | SetMolCharge           | charge:int                 |                  | Charge                  |
   +------------------------+----------------------------+------------------+-------------------------+
   | SetMolMultiplicity     | multiplicity:int           |                  | Multiplicity            |
   +------------------------+----------------------------+------------------+-------------------------+
   |.. centered::                             **Generals**                                            |
   +------------------------+----------------------------+------------------+-------------------------+
   | AddAtom                | atoms:list,                | Molecule         | Add atoms               |
   +------------------------+----------------------------+------------------+-------------------------+
   |                        | charge:int(999999),        |                  |                         |
   +------------------------+----------------------------+------------------+-------------------------+
   |                        | multiplicity:int(-1)       |                  |                         |
   +------------------------+----------------------------+------------------+-------------------------+
   | RemoveAtom             | atoms:int,                 | Molecule         | Remove atoms            |
   +------------------------+----------------------------+------------------+-------------------------+
   |                        | charge:int(999999),        |                  |                         |
   +------------------------+----------------------------+------------------+-------------------------+
   |                        | multiplicity:int(-1)       |                  |                         |
   +------------------------+----------------------------+------------------+-------------------------+