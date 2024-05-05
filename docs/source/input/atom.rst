Atom class
----------

The Atom class is basically the Atom class from the RDKit program, 
because it has the properties that are needed within the AMCESS program.

The web page where the documentation for the RDKit's Atom class is the following:

.. centered:: `<https://www.rdkit.org/docs/source/rdkit.Chem.rdchem.html>`_

The RDKit's Atom class is inside the Chem package (`rdkit.Chem <https://www.rdkit.org/docs/source/rdkit.html>`_), 
because it belongs to the rdchem module (`rdkit.Chem.rdchem <https://www.rdkit.org/docs/source/rdkit.Chem.html>`_).

Constructor
^^^^^^^^^^^

The Atom class constructor is presented in the following code lines
(note that the line number may not match the line number in the .py file):

.. _python-a1:
.. code-block:: python
    :caption: : Constructor of the Atom class.
    :linenos:
    :emphasize-lines: 2, 7, 10

     def __init__(
        self, element: str = "", x: float = 0.0, y: float = 0.0, z: float = 0.0
     ):
     """ ... """

        # Build: Atom class from RDKit
        super().__init__(element)

        # Object coordinates
        self.SetCoord(x, y, z)

Constructor's arguments are in lines 2 which are element: *str* and, x: *float, int*, y: *float, int*, z: *float, int*. 
The element argument receives the atomic symbol, while x, y, z are responsible for storing the coordinates of the atom.

In the line 7 AMCESS's Atom class inheritances all propierties from RDKit's Atom class. In the line 10, the SetCoord 
method is used to verify and storing the arguments x, y, z.

To avoid errors, the RDKit's Atom class is instantiated before assigning x, y, z as properties of the AMCESS's Atom 
object, because in the contrary case the x, y, z have to be put as arguments in super().__init__(); which would throw 
an error, due to that the RDKit's Atom class only receives as argument the atomic symbol.

Properties
^^^^^^^^^^

RDKit's Atom class is composed for he next properties:

.. tabularcolumns:: p{0.12\linewidth}p{0.196\linewidth}p{0.30\linewidth}p{0.30\linewidth}
.. table:: Table 1 : Properties from RDKit's Atom class (* Property name).
   :name: tab_a1
   :widths: 30, 40, 20, 10
   :class: longtable
   :align: center
   :width: 66%

   +------------------------+--------------------------+------------------+-------------------------+
   |**Propety**             |  **Argument**            | **Return**       |  **Description**        |
   +========================+==========================+==================+=========================+
   |.. centered::                             **Getters**                                           |
   +------------------------+--------------------------+------------------+-------------------------+
   | GetAtomMapNum          |                          | int              | Atom's indexes          |
   +------------------------+--------------------------+------------------+-------------------------+
   | GetAtomicNum           |                          | int              | Atomic number           |
   +------------------------+--------------------------+------------------+-------------------------+
   | GetBonds               |                          | tuple            | Bonds                   |
   +------------------------+--------------------------+------------------+-------------------------+
   | GetBoolProp            | key:str*                 | bool             | Property with bool value|
   +------------------------+--------------------------+------------------+-------------------------+
   | GetChiralTag           |                          | ChiralType       |                         |
   +------------------------+--------------------------+------------------+-------------------------+
   | GetDegree              |                          | int              | Bonds number            |
   +------------------------+--------------------------+------------------+-------------------------+
   | GetDoubleProp          | key:str*                 | float            |Property with float value|
   +------------------------+--------------------------+------------------+-------------------------+
   |GetExplicitBitVectProp  | key:str*                 | ExplicitBitVect  | Property with           |
   |                        |                          |                  | ExplicitBitVect value   |
   +------------------------+--------------------------+------------------+-------------------------+
   | GetExplicitValence     |                          | int              | Valence                 |
   +------------------------+--------------------------+------------------+-------------------------+
   | GetFormalCharge        |                          | int              | Charge                  |
   +------------------------+--------------------------+------------------+-------------------------+
   | GetHybridization       |                          |HybridizationType | Hybridization           |
   +------------------------+--------------------------+------------------+-------------------------+
   | GetIdx                 |                          | int              | Atom's index            |
   +------------------------+--------------------------+------------------+-------------------------+
   | GetImplicitValence     |                          | int              | Implicit Hs             |
   +------------------------+--------------------------+------------------+-------------------------+
   | GetIntProp             | key:str*                 | int              |Property with int value  |
   +------------------------+--------------------------+------------------+-------------------------+
   | GetIsAromatic          |                          | bool             | Aromatic                |
   +------------------------+--------------------------+------------------+-------------------------+
   | GetIsotope             |                          | int              | Isotope                 |
   +------------------------+--------------------------+------------------+-------------------------+
   | GetMass                |                          | float            | Atomic mass             |
   +------------------------+--------------------------+------------------+-------------------------+
   | GetMonomerInfo         |                          | AtomMonomerInfo  | MonomerInfo object      |
   +------------------------+--------------------------+------------------+-------------------------+
   | GetNeighbors           |                          | tuple            | Neighbors atoms         |
   +------------------------+--------------------------+------------------+-------------------------+
   | GetNoImplicit          |                          | bool             | Implicit Hs             |
   +------------------------+--------------------------+------------------+-------------------------+
   | GetNumExplicitHs       |                          | int              | Number of explicit Hs   |
   +------------------------+--------------------------+------------------+-------------------------+
   | GetNumImplicitHs       |                          | int              | Number of Implicit Hs   |
   +------------------------+--------------------------+------------------+-------------------------+
   |GetNumRadicalElectrons  |                          | int              | Unpaired e-s            |
   +------------------------+--------------------------+------------------+-------------------------+
   | GetOwningMol           |                          | Mol              | Mol object              |
   +------------------------+--------------------------+------------------+-------------------------+
   | GetPDBResidueInfo      |                          |AtomPDBResidueInfo|MonomerInfo object (PDB) |
   +------------------------+--------------------------+------------------+-------------------------+
   | GetProp                | key:str*,                | object           | Property with object    |
   |                        | autoConvert:bool(F)      |                  | value                   |
   +------------------------+--------------------------+------------------+-------------------------+
   | GetPropNames           |includePrivate:bool(F),   | list             | list of properties      |
   +------------------------+--------------------------+------------------+-------------------------+
   |                        |includeComputed:bool(F)   |                  |                         |
   +------------------------+--------------------------+------------------+-------------------------+
   | GetPropsAsDict         |includePrivate:bool(T),   | dict             | dict of properties      |
   +------------------------+--------------------------+------------------+-------------------------+
   |                        |includeComputed:bool(T),  |                  |                         |
   +------------------------+--------------------------+------------------+-------------------------+
   |                        |autoConvertStrings:bool(T)|                  |                         |
   +------------------------+--------------------------+------------------+-------------------------+
   | GetQueryType           |                          | str              |                         |
   +------------------------+--------------------------+------------------+-------------------------+
   | GetSmarts              |doKekule:bool(F),         | str              | SMART or SMILES         |
   +------------------------+--------------------------+------------------+-------------------------+
   |                        |allHsExplicit:bool(F),    |                  |                         |
   +------------------------+--------------------------+------------------+-------------------------+
   |                        |isomericSmiles:bool(T)    |                  |                         |
   +------------------------+--------------------------+------------------+-------------------------+
   | GetSymbol              |                          | str              | Atomic symbol           |
   +------------------------+--------------------------+------------------+-------------------------+
   | GetTotalDegree         |                          | int              | All Bonds               |
   +------------------------+--------------------------+------------------+-------------------------+
   | GetTotalNumHs          | includeNeighbors:bool(F) | int              | Hs number               |
   +------------------------+--------------------------+------------------+-------------------------+
   | GetTotalValence        |                          | int              | Total valence           |
   +------------------------+--------------------------+------------------+-------------------------+
   | GetUnsignedProp        | key:str*                 | int              | Property with value > 0 |
   +------------------------+--------------------------+------------------+-------------------------+
   |.. centered:: **Setters**                                                                       |
   +------------------------+--------------------------+------------------+-------------------------+
   | SetAtomMapNum          | mapno:int, strict:bool(F)|                  | Atom's index or errase  |
   +------------------------+--------------------------+------------------+-------------------------+
   | SetAtomicNum           | newNum:int               |                  | Atomic number           |
   +------------------------+--------------------------+------------------+-------------------------+
   | SetBoolProp            | key:str, val:bool        |                  | Property value          |
   +------------------------+--------------------------+------------------+-------------------------+
   | SetChiralTag           | what:ChiralType          |                  | Chiral                  |
   +------------------------+--------------------------+------------------+-------------------------+
   | SetDoubleProp          | key:str, val:float       |                  | Property value          |
   +------------------------+--------------------------+------------------+-------------------------+
   |SetExplicitBitVectProp  | key:str,                 |                  | Property value          |
   +------------------------+--------------------------+------------------+-------------------------+
   |                        | val:ExplicitBitVect      |                  |                         |
   +------------------------+--------------------------+------------------+-------------------------+
   | SetFormalCharge        | what:int                 |                  | Charge                  |
   +------------------------+--------------------------+------------------+-------------------------+
   | SetHybridization       | what:HybridizationType   |                  | Hybridization           |
   +------------------------+--------------------------+------------------+-------------------------+
   | SetIntProp             | key:str, val:int         |                  | Property value          |
   +------------------------+--------------------------+------------------+-------------------------+
   | SetIsAromatic          | what:bool                |                  | Aromatic                |
   +------------------------+--------------------------+------------------+-------------------------+
   | SetIsotope             | what:int                 |                  | Isotope                 |
   +------------------------+--------------------------+------------------+-------------------------+
   | SetMonomerInfo         | info:AtomMonomerInfo     |                  |                         |
   +------------------------+--------------------------+------------------+-------------------------+
   | SetNoImplicit          | what:bool                |                  | Atom without implicit Hs| 
   +------------------------+--------------------------+------------------+-------------------------+
   | SetNumExplicitHs       | what:int                 |                  | Hs explicits            |
   +------------------------+--------------------------+------------------+-------------------------+
   |SetNumRadicalElectrons  | num:int                  |                  | Unpaired e-s            |
   +------------------------+--------------------------+------------------+-------------------------+
   | SetPDBResidueInfo      | info:AtomMonomerInfo     |                  | MonomerInfo Object      |
   +------------------------+--------------------------+------------------+-------------------------+
   | SetProp                | key:str, val:str         |                  | Property value          |
   +------------------------+--------------------------+------------------+-------------------------+
   | SetUnsignedProp        | key:str, val:int         |                  | Property value (> 0)    |
   +------------------------+--------------------------+------------------+-------------------------+
   |.. centered:: **Generals**                                                                      |
   +------------------------+--------------------------+------------------+-------------------------+
   | ClearProp              | key:str*                 |                  | Errase property data    |
   +------------------------+--------------------------+------------------+-------------------------+
   | DescribeQuery          |                          | str              | Debugging               |
   +------------------------+--------------------------+------------------+-------------------------+
   | HasOwningMol           |                          | bool             | Check atom belongs mol  |
   +------------------------+--------------------------+------------------+-------------------------+
   | HasProp                | key:str*                 | int              | Check property          |
   +------------------------+--------------------------+------------------+-------------------------+
   | HasQuery               |                          | bool             | Check if there's a query|
   +------------------------+--------------------------+------------------+-------------------------+
   | HasValenceViolation    |                          | bool             | Check valence           |
   +------------------------+--------------------------+------------------+-------------------------+
   | InvertChirality        |                          | bool             |                         |
   +------------------------+--------------------------+------------------+-------------------------+
   | IsInRing               |                          | bool             | Check atom belongs ring |
   +------------------------+--------------------------+------------------+-------------------------+
   | IsInRingSize           | size:int                 | bool             | Check atom belongs ring |
   +------------------------+--------------------------+------------------+-------------------------+
   | Match                  | whar:Atom                | bool             | Atoms compare           |
   +------------------------+--------------------------+------------------+-------------------------+
   |NeedsUpdatePropertyCache|                          | bool             | Check mol valence       |
   +------------------------+--------------------------+------------------+-------------------------+
   | UpdatePropertyCache    | strict:bool(T)           |                  | Regenerate properties   |
   +------------------------+--------------------------+------------------+-------------------------+

Because the RDKit's Atom class does not contain properties associated with the coordinates of the atoms, 
see :ref:`tab_a1`, it is created an AMCESS's Atom class. This class inherits everything from RDKit's Atom 
class but also stores the coordinates of an atom in the **object variables**: x, y, z. Then, the atom's 
coordinates can be consulted through its object in the following way:

.. code-block:: python
    :caption: : Atom's coordinates.
    :linenos:

     atom = Atom('C',1.0,2.0,3.0)
     print(f"Coordinates (x,y,z): ({atom.x}, {atom.y}, {atom.z})")
     Coordinates (x,y,z): (1.0, 2.0, 3.0)

Atom class objects can be created without defining the coordinates, in which case the coordinates are 0.0.

.. code-block:: python
    :caption: : Atom object with coordinates by default.
    :linenos: 
    
     atom = Atom('C')
     print(f"Coordinates (x,y,z): ({atom.x}, {atom.y}, {atom.z})")
     Coordinates (x,y,z): (0.0, 0.0, 0.0)

Due to the above, the AMCESS's Atom class also has the properties shown in the following table:

.. tabularcolumns:: p{0.12\linewidth}p{0.196\linewidth}p{0.30\linewidth}p{0.30\linewidth}
.. table:: : Properties from AMCCES's Atom class.
   :name: tab_a2
   :widths: 30, 40, 20, 10
   :class: longtable
   :align: center
   :width: 66%

   +------------------------+--------------------------+------------------+-------------------------+
   |**Propety**             |  **Argument**            | **Return**       |  **Description**        |     
   +------------------------+--------------------------+------------------+-------------------------+
   |.. centered:: **Magics**                                                                        |
   +------------------------+--------------------------+------------------+-------------------------+
   | __str__                |                          | str              | Object's dictionary     |
   +------------------------+--------------------------+------------------+-------------------------+
   |.. centered:: **Getters**                                                                       |
   +------------------------+--------------------------+------------------+-------------------------+
   | GetCoord               |                          | tuple            | Atomic coordinates      |
   +------------------------+--------------------------+------------------+-------------------------+
   |.. centered:: **Setters**                                                                       |
   +------------------------+--------------------------+------------------+-------------------------+
   | SetCoord               | x:float, y:float, z:float|                  | Atomic coordinates      |
   +------------------------+--------------------------+------------------+-------------------------+

AMCCES's Atom class is enriched with the __str__ magic method, which prints the atom objects on the screen like dictionaries.

.. code-block:: python
    :caption: Listing 5 : Print from __str__ magic method.
    :linenos: 
    
     atom = Atom('C')
     str(atom)
     "{'element': 'C', 'x': 1.0, 'y': 0.4, 'z': 0.9}"
     print(atom)
     {'element': 'C', 'x': 1.0, 'y': 0.4, 'z': 0.9}