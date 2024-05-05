Cluster class
-------------

The Cluster class inherits from the Molecule class and Atom class (:ref:`fig_i1`).

Constructor
^^^^^^^^^^^

The fundamental property of Cluster class is the object variable self.cluster_dict (line=9 in :ref:`python-c1`), which is a dictionary of list of Molecule objects. The *args argument receives the different formats of the molecules or atoms. The formats allowed are: Cluster, Molecule, Mol, dict and list (lines=36--50). **The formats of the different files that can contain the molecular information and the smiles notation need to be added**.

.. _python-c1:
.. code-block:: python
    :caption: Listing 11 : Cluster class' constructor.
    :linenos:
    :emphasize-lines: 9, 36-50

    def __init__(
        self,
        *args,
        freeze_molecule = None,
        sphere_radius = None,
        sphere_center: tuple[float, float, float] = (0., 0., 0.),
        seed: int = None,
    ):
        self._cluster_dict: dict = dict()
        self._multiplicity: int = 1
        self._charge: int = 0

        # fixing molecule to NOT move or rotate
        # initialize with an empty list
        self._freeze_molecule: list = (  # noqa
            [] if freeze_molecule is None else freeze_molecule
        )

        self._sphere_radius: float = sphere_radius
        self._sphere_center: tuple[float, float, float] = sphere_center

        # initialize random generator
        if not seed:
            self._seed: int = np.random.randint(0, 999999)
        else:
            self._seed = seed

        # ----------------------------------------------------
        # attrs post-initialization

        self.cluster_atoms: list[Molecule] = list()

        # for count, mol in enumerate(args):
        for mol in args:
            size: int = len(self._cluster_dict)
            if isinstance(mol, Cluster):
                for j in mol._cluster_dict:
                    self._cluster_dict[size + j] = mol._cluster_dict[j]
                self._charge += mol.GetMolCharge()
                self.cluster_atoms += mol.GetMolList()
                # restarting the loop
                continue
            elif isinstance(mol, Molecule):
                new_molecule: Molecule = mol
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

            self.cluster_atoms += new_molecule.GetMolList() 
            # ! how is computed the cluster total multiplicity?
            self._charge += new_molecule.GetMolCharge()
            self._cluster_dict[size] = new_molecule

        # ! initializing Cluster as a 'Molecule' (sum of all individual ones)
        super().__init__(
            atoms=self.cluster_atoms,
            charge=self._charge,
            multiplicity=self._multiplicity,
        )

The self._sphere_radius and self._sphere_center varaibles are used to define the spherical box where is the cluster. On the other hand, the self._freeze_molecule property is used to select the molecule(s) that will be left still during the study. Finally, the self._seed stores the seed related to the random movement.

Properties
^^^^^^^^^^

Cluster class properties are:

.. tabularcolumns:: p{0.12\linewidth}p{0.196\linewidth}p{0.30\linewidth}p{0.30\linewidth}
.. table:: Table 5 : Cluster class properties.
   :name: tab_c1
   :widths: 30, 40, 20, 10
   :class: longtable
   :align: center
   :width: 66%

   +------------------------+----------------------------+-------------------+-------------------------+
   |**Propety**             |  **Argument**              | **Return**        |  **Description**        |
   +========================+============================+===================+=========================+
   |.. centered::                             **Magics**                                               |
   +------------------------+----------------------------+-------------------+-------------------------+
   | __add__                | other:Cluster              | Cluster           | Sum objects             |
   +------------------------+----------------------------+-------------------+-------------------------+
   | __mul__                | value:int                  | Cluster           | __rmul__                |
   +------------------------+----------------------------+-------------------+-------------------------+
   | __rmul__               | value:int                  | Cluster           | Multiply object         |
   +------------------------+----------------------------+-------------------+-------------------------+ 
   | __str__                |                             | str              | Print objects           |
   +------------------------+----------------------------+-------------------+-------------------------+
   |.. centered :: **Getters**                                                                         |
   +------------------------+----------------------------+-------------------+-------------------------+
   | GetClusterDict         |                            | dict              | Cluster object to dict  |
   +------------------------+----------------------------+-------------------+-------------------------+
   | GetClusterList         |                            | list              | Cluster object to list  |
   +------------------------+----------------------------+-------------------+-------------------------+
   | GetFreezeMol           |                            | list or int       | ID: freeze molecule     |
   +------------------------+----------------------------+-------------------+-------------------------+
   | GetMol                 |                            | Molecule          | Select molecule         |
   +------------------------+----------------------------+-------------------+-------------------------+
   | GetMols                |                            | list              | Molecule list           |
   +------------------------+----------------------------+-------------------+-------------------------+
   | GetRandomGen           |                            |np.random.Generator| Random number generator |
   +------------------------+----------------------------+-------------------+-------------------------+
   | GetSphereCenter        |                            | tuple             | Sphere center           |
   +------------------------+----------------------------+-------------------+-------------------------+
   | GetSphereR             |                            | float             | Sphere radiu            |
   +------------------------+----------------------------+-------------------+-------------------------+
   | GetTotalMol            |                            | int               | Molecules number        |
   +------------------------+----------------------------+-------------------+-------------------------+
   | GetSeed                |                            | int               | Seed for generator      |
   +------------------------+----------------------------+-------------------+-------------------------+
   |.. centered:: **Setters**                                                                          |
   +------------------------+----------------------------+-------------------+-------------------------+
   | SetFreezeMol           | id:int or list[int]        |                   | Freeze molecule(s)      |
   +------------------------+----------------------------+-------------------+-------------------------+
   | SetSeed                | seed:int                   |                   | Seed for generator      |
   +------------------------+----------------------------+-------------------+-------------------------+
   | SetSphereCenter        | center:tuple               |                   | Sphere center           |
   +------------------------+----------------------------+-------------------+-------------------------+
   | SetSphereR             | radius:float               |                   | Sphere radiu            |
   +------------------------+----------------------------+-------------------+-------------------------+
   |.. centered:: **Generals**                                                                         |
   +------------------------+----------------------------+-------------------+-------------------------+
   | CalCentRSphere         | tolerance_radius:float(1)  | Cluster           |Compute: Sphere's values |
   +------------------------+----------------------------+-------------------+-------------------------+
   | InitializeCluster      | max_closeness:float(1)     | Cluster           | Avoid overlap           |
   +------------------------+----------------------------+-------------------+-------------------------+
   | Overlapping            | first_coord:list,          | bool              | Check overlap           |
   +------------------------+----------------------------+-------------------+-------------------------+
   |                        | second_coord:list,         |                   |                         |
   +------------------------+----------------------------+-------------------+-------------------------+
   |                        | max_closeness:float(1)     |                   |                         |
   +------------------------+----------------------------+-------------------+-------------------------+
   | MoveMol                | molecule:int(0),           | Cluster           |Move atomo/molecule      |
   +------------------------+----------------------------+-------------------+-------------------------+
   |                        | max_step:float,            |                   |                         |
   +------------------------+----------------------------+-------------------+-------------------------+
   |                        | max_rotation:float,        |                   |                         |
   +------------------------+----------------------------+-------------------+-------------------------+
   |                        | max_closeness:int(1)       |                   |                         |
   +------------------------+----------------------------+-------------------+-------------------------+
   | RemoveMol              | molecule:int               | Cluster           | Remove atom/molecule    |
   +------------------------+----------------------------+-------------------+-------------------------+
   | RotateMol              | molecule:int,              | Cluster           | Rotate molecule         |
   +------------------------+----------------------------+-------------------+-------------------------+
   |                        | x:float(0),                |                   |                         |
   +------------------------+----------------------------+-------------------+-------------------------+
   |                        | y:float(0),                |                   |                         |
   +------------------------+----------------------------+-------------------+-------------------------+
   |                        | z:float(0)                 |                   |                         |
   +------------------------+----------------------------+-------------------+-------------------------+
   | TranslateMol           | molecule:int,              | Cluster           | Translate atom/molecule |
   +------------------------+----------------------------+-------------------+-------------------------+ 
   |                        | x:float(0),                |                   |                         |
   +------------------------+----------------------------+-------------------+-------------------------+
   |                        | y:float(0),                |                   |                         |
   +------------------------+----------------------------+-------------------+-------------------------+
   |                        | z:float(0)                 |                   |                         |
   +------------------------+----------------------------+-------------------+-------------------------+