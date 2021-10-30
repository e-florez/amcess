import numpy as np
from pyscf import gto, scf

from amcess.base_molecule import Cluster


class ElectronicEnergy:
    def __init__(self, object_system, sphere_center, sphere_radius) -> None:
        """
        Class to calculate electronic energy

        Attributes
        ----------
            molecule_object : object
                Object initialized with Molecule or Cluster class
        """

        self._object_system_initial = object_system
        self._object_system_before = object_system
        self._object_system_current = object_system
        self._sphere_center = sphere_center
        self._sphere_radius = sphere_radius

    # ===============================================================
    # PROPERTIES
    # ===============================================================
    @property
    def object_system_initial(self):
        return self._object_system_initial

    @object_system_initial.setter
    def object_system_initial(self, new_object_system):
        self._object_system_inital = new_object_system
        self._object_system_before = new_object_system
        self._object_system_current = new_object_system

    @property
    def object_system_before(self):
        return self._object_system_before

    @property
    def object_system_current(self):
        return self._object_system_current

    @object_system_current.setter
    def object_system_current(self, new_object_system):
        self._object_system_before = self._object_system_current
        self._object_system_current = new_object_system

    # ===============================================================
    # Methods
    # ===============================================================
    def object_system_current(self, new_object_system):
        self._object_system_before = self._object_system_current
        self._object_system_current = new_object_system

    def beforeatoms(self, system):
        """
        Store attribute of before object made with Molecule

        Args:
            system ([type]): [description]

        Returns:
            [type]: [description]
        """
        self.beforeatoms.system = system or self.beforeatoms.system
        return self.beforeatoms.system


def build_input_pyscf(x_random, obj_ee, type_search):
    """
    Build input to pyscf

    Parameters
    ----------
        x : array 1D
            possible new positions and angles.
        system_object : object
            Object initialized with Molecule or Cluster class
        icall : integer
            number of call

    Returns
    -------
        [type]: [description]
    """

    system_object = obj_ee._object_system_current

    new_geom = dict()
    # Rotate and translate
    for i in range(system_object.total_molecules):
        new_geom[i] = {
            "atoms": system_object.move_molecules(
                i,
                (x_random[i * 3], x_random[i * 3 + 1], x_random[i * 3 + 2]),
                (
                    x_random[(i + system_object.total_molecules) * 3],
                    x_random[(i + system_object.total_molecules) * 3 + 1],
                    x_random[(i + system_object.total_molecules) * 3 + 2],
                ),
                obj_ee._sphere_radius,
                obj_ee._sphere_center,
            )
            .get_molecule(i)
            .atoms
        }

    system_object = Cluster(*new_geom.values()).initialize_cluster()

    obj_ee.object_system_current(system_object)

    # Build input to pyscf
    symbols = system_object.symbols
    input_mol = "'"
    for i in range(system_object.total_atoms):
        input_mol += str(symbols[i])
        for j in range(3):
            input_mol += "  " + str(system_object.coordinates[i][j])
        if i < system_object.total_atoms - 1:
            input_mol += "; "
        else:
            input_mol += " '"

    return input_mol, system_object


def hf_pyscf(x, *args, ncall=[0]):
    """
    Calculate of electronic energy with pyscf

    Parameters
    ----------
        x : array 1D
            Possible new positions and angles
        shgo #?
        args : list
            basis set, Object of Molecule, name output xyz,
            type optimizaiton

    Returns
    -------
        e : float
            electronic energy

    """

    input_pyscf, new_object = build_input_pyscf(x, args[1], args[3])

    mol = gto.M(
        atom=input_pyscf,
        basis=args[0],
    )

    e = scf.HF(mol).kernel()
    args[2].write(str(new_object.total_atoms) + "\n")
    args[2].write("Energy: " + str(e) + "\n")
    l: int = 0
    for symbols in new_object.symbols:
        args[2].write(
            str(symbols)
            + "  "
            +
            # 1 A = 1.88973 Bohr
            str(new_object.atoms[l][1] / 1.88973)
            + "  "
            + str(new_object.atoms[l][2] / 1.88973)
            + "  "
            + str(new_object.atoms[l][3] / 1.88973)
            + "\n"
        )
        l: int = l + 1

    ncall[0] += 1  # count calls
    return e
