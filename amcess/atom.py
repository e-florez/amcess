from rdkit.Chem.rdchem import Atom as RDKAtom


class Atom(RDKAtom):
    # ===============================================================
    # Building
    # ===============================================================
    def __init__(
        self, element: str = "", x: float = 0.0, y: float = 0.0, z: float = 0.0
    ):
        """
        This class inherits attributes of the Atom class from RDKit,
        for more information:
            *) https://www.rdkit.org/docs/cppapi/classRDKit_1_1Atom.html
            *) https://www.rdkit.org/docs/source/rdkit.Chem.rdchem.html

        The idea of this class is to achieve of standard of amcess and
        add atom coordinates to Atom class from RDKit.

        !Class description:
        Representation of an individual atoms (<element> <X> <Y> <Z>)

        .. rubric:: Examples

        >>> Atom(element='H')
        {'element': 'H', 'x': 0., 'y': 0., 'z': 0.}
        NOTE: The coordinate by default are (0.,0.,0.)

        >>> Atom('F', 0, 0, 1.97)
        {'element': 'F', 'x': 0., 'y': 0., 'z': 1.97}

        .. rubric:: Returns

        atom : object
            object like dict {
                              'element': str,
                                    'x': float,
                                    'y': float,
                                    'z': float
                             }

        .. rubric:: Raises

        ValueError
            format MUST be (str, float, float, float) with NOT empty filed
        """
        # Build: Atom class from RDKit
        super().__init__(element)

        # Object coordinates
        self.SetCoord(x, y, z)

    # ===============================================================
    # MAGIC METHODS
    # ===============================================================
    def __str__(self):
        """Magic method '__str__' to print the object as a dictionary"""
        atom = self.GetSymbol()
        return f"atom: {atom}, x: {self.x}, y: {self.y}, z: {self.z}"

    # ===============================================================
    # PROPERTIES
    # ===============================================================

    #################################################################
    # ! Getter
    #################################################################
    def GetCoord(self) -> tuple[float, float, float]:
        "Return atomic coordinates"
        return (self.x, self.y, self.z)

    #################################################################
    # ! Setter
    #################################################################
    def SetCoord(self, x: float, y: float, z: float) -> None:
        "Change atomic coordinates"
        # Coordinate: Must be valid float or int
        for value in [x, y, z]:
            if not isinstance(value, (int, float)):
                raise ValueError(
                    "\n\nMust be valid NOT empty float"
                    f"\nyou get --> '{value}' with type:"
                    f"'{type(value).__name__}'"
                )
        self.x = x
        self.y = y
        self.z = z
