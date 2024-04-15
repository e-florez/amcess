import attr


from rdkit.Chem.rdchem import Atom as RDKAtom

from .data.atomic_data import atomic_mass


@attr.s(frozen=True)
class Atom(RDKAtom):
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
        object like dict {'element': str, 'x': float, 'y': float, 'z': float}

    .. rubric:: Raises

    ValueError
        format MUST be (str, float, float, float) with NOT empty filed
    """

    element: str = attr.ib()
    x: int = attr.ib()
    y: int = attr.ib()
    z: int = attr.ib()

    # ===============================================================
    # VALIDATORS
    # ===============================================================
    # ! super.__init__() already check element
    # @element.validator
    # def _check_valid_element(self, element, value):
    #     """Element: Must be valid NOT empty alphanumeric character"""
    #     if not value.isalnum():
    #         raise ValueError(
    #             "\n\nMust be valid NOT empty alphanumeric character"
    #             f"\nyou get --> '{value}'\n"
    #         )

    @x.validator
    @y.validator
    @z.validator
    def _check_valid_point(self, coordinate, value):
        """Coordinate: Must be valid float"""
        if not isinstance(value, (int, float)):
            raise ValueError(
                "\n\nMust be valid NOT empty float"
                f"\nyou get --> '{value}' with type: '{type(value).__name__}'"
            )

    # ===============================================================
    # MAGIC METHODS
    # ===============================================================
    def __str__(self):
        """Magic method '__str__' to print the object as a dictionary"""
        return str(attr.asdict(self))
    
    # ===============================================================
    # Building
    # ===============================================================
    def __attrs_pre_init__(self, element, x: float = 0.0, y: float = 0.0, z: float = 0.0):
        """
        __attrs_pre_init__ is automatically detected and run before attrs
        starts initializing. If __attrs_pre_init__ takes more than the
        self argument, the attrs-generated __init__ will call it with the
        same arguments it received itself. This is useful if you need to
        inject a call to super().__init__()
        """
        super().__init__(element)