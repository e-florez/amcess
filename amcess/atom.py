import attr


from rdkit.Chem.rdchem import Atom as RDKAtom

from .data.atomic_data import atomic_mass


@attr.s(frozen=True)
class Atom(RDKAtom):
    # """
    # Representation of an individual atomas (<element> <X> <Y> <Z>)

    # .. rubric:: Examples

    # >>> Atom(element='H', x=0, y=0, z=0)
    # {'element': 'H', 'x': 0, 'y': 0, 'z': 0}

    # >>> Atom('F', 0, 0, 1.97)
    # {'element': 'F', 'x': 0, 'y': 0, 'z': 1.97}

    # .. rubric:: Returns

    # atom : object
    #     object like dict {'element': str, 'x': float, 'y': float, 'z': float}

    # .. rubric:: Raises

    # ValueError
    #     format MUST be (str, float, float, float) with NOT empty filed
    # """

    element: str = attr.ib()
    x: int = attr.ib()
    y: int = attr.ib()
    z: int = attr.ib()

    # ===============================================================
    # VALIDATORS
    # ===============================================================
    @element.validator
    def _check_valid_element(self, element, value):
        """Element: Must be valid NOT empty alphanumeric character"""
        if not value.isalnum():
            raise ValueError(
                "\n\nMust be valid NOT empty alphanumeric character"
                f"\nyou get --> '{value}'\n"
            )

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
    
    # # ===============================================================
    # # Building
    # # ===============================================================
    def __attrs_pre_init__(self, element, x: float = 0.0, y: float = 0.0, z: float = 0.0):
        super().__init__(element)

