import numpy as np
import pytest

from amcess.atom import Atom

# ===============================================================
# Atom class
# ===============================================================
@pytest.mark.parametrize(
    "atom, expected_result",
    [
        (("H", 0, 0, 0), {"element": "H", "x": 0, "y": 0, "z": 0}),
    ],
)
def test_atom_str_magic_function(atom, expected_result):
    """Testing str magic function of Atom class

    Examples
    --------
    >>>Atom('F', 0, 0, 1.97)
    {'element': 'F', 'x': 0, 'y': 0, 'z': 1.97}
    """
    input_atom = str(Atom(*atom))

    assert input_atom == str(expected_result)


@pytest.mark.parametrize(
    "wrong_atom",
    [
        (("H", 0, 0, "")),
        (("H", 0, 0, "z"))
    ],
)
def test_atom_coordinates_fails(wrong_atom):
    """Testing Atom's coordinates failures"""
    with pytest.raises(ValueError):
        Atom(*wrong_atom)

@pytest.mark.parametrize(
    "wrong_atom",
    [
        (("", 0, 0, 0)),
        (("'", 0, 0, 0)),
        (("@", 0, 0, 0)),
    ],
)
def test_atom_symbol_fails(wrong_atom):
    """Testing Atom's symbol failures"""
    with pytest.raises(RuntimeError):
        Atom(*wrong_atom)
