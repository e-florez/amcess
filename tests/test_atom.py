import numpy as np
import pytest

from amcess.atom import Atom

COORDINATES = {
    "dummy": [("A", 10, 20, 30), ("X0", -0.5, 0, -10)],
    "hf": [("H", 0, 0, 0), ("F", 0.917, 0, 0)],
    "water": [
        ("O", 0, 0, 0),
        ("H", 0.58708, 0.75754, 0),
        ("H", -0.58708, 0.75754, 0),
    ],
}


# ===============================================================
# Atom class
# ===============================================================
@pytest.mark.parametrize(
    "atom, expected_result",
    [
        (("H", 0, 0, 0), {"element": "H", "x": 0, "y": 0, "z": 0}),
    ],
)
def test_atom_class(atom, expected_result):
    """Testing Atom class

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
        (("H", 0, 0, "z")),
        (("", 0, 0, 0)),
        (("'", 0, 0, 0)),
        (("@", 0, 0, 0)),
    ],
)
def test_atom_class_fails(wrong_atom):
    """Testing Atom class failures"""
    with pytest.raises(ValueError):
        Atom(*wrong_atom)


@pytest.mark.parametrize(
    "atom, expected_result",
    [
        (("O", 0, 0, 0), 16),
    ],
)
def test_atom_mass(atom, expected_result):
    """Testing Atomic mass"""
    input_atom = Atom(*atom)

    assert input_atom.atomic_mass - expected_result < 0.1


@pytest.mark.parametrize(
    "atom, expected_result",
    [
        (("Xe20", 0, 0, 0), "Xe20"),
    ],
)
def test_atom_symbol(atom, expected_result):
    """Testing Atomic mass"""
    input_atom = Atom(*atom)

    assert input_atom.symbol == expected_result
