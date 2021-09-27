def atomic_mass(symbol):
    """
    Atomic Mass for any element in the Periodic Table.
    Taken from `PubChem (NIH)`_

    .. warning:: Returns ZERO if it is not defined

    Examples
    --------
    >>> print(atomic_mass("O"))
    15.999
    >>> print(atomic_mass("ni"))
    58.693
    >>> print(atomic_mass("XYZ"))
    0

    Parameters
    ----------
    symbol : str
        atomic symbol case-insensitive

    Returns
    -------
    atomic mass : float
        atomic mass from the periodic table if symbol is not defined
        it returns ONE

    .. _PubChem (NIH):
        https://pubchem.ncbi.nlm.nih.gov/ptable/

    """

    ATOMIC_MASS = {
        "H": 1.008,
        "He": 4.0026,
        "Li": 7,
        "Be": 9.012183,
        "B": 10.81,
        "C": 12.011,
        "N": 14.007,
        "O": 15.999,
        "F": 18.99840316,
        "Ne": 20.18,
        "Na": 22.9897693,
        "Mg": 24.305,
        "Al": 26.981538,
        "Si": 28.085,
        "P": 30.973762,
        "S": 32.07,
        "Cl": 35.45,
        "Ar": 39.9,
        "K": 39.0983,
        "Ca": 40.08,
        "Sc": 44.95591,
        "Ti": 47.867,
        "V": 50.9415,
        "Cr": 51.996,
        "Mn": 54.93804,
        "Fe": 55.84,
        "Co": 58.93319,
        "Ni": 58.693,
        "Cu": 63.55,
        "Zn": 65.4,
        "Ga": 69.723,
        "Ge": 72.63,
        "As": 74.92159,
        "Se": 78.97,
        "Br": 79.9,
        "Kr": 83.8,
        "Rb": 85.468,
        "Sr": 87.62,
        "Y": 88.90584,
        "Zr": 91.22,
        "Nb": 92.90637,
        "Mo": 95.95,
        "Tc": 96.90636,
        "Ru": 101.1,
        "Rh": 102.9055,
        "Pd": 106.42,
        "Ag": 107.868,
        "Cd": 112.41,
        "In": 114.818,
        "Sn": 118.71,
        "Sb": 121.76,
        "Te": 127.6,
        "I": 126.9045,
        "Xe": 131.29,
        "Cs": 132.905452,
        "Ba": 137.33,
        "La": 138.9055,
        "Ce": 140.116,
        "Pr": 140.90766,
        "Nd": 144.24,
        "Pm": 144.91276,
        "Sm": 150.4,
        "Eu": 151.964,
        "Gd": 157.2,
        "Tb": 158.92535,
        "Dy": 162.5,
        "Ho": 164.93033,
        "Er": 167.26,
        "Tm": 168.93422,
        "Yb": 173.05,
        "Lu": 174.9668,
        "Hf": 178.49,
        "Ta": 180.9479,
        "W": 183.84,
        "Re": 186.207,
        "Os": 190.2,
        "Ir": 192.22,
        "Pt": 195.08,
        "Au": 196.96657,
        "Hg": 200.59,
        "Tl": 204.383,
        "Pb": 207,
        "Bi": 208.9804,
        "Po": 208.98243,
        "At": 209.98715,
        "Rn": 222.01758,
        "Fr": 223.01973,
        "Ra": 226.02541,
        "Ac": 227.02775,
        "Th": 232.038,
        "Pa": 231.03588,
        "U": 238.0289,
        "Np": 237.048172,
        "Pu": 244.0642,
        "Am": 243.06138,
        "Cm": 247.07035,
        "Bk": 247.07031,
        "Cf": 251.07959,
        "Es": 252.083,
        "Fm": 257.09511,
        "Md": 258.09843,
        "No": 259.101,
        "Lr": 266.12,
        "Rf": 267.122,
        "Db": 268.126,
        "Sg": 269.128,
        "Bh": 270.133,
        "Hs": 269.1336,
        "Mt": 277.154,
        "Ds": 282.166,
        "Rg": 282.169,
        "Cn": 286.179,
        "Nh": 286.182,
        "Fl": 290.192,
        "Mc": 290.196,
        "Lv": 293.205,
        "Ts": 294.211,
        "Og": 295.216,
    }

    try:
        alpha_symbol = "".join([s for s in str(symbol) if s.isalpha()])
        value = ATOMIC_MASS[str(alpha_symbol[:2]).title()]
    except (KeyError, TypeError):
        # print(f"\nWarning! symbol '{symbol}' is not in the Periodic Table\n")
        value = 0

    return value


def from_atomic_number_to_symbols(value):
    """
    Returns the symbol for an atomic number (Z) in the periodic table
    (0 <= Z <= 118)

    Examples
    --------
    >>> print(from_number_to_symbols(20))
    'Ca'
    >>> print(from_number_to_symbols(50))
    'Sn'


    Parameters
    ----------
    value : str or int
        atomic number (0, 118]

    Returns
    -------
    symbol : str
        atomic symbol from the periodic table if 'value' is not (0, 118]
        it returns the same value
    """

    ATOMS_Z_SYM = {
        1: "H",
        2: "He",
        3: "Li",
        4: "Be",
        5: "B",
        6: "C",
        7: "N",
        8: "O",
        9: "F",
        10: "Ne",
        11: "Na",
        12: "Mg",
        13: "Al",
        14: "Si",
        15: "P",
        16: "S",
        17: "Cl",
        18: "Ar",
        19: "K",
        20: "Ca",
        21: "Sc",
        22: "Ti",
        23: "V",
        24: "Cr",
        25: "Mn",
        26: "Fe",
        27: "Co",
        28: "Ni",
        29: "Cu",
        30: "Zn",
        31: "Ga",
        32: "Ge",
        33: "As",
        34: "Se",
        35: "Br",
        36: "Kr",
        37: "Rb",
        38: "Sr",
        39: "Y",
        40: "Zr",
        41: "Nb",
        42: "Mo",
        43: "Tc",
        44: "Ru",
        45: "Rh",
        46: "Pd",
        47: "Ag",
        48: "Cd",
        49: "In",
        50: "Sn",
        51: "Sb",
        52: "Te",
        53: "I",
        54: "Xe",
        55: "Cs",
        56: "Ba",
        57: "La",
        58: "Ce",
        59: "Pr",
        60: "Nd",
        61: "Pm",
        62: "Sm",
        63: "Eu",
        64: "Gd",
        65: "Tb",
        66: "Dy",
        67: "Ho",
        68: "Er",
        69: "Tm",
        70: "Yb",
        71: "Lu",
        72: "Hf",
        73: "Ta",
        74: "W",
        75: "Re",
        76: "Os",
        77: "Ir",
        78: "Pt",
        79: "Au",
        80: "Hg",
        81: "Tl",
        82: "Pb",
        83: "Bi",
        84: "Po",
        85: "At",
        86: "Rn",
        87: "Fr",
        88: "Ra",
        89: "Ac",
        90: "Th",
        91: "Pa",
        92: "U",
        93: "Np",
        94: "Pu",
        95: "Am",
        96: "Cm",
        97: "Bk",
        98: "Cf",
        99: "Es",
        100: "Fm",
        101: "Md",
        102: "No",
        103: "Lr",
        104: "Rf",
        105: "Db",
        106: "Sg",
        107: "Bh",
        108: "Hs",
        109: "Mt",
        110: "Ds",
        111: "Rg",
        112: "Cn",
        113: "Nh",
        114: "Fl",
        115: "Mc",
        116: "Lv",
        117: "Ts",
        118: "Og",
    }

    try:
        symbol = ATOMS_Z_SYM[int(value)]
    except (KeyError, ValueError):
        print("\nWarning! value is not an 'int' in (0, 118]\n")
        symbol = value

    return symbol


def from_symbols_to_atomic_number(symbol):
    """
    Returns atomic number (Z) in the periodic table (0 <= Z <= 118)
    from any valid symbol

    Example
    -------
    >>> print(from_symbols_to_atomic_number("Ca"))
    20
    >>> print(from_symbols_to_atomic_number("sn"))
    50

    Parameters
    ----------
    symbol : str
        atomic symbol from the periodic table if 'value' is not (0, 118]
        it returns the same value

    Returns
    -------
    value : str or int
        atomic number (0, 118]

    """

    ATOMS_SYM_Z = {
        "H": 1,
        "He": 2,
        "Li": 3,
        "Be": 4,
        "B": 5,
        "C": 6,
        "N": 7,
        "O": 8,
        "F": 9,
        "Ne": 10,
        "Na": 11,
        "Mg": 12,
        "Al": 13,
        "Si": 14,
        "P": 15,
        "S": 16,
        "Cl": 17,
        "Ar": 18,
        "K": 19,
        "Ca": 20,
        "Sc": 21,
        "Ti": 22,
        "V": 23,
        "Cr": 24,
        "Mn": 25,
        "Fe": 26,
        "Co": 27,
        "Ni": 28,
        "Cu": 29,
        "Zn": 30,
        "Ga": 31,
        "Ge": 32,
        "As": 33,
        "Se": 34,
        "Br": 35,
        "Kr": 36,
        "Rb": 37,
        "Sr": 38,
        "Y": 39,
        "Zr": 40,
        "Nb": 41,
        "Mo": 42,
        "Tc": 43,
        "Ru": 44,
        "Rh": 45,
        "Pd": 46,
        "Ag": 47,
        "Cd": 48,
        "In": 49,
        "Sn": 50,
        "Sb": 51,
        "Te": 52,
        "I": 53,
        "Xe": 54,
        "Cs": 55,
        "Ba": 56,
        "La": 57,
        "Ce": 58,
        "Pr": 59,
        "Nd": 60,
        "Pm": 61,
        "Sm": 62,
        "Eu": 63,
        "Gd": 64,
        "Tb": 65,
        "Dy": 66,
        "Ho": 67,
        "Er": 68,
        "Tm": 66,
        "Ir": 77,
        "Pt": 78,
        "Au": 79,
        "Hg": 80,
        "Tl": 81,
        "Pb": 82,
        "Bi": 83,
        "Po": 84,
        "At": 85,
        "Rn": 86,
        "Fr": 87,
        "Ra": 88,
        "Ac": 89,
        "Th": 90,
        "Pa": 91,
        "U": 92,
        "Np": 93,
        "Pu": 94,
        "Am": 95,
        "Cm": 96,
        "Bk": 97,
        "Cf": 98,
        "Es": 99,
        "Fm": 100,
        "Md": 101,
        "No": 102,
        "Lr": 103,
        "Rf": 104,
        "Db": 105,
        "Sg": 106,
        "Bh": 107,
        "Hs": 108,
        "Mt": 109,
        "Ds": 110,
        "Rg": 111,
        "Cn": 112,
        "Nh": 113,
        "Fl": 114,
        "Mc": 115,
        "Lv": 116,
        "Ts": 117,
        "Og": 118,
    }

    try:
        value = ATOMS_SYM_Z[str(symbol).title()]
    except (KeyError, TypeError):
        print("\nWarning! symbol is not in according to the periodic table\n")
        value = symbol

    return value
