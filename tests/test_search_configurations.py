import scipy

import pytest
from amcess.base_molecule import Molecule, Cluster
from amcess.search_configuration import SearchConfig


@pytest.mark.parametrize(
    "molecule1, molecule2, expected_coordinates",
    [
        (
            [("H", 0, 0, 0), ("F", 0.917, 0, 0)],
            [("H", 0, 0, 0), ("H", 0.74, 0, 0)],
            [
                ("H", 0, 0, 0),
                ("F", 0.917, 0, 0),
                ("H", 0, 0, 0),
                ("H", 0.74, 0, 0),
            ],
        ),
    ],
)
def test_cluster_object_coordinates_into_search_conf(
    molecule1, molecule2, expected_coordinates
):
    """[summary]
    Test: Passed of object from Cluser to SearchConfig
    Args:
        molecule1 ([dict]): dictionary with symbols and coordinates
                            of the first molecule
        molecule2 ([dict]): dictionary with symbols and coordinates
                            of the second molecule
        expected_coordinates ([list]): expected coordinates plus symbols
    """
    search_config = SearchConfig(Cluster(molecule1, molecule2))

    assert search_config._system_object.atoms == expected_coordinates


@pytest.mark.parametrize(
    "molecule1, molecule2, expected_coordinates",
    [
        (
            [("H", 0, 0, 0), ("F", 0.917, 0, 0)],
            [("H", 0, 0, 0), ("H", 0.74, 0, 0)],
            ["F", "H"],
        ),
    ],
)
def test_cluster_object_elements_into_search_conf(
    molecule1, molecule2, expected_coordinates
):
    """[summary]
    Test: Passed of object from Cluser to SearchConfig
    Args:
        molecule1 ([dict]): dictionary with symbols and coordinates
                            of the first molecule
        molecule2 ([dict]): dictionary with symbols and coordinates
                            of the second molecule
        expected_coordinates ([list]): expected coordinates plus symbols
    """
    search_config = SearchConfig(Cluster(molecule1, molecule2))
    # it is necessary to sort to avoid errors because the elements attr
    # can return in different order the symbols
    assert (
        sorted(search_config._system_object.elements) == expected_coordinates
    )


@pytest.mark.parametrize(
    "molecule1, molecule2, expected_bases",
    [
        (
            [("H", 0, 0, 0), ("F", 0.917, 0, 0)],
            [("H", 0, 0, 0), ("H", 0.74, 0, 0)],
            "sto-3g",
        ),
    ],
)
def test_search_conf_basis_default(molecule1, molecule2, expected_bases):
    """[summary]
    Test: Basis set by default (sto-3g)
    Args:
        molecule1 ([dict]): dictionary with symbols and coordinates
                            of the first molecule
        molecule2 ([dict]): dictionary with symbols and coordinates
                            of the second molecule
        expected_basis ([string]): basis set
    """
    search_config = SearchConfig(Cluster(molecule1, molecule2))
    assert search_config._bases_set == expected_bases


@pytest.mark.parametrize(
    "molecule1, molecule2, tolerance_radius, expecteted_contour_value",
    [
        (
            [("H", 0.0, 0.0, 0.0), ("H", 0.78, 0.0, 0.0)],
            [("H", 0.0, 0.0, 1.0), ("H", 0.78, 0.0, 1.0)],
            1.0,
            ((0.7407005820230608, 0.0, 0.5), 1.8936651230798374),
        ),
    ],
)
def test_spherical_contour(
    molecule1, molecule2, tolerance_radius, expecteted_contour_value
):
    """
    Test spherical contour cluster method in Cluster clase
    """
    obj = SearchConfig(Cluster(molecule1, molecule2))
    obj.spherical_contour_cluster(tolerance_radius)

    assert obj.sphere_center, obj.sphere_radius == expecteted_contour_value


def test_SC_init_first_value_error():
    """
    Test first TypeError into __init__ method of SearchConfig class
    when cluster object is None
    """
    with pytest.raises(TypeError) as e:
        SearchConfig()
    assert (
        str(e.value) == "AttributeError system_object isn't difinite\n"
        "It's NoneType"
    )


@pytest.mark.parametrize(
    "molecule1",
    [
        ([("H", 0.0, 0.0, 0.0), ("H", 0.78, 0.0, 0.0)]),
    ],
)
def test_SC_init_second_value_error(molecule1):
    """
    Test second TypeError into __init__ method of SearchConfig class
    when cluster object is not object type
    """
    with pytest.raises(TypeError) as e:
        SearchConfig(1.0)
    assert (
        str(e.value) == "AttributeError system_object isn't difinite\n"
        "as an object Cluster"
    )
    with pytest.raises(TypeError) as e:
        SearchConfig(1)
    assert (
        str(e.value) == "AttributeError system_object isn't difinite\n"
        "as an object Cluster"
    )
    with pytest.raises(TypeError) as e:
        SearchConfig(1.0)
    assert (
        str(e.value) == "AttributeError system_object isn't difinite\n"
        "as an object Cluster"
    )
    with pytest.raises(TypeError) as e:
        SearchConfig([1.0])
    assert (
        str(e.value) == "AttributeError system_object isn't difinite\n"
        "as an object Cluster"
    )
    with pytest.raises(TypeError) as e:
        SearchConfig(Molecule(molecule1))
    assert (
        str(e.value) == "AttributeError system_object isn't difinite\n"
        "as an object Cluster"
    )


@pytest.mark.parametrize(
    "molecule1, molecule2, expected_bounds",
    [
        (
            [("H", 0.0, 0.0, 0.0), ("H", 0.78, 0.0, 0.0)],
            [("H", 0.0, 0.0, 1.0), ("H", 0.78, 0.0, 1.0)],
            [
                (-1.6341135544995076, 1.6341135544995076),
                (-1.6341135544995076, 1.6341135544995076),
                (-1.6341135544995076, 1.6341135544995076),
                (-1.6341135544995076, 1.6341135544995076),
                (-1.6341135544995076, 1.6341135544995076),
                (-1.6341135544995076, 1.6341135544995076),
                (0, scipy.pi),
                (0, scipy.pi),
                (0, scipy.pi),
                (0, scipy.pi),
                (0, scipy.pi),
                (0, scipy.pi),
            ],
        ),
    ],
)
def test_SC_bounds_grep(molecule1, molecule2, expected_bounds):
    """
    Test @property associated with bounds variable
    """
    assert (
        SearchConfig(Cluster(molecule1, molecule2)).bounds == expected_bounds
    )


@pytest.mark.parametrize(
    "molecule1, molecule2, new_bounds",
    [
        (
            [("H", 0.0, 0.0, 0.0), ("H", 0.78, 0.0, 0.0)],
            [("H", 0.0, 0.0, 1.0), ("H", 0.78, 0.0, 1.0)],
            [
                (-1.6341135544995076, 1.6341135544995076),
                (0, scipy.pi),
                (0, scipy.pi),
                (0, scipy.pi),
                (0, scipy.pi),
                (0, scipy.pi),
                (0, scipy.pi),
            ],
        ),
    ],
)
def test_SC_bounds_set_TypeError(molecule1, molecule2, new_bounds):
    """
    Test @bounds.setter, ValueError associated with len(bounds)
    """
    with pytest.raises(ValueError) as e:
        SearchConfig(Cluster(molecule1, molecule2)).bounds = new_bounds
    assert (
        str(e.value) == "\n\nArray dimensions insufficient: "
        f"\ndimensions of old bounds: '{12}'\n"
        f"\ndimensions of new bounds: '{len(new_bounds)}'\n"
    )


@pytest.mark.parametrize(
    "molecule1, molecule2, expected_name",
    [
        (
            [("H", 0.0, 0.0, 0.0), ("H", 0.78, 0.0, 0.0)],
            [("H", 0.0, 0.0, 1.0), ("H", 0.78, 0.0, 1.0)],
            "configurations.xyz",
        ),
    ],
)
def test_SC_output_name_grep(molecule1, molecule2, expected_name):
    """
    Test @property associated with output_name variable
    """
    assert (
        SearchConfig(Cluster(molecule1, molecule2)).output_name
        == expected_name
    )


@pytest.mark.parametrize(
    "molecule1, molecule2, new_name",
    [
        (
            [("H", 0.0, 0.0, 0.0), ("H", 0.78, 0.0, 0.0)],
            [("H", 0.0, 0.0, 1.0), ("H", 0.78, 0.0, 1.0)],
            "output.xyz",
        ),
    ],
)
def test_SC_output_name_set(molecule1, molecule2, new_name):
    """
    Test @output_name.setter associated with output_name variable
    """
    obj_sc = SearchConfig(Cluster(molecule1, molecule2))
    obj_sc.output_name = new_name
    assert obj_sc.output_name == new_name


@pytest.mark.parametrize(
    "molecule1, molecule2",
    [
        (
            [("H", 0.0, 0.0, 0.0), ("H", 0.78, 0.0, 0.0)],
            [("H", 0.0, 0.0, 1.0), ("H", 0.78, 0.0, 1.0)],
        ),
    ],
)
def test_SC_output_name_set(molecule1, molecule2):
    """
    Test TypeError into @output_name.setter
    """
    with pytest.raises(TypeError) as e:
        obj_sc = SearchConfig(Cluster(molecule1, molecule2))
        obj_sc.output_name = 1
    assert (
        str(e.value) == "\n\nThe new name to output is not a string"
        f"\nplease, check: '{type(1)}'\n"
    )
    with pytest.raises(TypeError) as e:
        obj_sc = SearchConfig(Cluster(molecule1, molecule2))
        obj_sc.output_name = 1.0
    assert (
        str(e.value) == "\n\nThe new name to output is not a string"
        f"\nplease, check: '{type(1.0)}'\n"
    )
    with pytest.raises(TypeError) as e:
        obj_sc = SearchConfig(Cluster(molecule1, molecule2))
        obj_sc.output_name = (1,)
    assert (
        str(e.value) == "\n\nThe new name to output is not a string"
        f"\nplease, check: '{type((1,))}'\n"
    )


@pytest.mark.parametrize(
    "molecule1, molecule2, expected_search",
    [
        (
            [("H", 0.0, 0.0, 0.0), ("H", 0.78, 0.0, 0.0)],
            [("H", 0.0, 0.0, 1.0), ("H", 0.78, 0.0, 1.0)],
            1,
        ),
    ],
)
def test_SC_search_type_grep(molecule1, molecule2, expected_search):
    """
    Test @property associated with search methodology type
    """
    assert (
        SearchConfig(Cluster(molecule1, molecule2)).search_type
        == expected_search
    )


@pytest.mark.parametrize(
    "molecule1, molecule2",
    [
        (
            [("H", 0.0, 0.0, 0.0), ("H", 0.78, 0.0, 0.0)],
            [("H", 0.0, 0.0, 1.0), ("H", 0.78, 0.0, 1.0)],
        ),
    ],
)
def test_SC_search_type_TP_set(molecule1, molecule2):
    """
    Test Type Error @property associated with search methodology type
    """
    with pytest.raises(TypeError) as e:
        obj_sc = SearchConfig(Cluster(molecule1, molecule2))
        obj_sc.search_type = 1.0
    assert (
        str(e.value) == "\n\nThe new search methodology is not an integer"
        f"\nplease, check: '{type(1.0)}'\n"
    )
    with pytest.raises(TypeError) as e:
        obj_sc = SearchConfig(Cluster(molecule1, molecule2))
        obj_sc.search_type = (1.0,)
    assert (
        str(e.value) == "\n\nThe new search methodology is not an integer"
        f"\nplease, check: '{type((1.0,))}'\n"
    )
    with pytest.raises(TypeError) as e:
        obj_sc = SearchConfig(Cluster(molecule1, molecule2))
        obj_sc.search_type = "1.0"
    assert (
        str(e.value) == "\n\nThe new search methodology is not an integer"
        f"\nplease, check: '{type('1.0')}'\n"
    )


@pytest.mark.parametrize(
    "molecule1, molecule2, search_methodology",
    [
        (
            [("H", 0.0, 0.0, 0.0), ("H", 0.78, 0.0, 0.0)],
            [("H", 0.0, 0.0, 1.0), ("H", 0.78, 0.0, 1.0)],
            4,
        ),
    ],
)
def test_SC_search_type_VE_set(molecule1, molecule2, search_methodology):
    """
    Test ValueError @property associated with search methodology type
    when search metodology choose is higher than 3
    """
    with pytest.raises(TypeError) as e:
        SearchConfig(Cluster(molecule1, molecule2), search_methodology)
    assert (
        str(e.value)
        == "\n\nThe search methodology is associated with a integer \n"
        "1 -> Dual Annealing \n"
        "2 -> SHGO \n"
        "3 -> Bayessiana \n"
        f"\nplease, check: '{type(search_methodology)}'\n"
    )


@pytest.mark.parametrize(
    "molecule1, molecule2, new_search_methodology",
    [
        (
            [("H", 0.0, 0.0, 0.0), ("H", 0.78, 0.0, 0.0)],
            [("H", 0.0, 0.0, 1.0), ("H", 0.78, 0.0, 1.0)],
            2,
        ),
    ],
)
def test_SC_search_type_VE_set(molecule1, molecule2, new_search_methodology):
    """
    Test @search_type.setter when search metodology choose is lower than 3
    or equal to 3
    """
    obj_sc = SearchConfig(Cluster(molecule1, molecule2))
    obj_sc.search_type = new_search_methodology
    assert obj_sc.search_type == new_search_methodology


@pytest.mark.parametrize(
    "molecule1, molecule2, new_bases_set",
    [
        (
            [("H", 0.0, 0.0, 0.0), ("H", 0.78, 0.0, 0.0)],
            [("H", 0.0, 0.0, 1.0), ("H", 0.78, 0.0, 1.0)],
            "3-21g",
        ),
    ],
)
def test_SC_bases_set_grep(molecule1, molecule2, new_bases_set):
    """
    Test @bases_set.setter for new bases set
    """
    obj_sc = SearchConfig(Cluster(molecule1, molecule2))
    obj_sc.bases_set = new_bases_set
    assert obj_sc.bases_set == new_bases_set


@pytest.mark.parametrize(
    "molecule1, molecule2",
    [
        (
            [("H", 0.0, 0.0, 0.0), ("H", 0.78, 0.0, 0.0)],
            [("H", 0.0, 0.0, 1.0), ("H", 0.78, 0.0, 1.0)],
        ),
    ],
)
def test_SC_bases_set_set(molecule1, molecule2):
    """
    Test Type Error @property associated with search methodology type
    """
    with pytest.raises(TypeError) as e:
        obj_sc = SearchConfig(Cluster(molecule1, molecule2))
        obj_sc.bases_set = 1.0
    assert (
        str(e.value) == "\n\nThe new name to output is not a string"
        f"\nplease, check: '{type(1.0)}'\n"
    )
    with pytest.raises(TypeError) as e:
        obj_sc = SearchConfig(Cluster(molecule1, molecule2))
        obj_sc.bases_set = (1.0,)
    assert (
        str(e.value) == "\n\nThe new name to output is not a string"
        f"\nplease, check: '{type((1.0,))}'\n"
    )
    with pytest.raises(TypeError) as e:
        obj_sc = SearchConfig(Cluster(molecule1, molecule2))
        obj_sc.bases_set = 1
    assert (
        str(e.value) == "\n\nThe new name to output is not a string"
        f"\nplease, check: '{type(1)}'\n"
    )
