import os

import pytest
import scipy

from amcess.base_molecule import Cluster, Molecule
from amcess.electronic_energy import hf_pyscf
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
    "molecule1, molecule2, expected_basis",
    [
        (
            [("H", 0, 0, 0), ("F", 0.917, 0, 0)],
            [("H", 0, 0, 0), ("H", 0.74, 0, 0)],
            "sto-3g",
        ),
    ],
)
def test_search_conf_basis_default(molecule1, molecule2, expected_basis):
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
    assert search_config._basis_set == expected_basis


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
                (-2.6341135544995076, 1.6341135544995076),
                (-2.6341135544995076, 1.6341135544995076),
                (-2.6341135544995076, 1.6341135544995076),
                (-2.6341135544995076, 1.6341135544995076),
                (-2.6341135544995076, 1.6341135544995076),
                (-2.6341135544995076, 1.6341135544995076),
                (0, 3.16),
                (0, 3.16),
                (0, 3.16),
                (0, 3.16),
                (0, 3.16),
                (0, 3.16),
            ],
        ),
    ],
)
def test_SC_new_bounds_grep(molecule1, molecule2, new_bounds):
    """
    Test @property associated with bounds variable
    """
    obj_sc = SearchConfig(Cluster(molecule1, molecule2))
    obj_sc.bounds = new_bounds
    assert obj_sc.bounds == new_bounds


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
def test_SC_output_new_name_set(molecule1, molecule2, new_name):
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
def test_SC_output_name_TE_set(molecule1, molecule2):
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
def test_SC_search_type_TE(molecule1, molecule2, search_methodology):
    """
    Test ValueError associated with search methodology into __init__
    """
    with pytest.raises(ValueError) as e:
        SearchConfig(Cluster(molecule1, molecule2), search_methodology)
    assert (
        str(e.value) == "ValueError search_methodology is an integer\n"
        "between 1 and 3\n"
        f"\nplease, check: '{search_methodology}'\n"
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
def test_SC_search_type_TE_set(molecule1, molecule2, search_methodology):
    """
    Test ValueError @property associated with search methodology type
    when search metodology choose is higher than 3
    """
    with pytest.raises(ValueError) as e:
        obj_sc = SearchConfig(Cluster(molecule1, molecule2))
        obj_sc.search_type = search_methodology
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
def test_SC_new_search_type_set(molecule1, molecule2, new_search_methodology):
    """
    Test @search_type.setter when search metodology choose is lower than 3
    or equal to 3
    """
    obj_sc = SearchConfig(Cluster(molecule1, molecule2))
    obj_sc.search_type = new_search_methodology
    assert obj_sc.search_type == new_search_methodology


@pytest.mark.parametrize(
    "molecule1, molecule2, new_basis_set",
    [
        (
            [("H", 0.0, 0.0, 0.0), ("H", 0.78, 0.0, 0.0)],
            [("H", 0.0, 0.0, 1.0), ("H", 0.78, 0.0, 1.0)],
            "3-21g",
        ),
    ],
)
def test_SC_basis_set_grep(molecule1, molecule2, new_basis_set):
    """
    Test @basis_set.setter for new basis set
    """
    obj_sc = SearchConfig(Cluster(molecule1, molecule2))
    obj_sc.basis_set = new_basis_set
    assert obj_sc.basis_set == new_basis_set


@pytest.mark.parametrize(
    "molecule1, molecule2",
    [
        (
            [("H", 0.0, 0.0, 0.0), ("H", 0.78, 0.0, 0.0)],
            [("H", 0.0, 0.0, 1.0), ("H", 0.78, 0.0, 1.0)],
        ),
    ],
)
def test_SC_basis_set_set(molecule1, molecule2):
    """
    Test Type Error @property associated with basis set
    """
    with pytest.raises(TypeError) as e:
        obj_sc = SearchConfig(Cluster(molecule1, molecule2))
        obj_sc.basis_set = 1.0
    assert (
        str(e.value) == "\n\nThe new name to output is not a string"
        f"\nplease, check: '{type(1.0)}'\n"
    )
    with pytest.raises(TypeError) as e:
        obj_sc = SearchConfig(Cluster(molecule1, molecule2))
        obj_sc.basis_set = (1.0,)
    assert (
        str(e.value) == "\n\nThe new name to output is not a string"
        f"\nplease, check: '{type((1.0,))}'\n"
    )
    with pytest.raises(TypeError) as e:
        obj_sc = SearchConfig(Cluster(molecule1, molecule2))
        obj_sc.basis_set = 1
    assert (
        str(e.value) == "\n\nThe new name to output is not a string"
        f"\nplease, check: '{type(1)}'\n"
    )


@pytest.mark.parametrize(
    "molecule1, molecule2, new_sphere_center",
    [
        (
            [("H", 0.0, 0.0, 0.0), ("H", 0.78, 0.0, 0.0)],
            [("H", 0.0, 0.0, 1.0), ("H", 0.78, 0.0, 1.0)],
            (2.0, 0.0, 1.0),
        ),
    ],
)
def test_SC_new_sphere_center_set(molecule1, molecule2, new_sphere_center):
    """
    Test new sphere center @sphere_center.setter
    """
    obj_sc = SearchConfig(Cluster(molecule1, molecule2))
    obj_sc.sphere_center = new_sphere_center
    assert obj_sc.sphere_center == new_sphere_center


@pytest.mark.parametrize(
    "molecule1, molecule2",
    [
        (
            [("H", 0.0, 0.0, 0.0), ("H", 0.78, 0.0, 0.0)],
            [("H", 0.0, 0.0, 1.0), ("H", 0.78, 0.0, 1.0)],
        ),
    ],
)
def test_SC_sphere_center_TE_set(molecule1, molecule2):
    """
    Test TypeError @sphere_center.setter
    """
    with pytest.raises(TypeError) as e:
        obj_sc = SearchConfig(Cluster(molecule1, molecule2))
        obj_sc.sphere_center = 1.0
    assert (
        str(e.value)
        == "\n\nThe Sphere center must be a tuple with three elements: "
        "(float, float, float)"
        f"\nplease, check: '{type(1.0)}'\n"
    )
    with pytest.raises(TypeError) as e:
        obj_sc = SearchConfig(Cluster(molecule1, molecule2))
        obj_sc.sphere_center = 1
    assert (
        str(e.value)
        == "\n\nThe Sphere center must be a tuple with three elements: "
        "(float, float, float)"
        f"\nplease, check: '{type(1)}'\n"
    )
    with pytest.raises(TypeError) as e:
        obj_sc = SearchConfig(Cluster(molecule1, molecule2))
        obj_sc.sphere_center = "1.0"
    assert (
        str(e.value)
        == "\n\nThe Sphere center must be a tuple with three elements: "
        "(float, float, float)"
        f"\nplease, check: '{type('1.0')}'\n"
    )
    with pytest.raises(TypeError) as e:
        obj_sc = SearchConfig(Cluster(molecule1, molecule2))
        obj_sc.sphere_center = [1.0]
    assert (
        str(e.value)
        == "\n\nThe Sphere center must be a tuple with three elements: "
        "(float, float, float)"
        f"\nplease, check: '{type([1.0])}'\n"
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
def test_SC_sphere_center_VE_set(molecule1, molecule2):
    """
    Test ValueError @sphere_center.setter is a tuple with amount
    elements different to 3
    """
    with pytest.raises(ValueError) as e:
        obj_sc = SearchConfig(Cluster(molecule1, molecule2))
        obj_sc.sphere_center = (1.0,)
    assert (
        str(e.value)
        == "\n\nThe Sphere center must be a tuple with three elements: "
        "(float, float, float)"
        f"\nplease, check: '{(1.0,)}'\n"
    )
    with pytest.raises(ValueError) as e:
        obj_sc = SearchConfig(Cluster(molecule1, molecule2))
        obj_sc.sphere_center = (
            1.0,
            1.0,
        )
    assert (
        str(e.value)
        == "\n\nThe Sphere center must be a tuple with three elements: "
        "(float, float, float)"
        f"\nplease, check: '{(1.0,1.0,)}'\n"
    )
    with pytest.raises(ValueError) as e:
        obj_sc = SearchConfig(Cluster(molecule1, molecule2))
        obj_sc.sphere_center = (
            1.0,
            1.0,
            1.0,
            1.0,
        )
    assert (
        str(e.value)
        == "\n\nThe Sphere center must be a tuple with three elements: "
        "(float, float, float)"
        f"\nplease, check: '{(1.0,1.0,1.0,1.0,)}'\n"
    )


@pytest.mark.parametrize(
    "molecule1, molecule2, new_radius",
    [
        (
            [("H", 0.0, 0.0, 0.0), ("H", 0.78, 0.0, 0.0)],
            [("H", 0.0, 0.0, 1.0), ("H", 0.78, 0.0, 1.0)],
            2.5,
        ),
    ],
)
def test_SC_new_sphere_radius_set(molecule1, molecule2, new_radius):
    """
    Test new spehre radius  @sphere_radius.setter
    """
    obj_sc = SearchConfig(Cluster(molecule1, molecule2))
    obj_sc.sphere_radius = new_radius
    assert obj_sc.sphere_radius == new_radius


@pytest.mark.parametrize(
    "molecule1, molecule2",
    [
        (
            [("H", 0.0, 0.0, 0.0), ("H", 0.78, 0.0, 0.0)],
            [("H", 0.0, 0.0, 1.0), ("H", 0.78, 0.0, 1.0)],
        ),
    ],
)
def test_SC_sphere_radius_TE_set(molecule1, molecule2):
    """
    Test TypeError @sphere_radius.setter
    """
    with pytest.raises(TypeError) as e:
        obj_sc = SearchConfig(Cluster(molecule1, molecule2))
        obj_sc.sphere_radius = "1.0"
    assert (
        str(e.value) == "\n\nThe Sphere  Radius must be a float or int"
        f"\nplease, check: '{type('1.0')}'\n"
    )
    with pytest.raises(TypeError) as e:
        obj_sc = SearchConfig(Cluster(molecule1, molecule2))
        obj_sc.sphere_radius = [1.0]
    assert (
        str(e.value) == "\n\nThe Sphere  Radius must be a float or int"
        f"\nplease, check: '{type([1.0])}'\n"
    )
    with pytest.raises(TypeError) as e:
        obj_sc = SearchConfig(Cluster(molecule1, molecule2))
        obj_sc.sphere_radius = (1.0,)
    assert (
        str(e.value) == "\n\nThe Sphere  Radius must be a float or int"
        f"\nplease, check: '{type((1.0,))}'\n"
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
def test_SC_sphere_radius_VE_set(molecule1, molecule2):
    """
    Test ValueError @sphere_radius.setter
    """
    with pytest.raises(ValueError) as e:
        obj_sc = SearchConfig(Cluster(molecule1, molecule2))
        obj_sc.sphere_radius = 0.8
    assert (
        str(e.value) == "\n\nThe Sphere  Radius must be larger than 1 Angstrom"
        f"\nplease, check: '{0.8}'\n"
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
def test_SC_cost_function_TE_set(molecule1, molecule2):
    """
    Test TypeError @sphere_radius.setter
    """
    with pytest.raises(TypeError) as e:
        obj_sc = SearchConfig(Cluster(molecule1, molecule2))
        obj_sc.cost_function_number = 1.0
    assert (
        str(e.value) == "\n\nThe new cost function is not a integer"
        f"\nplease, check: '{type(1.0)}'\n"
    )
    with pytest.raises(TypeError) as e:
        obj_sc = SearchConfig(Cluster(molecule1, molecule2))
        obj_sc.cost_function_number = [1.0]
    assert (
        str(e.value) == "\n\nThe new cost function is not a integer"
        f"\nplease, check: '{type([1.0])}'\n"
    )
    with pytest.raises(TypeError) as e:
        obj_sc = SearchConfig(Cluster(molecule1, molecule2))
        obj_sc.cost_function_number = (1.0,)
    assert (
        str(e.value) == "\n\nThe new cost function is not a integer"
        f"\nplease, check: '{type((1.0,))}'\n"
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
def test_SC_cost_function_VE_set(molecule1, molecule2):
    """
    Test ValueError @sphere_radius.setter
    """
    with pytest.raises(ValueError) as e:
        obj_sc = SearchConfig(Cluster(molecule1, molecule2))
        obj_sc.cost_function_number = 2
    assert (
        str(e.value) == "\n\nThe new cost function is not implemeted "
        "\n 1 -> Hartree Fock into pyscf"
        f"\nplease, check: '{2}'\n"
    )


@pytest.mark.parametrize(
    "molecule1, molecule2, program_ee",
    [
        (
            [("H", 0.0, 0.0, 0.0), ("H", 0.78, 0.0, 0.0)],
            [("H", 0.0, 0.0, 1.0), ("H", 0.78, 0.0, 1.0)],
            2,
        ),
    ],
)
def test_SC_cost_function_VE_init(molecule1, molecule2, program_ee):
    """
    Test ValueError associated cost function into __init__
    """
    with pytest.raises(ValueError) as e:
        SearchConfig(
            Cluster(molecule1, molecule2),
            program_electronic_structure=program_ee,
        )
    assert (
        str(e.value) == "ValueError, only implemeted an option for\n"
        "electronic structure\n"
        f"\nplease, check: '{program_ee}'\n"
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
def test_SC_cost_function_grep(molecule1, molecule2):
    """
    Test default cost_function_ee @property
    """
    obj_sc = SearchConfig(Cluster(molecule1, molecule2))
    assert obj_sc.cost_function_ee == hf_pyscf


@pytest.mark.parametrize(
    "molecule1, molecule2",
    [
        (
            [("H", 0.0, 0.0, 0.0), ("H", 0.78, 0.0, 0.0)],
            [("H", 0.0, 0.0, 1.0), ("H", 0.78, 0.0, 1.0)],
        ),
    ],
)
def test_SC_cost_function_number_grep(molecule1, molecule2):
    """
    Test default cost_function_number @property
    """
    obj_sc = SearchConfig(Cluster(molecule1, molecule2))
    assert obj_sc.cost_function_number == 1


@pytest.mark.parametrize(
    "molecule1, molecule2, program_ee",
    [
        (
            [("H", 0.0, 0.0, 0.0), ("H", 0.78, 0.0, 0.0)],
            [("H", 0.0, 0.0, 1.0), ("H", 0.78, 0.0, 1.0)],
            1,
        ),
    ],
)
def test_SC_new_cost_function(molecule1, molecule2, program_ee):
    """
    Test ValueError @sphere_radius.setter
    """
    obj_sc = SearchConfig(Cluster(molecule1, molecule2))
    obj_sc.cost_function_number = program_ee
    assert obj_sc._func == hf_pyscf


@pytest.mark.parametrize(
    "molecule1, molecule2, default_tolerance_radius",
    [
        (
            [("H", 0.0, 0.0, 0.0), ("H", 0.78, 0.0, 0.0)],
            [("H", 0.0, 0.0, 1.0), ("H", 0.78, 0.0, 1.0)],
            1.0,
        ),
    ],
)
def test_SC_tolerance_radius_grep(
    molecule1, molecule2, default_tolerance_radius
):
    """
    Test ask tolerance radius
    """
    assert (
        SearchConfig(Cluster(molecule1, molecule2)).radius_contour
        == default_tolerance_radius
    )


@pytest.mark.parametrize(
    "molecule1, molecule2, new_tolerance_radius",
    [
        (
            [("H", 0.0, 0.0, 0.0), ("H", 0.78, 0.0, 0.0)],
            [("H", 0.0, 0.0, 1.0), ("H", 0.78, 0.0, 1.0)],
            2.0,
        ),
    ],
)
def test_SC_new_tolerance_radius_set(
    molecule1, molecule2, new_tolerance_radius
):
    """
    Test TypeError @radius_contour.setter
    """
    obj_sc = SearchConfig(Cluster(molecule1, molecule2))
    obj_sc.radius_contour = new_tolerance_radius
    assert obj_sc.radius_contour == new_tolerance_radius


@pytest.mark.parametrize(
    "molecule1, molecule2",
    [
        (
            [("H", 0.0, 0.0, 0.0), ("H", 0.78, 0.0, 0.0)],
            [("H", 0.0, 0.0, 1.0), ("H", 0.78, 0.0, 1.0)],
        ),
    ],
)
def test_SC_tolerance_radius_TE_set(molecule1, molecule2):
    """
    Test TypeError @radius_contour.setter
    """
    with pytest.raises(TypeError) as e:
        obj_sc = SearchConfig(Cluster(molecule1, molecule2))
        obj_sc.radius_contour = 1
    assert (
        str(e.value) == "\n\nThe new tolerance radius is not a float"
        f"\nplease, check: '{type(1)}'\n"
    )
    with pytest.raises(TypeError) as e:
        obj_sc = SearchConfig(Cluster(molecule1, molecule2))
        obj_sc.radius_contour = "1.0"
    assert (
        str(e.value) == "\n\nThe new tolerance radius is not a float"
        f"\nplease, check: '{type('1.0')}'\n"
    )
    with pytest.raises(TypeError) as e:
        obj_sc = SearchConfig(Cluster(molecule1, molecule2))
        obj_sc.radius_contour = [1.0]
    assert (
        str(e.value) == "\n\nThe new tolerance radius is not a float"
        f"\nplease, check: '{type([1.0])}'\n"
    )
    with pytest.raises(TypeError) as e:
        obj_sc = SearchConfig(Cluster(molecule1, molecule2))
        obj_sc.radius_contour = (1.0,)
    assert (
        str(e.value) == "\n\nThe new tolerance radius is not a float"
        f"\nplease, check: '{type((1.0,))}'\n"
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
def test_SC_run_method(molecule1, molecule2):
    """
    Test SC.run method
    """
    SearchConfig(Cluster(molecule1, molecule2)).run(NT=1, mxcycle=1)
    with open("configurations.xyz", "r") as f:
        readl = f.readline()
    os.remove("configurations.xyz")
    assert readl == "4\n"


@pytest.mark.parametrize(
    "molecule1, molecule2",
    [
        (
            [("H", 0.0, 0.0, 0.0), ("H", 0.78, 0.0, 0.0)],
            [("H", 0.0, 0.0, 1.0), ("H", 0.78, 0.0, 1.0)],
        ),
    ],
)
def test_SC_da_method(molecule1, molecule2):
    """
    Test SC.da method
    """
    SearchConfig(Cluster(molecule1, molecule2)).da(NT=1, mxcycle=1)
    with open("configurations.xyz", "r") as f:
        readl = f.readline()
    os.remove("configurations.xyz")
    assert readl == "4\n"
