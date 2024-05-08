import math
import os

import pytest

from amcess.ascec import Ascec
from amcess.molecule import Molecule
from amcess.cluster import Cluster
from amcess.search_engine import SearchConfig

METHODS = {
    "ASCEC": 0,
    "dual_annealing": 1,
}
available = list(METHODS.keys())


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
def test_CoordibatesClusterObject(molecule1, molecule2, expected_coordinates):
    """[summary]
    Test: Passed of object from Cluster to SearchConfig
    Args:
        molecule1 ([dict]): dictionary with symbols and coordinates
                            of the first molecule
        molecule2 ([dict]): dictionary with symbols and coordinates
                            of the second molecule
        expected_coordinates ([list]): expected coordinates plus symbols
    """
    search_config = SearchConfig(
        Cluster(molecule1, molecule2, sphere_center=(0.0, 0.0, 0.0), sphere_radius=3.0)
    )

    assert search_config.GetSystemObject().GetMolList() == expected_coordinates


@pytest.mark.parametrize(
    "molecule1, molecule2, expected_coordinates",
    [
        (
            [("H", 0, 0, 0), ("F", 0.917, 0, 0)],
            [("H", 0, 0, 0), ("H", 0.74, 0, 0)],
            ["H", "F", "H", "H"],
        ),
    ],
)
def test_SymbolsClusterObject(
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

    assert search_config.GetSystemObject().GetAtomicSymbols() == expected_coordinates


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
def test_BasisSetSearchConf(molecule1, molecule2, expected_basis):
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

    assert search_config.GetBasisSet() == expected_basis


@pytest.mark.parametrize(
    "molecule1, molecule2",
    [
        (
            [("H", 0.0, 0.0, 0.0), ("H", 0.78, 0.0, 0.0)],
            [("H", 0.0, 0.0, 1.0), ("H", 0.78, 0.0, 1.0)],
        ),
    ],
)
def test_ClusterObjectErrorOne(molecule1, molecule2):
    """
    Test first TypeError into __init__ method of SearchConfig class
    when cluster object is None
    """
    with pytest.raises(TypeError) as e:
        obj_sc = SearchConfig(Cluster(molecule1, molecule2))
        obj_sc.SetSystemObject(None)

    assert str(e.value) == "System_object isn't difinite\n" "It's NoneType"


@pytest.mark.parametrize(
    "molecule1, molecule2",
    [
        (
            [("H", 0.0, 0.0, 0.0), ("H", 0.78, 0.0, 0.0)],
            [("H", 0.0, 0.0, 1.0), ("H", 0.78, 0.0, 1.0)],
        ),
    ],
)
def test_ClusterObjectErrorTwo(molecule1, molecule2):
    """
    Test second TypeError into __init__ method of SearchConfig class
    when cluster object is not object type
    """
    obj_sc = SearchConfig(Cluster(molecule1, molecule2))
    with pytest.raises(TypeError) as e:
        obj_sc.SetSystemObject(1.0)
    assert (
        str(e.value) == "System_object isn't difinite as an object Cluster\n"
        f"please, check:\n'{1.0}'"
    )
    with pytest.raises(TypeError) as e:
        obj_sc.SetSystemObject(1)
    assert (
        str(e.value) == "System_object isn't difinite as an object Cluster\n"
        f"please, check:\n'{1}'"
    )
    with pytest.raises(TypeError) as e:
        obj_sc.SetSystemObject((1.0,))
    assert (
        str(e.value) == "System_object isn't difinite as an object Cluster\n"
        f"please, check:\n'{(1.0,)}'"
    )
    with pytest.raises(TypeError) as e:
        obj_sc.SetSystemObject([1.0])
    assert (
        str(e.value) == "System_object isn't difinite as an object Cluster\n"
        f"please, check:\n'{[1.0]}'"
    )
    with pytest.raises(TypeError) as e:
        obj_sc.SetSystemObject(Molecule(molecule1))
    assert (
        str(e.value) == "System_object isn't difinite as an object Cluster\n"
        f"please, check:\n'{Molecule(molecule1)}'"
    )


@pytest.mark.parametrize(
    "molecule1, molecule2, radius, expected_bounds",
    [
        (
            [("H", 0.0, 0.0, 0.0), ("H", 0.78, 0.0, 0.0)],
            [("H", 0.0, 0.0, 1.0), ("H", 0.78, 0.0, 1.0)],
            2.07335921293852,
            [
                (-2.07335921293852, 2.07335921293852),
                (-2.07335921293852, 2.07335921293852),
                (-2.07335921293852, 2.07335921293852),
                (-180.0, 180.0),
                (-180.0, 180.0),
                (-180.0, 180.0),
            ],
        ),
    ],
)
def test_ClusterObjectGrepBounds(molecule1, molecule2, radius, expected_bounds):
    """
    Test bounds of Cluster object
    """
    search_config = SearchConfig(Cluster(molecule1, molecule2, sphere_radius=radius))
    assert (
        search_config.GetBounds()
        == expected_bounds
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
                (0, 3.16),
                (0, 3.16),
                (0, 3.16),
            ],
        ),
    ],
)
def test_ClusterObjectSetBounds(molecule1, molecule2, new_bounds):
    """
    Test SetBounds
    """
    obj_sc = SearchConfig(Cluster(molecule1, molecule2))
    obj_sc.SetBounds(new_bounds)
    assert obj_sc.GetBounds() == new_bounds


@pytest.mark.parametrize(
    "molecule1, molecule2, new_bounds",
    [
        (
            [("H", 0.0, 0.0, 0.0), ("H", 0.78, 0.0, 0.0)],
            [("H", 0.0, 0.0, 1.0), ("H", 0.78, 0.0, 1.0)],
            [
                (-1.6341135544995076, 1.6341135544995076),
                (0, math.pi),
                (0, math.pi),
                (0, math.pi),
                (0, math.pi),
                (0, math.pi),
                (0, math.pi),
            ],
        ),
    ],
)
def test_ClusterObjectSetBoundsTypeError(molecule1, molecule2, new_bounds):
    """
    Test ValueError associated with len(bounds)
    """
    with pytest.raises(ValueError) as e:
        SearchConfig(Cluster(molecule1, molecule2)).SetBounds(new_bounds)

    assert (
        str(e.value) == "\n\nArray dimensions insufficient: "
        f"\ndimensions of old bounds: '{6}'\n"
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
def test_GetOutPutName(molecule1, molecule2, expected_name):
    """
    Test get out_put name associated with output_name variable
    """
    search_config = SearchConfig(Cluster(molecule1, molecule2))

    assert search_config.GetOutputName() == expected_name


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
def test_SetOutPutName(molecule1, molecule2, new_name):
    """
    Test set output_name associated with output_name variable
    """
    obj_sc = SearchConfig(Cluster(molecule1, molecule2))
    obj_sc.SetOutputName(new_name)
    assert obj_sc.GetOutputName() == new_name


@pytest.mark.parametrize(
    "molecule1, molecule2",
    [
        (
            [("H", 0.0, 0.0, 0.0), ("H", 0.78, 0.0, 0.0)],
            [("H", 0.0, 0.0, 1.0), ("H", 0.78, 0.0, 1.0)],
        ),
    ],
)
def test_SetOutPutNameTypeError(molecule1, molecule2):
    """
    Test TypeError into SetOutputName
    """
    with pytest.raises(TypeError) as e:
        obj_sc = SearchConfig(Cluster(molecule1, molecule2))
        obj_sc.SetOutputName(1)
    assert (
        str(e.value) == "\n\nThe new name to output is not a string"
        f"\nplease, check: '{type(1)}'\n"
    )
    with pytest.raises(TypeError) as e:
        obj_sc = SearchConfig(Cluster(molecule1, molecule2))
        obj_sc.SetOutputName(1.0)
    assert (
        str(e.value) == "\n\nThe new name to output is not a string"
        f"\nplease, check: '{type(1.0)}'\n"
    )
    with pytest.raises(TypeError) as e:
        obj_sc = SearchConfig(Cluster(molecule1, molecule2))
        obj_sc.SetOutputName((1,))
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
            "ASCEC",
        ),
    ],
)
def test_GetSearchType(molecule1, molecule2, expected_search):
    """
    Test associated with search methodology type
    """
    search_config = SearchConfig(Cluster(molecule1, molecule2))
    assert search_config.GetSearchType() == expected_search


@pytest.mark.parametrize(
    "molecule1, molecule2",
    [
        (
            [("H", 0.0, 0.0, 0.0), ("H", 0.78, 0.0, 0.0)],
            [("H", 0.0, 0.0, 1.0), ("H", 0.78, 0.0, 1.0)],
        ),
    ],
)
def test_SetSearchTypeTypeError(molecule1, molecule2):
    """
    Test Type Error associated with search methodology type
    """
    with pytest.raises(TypeError) as e:
        obj_sc = SearchConfig(Cluster(molecule1, molecule2))
        obj_sc.SetSearchType(1.0)
    assert (
        str(e.value) == "\n\nThe new search methodology is not a string"
        f"\nplease, check: '{type(1.0)}'\n"
    )
    with pytest.raises(TypeError) as e:
        obj_sc = SearchConfig(Cluster(molecule1, molecule2))
        obj_sc.SetSearchType((1.0,))
    assert (
        str(e.value) == "\n\nThe new search methodology is not a string"
        f"\nplease, check: '{type((1.0,))}'\n"
    )
    with pytest.raises(TypeError) as e:
        obj_sc = SearchConfig(Cluster(molecule1, molecule2))
        obj_sc.SetSearchType(1)
    assert (
        str(e.value) == "\n\nThe new search methodology is not a string"
        f"\nplease, check: '{type(1)}'\n"
    )


@pytest.mark.parametrize(
    "molecule1, molecule2, search_methodology",
    [
        (
            [("H", 0.0, 0.0, 0.0), ("H", 0.78, 0.0, 0.0)],
            [("H", 0.0, 0.0, 1.0), ("H", 0.78, 0.0, 1.0)],
            "SHGOS",
        ),
    ],
)
def test_SetSearchTypeAllow(molecule1, molecule2, search_methodology):
    """
    Test ValueError associated with search methodology type
    when search metodology choose is higher than 3
    """
    with pytest.raises(ValueError) as e:
        obj_sc = SearchConfig(Cluster(molecule1, molecule2))
        obj_sc.SetSearchType(search_methodology)

    assert str(e.value) == f"Invalid value. options are: {available}"


@pytest.mark.parametrize(
    "molecule1, molecule2, new_search_methodology",
    [
        (
            [("H", 0.0, 0.0, 0.0), ("H", 0.78, 0.0, 0.0)],
            [("H", 0.0, 0.0, 1.0), ("H", 0.78, 0.0, 1.0)],
            "dual_annealing",
        ),
    ],
)
def test_SetGetSearchType(molecule1, molecule2, new_search_methodology):
    """
    Test SetSearchType when search metodology choose is lower than 3
    or equal to 3
    """
    obj_sc = SearchConfig(Cluster(molecule1, molecule2))
    obj_sc.SetSearchType(new_search_methodology)

    assert obj_sc.GetSearchType() == new_search_methodology


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
def test_SetGetBasisSet(molecule1, molecule2, new_basis_set):
    """
    Test SetBasisSet for new basis set
    """
    obj_sc = SearchConfig(Cluster(molecule1, molecule2))
    obj_sc.SetBasisSet(new_basis_set)
    assert obj_sc.GetBasisSet() == new_basis_set


@pytest.mark.parametrize(
    "molecule1, molecule2",
    [
        (
            [("H", 0.0, 0.0, 0.0), ("H", 0.78, 0.0, 0.0)],
            [("H", 0.0, 0.0, 1.0), ("H", 0.78, 0.0, 1.0)],
        ),
    ],
)
def test_SetBasisSetTypeError(molecule1, molecule2):
    """
    Test Type Error associated with basis set
    """
    with pytest.raises(TypeError) as e:
        obj_sc = SearchConfig(Cluster(molecule1, molecule2))
        obj_sc.SetBasisSet(1.0)
    assert (
        str(e.value) == "\n\nThe new name to basis set is not a string"
        f"\nplease, check: '{type(1.0)}'\n"
    )
    with pytest.raises(TypeError) as e:
        obj_sc = SearchConfig(Cluster(molecule1, molecule2))
        obj_sc.SetBasisSet((1.0,))
    assert (
        str(e.value) == "\n\nThe new name to basis set is not a string"
        f"\nplease, check: '{type((1.0,))}'\n"
    )
    with pytest.raises(TypeError) as e:
        obj_sc = SearchConfig(Cluster(molecule1, molecule2))
        obj_sc.SetBasisSet(1)
    assert (
        str(e.value) == "\n\nThe new name to basis set is not a string"
        f"\nplease, check: '{type(1)}'\n"
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
def test_SetCostFunctionTypeError(molecule1, molecule2):
    """
    Test TypeError SetCostFunction
    """
    with pytest.raises(TypeError) as e:
        obj_sc = SearchConfig(Cluster(molecule1, molecule2))
        obj_sc.SetCostFun(1.0)
    assert (
        str(e.value) == "\n\nThe new cost function is not a string"
        f"\nplease, check: '{type(1.0)}'\n"
    )
    with pytest.raises(TypeError) as e:
        obj_sc = SearchConfig(Cluster(molecule1, molecule2))
        obj_sc.SetCostFun([1.0])
    assert (
        str(e.value) == "\n\nThe new cost function is not a string"
        f"\nplease, check: '{type([1.0])}'\n"
    )
    with pytest.raises(TypeError) as e:
        obj_sc = SearchConfig(Cluster(molecule1, molecule2))
        obj_sc.SetCostFun((1.0,))
    assert (
        str(e.value) == "\n\nThe new cost function is not a string"
        f"\nplease, check: '{type((1.0,))}'\n"
    )


# @pytest.mark.parametrize(
#     "molecule1, molecule2, new_methodology",
#     [
#         (
#             [("H", 0.0, 0.0, 0.0), ("H", 0.78, 0.0, 0.0)],
#             [("H", 0.0, 0.0, 1.0), ("H", 0.78, 0.0, 1.0)],
#             2,
#         ),
#     ],
# )
# def test_SC_methodology_VE_set(molecule1, molecule2, new_methodology):
#     """
#     Test TypeError @methodology.setter
#     """
#     with pytest.raises(TypeError) as e:
#         obj_sc = SearchConfig(Cluster(molecule1, molecule2))
#         obj_sc.methodology = new_methodology
#     assert (
#         str(e.value) == "\n\nThe new name to methodology is not a string"
#         f"\nplease, check: '{type(new_methodology)}'\n"
#     )


# @pytest.mark.parametrize(
#     "molecule1, molecule2, expected",
#     [
#         (
#             [("H", 0.0, 0.0, 0.0), ("H", 0.78, 0.0, 0.0)],
#             [("H", 0.0, 0.0, 1.0), ("H", 0.78, 0.0, 1.0)],
#             "HF",
#         ),
#     ],
# )
# def test_SC_methodology_number_grep(molecule1, molecule2, expected):
#     """
#     Test default methodology (Hartree-Fock)
#     """
#     obj_sc = SearchConfig(Cluster(molecule1, molecule2))
#     assert obj_sc.methodology == expected


# @pytest.mark.parametrize(
#     "molecule1, molecule2, expected",
#     [
#         (
#             [("H", 0.0, 0.0, 0.0), ("H", 0.78, 0.0, 0.0)],
#             [("H", 0.0, 0.0, 1.0), ("H", 0.78, 0.0, 1.0)],
#             "DFT",
#         ),
#     ],
# )
# def test_SC_methodology_number_set(molecule1, molecule2, expected):
#     """
#     Test default methodology (Hartree-Fock)
#     """
#     obj_sc = SearchConfig(Cluster(molecule1, molecule2))
#     obj_sc.methodology = expected
#     assert obj_sc.methodology == expected


# @pytest.mark.parametrize(
#     "molecule1, molecule2, expectation",
#     [
#         (
#             [("H", 0.0, 0.0, 0.0), ("H", 0.78, 0.0, 0.0)],
#             [("H", 0.0, 0.0, 1.0), ("H", 0.78, 0.0, 1.0)],
#             "pyscf",
#         ),
#     ],
# )
# def test_SC_new_cost_function(molecule1, molecule2, expectation):
#     """
#     Test ValueError @sphere_radius.setter
#     """
#     obj_sc = SearchConfig(Cluster(molecule1, molecule2))
#     assert obj_sc._func_cost == expectation


# @pytest.mark.parametrize(
#     "molecule1, molecule2",
#     [
#         (
#             [("H", 0.0, 0.0, 0.0), ("H", 0.78, 0.0, 0.0)],
#             [("H", 0.0, 0.0, 1.0), ("H", 0.78, 0.0, 1.0)],
#         ),
#     ],
# )
# def test_SC_run_ascec_method(molecule1, molecule2):
#     """
#     Test SC.run method for ascec
#     """
#     kwargs = {"number_temperatures": 1, "maxCycles": 10}
#     SearchConfig(Cluster(molecule1, molecule2)).run(**kwargs)
#     with open("configurations.xyz", "r") as f:
#         readl = f.readline()
#     os.remove("configurations.xyz")
#     assert readl == "4\n"


# @pytest.mark.parametrize(
#     "molecule1, molecule2, method_min",
#     [
#         (
#             [("H", 0.0, 0.0, 0.0), ("H", 0.78, 0.0, 0.0)],
#             [("H", 0.0, 0.0, 1.0), ("H", 0.78, 0.0, 1.0)],
#             "dual_annealing",
#         ),
#     ],
# )
# def test_SC_run_da_method(molecule1, molecule2, method_min):
#     """
#     Test SC.run method for dual annealing
#     """
#     SearchConfig(Cluster(molecule1, molecule2), search_methodology=method_min).run(
#         maxfun=1, maxiter=1
#     )
#     with open("configurations.xyz", "r") as f:
#         readl = f.readline()
#     os.remove("configurations.xyz")
#     assert readl == "4\n"


# # @pytest.mark.parametrize(
# #     "molecule1, molecule2, search_meth",
# #     [
# #         (
# #             [("H", 0.0, 0.0, 0.0), ("H", 0.78, 0.0, 0.0)],
# #             [("H", 0.0, 0.0, 1.0), ("H", 0.78, 0.0, 1.0)],
# #             "SHGO",
# #         ),
# #     ],
# # )
# # def test_SC_run_shgo_method(molecule1, molecule2, search_meth):
# #     """
# #     Test SC.run method for shgo
# #     """
# #     SearchConfig(Cluster(molecule1, molecule2), search_methodology=search_meth).run(
# #         sampling_method="sobol", n=1
# #     )
# #     with open("configurations.xyz", "r") as f:
# #         readl = f.readline()
# #     os.remove("configurations.xyz")
# #     assert readl == "4\n"


# @pytest.mark.parametrize(
#     "molecule1, molecule2, search_meth, initer, maxiter",
#     [
#         (
#             [("H", 0.0, 0.0, 0.0), ("H", 0.78, 0.0, 0.0)],
#             [("H", 0.0, 0.0, 3.0), ("H", 0.78, 0.0, 3.0)],
#             "Bayesian",
#             3,
#             3,
#         ),
#     ],
# )
# def test_SC_run_bayesian_method(molecule1, molecule2, search_meth, initer, maxiter):
#     """
#     Test SC.run method for bayesian
#     """
#     SearchConfig(Cluster(molecule1, molecule2), search_methodology=search_meth).run(
#         initer=initer, maxiter=maxiter
#     )
#     with open("configurations.xyz", "r") as f:
#         readl = f.readline()
#     os.remove("configurations.xyz")
#     assert readl == "4\n"


# @pytest.mark.parametrize(
#     "molecule1, molecule2, search_meth, initer, maxiter, num_cores",
#     [
#         (
#             [("H", 0.0, 0.0, 0.0), ("H", 0.78, 0.0, 0.0)],
#             [("H", 0.0, 0.0, 3.0), ("H", 0.78, 0.0, 3.0)],
#             "Bayesian",
#             3,
#             3,
#             2,
#         ),
#     ],
# )
# def test_SC_run_parallel_bayesian_method(
#     molecule1, molecule2, search_meth, initer, maxiter, num_cores
# ):
#     """
#     Test SC.run method for bayesian
#     """
#     SearchConfig(Cluster(molecule1, molecule2), search_methodology=search_meth).run(
#         initer=initer, maxiter=maxiter, num_cores=num_cores
#     )
#     with open("configurations.xyz", "r") as f:
#         readl = f.readline()
#     os.remove("configurations.xyz")
#     assert readl == "4\n"


# @pytest.mark.parametrize(
#     "molecule1, molecule2, search_meth, initer, maxiter, MCMC",
#     [
#         (
#             [("H", 0.0, 0.0, 0.0), ("H", 0.78, 0.0, 0.0)],
#             [("H", 0.0, 0.0, 3.0), ("H", 0.78, 0.0, 3.0)],
#             "Bayesian",
#             3,
#             3,
#             True,
#         ),
#     ],
# )
# def test_SC_run_MCMC_bayesian_method(
#     molecule1, molecule2, search_meth, initer, maxiter, MCMC
# ):
#     """
#     Test SC.run method for bayesian
#     """
#     SearchConfig(Cluster(molecule1, molecule2), search_methodology=search_meth).run(
#         initer=initer, maxiter=maxiter, MCMC=MCMC
#     )
#     with open("configurations.xyz", "r") as f:
#         readl = f.readline()
#     os.remove("configurations.xyz")
#     assert readl == "4\n"


# @pytest.mark.parametrize(
#     "molecule1, molecule2",
#     [
#         (
#             [("H", 0.0, 0.0, 0.0), ("H", 0.78, 0.0, 0.0)],
#             [("H", 0.0, 0.0, 3.0), ("H", 0.78, 0.0, 3.0)],
#         ),
#     ],
# )
# def test_system_object_grep(molecule1, molecule2):
#     """
#     Test Ascec.criterion else
#     """
#     obj_sc = SearchConfig(Cluster(molecule1, molecule2))
#     obj_sc.system_object


# @pytest.mark.parametrize(
#     "molecule1, molecule2, expected",
#     [
#         (
#             [("H", 0.0, 0.0, 0.0), ("H", 0.78, 0.0, 0.0)],
#             [("H", 0.0, 0.0, 3.0), ("H", 0.78, 0.0, 3.0)],
#             "pyscf",
#         ),
#     ],
# )
# def test_func_cost_grep(molecule1, molecule2, expected):
#     """
#     Test func_cost grep
#     """
#     assert SearchConfig(Cluster(molecule1, molecule2)).func_cost == expected


# # This test is for ascec criterion
# @pytest.mark.parametrize(
#     "molecule1, molecule2",
#     [
#         (
#             [("H", 0.0, 0.0, 0.0), ("H", 0.78, 0.0, 0.0)],
#             [("H", 0.0, 0.0, 3.0), ("H", 0.78, 0.0, 3.0)],
#         ),
#     ],
# )
# def test_Ascec_ascec_method(molecule1, molecule2):
#     """
#     Test Ascec.criterion else
#     """
#     obj_sc = SearchConfig(Cluster(molecule1, molecule2))
#     ascec = Ascec(
#         obj_sc._system_object,
#         obj_sc._search_methodology,
#         obj_sc._methodology,
#         obj_sc._basis_set,
#         program="pyscf",
#         bounds=obj_sc._bounds,
#     )
#     ascec.e_before = -0.0001
#     ascec.energy_current = 0.0001
#     ascec.ascec_criterion(100.0)
