# Copyright (c) Ryan Kingsbury
# This source code is licensed under the BSD-style license found in the
# LICENSE file in the root directory of this source tree.

"""
Tests for the Manning theory library
"""

import pytest
from numpy import allclose

from membrane_toolkit.core import (
    get_activity_coefficient_manning,
    diffusion_coefficient_manning,
    manning_eql,
)


def test_activity_against_lit_monovalent():
    """
    Test activity coefficient calculations for monovalent, condensed
    counter-ion case.

    Data reported in Fig. 5 & 10 of Kamcev, J.; Paul, D. R.; Freeman, B. D.
    Ion Activity Coefficients in Ion Exchange Polymers: Applicability of
    Manning’s Counterion Condensation Theory. Macromolecules 2015, 48 (21),
    8011–8024.
    """
    # CR61 cation exchange membrane
    xi = 1.83  # Manning parameter
    Cfix = -3.21  # fixed charge conc. in mol/L water sorbed
    # 0.01 M external NaCl concentration
    Cs = 3e-4  # mol/L water sorbed
    assert allclose(
        get_activity_coefficient_manning(xi, Cfix, Cs), 0.2 ** 0.5, atol=0.01
    )
    # 1 M external NaCl concentration
    Cs = 0.4  # mol/L water sorbed
    assert allclose(
        get_activity_coefficient_manning(xi, Cfix, Cs), 0.29 ** 0.5, atol=0.01
    )


def test_activity_against_lit_multivalent():
    """
    Test activity coefficient calculations for multivalent, condensed
    counter-ion case.

    Data reported in Fig. 6 of (1) Galizia, M.; Manning, G. S.; Paul,
    D. R.; Freeman, B. D. Ion partitioning between brines and ion exchange
    polymers. Polymer (Guildf). 2019, 165 (January), 91–100.
    """
    # CR61 cation exchange membrane
    xi = 1.83  # Manning parameter
    Cfix = -3.21  # fixed charge conc. in mol/L water sorbed
    # 0.1 M external NaCl concentration
    Cs = 0.035  # Mobile salt concentration is Cl- conc. / 2
    assert allclose(
        get_activity_coefficient_manning(
            xi, Cfix, Cs, z_counter=2, z_co=-1, nu_counter=1, nu_co=2
        ),
        0.15 ** 0.33,
        atol=0.01,
    )
    # 4 M external NaCl concentration
    Cs = 3.5  # Mobile salt concentration is Cl- conc. / 2
    assert allclose(
        get_activity_coefficient_manning(
            xi, Cfix, Cs, z_counter=2, z_co=-1, nu_counter=1, nu_co=2
        ),
        0.8 ** 0.33,
        atol=0.035,
    )


def test_activity_against_lit_monovalent_uncondensed():
    """
    Test activity coefficient calculations for monovalent, uncondensed
    counter-ion case.

    Data reported in Fig. 6 of Kamcev, J.; Paul, D. R.; Freeman, B. D.
    Equilibrium ion partitioning between aqueous salt solutions and inhomogeneous
    ion exchange membranes. Desalination 2018, 446 (September), 31–41.
    """
    # CA267 cation exchange membrane
    xi = 1.00  # Manning parameter
    Cfix = -2.66  # fixed charge conc. in mol/L water sorbed
    # 0.03 M external NaCl concentration
    Cs = 0.0007  # mol/L water sorbed
    # bulk solution activity coefficient is 0.85 at this concentration
    # so gamma_membrane^2 = C_bulk ^2 0.85 ^2 / C_+ C_-
    # C- = Cs and C+ = Cs + Cfix
    # so gamma_membrane =( 0.03^2 0.85^2 / 0.0007 / 2.6607 )^0.5 = 0.590
    assert allclose(
        get_activity_coefficient_manning(xi, Cfix, Cs), 0.590, atol=0.02
    )


def test_activity_continuity_monovalent():
    """
    Test that the activity coefficients are continuous across the critical
    value of manning parameter

    """
    # CA267 cation exchange membrane
    Cfix = -2.66  # fixed charge conc. in mol/L water sorbed
    # 0.03 M external NaCl concentration
    Cs = 0.0007  # mol/L water sorbed
    # the critical value of xi is 1 for a monovalent counter-ion
    assert allclose(
        get_activity_coefficient_manning(1.01, Cfix, Cs),
        get_activity_coefficient_manning(0.99, Cfix, Cs),
        atol=0.01
    )


def test_bad_input_type():
    with pytest.raises(Exception, match="Invalid"):
        get_activity_coefficient_manning(2, -3, 3e-4, type="blah")


def test_sign_mismatch_1():
    with pytest.raises(Exception, match="Mismatch"):
        get_activity_coefficient_manning(2, 3, 3e-4)


def test_sign_mismatch_2():
    with pytest.raises(Exception, match="Mismatch"):
        get_activity_coefficient_manning(2, -3, 3e-4, z_counter=-1)


def test_off_stoichiometry():
    with pytest.raises(Exception, match="stoichiometry"):
        get_activity_coefficient_manning(
            2, -3, 3e-4, z_counter=1, z_co=-2, nu_counter=1, nu_co=2
        )


def test_diffusion_against_lit():
    # Test diffusion coefficient calculations against
    # data reported in XXXXX
    pass


def test_concentration_against_lit():
    # Test diffusion coefficient calculations against
    # data reported in XXXXX
    pass


def test_A_factor():
    pass


def test_charge_mismatch():
    # test the warning when the sign of the counter-ion
    # and the sign of the fixed charges match
    pass


def test_symmetry():
    # verify that we get equal diffusion coefficients,
    # activity coefficients, and concentrations when
    # we reverse the charge of the membrane (i.e.,
    # reverse the co- and counter-ions)
    from pyEQL import Solution
    bulk_solution = Solution(
        [["Na+", "0.5 mol/L"], ["Cl-", "0.5 mol/L"]], temperature="25 degC"
    )
    # Manning parameter
    xi = 1.5
    # charge density
    CD = 1
    # volume fraction of water
    vol_frac = 0.3

    # CEM case
    Cfix = -1 * CD
    z_counter = 1
    z_co = -1
    nu_counter = 1
    nu_co = 1

    s_mem = manning_eql(bulk_solution, xi, str(Cfix) + "mol/L")

    Cc_ct = s_mem.get_amount("Na+", "mol/L").magnitude
    Cc_co = s_mem.get_amount("Cl-", "mol/L").magnitude

    Dc_ct = diffusion_coefficient_manning(
        xi, Cfix, Cc_co, vol_frac, "counter", nu_counter, nu_co, z_counter, z_co
    )
    Dc_co = diffusion_coefficient_manning(
        xi, Cfix, Cc_co, vol_frac, "co", nu_counter, nu_co, z_counter, z_co
    )

    Ac_ct = get_activity_coefficient_manning(
        xi, Cfix, Cc_co, "counter", nu_counter, nu_co, z_counter, z_co
    )
    Ac_co = get_activity_coefficient_manning(
        xi, Cfix, Cc_co, "co", nu_counter, nu_co, z_counter, z_co
    )

    # AEM case
    Cfix = CD
    z_counter = -1
    z_co = 1
    nu_counter = 1
    nu_co = 1

    s_mem = manning_eql(bulk_solution, xi, str(Cfix) + "mol/L")

    Ca_ct = s_mem.get_amount("Cl-", "mol/L").magnitude
    Ca_co = s_mem.get_amount("Na+", "mol/L").magnitude

    Da_ct = diffusion_coefficient_manning(
        xi, Cfix, Ca_co, vol_frac, "counter", nu_counter, nu_co, z_counter, z_co
    )
    Da_co = diffusion_coefficient_manning(
        xi, Cfix, Ca_co, vol_frac, "co", nu_counter, nu_co, z_counter, z_co
    )

    Aa_ct = get_activity_coefficient_manning(
        xi, Cfix, Ca_co, "counter", nu_counter, nu_co, z_counter, z_co
    )
    Aa_co = get_activity_coefficient_manning(
        xi, Cfix, Ca_co, "co", nu_counter, nu_co, z_counter, z_co
    )

    assert Cc_ct == Ca_ct
    assert Cc_co == Ca_co
    assert Dc_ct == Da_ct
    assert Dc_co == Da_co
    assert Ac_ct == Aa_ct
    assert Ac_co == Aa_co
