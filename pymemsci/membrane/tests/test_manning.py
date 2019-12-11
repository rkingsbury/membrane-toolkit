# Tests for the Manning theory library

import pytest

import pyEQL
import manning
from numpy import allclose


class TestActivity:
    def test_activity_against_lit_monovalent(self):
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
        Cfix = "-3.21 mol/L"  # fixed charge conc. in mol/L water sorbed
        # 0.01 M external NaCl concentration
        Cs = "3e-4 mol/L"
        assert allclose(
            manning.get_activity_coefficient_manning(xi, Cfix, Cs),
            0.2 ** 0.5,
            atol=0.01,
        )
        # 1 M external NaCl concentration
        Cs = "0.4 mol/L"
        assert allclose(
            manning.get_activity_coefficient_manning(xi, Cfix, Cs),
            0.29 ** 0.5,
            atol=0.01,
        )

    def test_activity_against_lit_multivalent(self):
        """
        Test activity coefficient calculations for multivalent, condensed
        counter-ion case.

        Data reported in Fig. 6 of (1) Galizia, M.; Manning, G. S.; Paul, 
        D. R.; Freeman, B. D. Ion partitioning between brines and ion exchange 
        polymers. Polymer (Guildf). 2019, 165 (January), 91–100.
        """
        # CR61 cation exchange membrane
        xi = 1.83  # Manning parameter
        Cfix = "-3.21 mol/L"  # fixed charge conc. in mol/L water sorbed
        # 0.1 M external NaCl concentration
        Cs = "0.035 mol/L"  # Mobile salt concentration is Cl- conc. / 2
        assert allclose(
            manning.get_activity_coefficient_manning(
                xi, Cfix, Cs, z_counter=2, z_co=-1, nu_counter=1, nu_co=2
            ),
            0.15 ** 0.33,
            atol=0.01,
        )
        # 4 M external NaCl concentration
        Cs = "3.5 mol/L"  # Mobile salt concentration is Cl- conc. / 2
        assert allclose(
            manning.get_activity_coefficient_manning(
                xi, Cfix, Cs, z_counter=2, z_co=-1, nu_counter=1, nu_co=2
            ),
            0.8 ** 0.33,
            atol=0.035,
        )

    def test_bad_input_type(self):
        with pytest.raises(Exception, match="Invalid"):
            manning.get_activity_coefficient_manning(
                2, "-3 mol/L", "3e-4 mol/L", type="blah"
            )

    def test_sign_mismatch_1(self):
        with pytest.raises(Exception, match="Mismatch"):
            manning.get_activity_coefficient_manning(2, "3 mol/L", "3e-4 mol/L")

    def test_sign_mismatch_2(self):
        with pytest.raises(Exception, match="Mismatch"):
            manning.get_activity_coefficient_manning(
                2, "-3 mol/L", "3e-4 mol/L", z_counter=-1
            )

    def test_diffusion_against_lit(self):
        # Test diffusion coefficient calculations against
        # data reported in XXXXX
        pass

    def test_concentration_against_lit(self):
        # Test diffusion coefficient calculations against
        # data reported in XXXXX
        pass

    def test_A_factor(self):
        pass

    def test_charge_mismatch(self):
        # test the warning when the sign of the counter-ion
        # and the sign of the fixed charges match
        pass


def test_symmetry():
    # verify that we get equal diffusion coefficients,
    # activity coefficients, and concentrations when
    # we reverse the charge of the membrane (i.e.,
    # reverse the co- and counter-ions)
    bulk_solution = pyEQL.Solution(
        [["Na+", "0.5 mol/L"], ["Cl-", "0.5 mol/L"]], temperature="25 degC"
    )
    # Manning parameter
    xi = 1.5
    # charge density
    CD = 1
    # volume fraction of water
    vol_frac = 0.3

    # CEM case
    Cfix = str(-1 * CD) + "mol/L"
    z_counter = 1
    z_co = -1
    nu_counter = 1
    nu_co = 1

    s_mem = manning.manning_eql(bulk_solution, xi, Cfix)

    Cc_ct = s_mem.get_amount("Na+", "mol/L")
    Cc_co = s_mem.get_amount("Cl-", "mol/L")

    Dc_ct = manning.get_diffusion_coefficient_manning(
        xi, Cfix, str(Cc_co), vol_frac, "counter", nu_counter, nu_co, z_counter, z_co
    )
    Dc_co = manning.get_diffusion_coefficient_manning(
        xi, Cfix, str(Cc_co), vol_frac, "co", nu_counter, nu_co, z_counter, z_co
    )

    Ac_ct = manning.get_activity_coefficient_manning(
        xi, Cfix, str(Cc_co), "counter", nu_counter, nu_co, z_counter, z_co
    )
    Ac_co = manning.get_activity_coefficient_manning(
        xi, Cfix, str(Cc_co), "co", nu_counter, nu_co, z_counter, z_co
    )

    # AEM case
    Cfix = str(CD) + "mol/L"
    z_counter = -1
    z_co = 1
    nu_counter = 1
    nu_co = 1

    s_mem = manning.manning_eql(bulk_solution, xi, Cfix)

    Ca_ct = s_mem.get_amount("Cl-", "mol/L")
    Ca_co = s_mem.get_amount("Na+", "mol/L")

    Da_ct = manning.get_diffusion_coefficient_manning(
        xi, Cfix, str(Ca_co), vol_frac, "counter", nu_counter, nu_co, z_counter, z_co
    )
    Da_co = manning.get_diffusion_coefficient_manning(
        xi, Cfix, str(Ca_co), vol_frac, "co", nu_counter, nu_co, z_counter, z_co
    )

    Aa_ct = manning.get_activity_coefficient_manning(
        xi, Cfix, str(Ca_co), "counter", nu_counter, nu_co, z_counter, z_co
    )
    Aa_co = manning.get_activity_coefficient_manning(
        xi, Cfix, str(Ca_co), "co", nu_counter, nu_co, z_counter, z_co
    )

    assert Cc_ct == Ca_ct
    assert Cc_co == Ca_co
    assert Dc_ct == Da_ct
    assert Dc_co == Da_co
    assert Ac_ct == Aa_ct
    assert Ac_co == Aa_co
