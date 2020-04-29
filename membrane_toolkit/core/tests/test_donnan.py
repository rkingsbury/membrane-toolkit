# Copyright (c) Ryan Kingsbury
# This source code is licensed under the BSD-style license found in the
# LICENSE file in the root directory of this source tree.

"""
Tests for Donnan exclusion methods
"""
import pytest

from membrane_toolkit.core import donnan_equilibrium


def test_donnan_equilibrium():
    """
    Cases to check:

    x same answer for positive or negative co-ion
    x invalid input stoichiometry (nu < 0)
    x invalid input stoichiometry (nu_co * z_co != nu_counter * z_counter)
    x invalid input stoichiometry (z_fix and z_counter same sign)
    x Cfix = 0
    - Cs = Cfix
    x monovalent, 1:1 case with gamma = 1 equals the result from the asinh formula
    - Numerical example with gamma not equal to one, monovalent
    - Numerical example with gamma not equal to one, multivalent

    Numerical examples drawn from data reported in Fig. 1 of (1) Kamcev, J.; Galizia, M.; Benedetti, F. M.; Jang,
    E.-S.; Paul, D. R.; Freeman, B.; Manning, G. S. Partitioning of Mobile Ions Between Ion Exchange Polymers and
    Aqueous Salt Solutions: Importance of Counter-ion Condensation. Phys. Chem. Chem. Phys. 2016, No. 8, 6021–6031.
    and Galizia, M.; Manning, G. S.; Paul, D. R.; Freeman, B. D. Ion partitioning between brines and ion exchange
    polymers. Polymer (Guildf). 2019, 165 (January), 91–100.
    """
    assert donnan_equilibrium(0.5, 4) == pytest.approx(0.061553, abs=1e-5)
    assert donnan_equilibrium(0.5, 4, z_counter=-1, z_co=1, z_fix=1) == pytest.approx(
        0.061553, abs=1e-5
    )
    assert donnan_equilibrium(
        0.5, 4, z_counter=2, z_co=-2, z_fix=-1
    ) == donnan_equilibrium(0.5, 4, z_counter=-2, z_co=2, z_fix=1)
    with pytest.raises(AssertionError):
        donnan_equilibrium(0.5, 4, z_counter=-1, z_fix=-1)
    with pytest.raises(AssertionError):
        donnan_equilibrium(0.5, 4, z_counter=-1, z_fix=-1)
    with pytest.raises(AssertionError):
        donnan_equilibrium(0.5, 4, z_counter=-1, z_co=1, z_fix=-1)
    with pytest.raises(AssertionError):
        donnan_equilibrium(0.5, 4, nu_co=-1)
    # when C_fix=0, co-ion concentration should equal bulk ion concentration times (z_co / z_ct) ** 1/(nu_co+nu_counter)
    assert donnan_equilibrium(0.5, 0) == pytest.approx(0.5, abs=1e-5)
    assert donnan_equilibrium(0.5, 0, z_counter=-1, z_co=1, z_fix=1) == pytest.approx(
        0.5, abs=1e-5
    )
    assert donnan_equilibrium(0.5, 0, nu_co=2, z_counter=2) == pytest.approx(
        0.5 * 2, abs=1e-5
    )
    assert donnan_equilibrium(
        0.5, 0, nu_co=2, z_co=1, z_counter=-2, z_fix=1
    ) == pytest.approx(0.5 * 2, abs=1e-5)
    assert donnan_equilibrium(0.5, 0, nu_co=3, z_counter=3) == pytest.approx(
        0.5 * 3, abs=1e-5
    )

    # CR61 membrane, 0.1 M NaCl, gamma=1
    assert donnan_equilibrium(0.1, 3.2) == pytest.approx(0.003, abs=1e-3)
    # CR61 membrane, 1 M NaCl, gamma = 0.42
    assert donnan_equilibrium(1.0, 3.2, gamma=0.42) == pytest.approx(0.1, abs=5e-2)
    # CR61 membrane, 0.3 M CaCl2, gamma=1
    assert donnan_equilibrium(
        0.3, 3.0, nu_counter=1, nu_co=2, z_counter=2
    ) == pytest.approx(0.25, abs=5e-2)
