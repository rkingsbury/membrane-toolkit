# Copyright (c) Ryan Kingsbury
# This source code is licensed under the BSD-style license found in the
# LICENSE file in the root directory of this source tree.

"""
Tests for unitized methods in core
"""
import pytest
from pint import DimensionalityError

from membrane_toolkit.core.unitized import ureg
from membrane_toolkit.core.unitized import (
    diffusion_coefficient_mackie_meares,
    apparent_permselectivity,
    nernst_potential,
    donnan_equilibrium
)


def test_diffusion_coefficient_mackie_meares():
    q1 = ureg.Quantity("1e-9 m**2/s")
    q2 = ureg.Quantity("0.4 dimensionless")
    q3 = ureg.Quantity("1e-5 m")
    assert diffusion_coefficient_mackie_meares(q1, q2) == ureg.Quantity(
        "6.25e-11 m**2/s"
    )
    with pytest.raises(ValueError):
        diffusion_coefficient_mackie_meares(1e-9, q2)
    with pytest.raises(DimensionalityError):
        diffusion_coefficient_mackie_meares(q1, q3)


def test_apparent_permselectivity():
    q1 = ureg.Quantity("-0.03 V")
    q2 = ureg.Quantity("-40 mV")
    q3 = ureg.Quantity("0.5 dimensionless")
    assert apparent_permselectivity(q1, q2, q3).magnitude == pytest.approx(0.75)
    assert apparent_permselectivity(q1, q2, q3).check(["dimensionless"])
    with pytest.raises(ValueError):
        apparent_permselectivity(30, q2)
    with pytest.raises(DimensionalityError):
        apparent_permselectivity(q3, q2)


def test_nernst_potential():
    q1 = ureg.Quantity("1 mol/L")
    q2 = ureg.Quantity("100 mmol/L")
    q3 = ureg.Quantity("0.5 dimensionless")
    q4 = ureg.Quantity("0.1 dimensionless")
    q5 = ureg.Quantity(30, ureg.degC)
    q6 = ureg.Quantity(25, ureg.degC)
    assert nernst_potential(q1, q2, 1, q6).magnitude == pytest.approx(0.05916, abs=1e-5)
    assert nernst_potential(q1, q2, 1, q6).units == "volt"
    assert nernst_potential(q3, q4, 1, q5).magnitude == pytest.approx(0.04204, abs=1e-5)
    assert nernst_potential(q3, q4, 1, q5).units == "volt"
    with pytest.raises(ValueError):
        nernst_potential(1, 0.1)
    with pytest.raises(DimensionalityError):
        nernst_potential(q3, q2)


def test_donnan_equilibrium():
    Cs = ureg.Quantity("0.5 mol/L")
    Cfix = ureg.Quantity("4 mol/L")
    q3 = ureg.Quantity("0.5 m")
    q4 = ureg.Quantity("-4 m")
    assert donnan_equilibrium(Cs, Cfix).magnitude == pytest.approx(0.061553, abs=1e-5)
    with pytest.raises(ValueError):
        donnan_equilibrium(0.5, Cfix)
    with pytest.raises(DimensionalityError):
        donnan_equilibrium(q3, q4)
