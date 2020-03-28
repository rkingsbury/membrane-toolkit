# Copyright (c) Ryan Kingsbury
# This source code is licensed under the BSD-style license found in the
# LICENSE file in the root directory of this source tree.

"""
Tests for unitized methods in core
"""
import pytest
from pint import DimensionalityError

from membrane_toolkit.core.unitized import ureg
from membrane_toolkit.core.unitized import *


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
    assert apparent_permselectivity(q1, q2, q3).check(['dimensionless'])
    with pytest.raises(ValueError):
        apparent_permselectivity(30, q2)
    with pytest.raises(DimensionalityError):
        apparent_permselectivity(q3, q2)