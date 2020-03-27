# Copyright (c) Ryan Kingsbury
# This source code is licensed under the BSD-style license found in the
# LICENSE file in the root directory of this source tree.

"""
Tests for unitized methods in core
"""
import pytest

from membrane_toolkit.core.unitized import ureg
from membrane_toolkit.core.unitized import *


def test_diffusion_coefficient_mackie_meares():
    q1 = ureg.Quantity('1e-9 m**2/s')
    q2 = ureg.Quantity('0.4 dimensionless')
    assert diffusion_coefficient_mackie_meares(q1, q2) == ureg.Quantity('6.25e-11 m**2/s')