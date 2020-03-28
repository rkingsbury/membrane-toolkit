# Copyright (c) Ryan Kingsbury
# This source code is licensed under the BSD-style license found in the
# LICENSE file in the root directory of this source tree.

"""
Tests for membrane potential methods
"""
import pytest

from membrane_toolkit.core.potential import *


def test_apparent_permselectivity():
    """
    - E_mem equals 0
    - E_mem equals E_ideal
    - t_counter is greater than 1
    - t_counter is less than 0
    - symmetry (flipping signs of potential)
    - 2 numerical examples
    """
    # input checking
    with pytest.raises(ValueError, match="must be between 0 and 1"):
        apparent_permselectivity(-40, -40, 2)
    with pytest.raises(ValueError, match="must be between 0 and 1"):
        apparent_permselectivity(-40, -40, -2)

    # basic calculations
    assert apparent_permselectivity(0, -40) == 0
    assert apparent_permselectivity(-40, -40) == 1
    assert apparent_permselectivity(-40, -40) == apparent_permselectivity(
        -40, -40)

    # numerical examples
    assert apparent_permselectivity(30, 40, 0.5) == pytest.approx(0.75)
    assert apparent_permselectivity(35, 37.8, 0.396) == pytest.approx(0.9386, abs=1e-4)
