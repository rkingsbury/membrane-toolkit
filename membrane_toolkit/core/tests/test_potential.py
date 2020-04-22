# Copyright (c) Ryan Kingsbury
# This source code is licensed under the BSD-style license found in the
# LICENSE file in the root directory of this source tree.

"""
Tests for membrane potential methods
"""
import pytest

from membrane_toolkit.core.potential import apparent_permselectivity, nernst_potential


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
    assert apparent_permselectivity(-40, -40) == apparent_permselectivity(-40, -40)

    # numerical examples
    assert apparent_permselectivity(30, 40, 0.5) == pytest.approx(0.75)
    assert apparent_permselectivity(35, 37.8, 0.396) == pytest.approx(0.9386, abs=1e-4)


def test_nernst_potential():
    """
    - a0 negative
    - aL negative
    - temperature negative
    - a0 = aL
    - flipping the sign of the counter ion
    - z_ct not equal to 1 or -1
    - 2 numerical examples
    """
    with pytest.raises(ValueError, match="invalid activity argument of a0"):
        nernst_potential(-0.2, 0.5)
    with pytest.raises(ValueError, match="invalid activity argument of aL"):
        nernst_potential(0.2, -0.5)
    with pytest.raises(ValueError, match="invalid temperature argument"):
        nernst_potential(0.2, 0.5, 1, -280)
    assert nernst_potential(0.5, 0.5, 2, 450) == 0
    assert nernst_potential(0.5, 0.1, -2) == -1 * nernst_potential(0.5, 0.1, 2)
    assert nernst_potential(1, 0.1) == pytest.approx(0.05916, abs=1e-5)
    assert nernst_potential(50, 2, 2, 76.85) == pytest.approx(0.048539, abs=1e-5)
