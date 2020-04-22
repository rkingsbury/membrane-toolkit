# Copyright (c) Ryan Kingsbury
# This source code is licensed under the BSD-style license found in the
# LICENSE file in the root directory of this source tree.

"""
Tests for diffusion coefficient methods
"""
import pytest

from membrane_toolkit.core.diffusion import diffusion_coefficient_mackie_meares


def test_diffusion_coefficient_mackie_meares():
    """
    Cases to check:
    - D_bulk is zero
    - phi_w is zero
    - both are zero
    - phi_w less than 0
    - phi_w greater than 1
    - phi_w equal to 1
    - numerical example
    """
    assert diffusion_coefficient_mackie_meares(0, 0.5) == 0
    assert diffusion_coefficient_mackie_meares(1e-9, 0) == 0
    assert diffusion_coefficient_mackie_meares(1e-9, 1) == pytest.approx(1e-9)
    assert diffusion_coefficient_mackie_meares(0, 0) == 0
    with pytest.raises(ValueError, match="between 0 and 1"):
        diffusion_coefficient_mackie_meares(0, -0.5)
    with pytest.raises(ValueError, match="between 0 and 1"):
        diffusion_coefficient_mackie_meares(0, 2)
    assert diffusion_coefficient_mackie_meares(1e-10, 0.4) == pytest.approx(6.25e-12)
