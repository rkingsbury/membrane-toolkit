# Copyright (c) Ryan Kingsbury
# This source code is licensed under the BSD-style license found in the
# LICENSE file in the root directory of this source tree.

"""
Unitized versions of all methods in membrane_toolkit.core

Unitized methods should have the same name as the base method and be wrapped with
appropriate units using pint's .wraps() decorator.
"""
from pint import UnitRegistry

ureg = UnitRegistry()

from membrane_toolkit.core.diffusion import (
    diffusion_coefficient_mackie_meares,
)


diffusion_coefficient_mackie_meares = ureg.wraps("=A", ("=A", ureg.dimensionless,))(
    diffusion_coefficient_mackie_meares
)