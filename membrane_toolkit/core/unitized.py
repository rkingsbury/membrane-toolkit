# Copyright (c) Ryan Kingsbury
# This source code is licensed under the BSD-style license found in the
# LICENSE file in the root directory of this source tree.

"""
Unitized versions of all methods in membrane_toolkit.core

Unitized methods should have the same name as the base method and be wrapped with
appropriate units using pint's .wraps() decorator. Unitized methods should perform
computations in base SI units unless there is a specific reason (performance,
numerical stability, etc.) to do otherwise.
"""
from pint import UnitRegistry
from membrane_toolkit.core import (
    diffusion_coefficient_mackie_meares,
    apparent_permselectivity,
    nernst_potential,
    donnan_equilibrium,
)

ureg = UnitRegistry()

# diffusion.py

diffusion_coefficient_mackie_meares = ureg.wraps(
    "m ** 2 / s", ("m ** 2 / s", ureg.dimensionless,)
)(diffusion_coefficient_mackie_meares)

# donnan.py

donnan_equilibrium = ureg.wraps(
    "mol/L", ("mol/L", "mol/L", None, None, None, None, None, None)
)(donnan_equilibrium)

# potential.py

apparent_permselectivity = ureg.wraps(
    ureg.dimensionless, ("V", "V", None)
)(apparent_permselectivity)

nernst_potential = ureg.wraps("V", ("=A", "=A", None, "degC"))(nernst_potential)
