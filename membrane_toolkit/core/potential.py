# Copyright (c) Ryan Kingsbury
# This source code is licensed under the BSD-style license found in the
# LICENSE file in the root directory of this source tree.

"""
Membrane potential related methods
"""


def apparent_permselectivity(
    E_mem: float,
    E_ideal: float,
    t_counter: float = 0.5,
) -> float:
    """
    Calculate the apparent permselectivity of a membrane from the membrane potential.

    Args:
        E_mem (float): Electrical potential across the membrane [V]
        E_ideal (float): Electrical potential across an ideally-selective membrane,
            usually calculated via the Nernst equation [V]
        t_counter (float): bulk solution transport number of the counter ion.
            [dimensionless]. Must be between 0 and 1.

    Returns:
        float: Apparent permselectivity of the membrane, [dimensionless]

    Notes:
        Apparent permselectivity ranges from 0 to 1, where 0 represents a non-selective
        membrane and 1 represents a perfectly-selective membrane. Apparent
        permselectivity is calculated from the membrane potential and differs from the
        "true" permselectivity. The equation is:

        .. math::

            \\alpha_{app} = \\frac{\\frac{E_{mem}}{E_{ideal}} + 1 - 2 t_{counter}}
                {2 t_{co}}

    References:
        Kingsbury, R. S.; Coronell, O. Modelling and validation of concentration
        dependence of ion exchange membrane permselectivity: significance of convection
        and Manning’s counter-ion condensation theory. J. Membrane Science XXXX.
    """
    # t_counter must be between 0 and 1
    if t_counter < 0 or t_counter > 1:
        raise ValueError(
            "Counter-ion transport number must be between 0 and 1. "
            "t_counter={:.3f}".format(t_counter)
        )

    return (E_mem / E_ideal + 1 - 2 * t_counter) / (2 - 2 * t_counter)
