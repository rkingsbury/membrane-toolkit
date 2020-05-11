# Copyright (c) Ryan Kingsbury
# This source code is licensed under the BSD-style license found in the
# LICENSE file in the root directory of this source tree.

"""
Membrane potential related methods
"""
import numpy as np


def apparent_permselectivity(
    E_mem: float, E_ideal: float, t_counter: float = 0.5,
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

        $$
            \\alpha_{app} = \\frac{\\frac{E_{mem}}{E_{ideal}} + 1 - 2 t_{counter}}
                {2 t_{co}}
        $$

    References:
        Kingsbury, R. S.; Coronell, O. Modelling and validation of concentration
        dependence of ion exchange membrane permselectivity: significance of convection
        and Manning’s counter-ion condensation theory. Submitted.
    """
    # t_counter must be between 0 and 1
    if t_counter < 0 or t_counter > 1:
        raise ValueError(
            "Counter-ion transport number must be between 0 and 1. "
            "t_counter={:.3f}".format(t_counter)
        )

    return (E_mem / E_ideal + 1 - 2 * t_counter) / (2 - 2 * t_counter)


def transport_number():
    """
    Calculate the transport number of an ionic species in the membrane
    """
    pass


def streaming_potential():
    pass


def membrane_potential():
    pass


def diffusion_potential():
    pass


def nernst_potential(
    a0: float, aL: float, z_ct: int = 1, temperature: float = 25
) -> float:
    """
    Calculate the ideal membrane or interface potential according to the Nernst
    equation.

    Args:
        a0 (float): Activity (or concentration) of the electrolyte on the left
            (x=0) side of the membrane. The units of a0 and aL must match.
        aL (float): Activity (or concentration) of the electrolyte on the left
            (x=0) side of the membrane. The units of a0 and aL must match.
        z_ct (int): signed charge of the counter-ion species. Default = +1
        temperature (float): Temperature of the electrolyte [C]. Default = 25 C

    Returns:
        The potential [V] of the left (x=0) interface with respect to the right (x=L)
        interface, according to the Nernst equation.

    Notes:
        The Nernst equation gives the electrical potential across an ideally-selective
        membrane separating two electrolyte solutions.

        $$
            E = \\frac{RT}{z_{ct} F}\\log{\\frac{a0}{aL}}
        $$

        where R (8.314 J/mol K) is the ideal gas constant and F (96485 C/mol) is the
        Faraday constant.

    References:
        Bard, A. J.; Faulkner, L. R. Electrochemical Methods: Fundamentals and
        Applications, 2nd ed.; John Wiley & Sons, 2001.

        Helfferich, F. Ion Exchange; McGraw-Hill: New York, 1962.

        Winger, A.; Bodamer, G.; Kunin, R. Some electrochemical properties of new
        synthetic ion exchange memebranes. J. Electrochem. Soc 1953, 100 (4), 178–184.
    """
    if a0 < 0:
        raise ValueError(
            "Received invalid activity argument of a0 = {}. Electrolyte"
            "activity must be positive.".format(a0)
        )

    if aL < 0:
        raise ValueError(
            "Received invalid activity argument of aL = {}. Electrolyte"
            "activity must be positive.".format(aL)
        )

    if temperature < -273.15:
        raise ValueError(
            "Received invalid temperature argument of {}. Temperature"
            "is below absolute zero!".format(temperature)
        )

    return 8.314 * (temperature + 273.15) / z_ct / 96485 * np.log(a0 / aL)
