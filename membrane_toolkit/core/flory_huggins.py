# Copyright (c) Ryan Kingsbury
# This source code is licensed under the BSD-style license found in the
# LICENSE file in the root directory of this source tree.

"""
Placeholder for code implementing the Flory-Huggins model of solvent
activity inside membranes.
"""


def flory_huggins_chi():
    """
    Calculate the Flory-Huggins chi parameter from water uptake.


# Calculate the Flory-Huggins interaction parameter by assuming that the water sorption data fit Flory-Huggins theory
#
# $$ \chi = \frac{\ln{\frac{a_w}{\phi_w} }-1+\phi_w}{(1-\phi_w)^2} $$


    def calculate_chi(row):
    chi = (math.log(row['Water Activity']/row['Vol_frac'])-1+row['Vol_frac'])/(1 - row['Vol_frac'])**2
    return chi

    """
    pass
