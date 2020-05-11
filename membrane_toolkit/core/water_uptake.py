# Copyright (c) Ryan Kingsbury
# This source code is licensed under the BSD-style license found in the
# LICENSE file in the root directory of this source tree.

"""
Methods related to water uptake and swelling
"""

# ### Calculate Water Volume  Fraction and Sorption Coefficient

# Water volume fraction is calculated from the uncorrected water uptake / swelling degree as
# 
# $$ \phi_w = \frac{SD}{(SD + \rho_w / \rho_p)} $$
# 
# where $\rho_w$ and $\rho_p$ are the density of water and the polymer, respectively
# 
# Water sorption (partition) coefficient is equal to
# 
# $$ K_w = \frac{\phi_w M_w}{C_w^s V_w} $$
# 
# Where $M_w$ and $V_w$ are the molar mass and molar volume of water (18.015 g/mol and 0.0182 L/mol) and $C_w^s$ is the mass concentration of water in the bulk solution (kg/L). Note that the water self-diffusion coefficient in pure water is $2.8x10^{-9} m^2/s$
# 
# See Geise et al. [doi: 10.1016/j.progpolymsci.2013.07.001](http://www.dx.doi.org/10.1016/j.progpolymsci.2013.07.001)
# 


def water_volume_fraction():
    #df_hydration['Vol_frac'] = df_hydration['SD'] / (df_hydration['SD'] + rho_water/rho_polymer)
    pass


def fixed_charge_concentration():
    # summary['Fixed Charge Conc. (mol/L)','mean'] = summary['IEC (meq/g)']['mean']/summary['SD']['mean']*0.998
    pass


def water_partition_coefficient():
    pass
