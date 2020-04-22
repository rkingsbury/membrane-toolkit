# Copyright (c) Ryan Kingsbury
# This source code is licensed under the BSD-style license found in the
# LICENSE file in the root directory of this source tree.

"""
Donnan exclusion module
"""


def equilibrate_donnan(
    C_bulk, C_fix, z_counter, z_co, nu_counter, nu_co, z_fix, a_bulk
):
    """
    Calculate the Donnan equilibrium at the interface between Phase1 and Phase2.

    Args:
        C_bulk: bulk salt concentration
        C_fix: fixed charge concentration
        z_counter: signed charge of the counter-ion
        z_co: signed charge of the co-ion
        nu_counter: stoichiometric coefficient of the counter-ion
        nu_co: stoichiometric coefficient of the co-ion
        z_fix: signed charge of the fixed groups
        a_bulk: activity of the salt in the bulk solution

    References:
        (1) Ohshima, H.; Ohki, S. Donnan potential and surface potential of a
        charged membrane. Biophys. J. 1985, 47 (5), 673–678.

    """
    # TODO - add argument checking (nu's can't be negative, z_counter and z_co must
    # have opposite signs)

    def solve_donnan(C_co):
        # calculate the counter-ion concentration
        # C_counter = (-z_co * C_co - z_fix * C_fix) / z_counter
        # calculate the activities
        # TODO make these functions of the respective concentrations
        a_co = 1
        a_counter = 1

        # calculate the gamma factor
        gamma = (
            nu_counter ** z_counter
            * nu_co ** z_co
            * a_counter ** (nu_counter + nu_co)
            / a_counter ** nu_counter
            / a_co ** nu_co
        )
        # return the error squared so we can utilize a scalar minimization routine
        return (
            C_co ** nu_co * ((z_co * C_co + z_fix * C_fix) / z_counter) ** nu_counter
            - gamma * C_bulk ** (nu_counter + nu_co)
        ) ** 2

    # solve the function above using one of scipy's nonlinear solvers
    from scipy.optimize import minimize_scalar

    # call a solver to solve for the co-ion concentration
    result = minimize_scalar(solve_donnan, method="bounded", bounds=[0, C_bulk / nu_co])

    # after solving for the co-ion concentration,
    # calculate the counter-ion concentraiton, Cg
    # result_counter = -(result_co * z_co + zfix * Cfix) / z_counter

    if result.success:
        return result.x
    else:
        raise ValueError("{} failed to find a solution".format(__name__))
