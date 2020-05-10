# Copyright (c) Ryan Kingsbury
# This source code is licensed under the BSD-style license found in the
# LICENSE file in the root directory of this source tree.

"""
Donnan exclusion module
"""


def donnan_equilibrium(
    C_bulk: float,
    C_fix: float,
    z_counter: int = 1,
    z_co: int = -1,
    nu_counter: int = 1,
    nu_co: int = 1,
    z_fix: int = -1,
    gamma: float = 1,
):
    """
    Calculate the Donnan equilibrium at the interface between Phase1 and Phase2.

    Args:
        C_bulk: bulk salt concentration [mol/L]
        C_fix: fixed charge concentration, without sign [mol/L]
        z_counter: signed charge of the counter-ion. Default = +1 [dimensionless]
        z_co: signed charge of the co-ion. By definition, this must have the same sign as the co-ion.
            Default = -1 [dimensionless]
        nu_counter: stoichiometric coefficient of the counter-ion.
            Default = 1 [dimensionless]
        nu_co: stoichiometric coefficient of the co-ion.
            Default = 1 [dimensionless]
        z_fix: signed charge of the fixed groups. Default = -1 [dimensionless]
        gamma: stoichiometrically-weighted ratio of salt activity coefficient in
            solution to that in the membrane Default = 1 [dimensionless].

    Returns:
        float: The co-ion concentration in the membrane [mol/L]. Note that for salts containing one or more
               multivalent ions, this is not the same as the mobile salt concentration.

    Raises:
        AssertionError: If the input stoichiometry is incorrect. Both \( \\nu \) must be positive, \( z_{counter} \)
            and \( z_{fix} \) must have opposite signs, and \( \\nu_{counter} * z_{counter} \) +
            \( \\nu_{co} * z_{co} \) must equal zero.

    Notes:
        The Donnan equilibrium between a membrane with fixed charged concentration \( \\bar C_{fix} \)
        (mol per L water sorbed) and a salt solution of bulk concentration \( C_s \) (mol/L) is given by:

        $$
        \\bar C_{co}^{\\nu_{co}} \\big ( \\frac{z_{co} \\bar C_{co} + z_{fix} \\bar C_{fix}}{z_{ct}} \
        \\big )^{\\nu_{ct}} = - \\Gamma \\nu_{ct}^{\\nu_{ct}} \\nu_{co}^{\\nu_{co}} C_s^{\\nu_{ct} + \\nu_{co}}
        $$

        where subscripts \( co \) and \( ct \) indicate the co-ion (same charge as the membrane) and counter-ion
        (opposite charge to the membrane), \( \\nu \) (dimensionless) are stoichiometric coefficients, and overbars
        indicate membrane-phase quantities, in units of moles per liter of water sorbed by the membrane. \( \\Gamma \)
        (dimensionless) is the ratio of activity coefficients in the bulk solution to those in the membrane, given by:

        $$
        \\Gamma = \\frac{\\gamma_{\\pm}^{\\nu_{ct} + \
            \\nu_{co}}}{\\bar \\gamma_{ct}^{\\nu_{ct}} \\bar \\gamma_{co}^{\\nu_{co}}}
        $$

        Traditionally, \( \\Gamma \) is either set to 1 (implying that ion activity coefficients are the same in the
        membrane and in bulk solution), or the ions are assumed to behave ideally in the membrane (activity coefficient
        in the membrane equal to 1), in which case $\Gamma$ equals the bulk solution activity coefficient. More
        recently, Manning theory has been used to compute the membrane-phase activity coefficients, making possible a
        direct calculation of \( \\Gamma \).

    References:
        Donnan, F. G. The theory of membrane equilibria. Chem. Rev. 1924, 1 (1), 73–90.

        Kamcev, J.; Galizia, M.; Benedetti, F. M.; Jang, E.-S.; Paul, D. R.;
        Freeman, B.; Manning, G. S. Partitioning of Mobile Ions Between Ion Exchange Polymers and Aqueous Salt
        Solutions: Importance of Counter-ion Condensation. Phys. Chem. Chem. Phys. 2016, No. 8, 6021–6031.

        Galizia, M.; Manning, G. S.; Paul, D. R.; Freeman, B. D. Ion partitioning between brines and ion exchange
        polymers. Polymer (Guildf). 2019, 165 (January), 91–100.

        Kingsbury, R. S.; Coronell, O. Modelling and validation of concentration dependence of ion exchange membrane
        permselectivity: significance of convection and Manning’s counter-ion condensation theory. Submitted.
    """
    # validate input arguments
    assert nu_counter > 0
    assert nu_co > 0
    assert nu_counter * z_counter == -1 * nu_co * z_co
    assert z_fix * z_counter < 0

    def _donnan_solver(C_co):
        # private function to interatively solve for co-ion concentration.

        # return the error squared so we can utilize a scalar minimization routine
        return C_co ** nu_co * (
            (z_co * C_co + z_fix * C_fix) / z_counter
        ) ** nu_counter + gamma * nu_co ** nu_co * nu_counter ** nu_counter * C_bulk ** (
            nu_counter + nu_co
        )

    # solve the function above using one of scipy's nonlinear solvers
    # from scipy.optimize import minimize
    from scipy.optimize import root_scalar

    # call a solver to solve for the co-ion concentration
    # result = minimize(_donnan_solver, 0.01, method="TNC", bounds=[(0, C_bulk * nu_co * 10)])
    result = root_scalar(
        _donnan_solver, x0=0.01, x1=C_bulk * nu_co, bracket=(0, C_bulk * nu_co * 2)
    )

    # after solving for the co-ion concentration,
    # calculate the counter-ion concentraiton, Cg
    # result_counter = -(result_co * z_co + zfix * Cfix) / z_counter

    if result.converged:
        return result.root
    else:
        raise ValueError("{} failed to find a solution".format(__name__))
