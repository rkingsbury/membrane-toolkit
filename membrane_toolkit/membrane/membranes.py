# Copyright (c) Ryan Kingsbury
# This source code is licensed under the BSD-style license found in the
# LICENSE file in the root directory of this source tree.

# Membrane class definitions


class Membrane:
    """
    Contain the properties of each IX membrane in an object

    Parameters
    ----------
    name : str
            A short identifier for the membrane (no spaces)

    IEC : str quantity, optional
            A string representing the ion exchange capacity of the membrane,
            mmol of charge per g of dry membrane, including sign.
            Defaults to -1.2 mmol/g if not specified.

    swelling : str quantity, optional
            A string representing the swelling degree of the membrane,
            g water per g membrane.
            Defaults to 0.20 if not specified

    # fixed_charge : str quantity, optional
    #         A string representing the fixed charge concentration in the
    #         membrane, including the unit and sign. e.g. '-4 mol/kg'
    #         Will be calculated from IEC and Swelling Degree if not specified.
    #
    description : str, optional
            A longer, descriptive name for the membrane, e.g. 'Neosepta ACS'

    permselectivity : number, optional
            A number between 0 and 1 representing the membrane permeselectivity.
            Defaults to 1 (ideal membrane) if not specified.

    thickness : str quantity , optional
            A string representing the thickness of the membrane, including
            the unit, e.g. '20 um'
            Defaults to '100 um' if not specified.

    area_resistance : str quantity, optional
            A string representing the area resistance of the membrane, including
            the unit, e.g. '3 ohm * cm **2'
            Defaults to '0 ohm * cm ** 2' (ideal membrane) if not specified

    K_ff : number, optional
            A "form factor" used to estimate ion diffusion coefficients inside
            the membrane based on diffusion coefficients in solution.
            D_membrane = D_solution * K_ff
            Defaults to 0.01 if not specified.

    D_water : str quantity, optional
            A string representing the diffusion coefficient of water, including
            the unit, e.g. '1.3e-9 m ** 2 /s'
            Defaults to 0 (ideal membrane) if not specified.

    t_water : number, optional
            The transport number of water, mol water per mol charge transported.

    """

    def __init__(self, name, FCD, **kwargs):

        self.id = name
        self.IEC = "1.2 mmol/g"
        self.swelling = "0.2 g/g"
        # self.fixed_charge_density = fixed_charge
        self.description = ""
        self.permselectivity = 1
        self.thickness = "100 um"
        self.area_resistance = "0 ohm * cm ** 2"
        self.K_ff = 0.01
        self.D_water = "0 m ** 2 / s"
        self.t_water = 5
        self.dielectric_constant = 30
        self.FCD = FCD

        for key in kwargs:
            if key == "description":
                self.description = kwargs[key]
            elif key == "IEC":
                self.IEC = kwargs[key]
            elif key == "swelling":
                self.swelling = kwargs[key]
            elif key == "permselectivity":
                self.permselectivity = kwargs[key]
            elif key == "thickness":
                self.thickness = kwargs[key]
            elif key == "area_resistance":
                self.area_resistance = kwargs[key]
            elif key == "K_ff":
                self.K_ff = kwargs[key]
            elif key == "D_water":
                self.D_water = kwargs[key]
            elif key == "t_water":
                self.t_water = kwargs[key]
            elif key == "dielectric_constant":
                self.dielectric_constant = kwargs[key]

    # simple methods to access the main properties

    def get_name(self):
        return self.id

    def get_description(self):
        return self.description

    def get_fixed_charge_density(self, swelling_degree=None):
        # if swelling_degree is None:
        #        SD = self.swelling
        #    else:
        #       SD = swelling_degree

        #    FCD = unit(self.IEC) / unit(SD) * unit('0.998 g/mL')

        #        return str(FCD.to('mol/L'))
        return self.FCD

    def get_permselectivity(self):
        return self.permselectivity

    def get_thickness(self):
        return self.thickness

    def get_area_resistance(self):
        return self.area_resistance

    def get_K_ff(self):
        return self.K_ff

    def get_D_water(self):
        return self.D_water

    def get_t_water(self):
        return self.t_water

    def get_dielectric_constant(self):
        # TODO calculate this as a volume fraction weighted average of
        # salt solution and dry polymer (~6)
        print(
            "warning:dielectric constant arbitrarily set at %s"
            % self.dielectric_constant
        )
        return self.dielectric_constant

    def get_manning_parameter(self, temperature="25 degC"):
        """
        Calculate the Manning parameter of the membrane

        Parameters
        ----------
        temperature : str quantity, optional
                    The temperature of the solution surrounding the membrane,
                    including the unit. Defaults to '25 degC' if omitted.
        Notes
        -----
        The Manning Parameter is defined as [#]:

        .. math:: \\psi = e^2 \\over {4 \\pi \\epsilon_0 \\epsilon k_b T b}

        where b is the average distance between fixed charge groups, calculated
        here based on the fixed charge concentration and the lattice distance.

        References
        ----------
        .. [#] J. Kamcev, D. R. Paul and B. D. Freeman, Macromolecules, 2015, acs.macromol.5b01654.
        .. [#] J. Kamcev, M. Galizia, F. M. Benedetti, E.-S. Jang, D. R. Paul, B. Freeman and G. S. Manning, Phys. Chem. Chem. Phys., 2016.

        See Also
        --------
        get_fixed_charge_density

        """
        # calculate the average volume per fixed charge groups as the reciprocal
        # of the molar concentration (times avogadro's number). Take the cube root
        # of the volume to get the average distance between molecules

        # TODO - should the units of fixed charge be mol/L membrane or per L of absorbed water?

        b = (abs(unit(self.get_fixed_charge_density()).to("mol/L")) * unit.N_A) ** (
            -1 / 3
        )
        temperature = unit(temperature)

        # calculate the manning parameters
        manning = unit.e ** 2 / (
            4
            * math.pi
            * unit.epsilon_0
            * self.get_dielectric_constant()
            * unit.k
            * temperature
            * b
        )

        return manning.to("dimensionless").magnitude

    # set output of the print() statement for the membrane
    def __str__(self):
        pass


class IXMembrane(Membrane):
    pass


class ROMembrane(Membrane):
    pass


class GasMembrane(Membrane):
    pass


class MDMembrane(Membrane):
    pass
