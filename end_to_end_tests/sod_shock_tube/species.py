class Species:
    def __init__(self, charge = 0., mass = 1., charge_over_mass = None, specific_heat_ratio = 5. / 3.):
        self.charge = charge
        self.mass = mass
        self.specific_heat_ratio = specific_heat_ratio

        if charge_over_mass is None:
            charge_over_mass = charge / mass
        self.charge_over_mass = charge_over_mass
