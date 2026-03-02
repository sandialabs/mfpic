class Species:
    def __init__(self, charge=0.0, mass=1.0, charge_over_mass=None, specific_heat_ratio=5.0 / 3.0):
        self.charge = charge
        self.mass = mass
        self.specific_heat_ratio = specific_heat_ratio

        if charge_over_mass is None:
            charge_over_mass = charge / mass
        self.charge_over_mass = charge_over_mass
