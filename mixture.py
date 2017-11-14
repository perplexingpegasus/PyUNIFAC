import numpy as np
from groups import Group


class Compound():
    def __init__(self, name, x, antoine_coefficients, groups):
        self.name = name
        self.x = x
        self.A, self.B, self.C = antoine_coefficients
        self.groups = []
        for key, value in groups.items():
            for i in range(value):
                self.groups.append(Group(key))
        self.q = sum([group.Q for group in self.groups])
        self.r = sum([group.R for group in self.groups])
        self.gamma = 1

    def y(self, pressure, temp):
        return self.x * self.gamma * self. p_sat(temp) / pressure

    def p_sat(self, temp, units='mmHG'):
        pressure = 10 ** (self.A - self.B / (temp + self.C))
        if units == 'mmHG':
            return pressure

    def phi(self, mixture):
        return self.x * self.r / sum([compound.x * compound.r for compound in mixture.compounds])

    def theta(self, mixture):
        return self.x * self.q / sum([compound.x * compound.q for compound in mixture.compounds])

    def ln_gamma_comb(self, mixture):
        phi = self.phi(mixture)
        theta = self.theta(mixture)
        return np.log(phi / self.x) + (1 - phi / self.x) - 5 * self.q * (np.log(phi / theta) + (1 - phi / theta))

    def ln_gamma_res(self, mixture):
        return sum(
            [
                group.ln_gamma(mixture.compounds, mixture.tempK) - group.ln_gamma([self], mixture.tempK)
                for group in self.groups
            ]
        )

    def get_activity_coefficient(self, mixture):
        self.gamma = np.exp(self.ln_gamma_comb(mixture) + self.ln_gamma_res(mixture))
        return self.gamma



class Mixture():
    def __init__(self, compounds, temp):
        self.tempC = temp
        if not all([type(compound) == Compound for compound in compounds]):
            raise TypeError("Parameter compounds must be an iterable with elements of class Compound")
        else:
            self.compounds = compounds

    @property
    def tempK(self):
        return self.tempC + 273.15

    def normalize_mole_fractions(self):
        mole_frac = sum([compound.x for compound in self.compounds])
        if mole_frac != 1:
            for compound in self.compounds:
                compound.x = compound.x / mole_frac

    def __add__(self, other):
        if type(other) == Compound:
            self.normalize_mole_fractions()
            self.compounds += other

    def __radd__(self, other):
        if type(other) == Compound:
            self.normalize_mole_fractions()
            self.compounds += other

    def __iadd__(self, other):
        if type(other) == Compound:
            self.normalize_mole_fractions()
            self.compounds += other


class VLE(Mixture):
    def __init__(self, compounds, temp):
        super().__init__(compounds, temp)
        self.pressure = self.pressure_calc()

    def __str__(self):
        string = 'Pressure: {} mmHg, Temperature {} C' \
                 '\n-------------------------------------------------------' \
                 '\n'.format(self.pressure, self.tempC)
        count = 1
        for compound in self.compounds:
            string += '{}) {} - \n' \
                      'P_sat: {} mmHg, x: {}, y: {}, gamma: {}' \
                      '\n\n'.format(
                count,
                compound.name,
                compound.p_sat(self.tempC),
                compound.x,
                compound.y(self.pressure, self.tempC),
                compound.gamma
                )
            count += 1
        return string

    def set_temperature(self, value):
        self.tempC = value
        self.pressure = self.pressure_calc()
        print('Pressure: {} mmHg\nTemperature: {} C'.format(self.pressure, self.tempC))

    def set_pressure(self, value, increment=10):
        while np.abs(self.pressure_calc() - value) > 0.000000001:
            initial = self.pressure_calc() - value
            self.tempC += increment
            final =  self.pressure_calc() - value
            if initial * final < 0:
                increment *= -1
            if np.abs(initial) < np.abs(final):
                increment /= 10
        self.pressure = self.pressure_calc()
        print('Pressure: {} mmHg\nTemperature: {} C'.format(self.pressure, self.tempC))

    def pressure_calc(self):
        return sum(
            [
                compound.p_sat(self.tempC) * compound.get_activity_coefficient(self) * compound.x
                for compound in self.compounds
            ]
        )


class LLE(Mixture):
    def __init__(self, compounds, temp):
        super().__init__(compounds, temp)