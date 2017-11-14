import numpy as np
a_ij_matrix = np.load('aij.npy', mmap_mode='r')

group_names = ['CH3', 'CH2', 'CH', 'ACH', 'AC', 'ACCH3', 'ACCH2', 'OH',
       'CH3OH', 'H2O', 'ACOH', 'CH3CO', 'CHO', 'CH3COO', 'CH2NH2', 'ACNH2',
       'COOH', 'CCL2', 'CCL3', 'CON(CH2)2']


class Group():
    def __init__(self, name):
        self.name = name
        self.index = group_names.index(name)
        self.a_ij = np.array(a_ij_matrix[self.index][:20])
        self.R = np.float(a_ij_matrix[self.index][20])
        self.Q = np.float(a_ij_matrix[self.index][21])

    def __str__(self):
        return '{} -------------------------\n\n'.format(self.name) + \
               'R: {}\nQ: {}\n\n'.format(self.R, self.Q) + \
               ''.join(['a_ij {}: {}\n'.format(group, a_ij) for group, a_ij in zip(group_names, self.a_ij)])


    def a_mj(self, other):
        return self.a_ij[group_names.index(other)]

    def psi(self, other, T):
        return np.exp(-self.a_mj(other) / T)

    @staticmethod
    def X(x, compounds):
        return x / sum([compound.x * len(compound.groups) for compound in compounds])

    def theta(self, x, compounds):
        numerator = self.X(x, compounds) * self.Q
        denominator = 0
        for compound in compounds:
            for group in compound.groups:
                denominator += group.X(compound.x, compounds) * group.Q
        return numerator / denominator

    def ln_gamma(self, compounds, T):
        sigma_1 = 0
        for compound in compounds:
            for group in compound.groups:
                sigma_1 += group.theta(compound.x, compounds) * group.psi(self.name, T)

        sigma_2 = 0
        for compound_i in compounds:
            for group_i in compound_i.groups:
                numerator = group_i.theta(compound_i.x, compounds) * self.psi(group_i.name, T)
                denominator = 0
                for compound_j in compounds:
                    for group_j in compound_j.groups:
                        denominator += group_j.theta(compound_j.x, compounds) * group_j.psi(group_i.name, T)
                sigma_2 += numerator / denominator

        return self.Q * (1 - np.log(sigma_1) - sigma_2)
