import mixture


# etoh = mixture.Compound(
#     'Ethanol',
#     0.4,
#     (8.1122, 1592.864, 226.184),
#     {'CH3': 1, 'CH2': 1, 'OH': 1}
# )
#
# h2o = mixture.Compound(
#     'Water',
#     0.6,
#     (8.07131, 1730.63, 233.426),
#     {'H2O': 1}
# )
#
# vodka = mixture.VLE([etoh, h2o], 25)
#
# vodka.set_pressure(760)

benzene = mixture.Compound(
    name='Benzene',
    moles=1,
    antoine_coefficients=(6.87987, 1196.76, 219.161),
    groups={'ACH': 6}
)

hexane = mixture.Compound(
    name='Hexane',
    moles=1,
    antoine_coefficients=(6.91058, 1189.64, 226.28),
    groups={'CH3': 2, 'CH2': 4}
)

phase_alpha = mixture.VLE([benzene, hexane], 18.85)

print('hydrocarbon phase:')
print(phase_alpha)

water = mixture.Compound(
    name='Water',
    moles=1.0,
    antoine_coefficients=(8.07131, 1730.63, 233.426),
    groups={'H2O': 1}
)

phase_beta = mixture.VLE([water], 18.85)

print('water phase:')
print(phase_beta)

total_pressure = 0
for phase in phase_alpha, phase_beta:
    for compound in phase.compounds:
        total_pressure += compound.gamma * compound.x(phase) * compound.p_sat(18.85)
print('Pressure: {}'.format(total_pressure))
for phase in phase_alpha, phase_beta:
    for compound in phase.compounds:
        print('y {}: {}'.format(compound.name, compound.gamma * compound.x(phase) * compound.p_sat(18.85) / total_pressure))