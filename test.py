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

water = mixture.Compound(
    'Water',
    0.3333333333,
    (8.07131, 1730.63, 233.426),
    {'H2O': 1}
)

benzene = mixture.Compound(
    'Benzene',
    0.3333333333,
    (6.87987, 1196.76, 219.161),
    {'ACH': 6}
)

hexane = mixture.Compound(
    'Hexane',
    0.333333333333333,
    (6.91058, 1189.64, 226.28),
    {'CH3': 2, 'CH2': 4}
)

mix = mixture.VLE([water, benzene, hexane], 25)

mix.set_temperature(18.85)
mix.get_ys()