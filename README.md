# PyUNIFAC

Python implementation of UNIFAC activity model for non-ideal mixtures

Currently the script models vapor liquid equilibrium of any number of components. Components can be initialized with the class "Component", specifying the groups present in the molecule and the total moles of the substance. These can then be added to a "Mixture" object where the temperature and pressure can be adjusted, and temperature can be iteratively solved for at a given pressure.

In the future this script will include modelling for liquid-liquid equilibrium as well.
