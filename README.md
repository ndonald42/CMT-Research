# Condensed Matter Theory Research
Undergraduate research in condensed matter theory. The focus was on topological entanglement entropy in Kitaev honeycomb models and Ising models.
# File 1: Kitaev.py
This Python file constructs the Kitaev honeycomb model for a desired lattice size. The Hamiltonian is then constructed and used to find the topological entanglement entropy and the partial Von Neumann entropies.
# File 2: Kitaev.nb
This Mathematica notebook constructs the Kitaev honeycomb model for a desired lattic size. The Hamiltonian is then constructed and diagonalizaed. The partition function, internal energy, specific heat, and entropy are all constructed and some are plotted for different coupling values. The magnetization is then constructed and plotted for several values of the coupling as well as the magnetic susceptibility.
# File 3: Ising Model.nb
This Mathematica notebook constructs the Ising model for a desired chain length. The Hamiltonian is then constructed and diagonalizaed. The partition function, internal energy, specific heat, and entropy are all constructed and some are plotted for different coupling values. The magnetization is then constructed and plotted for several values of the coupling as well as the magnetic susceptibility.
# File 4: Kitaev_Contour.py
This Python file imports large amounts of locally stored data which was originally generated using C++ code that was sent to a cluster. It constructs a contour plot to describe the phase transitions in a 24-site Kitaev model for different coupling values and magnetic field strength by looking at the topological entanglement entropy. It then looks at several cuts accross this diagram for fixed coupling values to compare curves.
# File 5: Kitaev_Mag_Sus.py
This Python imports large amounts of locally stored data which was originally generated using C++ code that was sent to a cluster. It plots the magnetic susceptibility and magentization for a 24-site Kitaev honeycomb model for different fixed values of magnetic field strength and coupling.
