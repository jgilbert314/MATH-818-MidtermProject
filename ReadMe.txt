This code is meant to implement some of the models for quantum circuits described in the article:
	Blais, Alexandre, et al. "Circuit quantum electrodynamics." Reviews of Modern Physics 93.2 (2021): 025005.

Currently the code implements equation 70 for the paper, the master equation for a harmonic oscillator connected to a transmission line.

A brief overview of the code:
	- Workbench:
		- User interface script. Parameters are be defined here, solution is calculated using ode45, solution is visualized.
	- masterEq:
		- Defines ODE for density matrix (eq 70)
	- calcDissMat:
		- Defines dissipative terms in master equation (eq 71)
	- calcHamiltionianOsc:
		- Defines Hamiltonian for harmonic oscillator (below eq 68)