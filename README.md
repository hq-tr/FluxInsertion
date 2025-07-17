# Insert flux to a polynomial.

Written by Trung

## DESCRIPTION
This directory contains a routine for constructing quasihole states of a given fractional quantum Hall phase.

Given any state $\psi(\{z\})$ where $\{z\}\equiv(z_1,z_2,...,z_N)$ are the holomorphic variables of $N$ electrons, a flux can be added at position $a$ by

$$\prod_{i=1}^{N}(z_i-a)\psi(\{z\})$$

Each many-body state can be expanded in the Fock basis (also called the monomial basis since in first quantization each basis state is a monomial). 

## USAGE
FQH_quasihole.jl


	FQH_state_v2.jl
	Density.jl
	Misc.jl

0. Before running, prepare a file that contains the jack polynomial of the ground state
	- ***IMPORTANT*** this state must be a pure, unnormalized jack polynomial.

	- A sample file is provided, named "7e19". It is the ground state of the Laughlin state with 7 electrons.

	- Also make sure Julia is installed. Julia is available for all users on Workstation#1. For others, check with

	julia --version


1. Run the routine with

	julia FluxInsertion.jl

The program will prompt the user to input the parameters from keyboard:

	Working on the disk(1) or sphere(2)?  

→ Input a number, 1 or 2, accordingly

	File name of the jack polynomial ground state: 

→ Input the file name prepared in Step 0.

	How many fluxes to add?  

→ Number of added flux(es).

Next, the program will prompt for the position of each flux. The position should be input as two floats in the same line, separated by space.

- If disk (1) is chosen, the two numbers are X and Y coordinates of the flux

- If sphere (2) is chosen, the two numbers are theta (polar angle) and phi (azimuthal angle)


2. The state will be normalized on the appropriate geometry based on the choice disk(1) or sphere(2).

	Output state is named `<fname>_flux_<geom> where <fname>` is name of the ground state file and `<geom>` is either `dsk` or `sph` depending on the geometry.

