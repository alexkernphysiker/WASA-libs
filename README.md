WASA analysis libraries
=======================
Libraries used for analysis of data obtained from the experiment WASA-at-COSY on searching eta-mesic 3He in May 2014.
All files are distributed under GPL license


Required software
=================
	ROOT
framework for calulations

	Pluto
library for Monte Carlo simulation of reaction products. Uses ROOT.

	gnuplot
software for plotting. Is used by some applications performing final analysis.


Needed environment variables
============================
	ROOTSYS
path where ROOT is installed

	PLUTOSYS
path where pluto is installed

	PLUTO_OUTPUT
path where pluto files are stored


Directories and files
=====================
	FitGen
submodule with mathematical routines

	Kinematics
library with classes for calculation of reactions' kinematics

	Experiment
library with commont routines connected with data access and global constants used in current experiment

	ReconstructionFit
library that provides reconstruction fit

