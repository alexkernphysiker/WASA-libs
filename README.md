WASA libraries
=======================================
Libraries used for analysis of data obtained from the experiment WASA-at-COSY on searching eta-mesic 3He in May 2014.


Required software
=================
	ROOT - framework for calulations
	Pluto - library for Monte Carlo simulation of reaction products. Uses ROOT.
	gnuplot - software for plotting. Is used by some applications performing final analysis.


Needed environment variables
============================
	ROOTSYS - path where ROOT is installed
	PLUTOSYS - path where pluto is installed
	PLUTO_OUTPUT - path where pluto files are stored


Directories
===========
	FitGen - submodule with mathematical routines
	Kinematics - library with classes for calculation of reactions' kinematics
	Experiment - library with commont routines connected with data access and global constants used in current experiment
	ReconstructionFit - library that provides reconstruction fit
	RunPlutoApp - sources of an application than runs pluto for simulation of reactions. Requires following console arguments:
		[all/over] - specifying pbeam range (all in experiment or over eta-creation threshold)
		reaction products - separated with space