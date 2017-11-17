WASA analysis libraries
=======================
Libraries used for analysis of data obtained from the experiment WASA-at-COSY on searching eta-mesic 3He in May 2014.
This repository contains the part of code that is common for both repositories with data analysis software:

    https://github.com/alexkernphysiker/WASA-analysis
    https://github.com/alexkernphysiker/WASA-analysis


All files are distributed under GPL license.


Other required repositories
===========================

If you have your git repository with cmake project you should add the following submodules

    git submodule add https://github.com/alexkernphysiker/math_h.git
    git submodule add https://github.com/alexkernphysiker/FitGen.git
    git submodule add https://github.com/alexkernphysiker/WASA-libs.git
    git submodule update --init --recursive

Then add to CMakeLists.txt

    add_subdirectory(math_h)
    add_subdirectory(FitGen)
    add_subdirectory(WASA-libs)
    include_directories(${MATH_H_INC})
    include_directories(${FITGEN_INC})
    include_directories(${WASA_LIBS_INC})

Then commit your changes


Directories and files
=====================
	Kinematics
library with classes for calculation of reactions' kinematics

	Experiment
library with commont routines connected with data access and global constants used in current experiment

	ReconstructionFit
library that provides reconstruction fit

