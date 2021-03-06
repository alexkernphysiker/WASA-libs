// this file is distributed under 
// GPL license
#ifndef JUDIQVAJ
# define JUDIQVAJ
#include <iostream>
#include <string>
#include <sstream>
#include <math_h/tabledata.h>
#include <Genetic/fit.h>
#include <Genetic/initialconditions.h>
namespace SimulationDataProcess{
    void ForwardEkinReconstructionFit(
	const std::string&&reconstructionname,
	const std::shared_ptr<Genetic::IParamFunc>func,
	const MathTemplates::SortedChain<MathTemplates::value<double>>&&E_d_bins,
	const MathTemplates::SortedChain<MathTemplates::value<double>>&&E_k_bins,
	const std::shared_ptr<Genetic::IInitialConditions>init
    );
};
#endif
