// this file is distributed under 
// MIT license
#ifndef JUDIQVAJ
# define JUDIQVAJ
#include <iostream>
#include <string>
#include <sstream>
#include <Genetic/fit.h>
#include <Genetic/initialconditions.h>
namespace SimulationDataProcess{
	void He3ForEtaFit(
		const std::string&&reconstructionname,
		const std::shared_ptr<Genetic::IParamFunc>func,
		const std::vector<MathTemplates::value<double>>&&E_d_bins,
		const std::vector<MathTemplates::value<double>>&&E_k_bins,
		const std::shared_ptr<Genetic::IInitialConditions>init,
		Genetic::RANDOM&R
	);
};
#endif