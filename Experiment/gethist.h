// this file is distributed under 
// GPL license
#ifndef DAIHADDG
#define DAIHADDG
#include <string>
#include <functional>
#include <vector>
#include <math_h/tabledata.h>
namespace ROOT_data{
	std::pair<double,double> PresentRuns(const std::string&reaction,const std::string&param_suffix="_");
	enum histsource{MC,DATA};
	MathTemplates::hist<double> Hist(
            histsource src, const std::string&reaction, 
            const std::vector<std::string>&path,const std::string&histname,
            const std::string&param_suffix="_"
        );
	MathTemplates::hist2d<double> Hist2d(
            histsource src, const std::string&reaction, 
            const std::vector<std::string>&path,const std::string&histname,
            const std::string&param_suffix="_"
        );
};
#endif
