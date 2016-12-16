// this file is distributed under 
// MIT license
#ifndef DAIHADDG
#define DAIHADDG
#include <string>
#include <functional>
#include <vector>
#include <math_h/hists.h>
namespace ROOT_data{
	std::pair<double,double> PresentRuns(std::string&&reaction);
	enum histsource{MC,DATA};
	MathTemplates::hist<double> Hist(histsource src, const std::string&reaction, const std::vector<std::string>&path,const std::string&histname);
	MathTemplates::hist2d<double> Hist2d(histsource src, const std::string&reaction, const std::vector<std::string>&path,const std::string&histname);
};
#endif
