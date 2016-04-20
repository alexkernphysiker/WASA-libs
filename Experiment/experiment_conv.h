// this file is distributed under 
// MIT license
#ifndef AXWBNBYL
#	define AXWBNBYL
#include <math_h/functions.h>
inline double NormPhi(double p){
	const double twopi=2.0*MathTemplates::PI();
	double phi=p;
	while(phi<0)phi+=twopi;
	while(phi>=twopi)phi-=twopi;
	return phi;
}
#define ALLRUNS int runindex=45873;runindex<=46884;runindex++
const double p_beam_low=1.426;
const double p_beam_hi=1.635;
struct trigger{unsigned char number; unsigned long scaling;};
const trigger trigger_he3_forward={.number=10,.scaling=1};
const double he3_forward_correct_energy=0.015;
//ToDo: add other triggers used in analysis



#endif