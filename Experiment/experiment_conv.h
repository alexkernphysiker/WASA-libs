// this file is distributed under 
// GPL license
#ifndef AXWBNBYL
#	define AXWBNBYL
#include <math_h/functions.h>
#define CSTR(A) (const_cast<char*>(string(A).c_str()))
const double p_beam_low=1.426;
const double p_beam_hi=1.635;
struct trigger{unsigned char number; unsigned long scaling;};
const trigger trigger_he3_forward={.number=10,.scaling=1};
const trigger trigger_gammas_central={.number=7,.scaling=10};
const trigger trigger_elastic1={.number=17,.scaling=4000};
const trigger trigger_elastic2={.number=21,.scaling=4000};
#endif
