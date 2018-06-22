// this file is distributed under 
// GPL license
#ifndef _________PARAMETERS______H________
#define _________PARAMETERS______H________
#include <string>
enum ParameterMode{normal,up,down};
void ChangedParameter(size_t index,ParameterMode mode);
const double&getParameter(size_t index);
#endif