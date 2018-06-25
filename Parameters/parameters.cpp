// this file is distributed under 
// GPL license
#include <vector>
#include <math_h/sigma.h>
#include "parameters.h"
using namespace std;
using namespace MathTemplates;
const vector<value<>> m_parameters{
    //beam momentum correction
    {0.0040,0.0004},
    //3Heeta analysis
    {0.01,0.0005},
    //ppn analysis
    {35.,1.},{42.,1.},//angles
    {-18.,1.},{-9.,1.},//time
    //3He ngamma analysis
    {0.005,0.001},//energy threshold
    {25.,2.},{-5.,2.},{40.,2.},//time
    {60.,5.},//eta theta
    {0.51,0.03},//he missing mass
    {2.70,0.03},{3.00,0.03},//gamma missing mass
    {0.45,0.04},{0.65,0.03},//gamma invariant mass
    //2gamma specific
    {-0.10,0.1},//cos theta
    //6gamma specific
    {0.1,0.01}//3pi0
};
vector<ParameterMode> m_correction;
void init_corrections(){
    if(m_correction.size()==0){
        for(size_t i=0;i<m_parameters.size();i++)
            m_correction.push_back(param_normal);
    }
}
void ChangedParameter(size_t index,ParameterMode mode){
    if(index>=m_parameters.size())
        throw Exception<vector<value<>>,0>("Parameter index out of range");
    init_corrections();
    m_correction[index]=mode;
}
double getParameter(size_t index){
    if(index>=m_parameters.size())
        throw Exception<vector<value<>>,0>("Parameter index out of range");
    init_corrections();
    const auto&P=m_parameters[index];
    switch(m_correction[index]){
        case param_normal:
            return P.val();
        case param_up:
            return P.max();
        case param_down:
            return P.min();
        default:
            throw Exception<vector<value<>>,1>("Parameter correction mode out of range");
    }
}
size_t ParametersCount(){
    return m_parameters.size();
}
