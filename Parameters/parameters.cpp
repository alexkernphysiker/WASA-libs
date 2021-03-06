// this file is distributed under 
// GPL license
#include <math.h>
#include <vector>
#include <math_h/sigma.h>
#include "parameters.h"
#include "systematic.h"
using namespace std;
using namespace MathTemplates;
const vector<value<>> m_parameters{
    //beam momentum correction
    {0.0040,0.0002},
    //3Heeta cut height
    {0.01,0.0005},
    //ppn analysis
    {35.,1.},{42.,1.},//angles
    {-18.,1.},{-9.,1.},//time
    //3He ngamma analysis
    {0.025,0.002},//energy threshold
    {15.,2.},{0.,2.},{30.,2.},//time
    {60.,2.},//eta theta
    {0.51,0.002},//he missing mass
    {2.70,0.01},{3.00,0.01},//gamma missing mass
    //2gamma specific
    {0.45,0.01},{0.65,0.01},//gamma invariant mass
    //6gamma specific
    {0.05,0.005},//3pi0
    {0.35,0.01},{0.65,0.01},//gamma invariant mass
    {4.5,0.2}//he3 theta cut
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
const vector<value<>> m_parameters_local{
    {0.543,0.001},{0.553,0.001},//he3eta fit parameters
    {60,5},//ppn fit
    {1,1},//bound state reaction number )
    {1,1},//upper limit fit:background power
    {0,2.5},//upper limit left
    {0,2.5}//upper limit right
};
const value<>&Parameter(size_t i){return m_parameters_local[i];}

RawSystematicError::~RawSystematicError(){}
Uncertainties<2>RawSystematicError::operator()(bool filter)const{
        MathTemplates::Uncertainties<2> X=m_func("_");
        for(const size_t index:m_parameters){
            const std::string I=(index<10)?"0"+std::to_string(index):std::to_string(index);
            const auto xm=m_func(I+"-"),xp=m_func(I+"+");
            const auto d=(sqrt(pow(xp.val()-X.val(),2))+sqrt(pow(xm.val()-X.val(),2)))/2.0;
	    SystematicParamRec M,P;
	    M.changed_value=take_uncertainty_component<1>(xm);
	    M.delta_value=sqrt(pow(xm.val()-X.val(),2));
	    M.delta_sigma=sqrt(sqrt(pow(pow(X.uncertainty<1>(),2)-pow(xm.uncertainty<1>(),2),2)));
	    M.delta_sigma_sq=pow(pow(X.uncertainty<1>(),2)-pow(xm.uncertainty<1>(),2),2);
	    P.changed_value=take_uncertainty_component<1>(xp);
	    P.delta_value=sqrt(pow(xp.val()-X.val(),2));
	    P.delta_sigma=sqrt(sqrt(pow(pow(X.uncertainty<1>(),2)-pow(xp.uncertainty<1>(),2),2)));
	    P.delta_sigma_sq=pow(pow(X.uncertainty<1>(),2)-pow(xp.uncertainty<1>(),2),2);
	    m_param_rec[index]=pair<SystematicParamRec,SystematicParamRec>(M,P);
	    if((!filter)||(M.delta_sigma==0)||(P.delta_sigma==0)||((M.delta_value/M.delta_sigma)>1)||((P.delta_value/P.delta_sigma)>1))
		X+=MathTemplates::uncertainties(0.0,0.0,d);
        }
        return X;
}
RawSystematicError2::~RawSystematicError2(){}
Uncertainties<3>RawSystematicError2::operator()(bool filter)const{
        MathTemplates::Uncertainties<3> X=m_func("_");
        for(const size_t index:m_parameters){
            const std::string I=(index<10)?"0"+std::to_string(index):std::to_string(index);
            const auto xm=m_func(I+"-"),xp=m_func(I+"+");
            const auto d=(sqrt(pow(xp.val()-X.val(),2))+sqrt(pow(xm.val()-X.val(),2)))/2.0;
	    SystematicParamRec M,P;
	    M.changed_value=take_uncertainty_component<1>(xm);
	    M.delta_value=sqrt(pow(xm.val()-X.val(),2));
	    M.delta_sigma=sqrt(sqrt(pow(pow(X.uncertainty<1>(),2)-pow(xm.uncertainty<1>(),2),2)));
	    M.delta_sigma_sq=pow(pow(X.uncertainty<1>(),2)-pow(xm.uncertainty<1>(),2),2);
	    P.changed_value=take_uncertainty_component<1>(xp);
	    P.delta_value=sqrt(pow(xp.val()-X.val(),2));
	    P.delta_sigma=sqrt(sqrt(pow(pow(X.uncertainty<1>(),2)-pow(xp.uncertainty<1>(),2),2)));
	    P.delta_sigma_sq=pow(pow(X.uncertainty<1>(),2)-pow(xp.uncertainty<1>(),2),2);
	    if(xp>0){
		    m_param_rec[index]=pair<SystematicParamRec,SystematicParamRec>(M,P);
		    if((!filter)||(M.delta_sigma==0)||(P.delta_sigma==0)||((M.delta_value/M.delta_sigma)>1)||((P.delta_value/P.delta_sigma)>1))
			X+=MathTemplates::uncertainties(0.0,0.0,0.0,P.delta_value);
	    }else{
		    m_param_rec[index]=pair<SystematicParamRec,SystematicParamRec>(P,M);
		    if((!filter)||(M.delta_sigma==0)||(P.delta_sigma==0)||((M.delta_value/M.delta_sigma)>1)||((P.delta_value/P.delta_sigma)>1))
			X+=MathTemplates::uncertainties(0.0,0.0,P.delta_value,0.0);
	    }
	    if(xm>X)
		X+=MathTemplates::uncertainties(0.0,0.0,0.0,M.delta_value);
	    else
		X+=MathTemplates::uncertainties(0.0,0.0,M.delta_value,0.0);
        }
        return X;
}
