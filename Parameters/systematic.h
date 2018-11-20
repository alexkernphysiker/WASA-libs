// this file is distributed under 
// GPL license
#ifndef _________SYSTEMATIC______H________
#define _________SYSTEMATIC______H________
#include <memory>
#include <functional>
#include <map>
#include <list>
#include <utility>
#include <vector>
#include <math_h/sigma3.h>
const size_t
he3eta_cut_left=0,
he3eta_cut_right=he3eta_cut_left+1,
ppn_fit_range=he3eta_cut_right+1,
bound_state_reaction_index=ppn_fit_range+1,
upper_limit_fit_power=bound_state_reaction_index+1,
upper_limit_left=upper_limit_fit_power+1,
upper_limit_right=upper_limit_left+1
;
const MathTemplates::value<>&Parameter(size_t i);
struct SystematicParamRec{
	MathTemplates::value<> changed_value;
	double delta_value,delta_sigma;
};

template<size_t index, size_t... indices>class SystematicError;
template<size_t index>
class SystematicError<index>{
    template<size_t i,size_t... ii>friend class SystematicError;
private:
    std::function<MathTemplates::Uncertainties<2>(const double&)> m_func;
    mutable std::vector<std::pair<SystematicParamRec,SystematicParamRec>> m_param_rec;
protected:
    inline double get()const{
        const auto&P=Parameter(index);
        return m_func(P.val()).val();
    }
public:
    inline MathTemplates::Uncertainties<2>operator()()const{
        const auto&P=Parameter(index);
        const auto x=m_func(P.val());
        const auto xm=m_func(P.min()),xp=m_func(P.max());
        const auto d=(abs(xm.val()-x.val())+abs(xp.val()-x.val()))/2.0;
	m_param_rec.clear();
	m_param_rec.push_back(std::make_pair<SystematicParamRec,SystematicParamRec>(
		{
		    .changed_value=MathTemplates::take_uncertainty_component<1>(xm),
		    .delta_value=sqrt(pow(xm.val()-x.val(),2)),
		    .delta_sigma=sqrt(abs(pow(x.template uncertainty<1>(),2)-pow(xm.template uncertainty<1>(),2)))
		},
		{
		    .changed_value=MathTemplates::take_uncertainty_component<1>(xp),
		    .delta_value=sqrt(pow(xp.val()-x.val(),2)),
		    .delta_sigma=sqrt(abs(pow(x.template uncertainty<1>(),2)-pow(xp.template uncertainty<1>(),2)))
		}
	));
        return x+MathTemplates::uncertainties(0.0,0.0,d);
    }
    inline const std::vector<std::pair<SystematicParamRec,SystematicParamRec>>&details()const{return m_param_rec;}
    inline SystematicError(std::function<MathTemplates::Uncertainties<2>(const double&)> func):
            m_func([func](const double&x){return func(x);}){}
    inline ~SystematicError(){}
};
template<size_t index, size_t... indices>
class SystematicError{
    template<size_t i,size_t... ii>friend class SystematicError;
private:
    mutable std::vector<std::pair<SystematicParamRec,SystematicParamRec>> m_param_rec;
    std::function<MathTemplates::Uncertainties<2>(const double&)> m_func;
protected:
public:
    inline MathTemplates::Uncertainties<2>operator()()const{
        const auto&P=Parameter(index);
        const auto x=m_func(P.val());
        const auto xm=m_func(P.min()),xp=m_func(P.max());
        const auto d=(sqrt(pow(xm.val()-x.val(),2))+sqrt(pow(xp.val()-x.val(),2)))/2.0;
	m_param_rec.insert(m_param_rec.begin(),std::make_pair<SystematicParamRec,SystematicParamRec>(
		{
		    .changed_value=MathTemplates::take_uncertainty_component<1>(xm),
		    .delta_value=sqrt(pow(xm.val()-x.val(),2)),
		    .delta_sigma=sqrt(abs(pow(x.template uncertainty<1>(),2)-pow(xm.template uncertainty<1>(),2)))
		},
		{
		    .changed_value=MathTemplates::take_uncertainty_component<1>(xp),
		    .delta_value=sqrt(pow(xp.val()-x.val(),2)),
		    .delta_sigma=sqrt(abs(pow(x.template uncertainty<1>(),2)-pow(xp.template uncertainty<1>(),2)))
		}
	));
        return x+MathTemplates::uncertainties(0.0,0.0,d);
    }
    inline const std::vector<std::pair<SystematicParamRec,SystematicParamRec>>&details()const{return m_param_rec;}
    template<typename... Args>
    inline SystematicError(std::function<MathTemplates::Uncertainties<2>(const double&,Args...)> func):
        m_func([func,this](const double&x){
	    SystematicError<indices...> calc([&x,func](Args... a){return func(x,a...);});
	    const auto res=calc();
	    m_param_rec.clear();
	    for(const auto&P:calc.details())m_param_rec.push_back(P);
	    return res;
	}){}
    inline ~SystematicError(){}
};

class RawSystematicError{
private:
    std::list<size_t> m_parameters;
    std::function<MathTemplates::Uncertainties<2>(const std::string&)> m_func;
    mutable std::map<size_t,std::pair<SystematicParamRec,SystematicParamRec>> m_param_rec;
public:
    template<class FUNC>
    RawSystematicError(const std::list<size_t>&params,FUNC func):m_parameters(params),m_func([func](const std::string&suffix){return func(suffix);}){}
    ~RawSystematicError();
    MathTemplates::Uncertainties<2>operator()()const;
    inline const std::map<size_t,std::pair<SystematicParamRec,SystematicParamRec>>&details()const{return m_param_rec;}
};

#endif
