// this file is distributed under 
// GPL license
#ifndef _________SYSTEMATIC______H________
#define _________SYSTEMATIC______H________
#include <memory>
#include <functional>
#include <map>
#include <list>
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

template<size_t index, size_t... indices>class SystematicError;
template<size_t index>
class SystematicError<index>{
    template<size_t i,size_t... ii>friend class SystematicError;
private:
    std::function<MathTemplates::Uncertainties<2>(const double&)> m_func;
    mutable std::vector<MathTemplates::Uncertainties<2>> m_upper;
    mutable std::vector<MathTemplates::Uncertainties<2>> m_lower;
    mutable std::vector<double> m_contrib;
    
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
        const auto d=abs(xm.val()-xp.val())/2.0;
	m_contrib.clear();m_upper.clear();m_lower.clear();
	m_contrib.push_back(d);
	m_upper.push_back(xp);
	m_lower.push_back(xm);
        return x+MathTemplates::uncertainties(0.0,0.0,d);
    }
    const std::vector<double>&contrib()const{return m_contrib;}
    const std::vector<MathTemplates::Uncertainties<2>>&upper()const{return m_upper;}
    const std::vector<MathTemplates::Uncertainties<2>>&lower()const{return m_lower;}
    inline SystematicError(std::function<MathTemplates::Uncertainties<2>(const double&)> func):
            m_func([func](const double&x){return func(x);}){}
    inline ~SystematicError(){}
};
template<size_t index, size_t... indices>
class SystematicError{
    template<size_t i,size_t... ii>friend class SystematicError;
private:
    mutable std::vector<double> m_contrib;
    mutable std::vector<MathTemplates::Uncertainties<2>> m_upper;
    mutable std::vector<MathTemplates::Uncertainties<2>> m_lower;
    std::function<MathTemplates::Uncertainties<2>(const double&)> m_func;
protected:
public:
    inline MathTemplates::Uncertainties<2>operator()()const{
        const auto&P=Parameter(index);
        const auto x=m_func(P.val());
        const auto xm=m_func(P.min()),xp=m_func(P.max());
        const auto d=abs(xm.val()-xp.val())/2.0;
	m_contrib.insert(m_contrib.begin(),d);
	m_upper.insert(m_upper.begin(),xp);
	m_lower.insert(m_lower.begin(),xm);
        return x+MathTemplates::uncertainties(0.0,0.0,d);
    }
    const std::vector<double>&contrib()const{return m_contrib;}
    const std::vector<MathTemplates::Uncertainties<2>>&upper()const{return m_upper;}
    const std::vector<MathTemplates::Uncertainties<2>>&lower()const{return m_lower;}
    template<typename... Args>
    inline SystematicError(std::function<MathTemplates::Uncertainties<2>(const double&,Args...)> func):
        m_func([func,this](const double&x){
	    SystematicError<indices...> calc([&x,func](Args... a){return func(x,a...);});
	    const auto res=calc();
	    m_contrib.clear();m_upper.clear();m_lower.clear();
	    for(const auto&p:calc.contrib())m_contrib.push_back(p);
	    for(const auto&p:calc.upper())m_upper.push_back(p);
	    for(const auto&p:calc.lower())m_lower.push_back(p);
	    return res;
	}){}
    inline ~SystematicError(){}
};

class RawSystematicError{
private:
    std::list<size_t> m_parameters;
    std::function<MathTemplates::Uncertainties<2>(const std::string&)> m_func;
    std::map<size_t,double> m_contrib;
    std::map<size_t,MathTemplates::Uncertainties<2>> m_values_up,m_values_down;
public:
    template<class FUNC>
    RawSystematicError(const std::list<size_t>&params,FUNC func):m_parameters(params),m_func([func](const std::string&suffix){return func(suffix);}){}
    ~RawSystematicError();
    MathTemplates::Uncertainties<2>operator()();
    const double&contrib(size_t p)const;
    const MathTemplates::Uncertainties<2>&upper(size_t p)const;
    const MathTemplates::Uncertainties<2>&lower(size_t p)const;
};

#endif
