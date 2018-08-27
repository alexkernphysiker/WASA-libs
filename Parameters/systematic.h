// this file is distributed under 
// GPL license
#ifndef _________SYSTEMATIC______H________
#define _________SYSTEMATIC______H________
#include <memory>
#include <functional>
#include <map>
#include <list>
#include <math_h/sigma3.h>
const size_t
he3eta_cut_left=0,
he3eta_cut_right=he3eta_cut_left+1,
ppn_fit_range=he3eta_cut_right+1,
bound_state_reaction_index=ppn_fit_range+1
;
const MathTemplates::value<>&Parameter(size_t i);

template<size_t index, size_t... indices>class SystematicError;
template<size_t index>
class SystematicError<index>{
    template<size_t i,size_t... ii>friend class SystematicError;
private:
    std::function<MathTemplates::Uncertainties<2>(const double&)> m_func;
protected:
    inline double get()const{
        const auto&P=Parameter(index);
        return m_func(P.val()).val();
    }
public:
    inline MathTemplates::Uncertainties<2>operator()()const{
        const auto&P=Parameter(index);
        const auto x=m_func(P.val());
        const double xm=m_func(P.min()).val(),xp=m_func(P.max()).val();
        const auto d=MathTemplates::uncertainties(0.0,0.0,abs(xm-xp)/2.0);
        return (d.template uncertainty<2>()>x.template uncertainty<1>())?x+d:x;
    }
    inline SystematicError(std::function<MathTemplates::Uncertainties<2>(const double&)> func):
            m_func([func](const double&x){return func(x);}){}
    inline ~SystematicError(){}
};
template<size_t index, size_t... indices>
class SystematicError{
    template<size_t i,size_t... ii>friend class SystematicError;
private:
    std::function<double(const double&)> ___m_func;
    std::function<MathTemplates::Uncertainties<2>(const double&)> m_func;
protected:
    inline double get()const{
        const auto&P=Parameter(index);
        return ___m_func(P.val());
    }
public:
    inline MathTemplates::Uncertainties<2>operator()()const{
        const auto&P=Parameter(index);
        const auto x=m_func(P.val());
        const double xm=___m_func(P.min()),xp=___m_func(P.max());
        const auto d=MathTemplates::uncertainties(0.0,0.0,abs(xm-xp)/2.0);
        return (d.template uncertainty<2>()>x.template uncertainty<1>())?x+d:x;
    }
    template<typename... Args>
    inline SystematicError(std::function<MathTemplates::Uncertainties<2>(const double&,Args...)> func):
            ___m_func([func](const double&x){return SystematicError<indices...>([&x,func](Args... a){return func(x,a...);}).get();}),
               m_func([func](const double&x){return SystematicError<indices...>([&x,func](Args... a){return func(x,a...);})();}){}
    inline ~SystematicError(){}
};

class RawSystematicError{
private:
    std::list<size_t> m_parameters;
    std::function<MathTemplates::Uncertainties<2>(const std::string&)> m_func;
    std::map<size_t,double> m_contrib;
public:
    template<class FUNC>
    RawSystematicError(const std::list<size_t>&params,FUNC func):m_parameters(params),m_func([func](const std::string&suffix){return func(suffix);}){}
    ~RawSystematicError();
    MathTemplates::Uncertainties<2>operator()();
    const double&contrib(size_t p)const;
};

#endif
