// this file is distributed under 
// GPL license
#ifndef _________SYSTEMATIC______H________
#define _________SYSTEMATIC______H________
#include <memory>
#include <functional>
#include <math_h/sigma3.h>
const size_t
he3eta_cut_left=0,
he3eta_cut_right=he3eta_cut_left+1,
ppn_fit_range=he3eta_cut_right+1
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
        return x+d;
    }
    template<class FUNC>
    inline SystematicError(FUNC func):m_func([func](const double&x){return func(x);}){}
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
        return x+d;
    }
    template<class FUNC>
    inline SystematicError(FUNC func):___m_func([func](const double&x){return SystematicError<indices...>([&x,func](auto... a){return func(x,a...);}).get();}),
                                         m_func([func](const double&x){return SystematicError<indices...>([&x,func](auto... a){return func(x,a...);})();}){}
    inline ~SystematicError(){}
};

#endif
