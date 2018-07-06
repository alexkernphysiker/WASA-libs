// this file is distributed under 
// GPL license
#ifndef _________PARAMETERS______H________
#define _________PARAMETERS______H________
#include <string>
enum ParameterMode{param_normal,param_up,param_down};
void ChangedParameter(size_t index,ParameterMode mode);
double getParameter(size_t index);
size_t ParametersCount();
//Parameter index constants
const size_t
pbeam_corr=0,
he3_cut_h=pbeam_corr+1,
ppn_th1=he3_cut_h+1,
ppn_th2=ppn_th1+1,
ppn_t1=ppn_th2+1,
ppn_t2=ppn_t1+1,
gamma_E_thr=ppn_t2+1,
time_dt=gamma_E_thr+1,
time_t1=time_dt+1,
time_t2=time_t1+1,
eta_theta_thr=time_t2+1,
he3mm_cut=eta_theta_thr+1,
gamma_mm_lo=he3mm_cut+1,
gamma_mm_hi=gamma_mm_lo+1,
gamma_im_lo=gamma_mm_hi+1,
gamma_im_hi=gamma_im_lo+1,
three_pi0=gamma_im_hi+1,
gamma_im_lo6=three_pi0+1,
gamma_im_hi6=gamma_im_lo6+1,
he3_theta_cut=gamma_im_hi6+1;
#endif
