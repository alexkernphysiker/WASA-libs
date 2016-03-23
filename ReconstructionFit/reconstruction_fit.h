// this file is distributed under 
// MIT license
#ifndef JUDIQVAJ
# define JUDIQVAJ
#include <iostream>
#include <string>
#include <sstream>
#include <gnuplot_wrap.h>
#include <math_h/error.h>
#include <math_h/structures.h>
#include <Genetic/fit.h>
namespace SimulationDataProcess{
	using namespace std;
	using namespace Genetic;
	using namespace MathTemplates;
	using namespace GnuplotWrap;
	string SimulationDataPath();
	template<class FITFUNC>
	void He3ForEtaFit(const string&&reconstructionname,const vector<value<double>>&&E_d_bins,const vector<value<double>>&&E_k_bins,const shared_ptr<IInitialConditions>init,RANDOM&R){
		auto params_shown=make_pair(0,2);
		vector<value<double>> theta_bins=BinsByStep(0.10,0.002,0.13);
		vector<Distribution2D<double>> E_sp2;
		for(const value<double>&bin:theta_bins)
			E_sp2.push_back(Distribution2D<double>(E_d_bins,E_k_bins));
		cout<<theta_bins.size()<<" theta bins"<<endl;
		{
			ifstream file;
			file.open(SimulationDataPath()+reconstructionname+".simulation.txt");
			if(file){
				cout<<"reading..."<<endl;
				string line;
				while(getline(file,line)){
					istringstream str(line);
					double Edep,theta,Ekin;
					str>>Edep>>theta>>Ekin;
					for(size_t i=0,n=theta_bins.size();i<n;i++)
						if(theta_bins[i].contains(theta))
							E_sp2[i].Fill(Edep,Ekin);
				}
				file.close();
				cout<<"done."<<endl;
			}else
				throw Exception<ifstream>("No input data");
		}
		cout<<"Init1"<<endl;
		auto points=make_shared<FitPoints>();
		for(size_t i=0,n=theta_bins.size();i<n;i++)
			E_sp2[i].FullCycle([&theta_bins,&E_sp2,&points,i](const point3d<value<double>>&P){
				if(!P.Z().contains(0))
					points<<Point({P.X().val(),theta_bins[i].val()},P.Y().val(),P.Z().val());
			});
		cout<<points->size()<<" points"<<endl;
		FitFunction<DifferentialMutations<>,FITFUNC,SumWeightedSquareDiff> fit(points);
		cout<<"Init2"<<endl;
		fit.Init(40*FITFUNC::ParamCount,init,R);
		cout<<"population "<<fit.PopulationSize()<<endl;
		cout<<"Fitting"<<endl;
		while(
			(!fit.AbsoluteOptimalityExitCondition(0.000001))&&
			(!fit.RelativeOptimalityExitCondition(0.0001))
		){
			fit.Iterate(R);
			cout<<fit.Optimality()<<"<S<"<<fit.Optimality(fit.PopulationSize()-1)<<"     \r";
		}
		cout<<"done.                                                                            "<<endl;
		{
			ofstream out;
			out.open(SimulationDataPath()+reconstructionname+".fit.txt");
			if(out){
				out<<fit;
				out.close();
			}else
				throw Exception<ofstream>("Cannot write output");
		}
		for(size_t i=0,n=theta_bins.size();i<n;i++){
			PlotHist2d<double>(sp2).Distr(E_sp2[i],"Theta=["+to_string(theta_bins[i].min())+":"+to_string(theta_bins[i].max())+"]");
			PlotHist2d<double>(normal).Distr(E_sp2[i],"Theta=["+to_string(theta_bins[i].min())+":"+to_string(theta_bins[i].max())+"]");
			double max=0;E_sp2[i].FullCycle([&max](const point3d<value<double>>&P){if(P.Z().val()>max)max=P.Z().val();});
			vector<point<double>> lo,hi;
			E_sp2[i].FullCycle([max,&lo,&hi](const point3d<value<double>>&P){
				if(!P.Z().contains(0)){
					auto p=point<double>(P.X().val(),P.Y().val());
					if(P.Z().val()>(max/2.0))
						hi.push_back(p);
					else
						lo.push_back(p);
				}
			});
			Plot<double>().Points(lo).Points(hi)
			.Line(SortedPoints<double>([&fit,i,&theta_bins](double x)->double{return fit({x,theta_bins[i].val()});},ChainWithStep(E_d_bins[0].val(),0.001,E_d_bins[E_d_bins.size()-1].val())));
		}
	}
};
#endif