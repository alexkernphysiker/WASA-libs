// this file is distributed under 
// GPL license
#include <math.h>
#include <TObject.h>
#include <TH1F.h>
#include <TH1.h>
#include <TH2.h>
#include <TFile.h>
#include <TTree.h>
#include <TDirectoryFile.h>
#include <math_h/error.h>
#include "experiment_conv.h"
#include "gethist.h"
#include "str_get.h"
namespace ROOT_data{
	using namespace std;
	using namespace MathTemplates;
	string inputpath=ENV(PRESEL_DATA);
	string outpath=ENV(OUTPUT_PLOTS);
	
	hist<double> ReadHist(const string&filename,const vector<string>&path,const string&histname){
	    hist<double> points;
	    TFile* file=TFile::Open(filename.c_str());
	    if(file){
		TDirectoryFile* dir1=file;
		for(string name:path){
		    TDirectoryFile* dir2=dynamic_cast<TDirectoryFile*>(dir1->Get(name.c_str()));
		    if(dir2)dir1=dir2;
		    else throw Exception<TDirectoryFile>("No directory "+name);
		}
		TH1F* histogram=dynamic_cast<TH1F*>(dir1->Get(histname.c_str()));
		if(histogram){
		    for(int i=1,N=histogram->GetNbinsX();i<=N;i++){
			double y=histogram->GetBinContent(i);
			double dy=sqrt(y);
			if(dy<1.0)
			    dy=1.0;
			double x=histogram->GetBinCenter(i);
			double dx=histogram->GetBinWidth(i)/2.0;
			points<<point<value<double>>({x,dx},{y,dy});
		    }
		}else throw Exception<TH1F>("No histogram "+histname+" in file "+filename);
		file->Close();
		delete file;
	    }
	    return points;
	}
	hist2d<double> ReadHist2D(const string&filename,const vector<string>&path,const string&histname){
	    TFile* file=TFile::Open(filename.c_str());
	    if(file){
		TDirectoryFile* dir1=file;
		for(string name:path){
		    TDirectoryFile* dir2=dynamic_cast<TDirectoryFile*>(dir1->Get(name.c_str()));
		    if(dir2)dir1=dir2;
		    else throw Exception<TDirectoryFile>("No directory "+name);
		}
		TH2F* histogram=dynamic_cast<TH2F*>(dir1->Get(histname.c_str()));
		if(histogram){
		    hist2d<double> res(
			BinsByCount(histogram->GetXaxis()->GetNbins(),
				    histogram->GetXaxis()->GetXmin(),
				    histogram->GetXaxis()->GetXmax()
			),
			BinsByCount(histogram->GetXaxis()->GetNbins(),
				    histogram->GetYaxis()->GetXmin(),
				    histogram->GetYaxis()->GetXmax()
			)
		    );
		    for(int i=1,N=res.X().size();i<N;i++){
			for(int j=1,M=res.Y().size();j<M;j++){
			    double y=histogram->GetBinContent(i,j);
			    double dy=sqrt(y);
			    if(dy<1.0)dy=1.0;
			    res.Bin(i,j)={y,dy};
			}
		    }
		    file->Close();
		    return res;
		}else throw Exception<TH1F>("No histogram "+histname+" in file "+filename);
		file->Close();
		delete file;
	    }
	    return hist2d<double>({},{});
	}

	
	hist<double> Hist(histsource src, const string&reaction, const vector<string>&path,const string&histname){
	    hist<double> res;
	    switch(src){
		case MC:{
		    for(ALLMC){
			hist<double> tmp=ReadHist(inputpath+"/MC"+reaction+"-"+to_string(runindex)+".root",path,histname);
			if(tmp.size()>0){
			    if(res.size()==0)
				res=tmp;
			    else
				res.imbibe(tmp);
			}
		    }
		}break;
		case DATA:{
		    for(ALLRUNS){
			hist<double> tmp=ReadHist(inputpath+"/Data"+reaction+"_run_"+to_string(runindex)+".root",path,histname);
			if(tmp.size()>0){
			    if(res.size()==0)
				res=tmp;
			    else
				res.imbibe(tmp);
			}
		    }
		}break;
	    };
	    return res;
	}
	hist2d< double > Hist2d(histsource src, const string& reaction, const vector< string >& path, const string& histname){
	    hist2d<double> res;
	    switch(src){
		case MC:{
		    res=ReadHist2D(inputpath+"/MC"+reaction+".root",path,histname);
		    for(ALLMC){
			hist2d<double> tmp=ReadHist2D(inputpath+"/MC"+reaction+"-"+to_string(runindex)+".root",path,histname);
			if(tmp.size()>0){
			    if(res.size()==0)
				res=tmp;
			    else
				res.imbibe(tmp);
			}
		    }
		}break;
		case DATA:{
		    for(ALLRUNS){
			hist2d<double> tmp=ReadHist2D(inputpath+"/Data"+reaction+"_run_"+to_string(runindex)+".root",path,histname);
			if(tmp.size()>0){
			    if(res.size()==0)
				res=tmp;
			    else
				res.imbibe(tmp);
			}
		    }
		}break;
	    };
	    return res;
	}
	pair<double,double> PresentRuns(string&&reaction){
	    size_t allruns=0,present_runs=0;
	    for(ALLRUNS){
		allruns++;
		TFile *file=TFile::Open((inputpath+"/Data"+reaction+"_run_"+to_string(runindex)+".root").c_str());
		if(file){
		    present_runs++;
		    file->Close();
		    delete file;
		}
	    }
	    return make_pair(double(present_runs),double(allruns));
	}

};
