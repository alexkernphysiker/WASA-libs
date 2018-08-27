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
			const auto Y=std_error(histogram->GetBinContent(i));
			const double x=histogram->GetBinCenter(i);
			const double dx=histogram->GetBinWidth(i)/2.0;
			points<<make_point(value<>(x,dx),Y);
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
			    res.Bin(i,j)=std_error(histogram->GetBinContent(i,j));
			}
		    }
		    file->Close();
                    delete file;
		    return res;
		}else throw Exception<TH1F>("No histogram "+histname+" in file "+filename);
		file->Close();
		delete file;
	    }
	    return hist2d<>({},{});
	}
const vector<size_t> valid_runs_list{
                                                                                          45935,45936,45937,45938,45939,
45940,45941,45942,45943,45944,45945,45946,45947,45948,45949,45950,45951,45952,45953,      45955,45956,45957,45958,45959,
45960,45961,45962,45963,45964,45965,45966,45967,45968,      45970,45971,45972,45973,45974,45975,45976,45977,45978,45979,
45980,45981,45982,45983,45984,45985,45986,45987,45988,45989,45990,45991,45992,45993,45994,45995,45996,45997,45998,45999,
46000,46001,46002,46003,      46005,46006,46007,46008,46009,46010,46011,      46013,46014,46015,46016,46017,46018,46019,
46020,46021,46022,46023,46024,46025,46026,46027,      46029,46030,46031,46032,46033,46034,46035,46036,46037,46038,46039,
46040,46041,46042,46043,46044,46045,46046,46047,46048,46049,46050,46051,46052,46053,46054,46055,46056,46057,46058,46059,
46060,46061,46062,46063,46064,46065,46066,46067,46068,46069,46070,46071,46072,46073,46074,46075,46076,      46078,46079,
46080,46081,46082,46083,46084,46085,46086,46087,46088,46089,46090,46091,46092,46093,46094,46095,46096,46097,46098,46099,
46100,46101,46102,46103,46104,            46107,46108,46109,                                    46116,46117,46118,46119,
46120,46121,46122,46123,46124,46125,46126,46127,46128,46129,46130,46131,46132,46133,46134,46135,46136,46137,46138,46139,
      46141,46142,46143,46144,      46146,46147,46148,46149,46150,46151,46152,46153,46154,46155,46156,46157,46158,46159,
46160,46161,46162,46163,46164,46165,46166,46167,      46169,46170,      46172,46173,46174,46175,46176,46177,46178,46179,
46180,46181,      46183,46184,46185,46186,46187,      46189,46190,46191,46192,46193,46194,46195,46196,46197,46198,46199,
46200,46201,46202,46203,46204,46205,46206,46207,46208,46209,46210,46211,46212,46213,46214,46215,46216,46217,46218,46219,
46220,46221,46222,46223,46224,46225,46226,46227,46228,46229,46230,46231,46232,46233,46234,46235,46236,46237,46238,46239,
46240,46241,46242,46243,46244,46245,46246,46247,46248,46249,46250,46251,46252,46253,46254,46255,46256,46257,46258,46259,
46260,46261,46262,46263,46264,46265,46266,46267,                  46271,46272,46273,46274,46275,46276,46277,46278,46279,
46280,46281,46282,46283,46284,46285,46286,46287,46288,46289,46290,46291,46292,46293,46294,46295,46296,46297,46298,46299,
46300,46301,46302,46303,46304,46305,46306,46307,46308,46309,46310,46311,46312,46313,46314,46315,46316,46317,46318,46319,
46320,46321,46322,46323,46324,46325,46326,46327,46328,46329,46330,46331,46332,46333,46334,46335,46336,46337,46338,46339,
46340,46341,46342,46343,46344,46345,46346,46347,46348,46349,46350,46351,46352,46353,46354,46355,46356,46357,46358,46359,
46360,46361,46362,46363,46364,46365,46366,46367,46368,46369,46370,46371,46372,46373,46374,46375,46376,46377,46378,46379,
46380,46381,46382,46383,46384,46385,46386,46387,46388,46389,46390,46391,46392,46393,46394,46395,46396,46397,46398,46399,
46400,46401,46402,46403,46404,46405,46406,46407,46408,46409,46410,46411,46412,46413,46414,46415,46416,46417,46418,46419,
46420,46421,46422,46423,46424,46425,46426,46427,46428,46429,46430,46431,46432,46433,46434,46435,46436,46437,46438,46439,
46440,46441,46442,46443,46444,46445,46446,46447,46448,      46450,46451,      46453,46454,46455,46456,46457,46458,46459,
46460,46461,46462,46463,46464,46465,46466,46467,46468,46469,46470,46471,46472,      46474,46475,46476,46477,46478,46479,
46480,46481,46482,46483,46484,46485,46486,46487,46488,46489,46490,46491,46492,46493,46494,46495,46496,46497,46498,46499,
46500,46501,46502,46503,46504,46505,46506,46507,46508,46509,46510,46511,46512,46513,46514,46515,46516,46517,46518,46519,
46520,46521,46522,46523,46524,46525,46526,46527,46528,46529,46530,46531,46532,46533,46534,46535,46536,46537,46538,46539,
46540,46541,46542,46543,46544,46545,46546,46547,46548,46549,46550,46551,      46553,46554,46555,46556,46557,46558,46559,
46560,46561,46562,46563,46564,46565,46566,46567,46568,46569,46570,46571,46572,46573,46574,      46576,46577,46578,46579,
46580,46581,46582,46583,46584,46585,46586,46587,46588,46589,46590,46591,46592,46593,46594,46595,46596,46597,46598,46599,
46600,46601,46602,46603,46604,46605,46606,46607,46608,46609,46610,46611,46612,46613,46614,46615,46616,46617,46618,46619,
46620,46621,46622,46623,46624,46625,46626,46627,46628,46629,      46631,46632,46633,46634,46635,46636,46637,46638,46639,
46640,46641,46642,46643,46644,46645,46646,46647,46648,46649,46650,46651,46652,46653,46654,46655,46656,46657,46658,46659,
46660,46661,46662,46663,46664,46665,46666,46667,46668,46669,46670,46671,46672,46673,46674,46675,      46677,46678,46679,
46680,46681,      46683,46684,46685,      46687,      46689,            46692,46693,46694,46695,46696,46697,46698,46699,
46700,46701,46702,46703,46704,46705,46706,46707,46708,46709,46710,46711,46712,46713,46714,46715,46716,46717,46718,46719,
46720,46721,46722,46723,46724,46725,46726,46727,46728,46729,46730,46731,46732,46733,46734,46735,46736,46737,
                                                46748,46749,46750,46751,46752,46753,46754,46755,46756,46757,46758,46759,
46760,46761,46762,46763,46764,46765,46766,46767,46768,46769,46770,46771,46772,46773,46774,46775,46776,46777,46778,46779,
46780,46781,46782,46783,46784,46785,46786,46787,46788,46789,46790,46791,46792,      46794,46795,46796,46797,46798,46799,
46800,46801,46802,46803,      46805,                              46811,46812,46813,46814,46815,46816,46817,46818,46819,
46820,46821,46822,46823,46824,46825,46826,46827,46828,46829,46830,46831,46832,46833,46834,46835,46836,46837,46838,46839,
46840,46841,46842,46843,46844,46845,46846,46847,46848,46849,46850,46851,46852,46853,46854,46855,46856,46857,46858,46859,
46860,46861,46862,46863,46864,46865,46866,46867,46868,46869,46870,46871,46872,46873,46874,46875,46876,46877,46878,46879,
46880,46881,46882,46883,46884
};
const vector<string> analyses{"All"};
vector<size_t> present_runs;
const string suffix_for_test="19+";
void create_runs_list(){
    if(present_runs.size()>0)return;
    for(const size_t runindex:valid_runs_list){
        bool accept=(runindex<=46492);
        if(accept)for(const string&reaction:analyses){
            TFile *file=TFile::Open((inputpath+"/Data"+reaction+to_string(runindex)+suffix_for_test+".root").c_str());
            if(file){
                file->Close();
                delete file;
            }else{
                accept=false;
            }
        }
        if(accept)present_runs.push_back(runindex);
    }
}
hist<double> Hist(histsource src, const string&reaction, const vector<string>&path,const string&histname,const string&param_suffix){
    create_runs_list();
    hist<double> res;
    switch(src){ 
        case MC:{
            for(int runindex=1;runindex<=10;runindex++){
                hist<double> tmp=ReadHist(inputpath+"/MC"+reaction+to_string(runindex)+param_suffix+".root",path,histname);
                if(tmp.size()>0){
                    if(res.size()==0)
                        res=tmp;
                    else
                        res.imbibe(tmp);
                }
            }
        }break;
        case DATA:{
            for(const size_t runindex:present_runs){
                hist<double> tmp=ReadHist(inputpath+"/Data"+reaction+to_string(runindex)+param_suffix+".root",path,histname);
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
hist2d< double > Hist2d(histsource src, const string& reaction, const vector< string >& path, const string& histname,const string&param_suffix){
    create_runs_list();
    hist2d<double> res;
    switch(src){
        case MC:{
            res=ReadHist2D(inputpath+"/MC"+reaction+".root",path,histname);
            for(int runindex=1;runindex<=10;runindex++){
                hist2d<double> tmp=ReadHist2D(inputpath+"/MC"+reaction+to_string(runindex)+param_suffix+".root",path,histname);
                if(tmp.size()>0){
                    if(res.size()==0)
                        res=tmp;
                    else
                        res.imbibe(tmp);
                }
            }
        }break;
        case DATA:{
            for(const size_t runindex:present_runs){
                hist2d<double> tmp=ReadHist2D(inputpath+"/Data"+reaction+to_string(runindex)+param_suffix+".root",path,histname);
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
pair<double,double> PresentRuns(const string&reaction,const string&param_suffix){
    create_runs_list();
    size_t prr=0;
        for(const size_t runindex:present_runs){
        TFile *file=TFile::Open((inputpath+"/Data"+reaction+to_string(runindex)+param_suffix+".root").c_str());
        if(file){
            prr++;
            file->Close();
            delete file;
        }
    }
    return make_pair(double(prr),double(valid_runs_list.size()));
}

};
