// this file is distributed under 
// GPL license
#include <list>
#include <string>
#include <sstream>
#include <random>
#include <PBeamSmearing.h>
#include <PReaction.h>
#include <Experiment/str_get.h>
#include <Experiment/experiment_conv.h>
#include <Kinematics/reactions.h>
#include "config.h"
using namespace std;
string ReplaceAll(const string&source, const string& from, const string& to) {
	string str=source;
	size_t start_pos = 0;
	while((start_pos = str.find(from, start_pos)) != string::npos) {
		str.replace(start_pos, from.length(), to);
		start_pos += to.length();
	}
	return str;
}
string ReplaceAll(string&& str, const string& from, const string& to){
	return ReplaceAll(str,from,to);
}
int main(int argc, char **arg){
	if(argc<3){
		cout<<"reaction expected"<<endl;
		return -1;
	}
	string react=arg[2];
	for(int i=3;i<argc;i++)
		react+=" "+string(arg[i]);
	PUSHD();
	CD(PLUTO);
	std::mt19937 gen;
	std::uniform_int_distribution<int> d(1,254);
	for(ALLMC){
	    Reaction he3eta(Particle::p(),Particle::d(),{Particle::he3(),Particle::eta()});
	    TF1 *mf=nullptr;
	    const double thr=1.573;
	    if(string(arg[1])=="all")
		mf=new TF1(CSTR("Uniform"),CSTR("1"),p_beam_low,p_beam_hi);
	    else
		mf=new TF1(CSTR("Uniform"),CSTR("1"),thr,p_beam_hi);
	    cout<<he3eta.P2Q(p_beam_low)<<endl;
	    cout<<he3eta.P2Q(p_beam_hi)<<endl;
	    cout<<he3eta.P2Q(thr)<<endl;
	    PBeamSmearing *smear = new PBeamSmearing(CSTR("beam_smear"),CSTR("Beam smearing"));
	    smear->SetReaction(CSTR("p + d"));
	    smear->SetMomentumFunction(mf);
	    makeDistributionManager()->Add(smear);
	    PUtils::SetSeed(d(gen));
	    cout<<react<<endl;
	    PReaction my_reaction(CSTR(to_string(p_beam_hi)),CSTR("p"),CSTR("d"),CSTR(react),
		CSTR(ReplaceAll(ReplaceAll(ReplaceAll(react," ",""),"[","_"),"]","_")+"-"+to_string(runindex))
		,1,0,0,0
	    );
	    my_reaction.Print();
	    my_reaction.Loop(1000000);
	}
	POPD();
	return 0;
}
