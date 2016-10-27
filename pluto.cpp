// this file is distributed under 
// MIT license
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
	Reaction he3eta(Particle::p(),Particle::d(),{Particle::he3(),Particle::eta()});
	TF1 *mf=nullptr;
	const double thr=1.573;
	if(string(arg[1])=="all")
		mf=new TF1("Uniform","1",p_beam_low,p_beam_hi);
	else
		mf=new TF1("Uniform","1",thr,p_beam_hi);
	cout<<he3eta.P2Q(p_beam_low)<<endl;
	cout<<he3eta.P2Q(p_beam_hi)<<endl;
	cout<<he3eta.P2Q(thr)<<endl;
	PBeamSmearing *smear = new PBeamSmearing("beam_smear", "Beam smearing");
	smear->SetReaction("p + d");
	smear->SetMomentumFunction(mf);
	makeDistributionManager()->Add(smear);
	std::mt19937 gen;
	std::uniform_int_distribution<int> d(1,254);
	PUtils::SetSeed(d(gen));
	cout<<react<<endl;
	PReaction my_reaction(const_cast<char*>(to_string(p_beam_hi).c_str()),"p","d",
		const_cast<char*>(react.c_str()),
		const_cast<char*>(ReplaceAll(ReplaceAll(ReplaceAll(react," ",""),"[","_"),"]","_").c_str())
		,1,0,0,0);
	my_reaction.Print();
	if(string(arg[1])=="all")
		my_reaction.Loop(10000000);
	else
		my_reaction.Loop(5000000);
	POPD();
	return 0;
}