// this file is distributed under 
// GPL license
#include <unistd.h>
#include <math.h>
#include <math_h/error.h>
#include <math_h/vectors.h>
#include <Experiment/experiment_conv.h>
#include "particles.h"
#include "reactions.h"
using namespace std;
using namespace MathTemplates;
Reaction::Reaction(const Particle& p, const Particle& t, const initializer_list<Particle>& products)
:m_projectile(p),m_target(t){
	for(const Particle&item:products)
		m_products.push_back(item);
}
Reaction::Reaction(const Reaction& source)
:m_projectile(source.m_projectile),m_target(source.m_target){
	for(const Particle&item:source.m_products)
		m_products.push_back(item);
}
Reaction::~Reaction(){}
const Particle& Reaction::projectile()const{return m_projectile;}
const Particle& Reaction::target()const{return m_target;}
const vector<Particle>& Reaction::products()const{return m_products;}
const double Reaction::M_before() const{return m_projectile.mass()+m_target.mass();}
const double Reaction::M_after() const{
	double res=0;
	for(const Particle&item:m_products)res+=item.mass();
	return res;
}

const double Reaction::E2Q(const double E) const{
	return sqrt(pow(E+M_before(),2)-pow(m_projectile.E2P(E),2))-M_after();
}
const double Reaction::P2Q(const double P) const{
	return sqrt(pow(m_projectile.P2E(P)+M_before(),2)-pow(P,2))-M_after();
}
const double Reaction::MissingMass(const vector<registered_particle_parameters>& data,const double Pbeam) const{
	auto PTotal=lorentz_byPM(Z()*Pbeam,m_projectile.mass())+lorentz_byPM(Zero(),m_target.mass());
	auto PReg=LorentzVector<>::zero();
	for(const auto&pr:data){
		double m=products()[pr.index].mass();
		PReg+=lorentz_byEM(m+pr.E,m,direction(pr.phi,pr.theta));
	}
	return (PTotal-PReg).M();
}
