// this file is distributed under 
// MIT license
#include <math.h>
#include "particles.h"
Particle::Particle(double m, int c):m_mass(m),m_charge(c){}
Particle::Particle(const Particle& source):m_mass(source.m_mass),m_charge(source.m_charge){}
Particle::~Particle(){}
const double Particle::mass() const{return m_mass;}
const int Particle::charge() const{return m_charge;}
const double Particle::E2P(const double E) const{
	return sqrt(pow(E+mass(),2)-pow(mass(),2));
}
const double Particle::P2E(const double P) const{
	return sqrt(pow(P,2)+pow(mass(),2))-mass();
}

const Particle& Particle::n(){
	static Particle res(0.93956,0);
	return res;
}
const Particle& Particle::p(){
	static Particle res(0.938272,1);
	return res;
}
const Particle& Particle::d(){
	static Particle res(1.875613,1);
	return res;
}
const Particle& Particle::he3(){
	static Particle res(2.808950,2);
	return res;
}
const Particle& Particle::eta(){
	static Particle res(0.5478,0);
	return res;
}
const Particle& Particle::pi0(){
	static Particle res(0.1350,0);
	return res;
}
const Particle& Particle::pi_minus(){
	static Particle res(0.1396,-1);
	return res;
}
const Particle& Particle::pi_plus(){
	static Particle res(0.1396,1);
	return res;
}