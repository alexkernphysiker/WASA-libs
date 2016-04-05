// this file is distributed under 
// MIT license
#include <math.h>
#include "particles.h"
Particle::Particle():m_mass(0),m_charge(0){}
Particle::Particle(double m, int c):m_mass(m),m_charge(c){}
Particle::Particle(const Particle& source):m_mass(source.m_mass),m_charge(source.m_charge){}
Particle::Particle(const Particle&& source):m_mass(source.m_mass),m_charge(source.m_charge){}
Particle& Particle::operator=(const Particle& source){
	m_mass=source.m_mass;
	m_charge=source.m_charge;
	return *this;
}
Particle& Particle::operator=(const Particle&& source){
	return operator=(source);
}
Particle::~Particle(){}
const double Particle::mass() const{return m_mass;}
const int Particle::charge() const{return m_charge;}
const double Particle::E2P(const double E) const{
	return sqrt(pow(E+mass(),2)-pow(mass(),2));
}
const double Particle::P2E(const double P) const{
	return sqrt(pow(P,2)+pow(mass(),2))-mass();
}
const Particle Particle::gamma(){
	return Particle(0.0,0);
}
const Particle Particle::n(){
	return Particle(0.93956,0);
}
const Particle Particle::p(){
	return Particle(0.938272,1);
}
const Particle Particle::d(){
	return Particle(1.875613,1);
}
const Particle Particle::he3(){
	return Particle(2.808950,2);
}
const Particle Particle::eta(){
	return Particle(0.5478,0);
}
const Particle Particle::pi0(){
	return Particle(0.1350,0);
}
const Particle Particle::pi_minus(){
	return Particle(0.1396,-1);
}
const Particle Particle::pi_plus(){
	return Particle(0.1396,1);
}