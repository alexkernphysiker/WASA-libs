// this file is distributed under 
// GPL license
#include <math.h>
#include "particles.h"
using namespace std;
using namespace MathTemplates;
const Vector4<double> Get4Vector(const particle_kinematics&data){
    return Vector4<double>::TimeDirLength4(data.particle.mass()+data.E,data.theta,data.phi,data.particle.mass());
}
const Vector4<double> Get4Vector(const std::vector<particle_kinematics>&data){
    auto total=Vector4<double>::zero();
    for(const auto&pr:data)
	total+=Get4Vector(pr);
    return total;
}
const double InvariantMass(const vector<particle_kinematics>& data){
	return Get4Vector(data).length4();
}

Particle::Particle():m_mass(INFINITY),m_charge(0){}
Particle::Particle(double m, int c):m_mass(m),m_charge(c){}
Particle::Particle(const Particle& source):m_mass(source.m_mass),m_charge(source.m_charge){}
bool Particle::operator==(const Particle& source) const{
	return (m_mass==source.m_mass)&&(m_charge==source.m_charge);
}
Particle::~Particle(){}
const double Particle::mass() const{return m_mass;}
const int Particle::charge() const{return m_charge;}
const double Particle::E2P(const double&E) const{
	return sqrt(pow(E+mass(),2)-pow(mass(),2));
}
const double Particle::P2E(const double&P) const{
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
const Particle Particle::he4(){
	return Particle(3.7274225,2);
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
