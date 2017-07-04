// this file is distributed under 
// GPL license
#include <math.h>
#include "particles.h"
using namespace std;
using namespace MathTemplates;
const Vector4<double> Get4Vector(const particle_kine&data){
    return Vector4<double>::byTime_Dir_and_Length4(data.particle.mass()+data.E,data.theta,data.phi,data.particle.mass());
}
const Vector4<double> Get4Vector(const std::vector<particle_kine>&data){
    auto total=Vector4<double>::zero();
    for(const auto&pr:data)
	total+=Get4Vector(pr);
    return total;
}
const Vector4<double> Get4Vector(const particle_kinp&data){
    return Vector4<double>::bySpaceC_and_Length4(Vector3<double>::Polar(data.P,data.theta,data.phi),data.particle.mass());
}
const Vector4<double> Get4Vector(const std::vector<particle_kinp>&data){
    auto total=Vector4<double>::zero();
    for(const auto&pr:data)
	total+=Get4Vector(pr);
    return total;
}

Particle::Particle():m_mass(INFINITY),m_charge(0){}
Particle::Particle(double m, int c):m_mass(m),m_charge(c){}
Particle::Particle(const Particle& source):m_mass(source.m_mass),m_charge(source.m_charge){}
bool Particle::operator==(const Particle& source) const{
	return (m_mass==source.m_mass)&&(m_charge==source.m_charge);
}
Particle::~Particle(){}
const double&Particle::mass() const{return m_mass;}
const int&Particle::charge() const{return m_charge;}
const double Particle::E2P(const double&E) const{
	return sqrt(pow(E+mass(),2)-pow(mass(),2));
}
const double Particle::P2E(const double&P) const{
	return sqrt(pow(P,2)+pow(mass(),2))-mass();
}
const Particle& Particle::gamma(){
	static Particle res(0.0,0);
	return res;
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
const Particle& Particle::he4(){
	static Particle res(3.7274225,2);
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
