// this file is distributed under 
// GPL license
#ifndef ______PARTICLES_H_____
#	define ______PARTICLES_H_____
#include <list>
#include <vector>
#include <math_h/vectors.h>
class Particle{
protected:
	Particle();
public:
	Particle(const Particle&source);
	bool operator==(const Particle&source)const;
	virtual ~Particle();
	
	const double mass()const;
	const int charge()const;
	const double E2P(const double&E)const;
	const double P2E(const double&P)const;
	
	static const Particle gamma();
	static const Particle n();
	static const Particle p();
	static const Particle d();
	static const Particle he3();
	static const Particle eta();
	static const Particle pi0();
	static const Particle pi_plus();
	static const Particle pi_minus();
private:
	Particle(double m, int c);
	double m_mass;
	int m_charge;
};
// E - GeV, theta,phi - radians
struct particle_kinematics{Particle particle; double E,theta,phi;};
const MathTemplates::Vector4<double> Get4Vector(const particle_kinematics&data);
const MathTemplates::Vector4<double> Get4Vector(const std::vector<particle_kinematics>&data);
const double InvariantMass(const std::vector<particle_kinematics>&data);
#endif
