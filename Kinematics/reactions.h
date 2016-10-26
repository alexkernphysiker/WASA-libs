// this file is distributed under 
// MIT license
#ifndef xylgjnjy
#	define xylgjnjy
#include <list>
#include <vector>
#include "particles.h"
class Reaction{
public:
	Reaction(const Particle&p,const Particle&t,const std::initializer_list<Particle>&products);
	Reaction(const Particle&&p,const Particle&&t,const std::initializer_list<Particle>&products);
	Reaction(const Reaction&source);
	virtual ~Reaction();
	const Particle&target()const;
	const Particle&projectile()const;
	const std::vector<Particle>&products()const;
	const double M_before()const;
	const double M_after()const;
	const double E2Q(const double E)const;
	const double P2Q(const double P)const;
	const double PbEr2Theta(const double Pbeam,const double Ereg)const;
	struct registered_particle_parameters{unsigned int index;double E;double theta;double phi;};
	const double MissingMass(const std::initializer_list<registered_particle_parameters>&data,const double Pbeam)const;
private:
	Particle m_projectile,m_target;
	std::vector<Particle> m_products;
};
struct particle_kinematics{Particle particle; double E,theta,phi;};
const double InvariantMass(const std::initializer_list<particle_kinematics>&data);
#endif 