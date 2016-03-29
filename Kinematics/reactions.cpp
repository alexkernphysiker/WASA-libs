// this file is distributed under 
// MIT license
#include <unistd.h>
#include <math.h>
#include <math_h/interpolate.h>
#include <math_h/error.h>
#include <Experiment/experiment_conv.h>
#include "particles.h"
#include "reactions.h"
namespace MathTemplates {
	template<class numt>
	class Vector3{
	private:
		numt m_x,m_y,m_z;
		Vector3(const numt&x,const numt&y,const numt&z):m_x(x),m_y(y),m_z(z){}
	public:
		static const Vector3 zero(){return Vector3(numt(0),numt(0),numt(0));}
		static const Vector3 DesCartes(const numt&x,const numt&y,const numt&z){return Vector3(x,y,z);}
		static const Vector3 DesCartes(const numt&&x,const numt&&y,const numt&&z){return Vector3(x,y,z);}
		static const Vector3 Polar(const numt&mag,const numt&theta,const numt&phi){
			return Vector3(mag*cos(phi)*sin(theta),mag*sin(phi)*sin(theta),mag*cos(theta));
		}
		static inline const Vector3 Polar(const numt&mag,const numt&&theta,const numt&&phi){return Polar(mag,theta,phi);}
		static inline const Vector3 Polar(const numt&&mag,const numt&theta,const numt&phi){return Polar(mag,theta,phi);}
		static inline const Vector3 Polar(const numt&&mag,const numt&&theta,const numt&&phi){return Polar(mag,theta,phi);}
		static inline const Vector3 Direction(const numt&theta,const numt&phi){return Polar(numt(1),theta,phi);}
		static inline const Vector3 Direction(const numt&&theta,const numt&phi){return Polar(numt(1),theta,phi);}
		static inline const Vector3 Direction(const numt&theta,const numt&&phi){return Polar(numt(1),theta,phi);}
		static inline const Vector3 Direction(const numt&&theta,const numt&&phi){return Polar(numt(1),theta,phi);}
		const numt&x()const{return m_x;}
		const numt&y()const{return m_y;}
		const numt&z()const{return m_z;}
		const numt mag_sqr()const{
			return sqrt(pow(m_x,2)+pow(m_y,2)+pow(m_z,2));
		}
		const numt mag()const{
			return sqrt(pow(m_x,2)+pow(m_y,2)+pow(m_z,2));
		}
		const numt cos_theta()const{
			return m_z/mag();
		}
		const numt sin_theta()const{
			return sqrt(pow(m_x,2)+pow(m_y,2))/mag();
		}
		const numt cos_phi()const{
			return m_x/sqrt(pow(m_x,2)+pow(m_y,2));
		}
		const numt sin_phi()const{
			return m_x/sqrt(pow(m_x,2)+pow(m_y,2));
		}
		
		Vector3(const Vector3&source):m_x(source.m_x),m_y(source.m_y),m_z(source.m_z){}
		Vector3&operator=(const Vector3&source){
			m_x=source.m_x;m_y=source.m_y;m_z=source.m_z;
			return *this;
		}
		Vector3&operator+=(const Vector3&second){
			m_x+=second.m_x;m_y+=second.m_y;m_z+=second.m_z;
			return *this;
		}
		Vector3&operator+=(const Vector3&&second){
			return operator+=(second);
		}
		const Vector3 operator+(const Vector3&second)const{
			return Vector3(m_x+second.m_x,m_y+second.m_y,m_z+second.m_z);
		}
		const Vector3 operator+(const Vector3&&second)const{
			return operator+(second);
		}
		Vector3&operator-=(const Vector3&second){
			m_x-=second.m_x;m_y-=second.m_y;m_z-=second.m_z;
			return *this;
		}
		Vector3&operator-=(const Vector3&&second){
			return operator-=(second);
		}
		const Vector3 operator-(const Vector3&second)const{
			return Vector3(m_x-second.m_x,m_y-second.m_y,m_z-second.m_z);
		}
		const Vector3 operator-(const Vector3&&second)const{
			return operator-(second);
		}

		Vector3&operator*=(const numt&second){
			m_x*=second;m_y*=second;m_z*=second;
			return *this;
		}
		Vector3&operator*=(const numt&&second){
			return operator*=(second);
		}
		const Vector3 operator*(const numt&second)const{
			return Vector3(m_x*second,m_y*second,m_z*second);
		}
		const Vector3 operator*(const numt&&second)const{
			return operator*(second);
		}

		Vector3&operator/=(const numt&second){
			m_x/=second;m_y/=second;m_z/=second;
			return *this;
		}
		Vector3&operator/=(const numt&&second){
			return operator/=(second);
		}
		const Vector3 operator/(const numt&second)const{
			return Vector3(m_x/second,m_y/second,m_z/second);
		}
		const Vector3 operator/(const numt&&second)const{
			return operator/(second);
		}

		const numt operator*(const Vector3&second)const{
			return (m_x*second.m_x)+(m_y*second.m_y)+(m_z*second.m_z);
		}
		const numt operator*(const Vector3&&second)const{
			return operator*(second);
		}
		
		const Vector3 VecP(const Vector3&second)const{
			return Vector3<numt>::DesCartes(
				(y()*second.z())-(second.y()*z()),
				(z()*second.x())-(second.z()*x()),
				(x()*second.y())-(second.x()*y())
			);
		}
		const Vector3 VecP(const Vector3&&second)const{
			return VecP(second);
		}
	};
	
	template<class numt>
	class FourMomentum{
	private:
		numt m_E;
		Vector3<numt> m_P;
		FourMomentum(const numt&E,const Vector3<numt>&P):m_E(E),m_P(P){
			if(m_E<0)throw Exception<FourMomentum>("Cannot set negative full energy");
		}
		FourMomentum(const numt&&E,const Vector3<numt>&P):FourMomentum(E,P){}
		FourMomentum(const numt&E,const Vector3<numt>&&P):FourMomentum(E,P){}
		FourMomentum(const numt&&E,const Vector3<numt>&&P):FourMomentum(E,P){}
	public:
		static const FourMomentum zero(){return FourMomentum(numt(0),Vector3<numt>::zero());}
		static const FourMomentum MassMomentum(const numt&m,const Vector3<numt>&p){
			if(m<0)throw Exception<FourMomentum>("Cannot set negative mass");
			return FourMomentum(sqrt(p.mag_sqr()+m*m),p);
		}
		static inline const FourMomentum MassMomentum(const numt&&m,const Vector3<numt>&p){return MassMomentum(m,p);}
		static inline const FourMomentum MassMomentum(const numt&m,const Vector3<numt>&&p){return MassMomentum(m,p);}
		static inline const FourMomentum MassMomentum(const numt&&m,const Vector3<numt>&&p){return MassMomentum(m,p);}
		static const FourMomentum MassEkinDir(const numt&m,const numt&E_k,const numt&theta,const numt&phi){
			if(m<0)throw Exception<FourMomentum>("Cannot set negative mass");
			if(E_k<0)throw Exception<FourMomentum>("Cannot set negative kinetic energy");
			numt E=m+E_k;
			numt p=sqrt((E*E)-(m*m));
			return FourMomentum(E,Vector3<numt>::Direction(theta,phi)*p);
		}
		static inline const FourMomentum MassEkinDir(const numt&m,const numt&E_k,const numt&&theta,const numt&&phi){return MassEkinDir(m,E_k,theta,phi);}
		static inline const FourMomentum MassEkinDir(const numt&m,const numt&&E_k,const numt&theta,const numt&phi){return MassEkinDir(m,E_k,theta,phi);}
		static inline const FourMomentum MassEkinDir(const numt&m,const numt&&E_k,const numt&&theta,const numt&&phi){return MassEkinDir(m,E_k,theta,phi);}
		static inline const FourMomentum MassEkinDir(const numt&&m,const numt&&E_k,const numt&&theta,const numt&&phi){return MassEkinDir(m,E_k,theta,phi);}
		
		const numt&E()const{return m_E;}
		const Vector3<numt>&P()const{return m_P;}
		const numt p()const{return m_P.mag();}
		const numt m()const{return sqrt((m_E*m_E)-m_P.mag_sqr());}
		const numt Ekin()const{return m_E-m();}
		
		FourMomentum(const FourMomentum&source):FourMomentum(source.m_E,source.m_P){}
		FourMomentum&operator=(const FourMomentum&source){m_E=source.m_E;m_P=source.m_P;}
		
		FourMomentum&operator+=(const FourMomentum&source){
			m_E+=source.m_E;
			m_P+=source.m_P;
			return *this;
		}
		FourMomentum&operator+=(const FourMomentum&&source){
			return operator+=(source);
		}
		const FourMomentum operator+(const FourMomentum&source)const{
			return FourMomentum(m_E+source.m_E,m_P+source.m_P);
		}
		const FourMomentum operator+(const FourMomentum&&source)const{
			return operator+(source);
		}

		FourMomentum&operator-=(const FourMomentum&source){
			m_E-=source.m_E;
			m_P-=source.m_P;
			return *this;
		}
		FourMomentum&operator-=(const FourMomentum&&source){
			return operator-=(source);
		}
		const FourMomentum operator-(const FourMomentum&source)const{
			return FourMomentum(m_E-source.m_E,m_P-source.m_P);
		}
		const FourMomentum operator-(const FourMomentum&&source)const{
			return operator-(source);
		}
	};
};
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

const double Reaction::EThreshold() const{
	return M_after()-M_before();
}
const double Reaction::PThreshold() const{
	return m_projectile.E2P(EThreshold());
}
const double Reaction::E2Q(const double E) const{
	return sqrt(pow(E+M_before(),2)-pow(m_projectile.E2P(E),2))-M_after();
}
const double Reaction::P2Q(const double P) const{
	return sqrt(pow(m_projectile.P2E(P)+M_before(),2)-pow(P,2))-M_after();
}
const double Reaction::PbEr2Theta(const double Pbeam, const double Ereg) const{
	if(m_products.size()!=2)
		throw Exception<Reaction>("Calculation of theta is available only for binary reactions");
	double E=m_projectile.P2E(Pbeam);
	double P=Pbeam;
	double Q=P2Q(P);
	double beta = P/(E+M_before());
	double gamma = 1./sqrt(1.-pow(beta,2));
	double T_reg_cm = Q/2.*(Q+2*m_products[1].mass())/(Q+M_after());
	double beta_reg_cm = sqrt(pow(T_reg_cm,2)+2*T_reg_cm*m_products[0].mass())/(T_reg_cm+m_products[0].mass());
	double p_reg_cm = beta_reg_cm*(m_products[0].mass()+T_reg_cm);
	double theta_reg_cm = acos((Ereg-(gamma-1)*m_products[0].mass()-gamma*T_reg_cm)/(gamma*beta*p_reg_cm));
	return atan(sin(theta_reg_cm)/(gamma*(cos(theta_reg_cm)+beta/beta_reg_cm)));
}

const double Reaction::MissingMass(const initializer_list<registered_particle_parameters>& data,const double Pbeam) const{
	auto PTotal=FourMomentum<double>::MassMomentum(m_projectile.mass(),Vector3<double>::Direction(0,0)*Pbeam)
		+FourMomentum<double>::MassMomentum(m_target.mass(),Vector3<double>::zero());
	auto PReg=FourMomentum<double>::zero();
	for(const auto&pr:data)
		PReg+=FourMomentum<double>::MassEkinDir(products()[pr.index].mass(),pr.E,pr.theta,pr.phi);
	return (PTotal-PReg).m();
}

const double InvariantMass(const initializer_list<particle_kinematics>& data){
	auto total=FourMomentum<double>::zero();
	for(const auto&pr:data)
		total+=FourMomentum<double>::MassEkinDir(pr.particle.mass(),pr.E,pr.theta,pr.phi);
	return total.m();
}
