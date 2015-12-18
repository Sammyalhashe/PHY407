#ifndef functions_hpp_IS_INCLUDED
#define functions_hpp_IS_INCLUDED

#include <cmath>
#include "QickArray.hpp"
#include "random.hpp"

//=================================================================
//package of functors, with a lot of overloads of valuePt() and
//() operator. The function valuePt() is allways the main 
//function of a functor, i.e. where the work is done.
//The () operators are merely for ease of use;
//output = functorname(&Func, double)
//=================================================================
using namespace std;

class Func{
  
public:

  virtual ~Func(){}

  virtual inline void setParams(QickArray& params_){};

  virtual bool getImpossible(){ return false; }

  virtual double valuePt(double x, double y){ return 42; }

  virtual double valuePt(QickArray &pos){ return 42; }

  virtual double valuePt(QickArray &pos, int dummy){ return 42; }

  virtual QickArray& valuePt(QickArray& pos, int i, 
			     Random& ran){ return pos; }

  virtual QickArray& valuePt(QickArray& pos, int i,  
			     QickArray& x, Random& ran){ return pos; }

  virtual QickArray& valuePt(QickArray &pos, bool dummy){ return pos; }

  virtual double valuePt(Func &func, QickArray &pos){ return 42; }

  virtual QickArray& valuePt(QickArray &pos, Func &func){ 
    return pos; 
  }

  virtual QickArray& valuePt(QickArray &pos, double x){ 
    return pos; 
  }

  virtual inline double operator()(double x, double y){ return 42; }

  virtual inline double operator()(QickArray &pos){ return 42; }

  virtual inline QickArray& operator()(QickArray& pos, int i,  
				   Random& ran){ return pos; }

  virtual inline QickArray& operator()(QickArray& pos, int i,  
				   QickArray& x, Random& ran){ return pos; }

  virtual inline QickArray& operator()(QickArray &pos, bool dummy)
  { return pos; }

  virtual inline double operator()(Func &func, QickArray &pos){ return 42; }

  virtual inline QickArray& operator()(QickArray &pos, Func &func){ 
    return pos; 
  }

  virtual inline QickArray& operator()(QickArray &pos, double x){ 
    return pos; 
  }

  virtual inline double sqr(double x){ return x*x; }

};


//=================================================================
// functor evaluating first derivate returning positions x dim array
//=================================================================
class Nabla : public virtual Func{

private:
  QickArray dummy;

public:

  virtual QickArray& valuePt(QickArray &pos, Func &func);
  
  virtual inline QickArray& operator()(QickArray &pos, Func &func){ 
    
    return valuePt(pos, func);
  }

};


//=================================================================
// functor evaluating second derivate summing over all positions
//=================================================================
class Nabla2 : public virtual Func{

public:

  virtual double valuePt(Func &func, QickArray &pos);

  virtual inline double operator()(Func &func, QickArray &pos){ 

    return valuePt(func, pos);
  }

};


//=================================================================
//a super class giving the interface for all wave functions
//no other use what so ever
//all wave functions should be on the form u(r) where wf_T=exp(u(r))
//=================================================================
class Wave : public virtual Func{

public:

  virtual inline void setParams(QickArray& params_){}

  virtual double valuePt(QickArray &pos){return 42.;}

  virtual inline double operator()(QickArray &pos){ 
    return valuePt(pos);
  }

};

//=================================================================
//functor for 1s wave function
//=================================================================
class Wave1s : public virtual Wave{

private:
  double alpha;

public:

  virtual inline void setParams(QickArray& params_){ alpha=params_(0); }

  virtual double valuePt(QickArray &pos);

  virtual inline double operator()(QickArray &pos){ 
    return valuePt(pos);
  }
};

//=================================================================
// u(r)=-alpha*(r)
//=================================================================
inline double Wave1s:: valuePt(QickArray &pos){

  int dim=pos.no_of_elements();
  double ans=0;
  for(int i=0; i!=dim; i++)
    ans+=sqr(pos(i));

  return -alpha*sqrt(ans);

}

//=================================================================
//functor for a simple correlation function
//=================================================================
class Correlation : public virtual Wave{

private:
  double beta;
  
public:

  virtual inline void setParams(QickArray& params_){ beta=params_(1); }

  virtual double valuePt(QickArray &pos);

  virtual inline double operator()(QickArray &pos){ 
    return valuePt(pos);
  }
};

//=================================================================
// u(r_ij)=.5*r_ij/(1+beta*r_ij)
//=================================================================
inline double Correlation:: valuePt(QickArray &pos){

  int dim=pos.no_of_elements();
  double r_ij=0;
  for(int i=0; i!=dim; i++)
    r_ij+=sqr(pos(i));
  r_ij = sqrt(r_ij);

  return .5*r_ij/(1+beta*r_ij);

}

//=================================================================
// wave function for helium
//=================================================================
class Helium : public virtual Wave{

private:
  Wave1s wave1s;
  Correlation corr;
  QickArray dummy;
public:

  virtual inline void setParams(QickArray& params_){ 
    wave1s.setParams(params_);
    corr.setParams(params_);
    dummy.redim(2);
  }

  virtual double valuePt(QickArray &pos);

  virtual inline double operator()(QickArray &pos){ 
    return valuePt(pos);
  }
};

//=================================================================
// functor making sure distance is correct, used for homogenous
// systems
//=================================================================
class Dist : public virtual Func{

private:
  QickArray tmp_pos;

public:

  virtual QickArray& valuePt(QickArray &pos, double l);
  
  virtual inline QickArray& operator()(QickArray &pos, double l){ 
    
    return valuePt(pos, l);
  }

};

//=================================================================
//=================================================================
inline QickArray& Dist:: valuePt(QickArray &pos, double l){

  int dim=pos.no_of_elements();
  tmp_pos.redim(dim);
  tmp_pos = pos;
  for(int i=0; i!=dim; i++)
    if(pos(i) > l/2)
      tmp_pos(i) = l-pos(i);
  return tmp_pos;
}

//=================================================================
// fuctor for trial wave function for N bosons, a la DuBois;
// gaussian times sum_{i < j}{ 1-a/r_ij, r_ij > a; 0, r_ij < a} 
//=================================================================
class DuBoisWave : public virtual Wave{

protected:
  int moved_particle;
  double alpha, beta, a, d, lambda;
  QickArray tmp_pos;

public:

  virtual inline void setParams(QickArray& params_){ 
    alpha          = params_(0);
    beta           = params_(1);
    a              = params_(2);
    moved_particle = int(params_(3));
    d              = params_(4);
    lambda         = params_(5);
  }

  virtual double valuePt(QickArray &pos);

  virtual inline double operator()(QickArray &pos){ 
    return valuePt(pos);
  }
};

//=================================================================
//Function for Gaussian part of DuBois wave function, only 
//for testing
//=================================================================
class DuBoisGauss : public virtual DuBoisWave{

public:

  virtual double valuePt(QickArray &pos);

  virtual inline double operator()(QickArray &pos){ 
    return valuePt(pos);
  }
};

//=================================================================
//Function for Jastrow part of DuBois wave function, only 
//for testing
//=================================================================
class DuBoisJastrow : public virtual DuBoisWave{

public:

  virtual inline void setParams(QickArray& params_){ 
    a              = params_(2);
    d              = params_(4);
  }

  virtual double valuePt(QickArray &pos);

  virtual inline double operator()(QickArray &pos){ 
    return valuePt(pos);
  }
};

//=================================================================
// DuBois trial wave function for N bosons, returning all of wave
// function, not just one-particle part as usually done, only
// for testing
//=================================================================
class DuBoisWaveAll : public virtual Wave{

private:
  double alpha, beta, a, d, lambda;
  QickArray tmp_pos;

public:

  virtual inline void setParams(QickArray& params_){ 
    alpha          = params_(0);
    beta           = params_(1);
    a              = params_(2);
    d              = params_(4);
    lambda         = params_(5);
  }

  virtual double valuePt(QickArray &pos);

  virtual inline double operator()(QickArray &pos){ 
    return valuePt(pos);
  }
};


//=================================================================
// analytic quantum force for duBois wave function
//=================================================================
class DuBoisQForce : public virtual Func{

private:
  int moved_particle;
  double alpha, beta, a, d, lambda;
  double u(double r);
  double du(double r);
  QickArray pos_ij, dF, dG;

public:

  virtual inline void setParams(QickArray& params_){ 
    alpha          = params_(0);
    beta           = params_(1);
    a              = params_(2);
    d              = params_(4);
    lambda         = params_(5);
  }

  virtual QickArray& valuePt(QickArray &pos, bool dummy);

  virtual inline QickArray& operator()(QickArray &pos, bool dummy){ 
    return valuePt(pos, dummy);
  }
};

inline double DuBoisQForce::u(double r){

  return log(1-a/r);
}
inline double DuBoisQForce::du(double r){

  return a/(r*(r-a));
}

//=================================================================
// analytic local energy for duBois wave function and harm. osc.
// with distortion lambda in z-direction
//=================================================================
class DuBoisLocalEnergy : public virtual Func{

private:
  double alpha, beta, a, d;
  double lambda; //omega_z/omega_rho
  double jackFeen;
  DuBoisJastrow fun;
  Nabla nabla;
  Nabla2 nabla2;
  double u(double r);
  double du(double r);
  double ddu(double r);
  QickArray pos_ij, pos_i, pos_ik, tmp_pos;

public:

  virtual inline void setParams(QickArray& params_){ 
    alpha    = params_(0);
    beta     = params_(1);
    a        = params_(2);
    d        = params_(4);
    lambda   = params_(5);
    jackFeen = params_(7);
    fun.setParams(params_);
  }

  virtual double valuePt(QickArray &pos);

  virtual double valuePt(QickArray &pos, int dummy);

  virtual inline double operator()(QickArray &pos){ 
    if(jackFeen==1)
      return valuePt(pos,1);
    else
      return valuePt(pos);
  }
};

inline double DuBoisLocalEnergy::u(double r){

  return log(1-a/r);
}
inline double DuBoisLocalEnergy::du(double r){

  return a/(r*(r-a));
}

inline double DuBoisLocalEnergy::ddu(double r){

  return a*(a-2*r)/sqr(r*(r-a));
}

//=================================================================
// analytic potential energy for bosons in harm. osc. with
// distortion in z-direction
//=================================================================
class DuBoisPotentialEnergy : public virtual Func{

private:
  double lambda; //omega_z/omega_rho

public:

  virtual inline void setParams(QickArray& params_){ 
    lambda = params_(5);
  }

  virtual double valuePt(QickArray &pos);

  virtual inline double operator()(QickArray &pos){ 
    return valuePt(pos);
  }
};

//=================================================================
// fuctor for trial wave function for N bosons, a la DuBois;
// gaussian times sum_{i < j}{ 1-a/r_ij, r_ij > a; 0, r_ij < a} 
//=================================================================
class DuBoisVortexWave : public virtual Wave{

protected:
  int moved_particle;
  double alpha, beta, a, d, lambda;
  QickArray tmp_pos;

public:

  virtual inline void setParams(QickArray& params_){ 
    alpha          = params_(0);
    beta           = params_(1);
    a              = params_(2);
    moved_particle = int(params_(3));
    d              = params_(4);
    lambda         = params_(5);
  }

  virtual double valuePt(QickArray &pos);

  virtual inline double operator()(QickArray &pos){ 
    return valuePt(pos);
  }
};

//=================================================================
// DuBois trial wave function for N bosons, returning all of wave
// function, not just one-particle part as usually done, only
// for testing
//=================================================================
class DuBoisVortexWaveAll : public virtual Wave{

private:
  double alpha, beta, a, d, lambda;
  QickArray tmp_pos;

public:

  virtual inline void setParams(QickArray& params_){ 
    alpha          = params_(0);
    beta           = params_(1);
    a              = params_(2);
    d              = params_(4);
    lambda         = params_(5);
  }

  virtual double valuePt(QickArray &pos);

  virtual inline double operator()(QickArray &pos){ 
    return valuePt(pos);
  }
};

//=================================================================
// analytic quantum force for duBoisVortex wave function
//=================================================================
class DuBoisVortexQForce : public virtual Func{

private:
  double alpha, beta, a, d, lambda;
  double u(double r);
  double du(double r);
  QickArray pos_ij, dF, dG;

public:

  virtual inline void setParams(QickArray& params_){ 
    alpha          = params_(0);
    beta           = params_(1);
    a              = params_(2);
    d              = params_(4);
    lambda         = params_(5);
  }

  virtual QickArray& valuePt(QickArray &pos, bool dummy);

  virtual inline QickArray& operator()(QickArray &pos, bool dummy){ 
    return valuePt(pos, dummy);
  }
};

inline double DuBoisVortexQForce::u(double r){

  return log(1-a/r);
}
inline double DuBoisVortexQForce::du(double r){

  return a/(r*(r-a));
}

//=================================================================
// analytic local energy for duBois wave function and harm. osc.
// with distortion lambda in z-direction
//=================================================================
class DuBoisVortexLocalEnergy : public virtual Func{

private:
  double alpha, beta, a, d;
  double lambda, kappa; //omega_z/omega_rho
  double jackFeen;
  DuBoisJastrow fun;
  Nabla nabla;
  Nabla2 nabla2;
  double u(double r);
  double du(double r);
  double ddu(double r);
  QickArray pos_ij, pos_i, pos_ik, tmp_pos;

public:

  virtual inline void setParams(QickArray& params_){ 
    alpha    = params_(0);
    beta     = params_(1);
    a        = params_(2);
    d        = params_(4);
    lambda   = params_(5);
    kappa    = params_(6);
    jackFeen = params_(7);
    fun.setParams(params_);
  }

  virtual double valuePt(QickArray &pos);

  virtual double valuePt(QickArray &pos, int dummy);

  virtual inline double operator()(QickArray &pos){ 
    if(jackFeen==1.)
      return valuePt(pos,1);
    else
      return valuePt(pos);
  }
};

inline double DuBoisVortexLocalEnergy::u(double r){

  return log(1-a/r);
}

inline double DuBoisVortexLocalEnergy::du(double r){

  return a/(r*(r-a));
}

inline double DuBoisVortexLocalEnergy::ddu(double r){

  return a*(a-2*r)/sqr(r*(a-r));
}

//=================================================================
// analytic potential energy for bosons in harm. osc. with
// distortion in z-direction
//=================================================================
class DuBoisVortexPotentialEnergy : public virtual Func{

private:
  double lambda, kappa; //omega_z/omega_rho

public:

  virtual inline void setParams(QickArray& params_){ 
    lambda = params_(5);
    kappa  = params_(6);
  }

  virtual double valuePt(QickArray &pos);

  virtual inline double operator()(QickArray &pos){ 
    return valuePt(pos);
  }
};

//=================================================================
// analytic local energy for duBois wave function and harm. osc.
// with distortion lambda in z-direction
//=================================================================
class DuBoisVortexExcitationEnergy : public virtual Func{

private:
  double alpha, beta, a, d;
  double lambda, kappa2; //omega_z/omega_rho
  double jackFeen;
  DuBoisJastrow fun;
  Nabla nabla;
  Nabla2 nabla2;

public:

  virtual inline void setParams(QickArray& params_){ 
    alpha    = params_(0);
    beta     = params_(1);
    a        = params_(2);
    d        = params_(4);
    lambda   = params_(5);
    kappa2   = sqr(params_(6));
    jackFeen = params_(7);
    fun.setParams(params_);
  }

  virtual double valuePt(QickArray &pos);

  virtual inline double operator()(QickArray &pos){ 
      return valuePt(pos);
  }
};


//=================================================================
// functor for trial wave function for N bosons, a la DuBois;
// f_vor times gaussian times 
// sum_{i < j}{ 1-a/r_ij, r_ij > a; 0, r_ij < a} 
//=================================================================
class ReattoVortexWave : public virtual Wave{

protected:
  int moved_particle;
  double alpha, beta2, a, d, lambda;
  QickArray tmp_pos;
  double u_vor(double rho2);

public:

  virtual inline void setParams(QickArray& params_){ 
    alpha          = params_(0);
    beta2          = sqr(params_(1));
    a              = params_(2);
    moved_particle = int(params_(3));
    d              = params_(4);
    lambda         = params_(5);
  }

  virtual double valuePt(QickArray &pos);

  virtual inline double operator()(QickArray &pos){ 
    return valuePt(pos);
  }
};

//u_vor = log(1-exp(-(rho/beta)^2))
inline double ReattoVortexWave::u_vor(double rho2){
  return log(1-exp(-rho2/beta2));
}

//=================================================================
// Reatto trial wave function for N bosons, returning all of wave
// function, not just one-particle part as usually done, only
// for testing
//=================================================================
class ReattoVortexWaveAll : public virtual Wave{

private:
  double alpha, beta2, a, d, lambda;
  QickArray tmp_pos;
  double u_vor(double rho2);

public:

  virtual inline void setParams(QickArray& params_){ 
    alpha          = params_(0);
    beta2          = sqr(params_(1));
    a              = params_(2);
    d              = params_(4);
    lambda         = params_(5);
  }

  virtual double valuePt(QickArray &pos);

  virtual inline double operator()(QickArray &pos){ 
    return valuePt(pos);
  }
};

//u_vor = log(1-exp(-(rho/beta)^2))
inline double ReattoVortexWaveAll::u_vor(double rho2){
  return log(1-exp(-rho2/beta2));
}

//=================================================================
// analytic quantum force for duBoisVortex wave function
//=================================================================
class ReattoVortexQForce : public virtual Func{

private:
  double alpha, beta2, a, d, lambda;
  double u(double r);
  double du(double r);
  double u_vor(double r2);
  double du_vor(double r2);
  QickArray pos_ij, dF, dG;

public:

  virtual inline void setParams(QickArray& params_){ 
    alpha          = params_(0);
    beta2          = sqr(params_(1));
    a              = params_(2);
    d              = params_(4);
    lambda         = params_(5);
  }

  virtual QickArray& valuePt(QickArray &pos, bool dummy);

  virtual inline QickArray& operator()(QickArray &pos, bool dummy){ 
    return valuePt(pos, dummy);
  }
};

inline double ReattoVortexQForce::u(double r){

  return log(1-a/r);
}
inline double ReattoVortexQForce::du(double r){

  return a/(r*(r-a));
}
inline double ReattoVortexQForce::u_vor(double r2){

  return log(1-exp(-r2/beta2));
}
inline double ReattoVortexQForce::du_vor(double r2){

  return 2*sqrt(r2)/beta2*(1/(1-exp(-r2/beta2))-1);
}

//=================================================================
// analytic local energy for duBois wave function and harm. osc.
// with distortion lambda in z-direction
//=================================================================
class ReattoVortexLocalEnergy : public virtual Func{

private:
  double alpha, beta2, a, d, x;
  double lambda, kappa; //omega_z/omega_rho
  double jackFeen;
  DuBoisJastrow fun;
  Nabla nabla;
  Nabla2 nabla2;
  double u(double r);
  double du(double r);
  double ddu(double r);
  double u_vor(double r2);
  double du_vor(double r2);
  double ddu_vor(double r2);
  QickArray pos_ij, pos_i, pos_ik, tmp_pos;

public:

  virtual inline void setParams(QickArray& params_){ 
    alpha    = params_(0);
    beta2    = sqr(params_(1));
    a        = params_(2);
    d        = params_(4);
    lambda   = params_(5);
    kappa    = params_(6);
    jackFeen = params_(7);
    fun.setParams(params_);
  }

  virtual double valuePt(QickArray &pos);

  virtual double valuePt(QickArray &pos, int dummy);

  virtual inline double operator()(QickArray &pos){ 
    if(jackFeen==1.)
      return valuePt(pos,1);
    else
      return valuePt(pos);
  }
};

inline double ReattoVortexLocalEnergy::u(double r){

  return log(1-a/r);
}

inline double ReattoVortexLocalEnergy::du(double r){

  return a/(r*(r-a));
}

inline double ReattoVortexLocalEnergy::ddu(double r){

  return a*(a-2*r)/sqr(r*(a-r));
}

inline double ReattoVortexLocalEnergy::u_vor(double r2){

  return log(1-exp(-r2/beta2));
}
inline double ReattoVortexLocalEnergy::du_vor(double r2){

  return 2*sqrt(r2)/beta2*(1/(1-exp(-r2/beta2))-1);
}
inline double ReattoVortexLocalEnergy::ddu_vor(double r2){

  x=1-exp(-r2/beta2);
  return (1/x-1)*(2/beta2-4*r2/sqr(beta2)/x);
}

//=================================================================
// analytic potential energy for bosons in harm. osc. with
// distortion in z-direction
//=================================================================
class ReattoVortexPotentialEnergy : public virtual Func{

private:
  double lambda, kappa; //omega_z/omega_rho

public:

  virtual inline void setParams(QickArray& params_){ 
    lambda = params_(5);
    kappa  = params_(6);
  }

  virtual double valuePt(QickArray &pos);

  virtual inline double operator()(QickArray &pos){ 
    return valuePt(pos);
  }
};

//=================================================================
// analytic local energy for duBois wave function and harm. osc.
// with distortion lambda in z-direction
//=================================================================
class ReattoVortexExcitationEnergy : public virtual Func{

private:
  double alpha, beta2, a, d, x;
  double lambda, kappa2; //lambda=omega_z/omega_rho
  double jackFeen;
  DuBoisJastrow fun;

public:

  virtual inline void setParams(QickArray& params_){ 
    alpha    = params_(0);
    beta2    = sqr(params_(1));
    a        = params_(2);
    d        = params_(4);
    lambda   = params_(5);
    kappa2   = sqr(params_(6));
    jackFeen = params_(7);
    fun.setParams(params_);
  }

  virtual double valuePt(QickArray &pos);

  virtual inline double operator()(QickArray &pos){ 
      return valuePt(pos);
  }
};


//=================================================================
// functor for trial wave function for N bosons, a la DuBois;
// f_vor times gaussian times 
// sum_{i < j}{ 1-a/r_ij, r_ij > a; 0, r_ij < a} 
//=================================================================
class Reatto2VortexWave : public virtual Wave{

protected:
  int moved_particle;
  double alpha, beta, a, d, lambda;
  QickArray tmp_pos;
  double u_vor(double rho2);

public:

  virtual inline void setParams(QickArray& params_){ 
    alpha          = params_(0);
    beta           = params_(1);
    a              = params_(2);
    moved_particle = int(params_(3));
    d              = params_(4);
    lambda         = params_(5);
  }

  virtual double valuePt(QickArray &pos);

  virtual inline double operator()(QickArray &pos){ 
    return valuePt(pos);
  }
};

//u_vor = log(1-exp(-(rho/beta)))
inline double Reatto2VortexWave::u_vor(double rho){
  return log(1-exp(-rho/beta));
}

//=================================================================
// Reatto2 trial wave function for N bosons, returning all of wave
// function, not just one-particle part as usually done, only
// for testing
//=================================================================
class Reatto2VortexWaveAll : public virtual Wave{

private:
  double alpha, beta, a, d, lambda;
  QickArray tmp_pos;
  double u_vor(double rho);

public:

  virtual inline void setParams(QickArray& params_){ 
    alpha          = params_(0);
    beta           = params_(1);
    a              = params_(2);
    d              = params_(4);
    lambda         = params_(5);
  }

  virtual double valuePt(QickArray &pos);

  virtual inline double operator()(QickArray &pos){ 
    return valuePt(pos);
  }
};

//u_vor = log(1-exp(-(rho/beta)^2))
inline double Reatto2VortexWaveAll::u_vor(double rho){
  return log(1-exp(-rho/beta));
}

//=================================================================
// analytic quantum force for duBoisVortex wave function
//=================================================================
class Reatto2VortexQForce : public virtual Func{

private:
  double alpha, beta, a, d, lambda;
  double u(double r);
  double du(double r);
  double u_vor(double r);
  double du_vor(double r);
  QickArray pos_ij, dF, dG;

public:

  virtual inline void setParams(QickArray& params_){ 
    alpha          = params_(0);
    beta           = params_(1);
    a              = params_(2);
    d              = params_(4);
    lambda         = params_(5);
  }

  virtual QickArray& valuePt(QickArray &pos, bool dummy);

  virtual inline QickArray& operator()(QickArray &pos, bool dummy){ 
    return valuePt(pos, dummy);
  }
};

inline double Reatto2VortexQForce::u(double r){

  return log(1-a/r);
}
inline double Reatto2VortexQForce::du(double r){

  return a/(r*(r-a));
}
inline double Reatto2VortexQForce::u_vor(double r){

  return log(1-exp(-r/beta));
}
inline double Reatto2VortexQForce::du_vor(double r){

  return 1/beta*(1/(1-exp(-r/beta))-1);
}

//=================================================================
// analytic local energy for duBois wave function and harm. osc.
// with distortion lambda in z-direction
//=================================================================
class Reatto2VortexLocalEnergy : public virtual Func{

private:
  double alpha, beta, a, d, x;
  double lambda, kappa; //omega_z/omega_rho
  double jackFeen;
  DuBoisJastrow fun;
  Nabla nabla;
  Nabla2 nabla2;
  double u(double r);
  double du(double r);
  double ddu(double r);
  double u_vor(double r);
  double du_vor(double r);
  double ddu_vor(double r);
  QickArray pos_ij, pos_i, pos_ik, tmp_pos;

public:

  virtual inline void setParams(QickArray& params_){ 
    alpha    = params_(0);
    beta     = params_(1);
    a        = params_(2);
    d        = params_(4);
    lambda   = params_(5);
    kappa    = params_(6);
    jackFeen = params_(7);
    fun.setParams(params_);
  }

  virtual double valuePt(QickArray &pos);

  virtual double valuePt(QickArray &pos, int dummy);

  virtual inline double operator()(QickArray &pos){ 
    if(jackFeen==1.)
      return valuePt(pos,1);
    else
      return valuePt(pos);
  }
};

inline double Reatto2VortexLocalEnergy::u(double r){

  return log(1-a/r);
}

inline double Reatto2VortexLocalEnergy::du(double r){

  return a/(r*(r-a));
}

inline double Reatto2VortexLocalEnergy::ddu(double r){

  return a*(a-2*r)/sqr(r*(a-r));
}

inline double Reatto2VortexLocalEnergy::u_vor(double r){

  return log(1-exp(-r/beta));
}
inline double Reatto2VortexLocalEnergy::du_vor(double r){

  return 1/beta*(1/(1-exp(-r/beta))-1);
}
inline double Reatto2VortexLocalEnergy::ddu_vor(double r){

  x=1-exp(-r/beta);
  return (1-1/x)/(sqr(beta)*x);
}

//=================================================================
// analytic potential energy for bosons in harm. osc. with
// distortion in z-direction
//=================================================================
class Reatto2VortexPotentialEnergy : public virtual Func{

private:
  double lambda, kappa; //omega_z/omega_rho

public:

  virtual inline void setParams(QickArray& params_){ 
    lambda = params_(5);
    kappa  = params_(6);
  }

  virtual double valuePt(QickArray &pos);

  virtual inline double operator()(QickArray &pos){ 
    return valuePt(pos);
  }
};

//=================================================================
// analytic local energy for duBois wave function and harm. osc.
// with distortion lambda in z-direction
//=================================================================
class Reatto2VortexExcitationEnergy : public virtual Func{

private:
  double alpha, beta, a, d, x;
  double lambda, kappa2; //lambda=omega_z/omega_rho
  double jackFeen;
  DuBoisJastrow fun;

public:

  virtual inline void setParams(QickArray& params_){ 
    alpha    = params_(0);
    beta     = params_(1);
    a        = params_(2);
    d        = params_(4);
    lambda   = params_(5);
    kappa2   = sqr(params_(6));
    jackFeen = params_(7);
    fun.setParams(params_);
  }

  virtual double valuePt(QickArray &pos);

  virtual inline double operator()(QickArray &pos){ 
      return valuePt(pos);
  }
};

//==================================================================
//Finds K given d to fullfill healing conditions for Polls Jastrow
//Function using Newton-Raphson
//==================================================================
class PollsKd : public virtual Func{

private:
  
  inline double F(double K, double d, double a){ 
    //return K*d - tan(K*(d-a)); }
    return 1/(K*d) - 1/tan(K*(d-a)); }
  
  inline double J(double K, double d, double a){ 
    //return d-(d-a)/sqr(cos(K*(d-a))); }
    return -1/(sqr(K)*d)+(d-a)*(1/sqr(cos(K*(d-a)))); }

public:

  virtual double valuePt(double d, double a);

  virtual inline double operator()(double d, double a){
    return valuePt(d,a);
  }

};


//=================================================================
// trial wave function for N bosons, a la Polls
// gaussian times sum_{i < j}{ d/rij*sin(K(r_ij-a))/sin(K(d-a)) }
//=================================================================
class PollsWave : public virtual Wave{

protected:
  int moved_particle;
  double alpha, beta, a, d, K, lambda;
  PollsKd Kd;  
  QickArray tmp_pos;

public:

  virtual inline void setParams(QickArray& params_){ 
    alpha = params_(0);
    beta  = params_(1);
    a     = params_(2);
    if(d!=params_(4)){
      d = params_(4);
      K = Kd(d,a);
      cerr << K << endl;
    }
    moved_particle = int(params_(3));
    lambda         = params_(5);
  }

  virtual double valuePt(QickArray &pos);

  virtual inline double operator()(QickArray &pos){ 
    return valuePt(pos);
  }
};

//=================================================================
//Function for Gaussian part of Polls wave function, only 
//for testing
//=================================================================
class PollsGauss : public virtual PollsWave{

public:

  virtual double valuePt(QickArray &pos);

  virtual inline double operator()(QickArray &pos){ 
    return valuePt(pos);
  }
};

//=================================================================
//Function for Correlation part of Polls wave function, only 
//for testing
//=================================================================
class PollsJastrow : public virtual PollsWave{
private:
  double K,a;
  PollsKd Kd;
public:

  virtual inline void setParams(QickArray& params_){ 
    a = params_(2);
    if(d!=params_(4)){
      d = params_(4);
      K = Kd(d,a);
    }
  }

  virtual double valuePt(QickArray &pos);

  virtual inline double operator()(QickArray &pos){ 
    return valuePt(pos);
  }
};

//=================================================================
// Polls trial wave function for N bosons, returning all of wave
// function, not just one-particle part as usually done, only
// for testing
//=================================================================
class PollsWaveAll : public virtual Wave{

private:
  double alpha, beta, d, a, K, lambda;
  PollsKd Kd;
  QickArray tmp_pos;

public:

  virtual inline void setParams(QickArray& params_){ 
    alpha = params_(0);
    beta  = params_(1);
    a = params_(2);
    if(d!=params_(4)){
      d = params_(4);
      K = Kd(d,a);
    }
    lambda         = params_(5);
  }

  virtual double valuePt(QickArray &pos);

  virtual inline double operator()(QickArray &pos){ 
    return valuePt(pos);
  }
};


//=================================================================
// analytic quantum force for Polls wave function
//=================================================================
class PollsQForce : public virtual Func{

private:
  int moved_particle;
  double alpha, beta, d, a, K, lambda;
  PollsKd Kd;
  double u(double r);
  double du(double r);
  QickArray pos_ij, dF, dG;

public:

  virtual inline void setParams(QickArray& params_){ 
    alpha = params_(0);
    beta  = params_(1);
    a     = params_(2);
    moved_particle = int(params_(3));
    if(d!=params_(4)){
      d = params_(4);
      K = Kd(d,a);
    }
    lambda         = params_(5);
  }

  virtual QickArray& valuePt(QickArray &pos, bool dummy);

  virtual inline QickArray& operator()(QickArray &pos, bool dummy){ 
    return valuePt(pos, dummy);
  }
};

inline double PollsQForce::u(double r){

  if(r<d)
    return log(d/r*sin(K*(r-a))/sin(K*(d-a)));
  else
    return 0.;

}

inline double PollsQForce::du(double r){

  if(r<d)
    return K/tan(K*(r-a))-1/r;
  else
    return 0.;

}

//=================================================================
// analytic local energy for duBois wave function and harm. osc.
// with distortion lambda in z-direction
//=================================================================
class PollsLocalEnergy : public virtual Func{

private:
  double alpha, beta, d, a, K;
  double lambda; //omega_z/omega_rho
  double jackFeen;
  DuBoisJastrow fun;
  Nabla nabla;
  Nabla2 nabla2;
  PollsKd Kd;
  double u(double r);
  double du(double r);
  double ddu(double r);
  QickArray pos_ij, pos_i, pos_ik, tmp_pos;
  
public:

  virtual inline void setParams(QickArray& params_){ 
    alpha = params_(0);
    beta  = params_(1);
    a     = params_(2);
    if(d!=params_(4)){
      d = params_(4);
      K = Kd(d,a);
    }
    lambda   = params_(5);
    jackFeen = params_(7);
    fun.setParams(params_);
  }

  virtual double valuePt(QickArray &pos);

  virtual double valuePt(QickArray &pos, int dummy);

  virtual inline double operator()(QickArray &pos){ 
    if(jackFeen==1.)
      return valuePt(pos,true);
    else
      return valuePt(pos);
  }
};

inline double PollsLocalEnergy::u(double r){

  if(r<d)
    return log(d/r*sin(K*(r-a))/sin(K*(d-a)));
  else
    return 0.;

}

inline double PollsLocalEnergy::du(double r){

  if(r<d)
    return K/tan(K*(r-a))-1/r;
  else
    return 0.;

}

inline double PollsLocalEnergy::ddu(double r){

  if(r<d)
    return 1/sqr(r) - sqr(K/sin(K*(r-a)));
  else
    return 0.;

}


//=================================================================
// trial wave function for N bosons, a la Polls with vortex in z
// 1p wave times sum_{i < j}{ d/rij*sin(K(r_ij-a))/sin(K(d-a)) }
//=================================================================
class PollsVortexWave : public virtual Wave{

protected:
  int moved_particle;
  double alpha, beta, a, d, K, lambda;
  PollsKd Kd;  
  QickArray tmp_pos;

public:

  virtual inline void setParams(QickArray& params_){ 
    alpha = params_(0);
    beta  = params_(1);
    a     = params_(2);
    moved_particle = int(params_(3));
    if(d!=params_(4)){
      d = params_(4);
      K = Kd(d,a);
    }
    lambda   = params_(5);
  }

  virtual double valuePt(QickArray &pos);

  virtual inline double operator()(QickArray &pos){ 
    return valuePt(pos);
  }
};


//=================================================================
// Polls trial wave function for N bosons, returning all of wave
// function, not just one-particle part as usually done, only
// for testing
//=================================================================
class PollsVortexWaveAll : public virtual Wave{

private:
  double alpha, beta, d, a, K, lambda;
  PollsKd Kd;
  QickArray tmp_pos;

public:

  virtual inline void setParams(QickArray& params_){ 
    alpha = params_(0);
    beta  = params_(1);
    a = params_(2);
    if(d!=params_(4)){
      d = params_(4);
      K = Kd(d,a);
    }
    lambda   = params_(5);
  }

  virtual double valuePt(QickArray &pos);

  virtual inline double operator()(QickArray &pos){ 
    return valuePt(pos);
  }
};


//=================================================================
// analytic quantum force for Polls wave function and vortex
//=================================================================
class PollsVortexQForce : public virtual Func{

private:
  int moved_particle;
  double alpha, beta, d, a, K, lambda;
  PollsKd Kd;
  double u(double r);
  double du(double r);
  QickArray pos_ij, dF, dG;

public:

  virtual inline void setParams(QickArray& params_){ 
    alpha = params_(0);
    beta  = params_(1);
    a     = params_(2);
    moved_particle = int(params_(3));
    if(d!=params_(4)){
      d = params_(4);
      K = Kd(d,a);
    }
    lambda   = params_(5);
  }

  virtual QickArray& valuePt(QickArray &pos, bool dummy);

  virtual inline QickArray& operator()(QickArray &pos, bool dummy){ 
    return valuePt(pos, dummy);
  }
};

inline double PollsVortexQForce::u(double r){

  if(r<d)
    return log(d/r*sin(K*(r-a))/sin(K*(d-a)));
  else
    return 0.;

}

inline double PollsVortexQForce::du(double r){

  if(r<d)
    return K/tan(K*(r-a))-1/r;
  else
    return 0.;

}

//=================================================================
// analytic local energy for duBois wave function and harm. osc.
// with distortion lambda in z-direction and vortex
//=================================================================
class PollsVortexLocalEnergy : public virtual Func{

private:
  double alpha, beta, d, a, K, kappa;
  double lambda; //omega_z/omega_rho
  double jackFeen;
  PollsKd Kd;
  double u(double r);
  double du(double r);
  double ddu(double r);
  QickArray pos_ij, pos_i, pos_ik, tmp_pos;

public:

  virtual inline void setParams(QickArray& params_){ 
    alpha = params_(0);
    beta  = params_(1);
    a     = params_(2);
    if(d!=params_(4)){
      d = params_(4);
      K = Kd(d,a);
    }
    lambda   = params_(5);
    kappa    = params_(6);
    jackFeen = params_(7);
  }

  virtual double valuePt(QickArray &pos);

  virtual double valuePt(QickArray &pos, int dummy);

  virtual inline double operator()(QickArray &pos){ 
    if(jackFeen==1.)
      return valuePt(pos,true);
    else
      return valuePt(pos);
  }
};

inline double PollsVortexLocalEnergy::u(double r){

  if(r<d)
    return log(d/r*sin(K*(r-a))/sin(K*(d-a)));
  else
    return 0.;

}

inline double PollsVortexLocalEnergy::du(double r){

  if(r<d)
    return K/tan(K*(r-a))-1/r;
  else
    return 0.;

}

inline double PollsVortexLocalEnergy::ddu(double r){

  if(r<d)
    return 1/sqr(r) - sqr(K/sin(K*(r-a)));
  else
    return 0.;

}

//=================================================================
//Movement of particle
//=================================================================
class Move : public virtual Func{
  
private:

  double step_size;
  QickArray pos_;

public:

  virtual inline void setParams(QickArray& params_){
    step_size = params_(1);
  }

  virtual bool getImpossible(){ return false; }

  virtual QickArray& valuePt(QickArray& pos, int i, Random& ran);

  virtual QickArray& valuePt(QickArray& pos, int i, QickArray& x, 
			 Random& ran);

  virtual inline QickArray& operator()(QickArray& pos, int i,  
				   Random& ran){
    return valuePt(pos, i, ran);
  }

  virtual inline QickArray& operator()(QickArray& pos, int i, 
				   QickArray& x, Random& ran){
    return valuePt(pos, i, x, ran);
  }

};

//=================================================================
//Movement of particle with uniform distribution; 
//step_length*(+/-.5)
//=================================================================
inline QickArray& Move::valuePt(QickArray& pos, int i, Random& ran){

  int cols, rows, dummy;
  dummy = pos.get_dimension_info(&cols,&rows);

  pos_.redim(rows);
  for(int j=0; j!=rows; j++)
    pos_(j) = pos(i,j)+(ran.ran1()-.5)*step_size;
  return pos_;
}

//=================================================================
//Movement of particle with uniform distribution; 
//step_length*(+/-.5)+x
//=================================================================
inline QickArray& Move::valuePt(QickArray& pos, int i, QickArray& x,
			    Random& ran){

  int cols, rows, dummy;
  dummy = pos.get_dimension_info(&cols,&rows);

  pos_.redim(rows);
  for(int j=0; j!=rows; j++)
    pos_(j) = pos(i,j)+x(i,j)+ran.gran(1,0)*step_size;
  return pos_;

}



//=================================================================
//DuBois movement of particle
//=================================================================
class DuBoisMove : public virtual Move{

private:

  double a, step_length;
  QickArray pos_, prepos_;
  bool impossible;

public:

  virtual inline void setParams(QickArray& params_){
    a = params_(0);
    step_length = params_(1);
  }

  virtual bool getImpossible(){ return impossible; }

  virtual QickArray& valuePt(QickArray& pos, int i, Random& ran);

  virtual QickArray& valuePt(QickArray& pos, int i, QickArray& x, 
			 Random& ran);

  virtual inline QickArray& operator()(QickArray& pos, int i, 
				   Random& ran){
    return valuePt(pos, i, ran);
  }

  virtual inline QickArray& operator()(QickArray& pos, int i,  
				   QickArray& x, Random& ran){
    return valuePt(pos, i, x, ran);
  }

};


//=================================================================
// trial wave function for liquid He4 given simulation 
// cube with sides of length l
//=================================================================
class LiqHe4Wave : public virtual Wave{

private:
  int moved_particle;
  double correlation5, l;
  Dist distance;

public:

  virtual inline void setParams(QickArray& params_){ 
    l=params_(0);
    correlation5=params_(1);
    moved_particle=int(params_(2));
  }

  virtual double valuePt(QickArray &pos);

  virtual inline double operator()(QickArray &pos){ 
    return valuePt(pos);
  }
};


//=================================================================
// trial wave function for liquid He4 given simulation 
// cube with sides of length l
//=================================================================
class LiqHe4WaveAll : public virtual Wave{

private:
  int moved_particle;
  double correlation5, l;
  Dist distance;

public:

  virtual inline void setParams(QickArray& params_){ 
    l=params_(0);
    correlation5=params_(1);
  }

  virtual double valuePt(QickArray &pos);

  virtual inline double operator()(QickArray &pos){ 
    return valuePt(pos);
  }
};

//=================================================================
// Numeric quantum force for liquid He4 with
// lennard-jones potential
//=================================================================
class LiqHe4QForceNum : public virtual Func{

protected:
  int moved_particle;
  double correlation5, l;
  Dist distance;
  QickArray q_force;

public:

  virtual inline void setParams(QickArray& params_){ 
    l=params_(0);
    correlation5=params_(1);
    moved_particle=int(params_(2));
  }

  virtual QickArray& valuePt(QickArray &pos, bool dummy);

  virtual inline QickArray& operator()(QickArray &pos, bool dummy){ 
    return valuePt(pos, dummy);
  }
};

//=================================================================
// Numeric quantum force for liquid He4 with
// lennard-jones potential
//=================================================================
class LiqHe4WfExp : public virtual LiqHe4QForceNum{

public:

  virtual inline void setParams(QickArray& params_){ 
    l=params_(0);
    correlation5=params_(1);
  }

  virtual double valuePt(QickArray &pos);

  virtual inline double operator()(QickArray &pos){ 
    return valuePt(pos);
  }
};

//=================================================================
// analytic quantum force for liquid He4 with
// lennard-jones potential
//=================================================================
class LiqHe4QForce : public virtual Func{

private:
  int moved_particle;
  double correlation5, l;
  Dist distance;
  QickArray q_force;

public:

  virtual inline void setParams(QickArray& params_){ 
    l=params_(0);
    correlation5=params_(1);
    moved_particle=int(params_(2));
  }

  virtual QickArray& valuePt(QickArray &pos, bool dummy);

  virtual inline QickArray& operator()(QickArray &pos, bool dummy){ 
    return valuePt(pos, dummy);
  }
};

//=================================================================
// analytic local energy for liquid He4 with
// lennard-jones potential
//=================================================================
class LiqHe4LocalEnergy : public virtual Func{

private:
  double correlation5, l;
  Dist distance;

public:

  virtual inline void setParams(QickArray& params_){ 
    l=params_(0);
    correlation5=params_(1);
  }

  virtual double valuePt(QickArray &pos);

  virtual inline double operator()(QickArray &pos){ 
    return valuePt(pos);
  }
};

//=================================================================
//functor returning potential, obsolete
//=================================================================
class Potential : public virtual Func{

public:

  virtual double valuePt(QickArray &pos);

  virtual inline double operator()(QickArray &pos){ 
    return valuePt(pos);
  }

};

//=================================================================
// fuctor for trial wave function for N bosons, a la Mateusz;
// gaussian times sum_{i < j}{ 1-a/r_ij, r_ij > a; 0, r_ij < a} 
//=================================================================
class MateuszWave : public virtual Wave{

protected:
  int moved_particle;
  double alpha, beta, a, a_1, a_2, b_1, b_2;

public:

  virtual inline void setParams(QickArray& params_){ 
    alpha          = params_(0);
    beta           = params_(1);
    a              = params_(2);
    moved_particle = int(params_(3));
    a_1            = params_(11);
    a_2            = params_(12);
    b_1            = params_(13);
    b_2            = params_(14);
  }

  virtual double valuePt(QickArray &pos);

  virtual inline double operator()(QickArray &pos){ 
    return valuePt(pos);
  }
};

//=================================================================
// Polls trial wave function for N bosons, returning all of wave
// function, not just one-particle part as usually done, only
// for testing
//=================================================================
class MateuszWaveAll : public virtual Wave{

protected:
  int moved_particle;
  double alpha, beta, a, a_1, a_2, b_1, b_2;

public:

  virtual inline void setParams(QickArray& params_){ 
    alpha          = params_(0);
    beta           = params_(1);
    a              = params_(2);
    moved_particle = int(params_(3));
    a_1            = params_(11);
    a_2            = params_(12);
    b_1            = params_(13);
    b_2            = params_(14);
  }

  virtual double valuePt(QickArray &pos);

  virtual inline double operator()(QickArray &pos){ 
    return valuePt(pos);
  }

};

//=================================================================
// analytic quantum force for duBois wave function
//=================================================================
class MateuszQForce : public virtual Func{

private:
  int moved_particle;
  double alpha, beta, a;
  double u(double r);
  double du(double r);
  QickArray pos_ij, dF, dG;

public:

  virtual inline void setParams(QickArray& params_){ 
    alpha          = params_(0);
    beta           = params_(1);
    a              = params_(2);
  }

  virtual QickArray& valuePt(QickArray &pos, bool dummy);

  virtual inline QickArray& operator()(QickArray &pos, bool dummy){ 
    return valuePt(pos, dummy);
  }
};

inline double MateuszQForce::u(double r){

  return log(1-a/r);
}
inline double MateuszQForce::du(double r){

  return a/(r*(r-a));
}

//=================================================================
// analytic local energy for duBois wave function and harm. osc.
// with distortion lambda in z-direction
//=================================================================
class MateuszLocalEnergy : public virtual Func{

private:
  double alpha, beta, a, V_a, V_r, mu_a, mu_r;
  double lambda; //omega_z/omega_rho
  Nabla nabla;
  Nabla2 nabla2;
  double u(double r);
  double du(double r);
  double ddu(double r);

public:

  virtual inline void setParams(QickArray& params_){ 
    alpha  = params_(0);
    beta   = params_(1);
    a      = params_(2);
    lambda = params_(4);
    V_a    = params_(7);
    V_r    = params_(8);
    mu_a   = params_(9);
    mu_r   = params_(10);
  }

  virtual double valuePt(QickArray &pos);

  virtual inline double operator()(QickArray &pos){ 
    return valuePt(pos);
  }
};

inline double MateuszLocalEnergy::u(double r){

  return log(1-a/r);
}
inline double MateuszLocalEnergy::du(double r){

  return a/(r*(r-a));
}

inline double MateuszLocalEnergy::ddu(double r){

  return a*(a-2*r)/sqr(r*(a-r));
}


//=================================================================
// analytic potential energy for bosons in harm. osc. with
// distortion in z-direction
//=================================================================
class MateuszPotentialEnergy : public virtual Func{

private:
  double lambda, V_a, V_r, mu_a, mu_r; //omega_z/omega_rho

public:

  virtual inline void setParams(QickArray& params_){ 
    lambda = params_(1);
    V_a    = params_(7);
    V_r    = params_(8);
    mu_a   = params_(9);
    mu_r   = params_(10);
  }

  virtual double valuePt(QickArray &pos);

  virtual inline double operator()(QickArray &pos){ 
    return valuePt(pos);
  }
};


#endif
