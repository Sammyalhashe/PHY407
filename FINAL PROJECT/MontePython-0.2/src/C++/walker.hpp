#ifndef walker_hpp_IS_INCLUDED
#define walker_hpp_IS_INCLUDED

#include "QickArray.hpp"
#include "particle.hpp"
#include "functions.hpp"

//=================================================================
// class containing functions with regards to a walker
//=================================================================
class Walker{

private:

  double *particles_old;
  double *particles_new;
  double *qf;
  double *q_f;
  double *nabla_u;
  double *o_pure;
  double *o_pure_new;
  double *o_pure_further;
  double *o_pure_further_new;
  double *e_l_ptr;
  double *wave_func_ptr;
  double *no_of_particles_ptr;
  double *dim_ptr;
  double *particles_moved_ptr;
  double *size_ptr;
  double *replicate_ptr;
  double *dead_ptr;
  double *old_eq_new_ptr;

  QickArray& posVec();
  QickArray& newPosVec();

  QickArray tmp_2d;
  QickArray tmp2_2d;
  QickArray pos;

  inline double sqr(double x){ return x*x; }

  int packWalker();

  void unpackWalker();

public:

  Walker(){
    //std::cerr << "I'm alive!" << std::endl;
  }

  ~Walker(){
    //std::cerr << "I'm a dead man walkin'!" << std::endl;
  }

  //must be called right after constructor:
  void initialize(int no_of_particles_, int dim_, 
		  double *buffer, int buffer_size_);

  double getParticlePosition(int i, int j); // returns j'th coord position 
                                    // of i'th particle
  double getNewParticlePosition(int i, int j); // returns j'th coord position 
                                    // of i'th particle's new pos

  double getOPureFurther(int i, int j);

  void setParticlePosition(int i, int j, double x); // set's value of 
                                                    //j'th coord,
                                                    // i'th particle
  void updateParticlePosition(int i); // pos=new_pos
  void resetParticlePosition(int i); // new_pos=pos

  QickArray& getAllPositions(){ return newPosVec(); }

  double getLocalEnergy(Func& wf, Func &nabla, Func &nabla2, 
			Func &potential, bool update=true);

  // function to use if analytic E_L is known
  double getLocalEnergy(Func& e_local, bool update=true);

  // function to use if analytic E_L is known
  double getOtherEnergy(Func& e_other);

  double getCoM();

  double getWaveFunction(Func &wf, bool update=true);

  QickArray& getQuantumForce(Func &nabla, Func &wf, 
			     bool update=true);
  QickArray& getNewQuantumForce(Func &nabla, Func &wf, 
			     bool update=true);
  
  // function to use if analytic F_q is known
  QickArray& getQuantumForce(Func &q_force, bool update=true);
  // function to use if analytic F_q is known
  QickArray& getNewQuantumForce(Func &q_force, bool update=true);
  
  double getGreensFunction(Func &nabla, Func &wf, double D, double tau);
  
  double getNewGreensFunction(Func &nabla, Func &wf, double D, double tau);
  
  // function to use if analytic F_q is known
  double getGreensFunction(Func &q_force, double D, double tau);

  // function to use if analytic F_q is known
  double getNewGreensFunction(Func &q_force, double D, double tau);
  
  inline void killWalker(){ 
    dead_ptr[0]=true; 
  }

  inline void makeWalker(){ 
    replicate_ptr[0]++; 
  }

  inline void madeWalker(){
    replicate_ptr[0]--; 
  }

  inline void calmWalker(){ 
    replicate_ptr[0] = 0; 
  }

  inline bool isDead(){ return *dead_ptr; }

  inline bool tooAlive(){ 
    double& replicate = *replicate_ptr;
    if(replicate > 0)
      return true;
    else
      return false;
  }

  inline int getNoOfParticles(){ return int(*no_of_particles_ptr); }

  inline int getDimension(){ return int(*dim_ptr); }


};

//=================================================================
// function to get position of j-th dimension of particle 
// number i
// this function returns the saved position
//=================================================================
inline double Walker:: getParticlePosition(int i, int j){
  double& no_of_particles = *no_of_particles_ptr;
  return particles_old[i+j*int(no_of_particles)];
}

//=================================================================
// function to get position of j-th dimension of particle 
// number i
// this function returns the trial position
//=================================================================
inline double Walker:: getNewParticlePosition(int i, int j){
  double& no_of_particles = *no_of_particles_ptr;
  return particles_new[i+j*int(no_of_particles)];
}

inline double Walker:: getOPureFurther(int i, int j){
  double& no_of_particles = *no_of_particles_ptr;
  return o_pure_further[i+j*int(no_of_particles)];
}

//=================================================================
// this function is used to set the trial position x in 
// coord (i,j)
//=================================================================
inline void Walker:: setParticlePosition(int i, int j, double x){
  double& no_of_particles = *no_of_particles_ptr;
  particles_new[i+j*int(no_of_particles)] = x;
  old_eq_new_ptr[0] = false;
}

//=================================================================
// function to save accepted position for particle i
//=================================================================
inline void Walker:: updateParticlePosition(int i){
  double& no_of_particles = *no_of_particles_ptr;
  double& dim = *dim_ptr;
  for(int j = 0; j!=dim; j++)
    particles_old[i+j*int(no_of_particles)] = 
      particles_new[i+j*int(no_of_particles)];
  old_eq_new_ptr[0] = true;
}

//=================================================================
// function to reset to saved position when move not accepted
//=================================================================
inline void Walker:: resetParticlePosition(int i){
  double& no_of_particles = *no_of_particles_ptr;
  double& dim = *dim_ptr;
  for(int j = 0; j!=dim; j++)
    particles_new[i+j*int(no_of_particles)] = 
      particles_old[i+j*int(no_of_particles)];
  old_eq_new_ptr[0] = true;
}

 //=================================================================
// function returning numerical local energy, computing new
// energy if update is true, else returning stored energy for
// this walker
// wf must be (given wf_T=exp(u(R))) on the form u(R)
//=================================================================
inline double Walker:: getLocalEnergy(Func &wf, Func &nabla, Func &nabla2,
				      Func &potential, bool update){
  double& e_l = *e_l_ptr;
  double& no_of_particles = *no_of_particles_ptr;
  double& dim = *dim_ptr;

  if(update){
    e_l = -nabla2(wf,newPosVec());
    tmp_2d = nabla(newPosVec(),wf);
    for(int i=0; i!=no_of_particles; i++)
      for(int x=0; x!=dim; x++){
	e_l -= sqr(tmp_2d(i,x));
      }
    e_l += potential(newPosVec());
  }

  e_l_ptr[0] = e_l;
  return e_l;

} // end getLocalEnergy


//=================================================================
// function returning analytical local energy, computing new
// energy if update is true
//=================================================================
inline double Walker:: getLocalEnergy(Func &e_local, bool update){

  double& e_l = *e_l_ptr;

  if(update)
    e_l = e_local(newPosVec());
  e_l_ptr[0] = e_l;
  return e_l;
} // end getLocalEnergy

//=================================================================
// function returning analytical other energy, computing new
// energy if update is true
//=================================================================
inline double Walker:: getOtherEnergy(Func &e_other){

  return e_other(newPosVec());
} // end getLocalEnergy



//=================================================================
// function returning wave function, recalculating if update
//=================================================================
inline double Walker:: getWaveFunction(Func &wf, bool update){

  double& wave_func = *wave_func_ptr;

  if(update)
    wave_func = wf(newPosVec());
  wave_func_ptr[0] = wave_func;
  return wave_func;

} // end getWaveFunction


//=================================================================
// numerical q-force, recalculated if update
//=================================================================
inline QickArray& Walker:: getQuantumForce(Func &nabla, Func &wf, 
					   bool update){

  double& no_of_particles = *no_of_particles_ptr;
  double& dim = *dim_ptr;

  if(update){
    tmp_2d = nabla(posVec(), wf);
    tmp_2d *= 2;
    for(int i=0; i!=no_of_particles; i++)
      for(int j=0; j!=dim; j++)
	qf[i+j*int(no_of_particles)] = tmp_2d(i,j);
  }
  else
    for(int i=0; i!=no_of_particles; i++)
      for(int j=0; j!=dim; j++)
	tmp_2d(i,j) = qf[i+j*int(no_of_particles)];

  return tmp_2d;
} // end getQuantumForce

//=================================================================
// same as getQuantumForce()
//=================================================================
inline QickArray& Walker:: getNewQuantumForce(Func &nabla, Func &wf, 
					      bool update){

  double& no_of_particles = *no_of_particles_ptr;
  double& dim = *dim_ptr;

  if(update){
    tmp_2d = nabla(newPosVec(), wf);
    tmp_2d *= 2;
    for(int i=0; i!=no_of_particles; i++)
      for(int j=0; j!=dim; j++)
	qf[i+j*int(no_of_particles)] = tmp_2d(i,j);
  }
  else
    for(int i=0; i!=no_of_particles; i++)
      for(int j=0; j!=dim; j++)
	tmp_2d(i,j) = qf[i+j*int(no_of_particles)];

  return tmp_2d;

} // end getQuantumForce


//=================================================================
// analytical q-force, for saved position
//=================================================================
inline QickArray& Walker:: getQuantumForce(Func &q_force, bool update){

  double& no_of_particles = *no_of_particles_ptr;
  double& dim = *dim_ptr;

  if(update){
    tmp_2d = q_force(posVec(), false);
    for(int i=0; i!=no_of_particles; i++)
      for(int j=0; j!=dim; j++)
	qf[i+j*int(no_of_particles)] = tmp_2d(i,j);
  }
  else
    for(int i=0; i!=no_of_particles; i++)
      for(int j=0; j!=dim; j++)
	tmp_2d(i,j) = qf[i+j*int(no_of_particles)];

  return tmp_2d;

} // end getQuantumForce

//=================================================================
// analytical q-force, for trial position
//=================================================================
inline QickArray& Walker:: getNewQuantumForce(Func &q_force, bool update){

  double& no_of_particles = *no_of_particles_ptr;
  double& dim = *dim_ptr;

  if(update){
    tmp_2d = q_force(newPosVec(), false);
    for(int i=0; i!=no_of_particles; i++)
      for(int j=0; j!=dim; j++)
	qf[i+j*int(no_of_particles)] = tmp_2d(i,j);
  }
  else
    for(int i=0; i!=no_of_particles; i++)
      for(int j=0; j!=dim; j++)
	tmp_2d(i,j) = qf[i+j*int(no_of_particles)];

  return tmp_2d;

} // end getQuantumForce

//=================================================================
// analytical greens function, probality of moving from old to
// new position
//=================================================================
inline double Walker:: getGreensFunction(Func &q_force, 
					 double D, double tau){

  double& no_of_particles = *no_of_particles_ptr;
  double& dim = *dim_ptr;

  tmp2_2d = getQuantumForce(q_force, false);

  double G = 0;
  for(int i=0; i!=no_of_particles; i++)
    for(int j=0; j!=dim; j++)
      G += sqr(getNewParticlePosition(i,j) - 
	       getParticlePosition(i,j) - D*tau*tmp2_2d(i,j));
  
  return (-G/4/D/tau);
} // end getQuantumForce

//=================================================================
// analytical greens function, probality of moving from new to
// old position
//=================================================================
inline double Walker:: getNewGreensFunction(Func &q_force, 
					    double D, double tau){

  double& no_of_particles = *no_of_particles_ptr;
  double& dim = *dim_ptr;

  tmp2_2d = getNewQuantumForce(q_force, true);

  double G = 0;
  for(int i=0; i!=no_of_particles; i++)
    for(int j=0; j!=dim; j++)
      G += sqr(getParticlePosition(i,j) - 
	       getNewParticlePosition(i,j) - D*tau*tmp2_2d(i,j));

  return (-G/4/D/tau);
} // end getQuantumForce


//=================================================================
// analytical greens function, probality of moving from old to
// new position
//=================================================================
inline double Walker:: getGreensFunction(Func &nabla, Func &wf, 
					 double D, double tau){

  double& no_of_particles = *no_of_particles_ptr;
  double& dim = *dim_ptr;
  tmp2_2d = getQuantumForce(nabla, wf, false);

  double G = 0;
  for(int i=0; i!=no_of_particles; i++)
    for(int j=0; j!=dim; j++)
      G += sqr(getNewParticlePosition(i,j) - 
	       getParticlePosition(i,j) - D*tau*tmp2_2d(i,j));

  return (-G/4/D/tau);
} // end getQuantumForce

//=================================================================
// analytical greens function, probality of moving from new to
// old position
//=================================================================
inline double Walker:: getNewGreensFunction(Func &nabla, Func &wf, 
					    double D, double tau){

  double& no_of_particles = *no_of_particles_ptr;
  double& dim = *dim_ptr;

  tmp2_2d = getNewQuantumForce(nabla, wf, true);

  double G = 0;
  for(int i=0; i!=no_of_particles; i++)
    for(int j=0; j!=dim; j++)
      G += sqr(getParticlePosition(i,j) - 
	       getNewParticlePosition(i,j) - D*tau*tmp2_2d(i,j));

  return (-G/4/D/tau);
} // end getQuantumForce

//=================================================================
// function returning Center-of-Mass position of this walker
//=================================================================
inline double Walker:: getCoM(){

  double& no_of_particles = *no_of_particles_ptr;
  double& dim = *dim_ptr;

  double com=0;
  for(int i=0; i!=no_of_particles; i++)
    for(int j=0; j!=dim; j++)
      com += sqr(getParticlePosition(i,j));
  return com;
}

//=================================================================
// function returning a QickArray of saved positions 
// of all particles
//=================================================================
inline QickArray& Walker:: posVec(){
  double& no_of_particles = *no_of_particles_ptr;
  double& dim = *dim_ptr;
  for(int i=0; i!=no_of_particles; i++)
    for(int j=0; j!=dim; j++)
      pos(i,j) = getParticlePosition(i,j);
  return pos;
}

//=================================================================
// function returning a QickArray of trial positions 
// of all particles
//=================================================================
inline QickArray& Walker:: newPosVec(){
  double& no_of_particles = *no_of_particles_ptr;
  double& dim = *dim_ptr;
  for(int i=0; i!=no_of_particles; i++)
    for(int j=0; j!=dim; j++)
      pos(i,j) = getNewParticlePosition(i,j);
  return pos;
}

#endif
