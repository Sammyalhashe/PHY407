#include "mc.hpp"

using namespace std;

void MCpp:: MetropolisStep(Walker& walker, Random& ran, Func& move,
			   Func& wf, QickArray& params, 
			   QickArray& counters, double step_length){

  int trials = int(counters(0));
  int accept = int(counters(1));

  QickArray move_params(2);
  move_params(0) = params(2);
  move_params(1) = step_length;

  int particles  = walker.getNoOfParticles();
  int dimensions = walker.getDimension();

  move.setParams(move_params);

  QickArray pos(dimensions);

  for(int i=0; i!=particles; i++){

    // update parameters
    params(3) = i;
    wf.setParams(params);

    double wf_old = walker.getWaveFunction(wf, true);

    pos=move(walker.getAllPositions(), i, ran);
      
    for(int j=0; j!=dimensions; j++){
      // trial move in steps of a
      walker.setParticlePosition(i,j,pos(j));
      
    }
    
    double wf_new = walker.getWaveFunction(wf, true);

    double w = exp(2*(wf_new-wf_old));

    trials++;
    // metropolis test
    if(ran.ran1() <= w){
      walker.updateParticlePosition(i);
      accept++;
    }
    else
      walker.resetParticlePosition(i);

  }

  counters(0) = double(trials);
  counters(1) = double(accept);

}

bool MCpp::diffuse(Walker& walker, Random& ran, Func& q_force, Func& wf_all,
		   Func& nabla, QickArray& params, QickArray& counters,
		   int numerical, double wf_x, double wf_y, double D,
		   double tau){


  double greens_xy, greens_yx;

  int trials        = int(counters(0));
  int accept        = int(counters(1));
  int accept_walker = int(counters(2));
  
  if(numerical==0)
    greens_xy = walker.getGreensFunction(q_force, D, tau);
  else
    greens_xy = walker.getGreensFunction(nabla, wf_all, D, tau);

  if(numerical==0)
    greens_yx = walker.getNewGreensFunction(q_force, D, tau);
  else
    greens_yx = walker.getNewGreensFunction(nabla, wf_all, D, tau);

  double w = exp(greens_yx-greens_xy+2*(wf_y-wf_x));

  trials++;
  // metropolis:
  if(w>ran.ran1()){
    walker.updateParticlePosition(int(params(3)));
    accept++;
    accept_walker++;
    counters(0) = trials;
    counters(1) = accept;
    counters(2) = accept_walker;
    return true; // do branching if accepted
  }
  else{
    walker.resetParticlePosition(int(params(3)));
    counters(0) = trials;
    return false; // skip branching if not accepted
  }

}

void MCpp::branch(Walker& walker, Random& ran, Func& local_e, Func& pot_e,
		  Func& wf, Func& wf_all, Func& nabla, Func& nabla2, 
		  int numerical, double tau, double e_local_x, double e_trial){

  double e_local_y, branching; 

  if(numerical==0)
    e_local_y = walker.getLocalEnergy(local_e);
  else
    e_local_y = walker.getLocalEnergy(wf_all, nabla, nabla2, pot_e);
  
  branching = -(.5*(e_local_x+e_local_y)-e_trial)*tau;
  branching = exp(branching);
  
  // random integer with average value equal to the branching
  int MB = int(branching);//+ran.ran1());
  if(branching-MB > ran.ran1())
    ++MB;

  // add MB-1 copies at end of list
  // if MB=0, mark this config as dead
  if(MB<20 && MB>0) for(int n=0; n!=MB-1; n++){
    walker.makeWalker();
  }
  else if(MB == 0){
    walker.killWalker();
  }

}

void MCpp::MonteCarloStep(Walker& walker, Random& ran, Func& wf, Func& wf_all,
			  Func& local_e, Func& q_force, Func& move, 
			  Func& pot_e, Func& nabla, Func& nabla2, 
			  QickArray& params, QickArray& counters,
			  int metropolis, double D, double tau, 
			  double e_trial, bool update, int numerical){

  double e_local_x, wf_x=1, wf_y=1;

  int particles  = walker.getNoOfParticles();
  int dimensions = walker.getDimension();

  QickArray move_params(2);
  move_params(0) = params(2);
  move_params(1) = sqrt(2.*D*tau);
  move.setParams(move_params);

  static QickArray fq_x(particles, dimensions);
  static QickArray fq_y(particles, dimensions);

  static QickArray pos(dimensions);
  static QickArray all_pos(dimensions);

  // initial info for this walker if update
  // update allways true first time
  if(numerical==0)
    e_local_x = walker.getLocalEnergy(local_e,update);
  else
    e_local_x = walker.getLocalEnergy(wf_all, nabla, nabla2, pot_e, 
				      update);

  if(metropolis==2)
    fq_x *= 0;
  else{
    if(numerical==0)
      fq_x = walker.getQuantumForce(q_force);
    else
      fq_x = walker.getQuantumForce(nabla, wf_all);
    
    fq_x *= D*tau;
  }

  for(int i=0; i!=particles; i++){
    
    params(3) = i;
    if(metropolis==1){
      wf.setParams(params);
      wf_x = walker.getWaveFunction(wf);
    }

    q_force.setParams(params);

    pos = move(walker.getAllPositions(), i, fq_x, ran);

    for(int j=0; j!=dimensions; j++){
      // trial move
      walker.setParticlePosition(i,j,pos(j));
    }

    if(metropolis!=1)
      walker.updateParticlePosition(i);
    else{
      wf_y = walker.getWaveFunction(wf);
      diffuse(walker, ran, q_force, wf_all, nabla, params, counters,
		     numerical, wf_x, wf_y, D, tau);
    }

  }

  if(metropolis==2){
    // drift
    move_params(1) = 0.;
    move.setParams(move_params);
    all_pos = walker.getAllPositions();
    if(numerical==0)
      fq_x = walker.getQuantumForce(q_force);
    else
      fq_x = walker.getQuantumForce(nabla, wf_all);
    
    fq_x *= 0.5*D*tau;
      
    for(int i=0; i!=particles; i++){
      
      pos = move(all_pos, i, fq_x, ran);
      
      for(int j=0; j!=dimensions; j++){
	// trial move
	walker.setParticlePosition(i,j,pos(j));
      }
      
      walker.updateParticlePosition(i);
    }
    
    // drift
    if(numerical==0)
      fq_y = walker.getNewQuantumForce(q_force);
    else
      fq_y = walker.getNewQuantumForce(nabla, wf_all);
    
    fq_y *= 0.5*D*tau;
    
    fq_y += fq_x;
    
    fq_y *= .5;
      
    for(int i=0; i!=particles; i++){
      
      pos = move(all_pos, i, fq_y, ran);
      
      for(int j=0; j!=dimensions; j++){
	// trial move
	walker.setParticlePosition(i,j,pos(j));
      }
      
      walker.updateParticlePosition(i);
    }
    
    branch(walker, ran, local_e, pot_e, wf, wf_all, nabla, nabla2, 
	   numerical, tau, e_local_x, e_trial);

    if(numerical==0)
      fq_x = walker.getNewQuantumForce(q_force);
    else
      fq_x = walker.getNewQuantumForce(nabla, wf_all);
    
    fq_x *= D*tau;

    for(int i=0; i!=particles; i++){
      
      pos = move(all_pos, i, fq_x, ran);
      
      for(int j=0; j!=dimensions; j++){
	walker.setParticlePosition(i,j,pos(j));
      }
      
      walker.updateParticlePosition(i);
      
    }

    move_params(1) = sqrt(2.*D*tau);
    move.setParams(move_params);
  }
  
  if(move.getImpossible())
    walker.killWalker();
  else if(e_trial!=0 && metropolis!=2)
    branch(walker, ran, local_e, pot_e, wf, wf_all, nabla, nabla2, 
	   numerical, tau, e_local_x, e_trial);

}
