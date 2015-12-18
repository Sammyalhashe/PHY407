#include "warray.hpp"
#include "walker.hpp"

using namespace std;

void Warray::initialize(int no_of_walkers_, int particles_, int dim_,
			double* a, int a_size){
  
  no_of_walkers = no_of_walkers_;
  walker_array_length = no_of_walkers;
  particles = particles_;
  dim = dim_;

  delete mc;
  mc = new MCpp();

  delete[] walkers;
  walkers = new Walker[no_of_walkers];

  int stride = particles*dim*9+9;
//   if(a_size != stride*no_of_walkers)
//     cerr << "Inconsistent size of walker_buffer\n" 
// 	 << "Should be " << stride*no_of_walkers 
// 	 << " but is " << a_size << ". This is bad!" << endl;


  double* walker_ptr = a;
  for(int i=0; i!=no_of_walkers; i++){
    walkers[i].initialize(particles,dim,walker_ptr,
			  stride);
    walker_ptr += stride;
  }

}


void Warray::UniDist(Random& ran){

  for(int k=0; k!=no_of_walkers; k++)
    for(int i=0; i!=particles; i++){
      for(int j=0; j!=dim; j++)
	walkers[k].setParticlePosition(i,j,ran.ran1()-.5);
      walkers[k].updateParticlePosition(i);
    }
  
}

void Warray::MetropolisSteps(Random& ran, Func& move,
			     Func& wf, QickArray& params, 
			     QickArray& counters, double step_length){

  for(int i=0; i!=no_of_walkers; i++){
    mc->MetropolisStep(walkers[i], ran, move, wf, params, 
		       counters, step_length);
  }

}
    
void Warray::MonteCarloSteps(Random& ran, Func& wf, Func& wf_all,
			     Func& local_e, Func& q_force, Func& move, 
			     Func& pot_e, Func& nabla, Func& nabla2, 
			     QickArray& params, QickArray& counters,
			     int metropolis, double D, double tau, 
			     double e_trial, bool update, int numerical){

  for(int i=0; i!=no_of_walkers; i++)
    mc->MonteCarloStep(walkers[i], ran, wf, wf_all, local_e, q_force,
		       move, pot_e, nabla, nabla2, params,
		       counters, metropolis, D, tau, e_trial, 
		       update, numerical);

}

double Warray::getLocalEnergies(Func& wf, Func &nabla, Func &nabla2, 
				Func &potential, bool update){
  double res=0;
  for(int i=0; i!=no_of_walkers; i++)
    res += walkers[i].getLocalEnergy(wf, nabla, nabla2,
				     potential, update);
  return res;

}

// function to use if analytic E_L is known
double Warray::getLocalEnergies(Func& e_local, bool update){

  double res=0;
  for(int i=0; i!=no_of_walkers; i++)
    res += walkers[i].getLocalEnergy(e_local, update);
  return res;

}

// function to use to get potential and vortex energies
double Warray::getOtherEnergies(Func& e_other){

  double res=0;
  for(int i=0; i!=no_of_walkers; i++)
    res += walkers[i].getOtherEnergy(e_other);
  return res;

}

QickArray* Warray::walkers2grid3D(QickArray& grid, double dx, double dy, 
				  double dz){

  // Function putting walkers in a 3d histogram, i.e. wf distribution

  int lx, ly, lz, dummy, xpos=0, ypos=0, zpos=0;
  dummy = grid.get_dimension_info(&lx,&ly,&lz);

  QickArray pos(particles,dim);

  double tmp_pos;

  for(int k=0; k!=no_of_walkers; k++){
    pos *= 0;
    for(int i=0; i!=particles; i++)
      for(int j=0; j!=dim; j++){
	tmp_pos = walkers[k].getParticlePosition(i,j);
	if(tmp_pos!=0.)
	  pos(i,j) += tmp_pos;
      }

    // letting center of the histogram be in the center of the grid
    // allways rounding down
    for(int i=0; i!=particles; i++){
      if(pos(i,0)>=0)
        xpos = int(pos(i,0)/dx) + lx/2;
	if(xpos>=lx) xpos=lx-1;
      else{
        xpos = int(pos(i,0)/dx) + lx/2 - 1;
	if(xpos<0) xpos=0;
      }
      if(pos(i,1)>=0)
        ypos = int(pos(i,1)/dy) + ly/2;
	if(ypos>=ly) ypos=ly-1;
      else{
        ypos = int(pos(i,1)/dy) + ly/2 - 1;
	if(ypos<0) ypos=0;
      }
      if(pos(i,2)>=0)
        zpos = int(pos(i,2)/dz) + lz/2;
	if(zpos>=lz) zpos=lz-1;
      else{
        zpos = int(pos(i,2)/dz) + lz/2 - 1;
	if(zpos<0) zpos=0;
      }
      grid(xpos,ypos,zpos)++;
    }
  }
  return &grid;
}

QickArray* Warray::walkers2grid3DPure(QickArray& grid, double dx, double dy, 
				      double dz, double weight){

  // Function putting walkers in a 3d histogram, i.e. wf distribution

  int lx, ly, lz, dummy, xpos=0, ypos=0, zpos=0;
  dummy = grid.get_dimension_info(&lx,&ly,&lz);

  QickArray pos(particles,dim);

  double tmp_pos;

  for(int k=0; k!=no_of_walkers; k++){
    pos *= 0;
    for(int i=0; i!=particles; i++)
      for(int j=0; j!=dim; j++){
	tmp_pos = walkers[k].getOPureFurther(i,j)/weight;
	if(tmp_pos!=0.)
	  pos(i,j) += tmp_pos;
      }

    // letting center of the histogram be in the center of the grid
    // allways rounding down
    for(int i=0; i!=particles; i++){
      if(pos(i,0)>=0)
        xpos = int(pos(i,0)/dx) + lx/2;
      if(xpos>=lx) xpos=lx-1;
      else{
        xpos = int(pos(i,0)/dx) + lx/2 - 1;
	if(xpos<0) xpos=0;
      }
      if(pos(i,1)>=0)
        ypos = int(pos(i,1)/dy) + ly/2;
	if(ypos>=ly) ypos=ly-1;
      else{
        ypos = int(pos(i,1)/dy) + ly/2 - 1;
	if(ypos<0) ypos=0;
      }
      if(pos(i,2)>=0)
        zpos = int(pos(i,2)/dz) + lz/2;
	if(zpos>=lz) zpos=lz-1;
      else{
        zpos = int(pos(i,2)/dz) + lz/2 - 1;
	if(zpos<0) zpos=0;
      }
      grid(xpos,ypos,zpos)+=1.;
    }
  }
  return &grid;
}
