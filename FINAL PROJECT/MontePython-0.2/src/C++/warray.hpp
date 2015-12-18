#ifndef warray_hpp_IS_INCLUDED
#define warray_hpp_IS_INCLUDED


#include <iostream>

#include "QickArray.hpp"
#include "random.hpp"
#include "walker.hpp"
#include "mc.hpp"

class Warray{


protected:
  
  Walker *walkers;

  MCpp *mc;

  int no_of_walkers, particles, dim;

  int walker_array_length;


public:

  Warray(){ 
    walkers = new Walker[0];
    mc = new MCpp();
  };

  ~Warray(){ 
    delete[] walkers;
    delete mc;
  }

  void initialize(int no_of_walkers_, int particles_, int dim_, 
		  double* a, int a_size);

  void UniDist(Random& ran);

  double getOFurther(int i,int j, int k){ return walkers[i].getOPureFurther(j,k); }

  void MetropolisSteps(Random& ran, Func& move,
		      Func& wf, QickArray& params, 
		      QickArray& counters, double step_length);
    
  void MonteCarloSteps(Random& ran, Func& wf, Func& wf_all,
		       Func& local_e, Func& q_force, Func& move, 
		       Func& pot_e, Func& nabla, Func& nabla2, 
		       QickArray& params, QickArray& counters,
		       int metropolis, double D, double tau, 
		       double e_trial, bool update, int numerical);

  double getLocalEnergies(Func& wf, Func &nabla, Func &nabla2, 
			Func &potential, bool update=true);

  // function to use if analytic E_L is known
  double getLocalEnergies(Func& e_local, bool update=true);

  // function to use to get potential and vortex energies
  double getOtherEnergies(Func& e_other);

  QickArray* walkers2grid3D(QickArray& grid, double dx, double dy, double dz); 

  QickArray* walkers2grid3DPure(QickArray& grid, 
				double dx, double dy, double dz,
				double weight); 
};

#endif
