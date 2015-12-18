#ifndef MC_HPP
#define MC_HPP
#include <cmath>
#include "QickArray.hpp"
#include "random.hpp"
#include "walker.hpp"
#include "functions.hpp"

class MCpp{

public:
  void MetropolisStep(Walker& walker, Random& ran, Func& move,
		      Func& wf, QickArray& params, 
		      QickArray& counters, double step_length);

  bool diffuse(Walker& walker, Random& ran, Func& q_force, Func& wf_all,
	       Func& nabla, QickArray& params, QickArray& counters,
	       int numerical, double wf_x, double wf_y, double D,
	       double tau);


  void branch(Walker& walker, Random& ran, Func& local_e, Func& pot_e,
	      Func& wf, Func& wf_all, Func& nabla, Func& nabla2, 
	      int numerical, double tau, double e_local_x, double e_trial);

  void MonteCarloStep(Walker& walker, Random& ran, Func& wf, Func& wf_all,
		      Func& local_e, Func& q_force, Func& move, 
		      Func& pot_e, Func& nabla, Func& nabla2, 
		      QickArray& params, QickArray& counters,
		      int metropolis, double D, double tau, 
		      double e_trial, bool update, int numerical);

};

#endif
