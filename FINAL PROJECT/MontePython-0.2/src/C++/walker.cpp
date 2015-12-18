#include "walker.hpp"

using namespace std;

//=================================================================
// initialization of walker class taking as input an int array
// containing number of particles and dimension, making an
// array of particles according to these
//=================================================================
void Walker:: initialize(int no_of_particles_, int dim_, 
			 double *buffer, int buffer_size_){

  int array_size = no_of_particles_*dim_;

  particles_old = (double*) buffer;
  buffer += array_size;
  particles_new = (double*) buffer;
  buffer += array_size;
  qf = (double*) buffer;
  buffer += array_size;
  q_f = (double*) buffer;
  buffer += array_size;
  nabla_u = (double*) buffer;
  buffer += array_size;

  o_pure  = (double*) buffer;
  buffer += array_size;
  o_pure_new = (double*) buffer;
  buffer += array_size;
  o_pure_further  = (double*) buffer;
  buffer += array_size;
  o_pure_further_new  = (double*) buffer;
  buffer += array_size;

  e_l_ptr = (double*) buffer;
  buffer++;
  wave_func_ptr = (double*) buffer;
  buffer++;
  no_of_particles_ptr = (double*) buffer;
  buffer++;
  dim_ptr = (double*) buffer;
  buffer++;
  particles_moved_ptr = (double*) buffer;
  buffer++;
  replicate_ptr = (double*) buffer;
  buffer++;
  size_ptr = (double*) buffer;
  buffer++;
  dead_ptr = (double*) buffer;
  buffer++;
  old_eq_new_ptr = (double*) buffer;
  
  tmp_2d.redim(int(*no_of_particles_ptr),int(*dim_ptr));
  tmp2_2d.redim(int(*no_of_particles_ptr),int(*dim_ptr));
  pos.redim(int(*no_of_particles_ptr),int(*dim_ptr));

}
