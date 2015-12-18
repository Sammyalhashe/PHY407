#include "particle.hpp"

//==========================================================
// Particle constructor
//==========================================================
void Particle:: initialize(int dim_){
  dim=dim_;
  pos = new double[dim];
  new_pos = new double[dim];
}

//==========================================================
// Particle destructor
//==========================================================
Particle:: ~Particle(){
  delete[] pos;
  delete[] new_pos;
}
