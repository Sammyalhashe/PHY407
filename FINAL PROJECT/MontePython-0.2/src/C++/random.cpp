#include "mpi.h"
#include "random.hpp"

Random::Random(){
  
  int size,rank;
  
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  
  for(int i=0; i!=rank+1; i++)
    s.ResetNextSubstream();
  
}

Random::Random(long idum_){
  idum=-(long)(abs(idum_));
  unsigned long seed[6];
  for(int i=0; i!=6; i++)
    seed[i]=idum;

  int size,rank;
  
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Comm_size(MPI_COMM_WORLD,&size);

  RngStream::SetPackageSeed(seed);

  for(int i=0; i!=rank+1; i++)
    s.ResetNextSubstream();

}

long Random::get_idum(){ 
  unsigned long seed[6];
  s.GetState(seed);
  return seed[0];
}


// normal random variate generator using box-muller
double Random::gran(double s, double m)
{        /* mean m, standard deviation s */
  double x1, x2, w, y1;
  static double y2;
  static int use_last = 0;

  if (use_last)        /* use value from previous call */
    {
      y1 = y2;
      use_last = 0;
    }
  else
    {
      do {
	x1 = 2.0 * ran1() - 1.0;
	x2 = 2.0 * ran1() - 1.0;
	w = x1 * x1 + x2 * x2;
      } while ( w >= 1.0 );
      w = sqrt( (-2.0 * log( w ) ) / w );
      y1 = x1 * w;
      y2 = x2 * w;
      use_last = 1;
    }
  return( m + y1 * s );
}


double Random::ran1(){return s.RandU01();}


