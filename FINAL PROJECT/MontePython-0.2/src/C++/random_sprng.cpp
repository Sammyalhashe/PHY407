#include "random.hpp"

#define SIMPLE_SPRNG            /* simple interface                        */
#define USE_MPI                 /* MPI version of SPRNG                    */
#include <sprng.h>

Random::Random(){
  idum=-1;
  gtype=2;
  init_sprng(gtype,idum+1,SPRNG_DEFAULT);
}
Random::Random(long idum_){
  idum=-(long)(abs(idum_));
  gtype=2;
  init_sprng(gtype,idum,SPRNG_DEFAULT);
}

long Random::get_idum(){ return isprng(); }


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

double Random::ran1(){return sprng();}

