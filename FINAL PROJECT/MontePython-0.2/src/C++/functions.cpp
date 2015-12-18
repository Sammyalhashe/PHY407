#include "functions.hpp"

double Helium::valuePt(QickArray &pos){

  int cols, rows, dim;
  dim = pos.get_dimension_info(&cols,&rows);
  if(dim!=2 || cols!=2)
    std::cerr << "pos must be a 2d 2 particles x pos array" << std::endl;

  double ans=0;
  for(int i=0; i!=cols; i++){
    for(int j=0; j!=rows; j++)
      dummy(j)=pos(i,j);
    ans+=wave1s(dummy);
  }
  
  for(int j=0; j!=rows; j++)
    dummy(j)=pos(1,j)-pos(0,j);

  return (ans+corr(dummy));

}


//=================================================================
//Movement of particle ensuring hard spere potential; 
//step_length*(+/-.5)
//=================================================================
QickArray& DuBoisMove::valuePt(QickArray& pos, int particle, 
			   Random& ran){

  int cols, rows, dummy;
  dummy = pos.get_dimension_info(&cols,&rows);

  if(!pos_.is_initialized())
    pos_.redim(rows);

  double dist=0;

  bool too_close = true;
  //hard sphere precautions
  while(too_close){
    for(int j=0; j!=rows; j++)
      pos_(j) = pos(particle,j)+step_length*(ran.ran1()-.5);
    //was the move allowed?
    too_close = false;
    for(int i=0; i!=cols; i++) if(i!=particle){
      dist = 0;
      for(int j=0; j!=rows; j++)
	dist += sqr(pos_(j)-pos(i,j));
      if(dist<=sqr(a)){
	too_close = true;
	break; // no point in checking the rest of the particles
      }
    }
  }

  return pos_;
}

//=================================================================
//Movement of particle ensuring hard spere potential; 
//step_length*(+/-.5)
//This overload gives posibility to move with a constant
//in adition to the random movement and has gaussian movement
//=================================================================
QickArray& DuBoisMove::valuePt(QickArray& pos, int particle, 
			       QickArray& x, Random& ran){

  int cols, rows, dummy;
  dummy = pos.get_dimension_info(&cols,&rows);

  if(!pos_.is_initialized()){
    pos_.redim(rows);
  }
  if(!prepos_.is_initialized()){
    prepos_.redim(rows);
  }

  double dist=0;
  double premove_dist=0;

  bool too_close = true;

  impossible = false;

  double count = 0.;
  //hard sphere precautions
  while(too_close){
    count+=1.;
    for(int j=0; j!=rows; j++){
      if(step_length==0)
	pos_(j) = pos(particle,j)+x(particle,j);
      else
	pos_(j) = pos(particle,j)+x(particle,j)+ran.gran(1,0)*step_length;
      prepos_(j) = pos(particle,j)+x(particle,j);
    }
    //was the move allowed?
    too_close = false;
    for(int i=0; i!=cols; i++) if(i!=particle){
      dist = 0;
      premove_dist = 0;
      for(int j=0; j!=rows; j++){
	dist += sqr(pos_(j)-pos(i,j));
	premove_dist += sqr(prepos_(j)-pos(i,j));
      }
      if(dist<=sqr(a)){
	// check if particles will "ever" get seperated:
	if((sqrt(premove_dist)+5.*step_length)>a) 
	  too_close = true;
	else
	  impossible = true;
	break; // no point in checking the rest of the particles
      }
    }
  }
  return pos_;
}

//==================================================================
//Finds K given d to fullfill healing conditions for Polls Jastrow
//Function using Newton-Raphson
//==================================================================
double PollsKd::valuePt(double d, double a){

  double K=0;

  double K_plus=1e-16;

  double conv=1e-25;

  while(fabs(K-K_plus)>conv){
    K=K_plus;
    K_plus = K-F(K,d,a)/J(K,d,a);
  }

  return K;

}

//=================================================================
// u(alpha,lambda,r_i) + v(a,r_ij)
//=================================================================
double PollsVortexWave:: valuePt(QickArray &pos){

  int cols, rows, dim;
  dim = pos.get_dimension_info(&cols,&rows);

  if(!tmp_pos.is_initialized())
    tmp_pos.redim(rows);
  
  double g=0, rho_i=0;
  // g(alpha,lambda,r_ij)
  for(int j=0; j!=rows; j++)
    if(j==2)
      g+=lambda*sqr(pos(moved_particle,j));
    else{
      rho_i+=sqr(pos(moved_particle,j));
    }

  g += rho_i;
  g  = log(sqrt(rho_i))-alpha*g;

  double f = 0;
  double r_ij;
  // f(a,r_ij)
  for(int i=0; i!=moved_particle; i++){
    r_ij = 0;
    for(int j=0; j!=rows; j++)
      r_ij += sqr(pos(moved_particle,j)-pos(i,j));
    r_ij = sqrt(r_ij);
    if(r_ij > a)
      if(r_ij < d)
        f += log(d/r_ij*sin(K*(r_ij-a))/sin(K*(d-a)));
      else
	f += 0.;
    else{
      f = log(0.);
      break;
    }
  }
  for(int i=moved_particle+1; i!=cols; i++){
    r_ij = 0;
    for(int j=0; j!=rows; j++)
      r_ij += sqr(pos(moved_particle,j)-pos(i,j));
    r_ij = sqrt(r_ij);
    if(r_ij > a)
      if(r_ij < d)
        f += log(d/r_ij*sin(K*(r_ij-a))/sin(K*(d-a)));
      else
	f += 0.;
    else{
      f = log(0.);
      break;
    }
  }

  return g+f;

}

//=================================================================
// g(alpha,lambda,r_i)f(a,r_ij)
//=================================================================
double PollsVortexWaveAll:: valuePt(QickArray &pos){

  int cols, rows, dim;
  dim = pos.get_dimension_info(&cols,&rows);

  if(!tmp_pos.is_initialized())
    tmp_pos.redim(rows);

  double f = 0;

  // f(a,r_ij)
  for(int moved_particle=0; moved_particle!=cols-1; moved_particle++)
    for(int i=moved_particle+1; i!=cols; i++){
      for(int j=0; j!=rows; j++)
	tmp_pos(j) = pos(moved_particle,j)-pos(i,j);
      double r_ij = 0;
      for(int j=0; j!=rows; j++){
	r_ij += sqr(tmp_pos(j));
      }
      r_ij = sqrt(r_ij);
      if(r_ij > a)
        if(r_ij < sqr(d))
	  f += log(d/r_ij*sin(K*(r_ij-a))/sin(K*(d-a)));
	else
	  f += 0;
      else{
        f = log(0.);
	break;
      }
    }

  // g(alpha,lambda,r_ij)
  double g=0, G=0, rho_i=0;
  // g(alpha,lambda,r_ij)
  for(int moved_particle=0; moved_particle!=cols; moved_particle++){
    g = rho_i = 0;
    for(int j=0; j!=rows; j++)
      if(j==2)
	g+=lambda*sqr(pos(moved_particle,j));
      else{
	rho_i+=sqr(pos(moved_particle,j));
	
      }	
    g += rho_i;
    G += log(sqrt(rho_i))-alpha*g;
  }

  return f+G;

}

//=================================================================
// 2*nabla psi / psi !?!?!?!?!?!?!?!?!?!?!?!?!?!?!?!?!?!?!?!?!?!?!?!?
//=================================================================
QickArray& PollsVortexQForce:: valuePt(QickArray &pos, bool dummy){

  int cols, rows, dim;
  dim = pos.get_dimension_info(&cols,&rows);

  if(!pos_ij.is_initialized())
    pos_ij.redim(rows);

  if(!dF.is_initialized()){
    dF.redim(cols,rows);
    dG.redim(cols,rows);
  }

  dF *= 0;

  double r_ij, xy_fac, rho_i;

  for(int i=0; i!=cols; i++){

    rho_i=0;
    for(int x=0; x!=2; x++)
      rho_i += sqr(pos(i,x));

    xy_fac = 1-1/(2*alpha*rho_i);

    for(int x=0; x!=rows; x++)
      if(x==2)
	dG(i,x) = lambda*pos(i,x);
      else
	dG(i,x) = xy_fac*pos(i,x);

    for(int j=0; j!=cols; j++) if(i!=j){

      r_ij=0;
      for(int x=0; x!=rows; x++){
	pos_ij(x) = pos(i,x)-pos(j,x);
	r_ij += sqr(pos_ij(x));
      }
      r_ij = sqrt(r_ij);
      
      pos_ij *= du(r_ij)/r_ij;

      for(int x=0; x!=rows; x++)
	dF(i,x) += pos_ij(x);

    }

  }

  dG *= -2*alpha;

  dF += dG;

  dF *= 2;

  return dF;

}

double PollsVortexLocalEnergy:: valuePt(QickArray &pos){

  int cols, rows, dim;
  dim = pos.get_dimension_info(&cols,&rows);

  if(!pos_ij.is_initialized()){
    pos_ij.redim(rows);
    pos_i.redim(rows);
    pos_ik.redim(rows);
  }

  double kinetic_e = 2*alpha*cols*(4+lambda);

  double e_lg2=0, e_lg3=0, e_lf1=0, e_lf2=0, e_lfg1=0, e_lfg2=0;
  double r_ij, r_ik, rho_i, du_ij;

  for(int i=0; i!=cols; i++){

    rho_i=0;
    for(int x=0; x!=rows; x++)
      if(x==2){
	pos_i(x) = lambda*pos(i,x);
	e_lg2   += sqr(pos_i(x));
      }
      else{
	pos_i(x) = pos(i,x);
	rho_i   += sqr(pos_i(x));
      }

    e_lg2 += rho_i;

    e_lg3 += 1/rho_i;

    if(a) for(int j=0; j!=cols; j++) if(j!=i){

      r_ij=0;
      for(int x=0; x!=rows; x++){
	pos_ij(x) = pos(i,x)-pos(j,x);
	r_ij     += sqr(pos_ij(x));
      }
      r_ij = sqrt(r_ij);
      
      if(r_ij<d){
	du_ij   = du(r_ij)/r_ij;
	pos_ij *= du_ij;
	
	for(int x=0; x!=rows; x++){
	  e_lfg1 += pos_ij(x)*pos_i(x);
	  if(x!=2)
	    e_lfg2 += pos_ij(x)*pos_i(x)/rho_i;
	}
	e_lf1 += ddu(r_ij) + 2*du_ij;
	
	for(int k=0; k!=cols; k++) if(k!=i){
	  
	  r_ik=0;
	  for(int x=0; x!=rows; x++){
	    pos_ik(x) = pos(i,x)-pos(k,x);
	    r_ik     += sqr(pos_ik(x));
	  }
	  r_ik = sqrt(r_ik);
	  if(r_ik<d){
	    pos_ik *= du(r_ik)/r_ik;
	    
	    for(int x=0; x!=rows; x++)
	      e_lf2 += pos_ij(x)*pos_ik(x);
	  }
	}
      }
    }
  }

  kinetic_e += -sqr(2*alpha)*e_lg2 - e_lg3 - e_lf1 - e_lf2 
    + 4*alpha*e_lfg1 - 2*e_lfg2;

  double potential_e = 0, tmp;

  for(int i=0; i!=cols; i++){
    rho_i=0;
    for(int j=0; j!=rows; j++){
      tmp = pos(i,j);
      if(j==2){
	potential_e += sqr(lambda*tmp);
      }
      else{
	rho_i += sqr(tmp);
      }
    }
    potential_e += rho_i + sqr(kappa)/rho_i;
  }

  return (kinetic_e+potential_e);

}

double PollsVortexLocalEnergy:: valuePt(QickArray &pos, int dummy){

  int cols, rows, dim;
  dim = pos.get_dimension_info(&cols,&rows);

  if(!tmp_pos.is_initialized()){
    tmp_pos.redim(rows);
    pos_i.redim(rows);
  }

  double r_ij;

  double e_lg2=0, e_lg3=0;
  double kinetic_e=0;
  double du_ij;
  double rho_i;
  for(int i=0; i!=cols-1; i++){
    rho_i=0;
    for(int x=0; x!=rows; x++)
      if(x==2){
	pos_i(x) = lambda*pos(i,x);
	e_lg2   += sqr(pos_i(x));
      }
      else{
	pos_i(x) = pos(i,x);
	rho_i   += sqr(pos_i(x));
      }

    e_lg3 += 1/rho_i;
    for(int j=i+1; j!=cols; j++){
      r_ij = 0;
      for(int x=0; x!=rows; x++)
	r_ij += sqr(pos(i,x)-pos(j,x));
      r_ij = sqrt(r_ij);
      if(r_ij<d){
	du_ij = du(r_ij);
	kinetic_e += -2*du_ij/r_ij+sqr(du_ij)-ddu(r_ij);
	//kinetic_e += sqr(du_ij);
      }
    }
  }
  kinetic_e *= .5;

  kinetic_e += alpha*cols*(2+lambda) - e_lg3;

  double potential_e = 0, tmp;

  for(int i=0; i!=cols; i++){
    rho_i=0;
    for(int j=0; j!=rows; j++){
      tmp = pos(i,j);
      if(j==2){
	potential_e += sqr(lambda*tmp);
      }
      else{
	rho_i += sqr(tmp);
      }
    }
    potential_e += rho_i + sqr(kappa)/rho_i;
  }

  return kinetic_e+potential_e;

}

//=================================================================
// g(alpha,lambda,r_i)f(a,r_ij)
//=================================================================
double PollsWave:: valuePt(QickArray &pos){

  int cols, rows, dim;
  dim = pos.get_dimension_info(&cols,&rows);

  if(!tmp_pos.is_initialized())
    tmp_pos.redim(rows);
  
  double f = 0;
  
  // f(a,r_ij)
  for(int i=moved_particle+1; i!=cols; i++){
    for(int j=0; j!=rows; j++)
      tmp_pos(j) = pos(moved_particle,j)-pos(i,j);
    double r_ij = 0;
    for(int j=0; j!=rows; j++){
      r_ij += sqr(tmp_pos(j));
    }
    r_ij = sqrt(r_ij);
    if(r_ij>a)
      if(r_ij<d)
	f += log(d/r_ij*sin(K*(r_ij-a))/sin(K*(d-a)));
      else
	f += 0;
    else{
      f=log(0.);
      break;
    }
  }
  for(int i=0; i!=moved_particle; i++){
    for(int j=0; j!=rows; j++)
      tmp_pos(j) = pos(moved_particle,j)-pos(i,j);
    double r_ij = 0;
    for(int j=0; j!=rows; j++){
      r_ij += sqr(tmp_pos(j));
    }
    r_ij = sqrt(r_ij);
    if(r_ij>a)
      if(r_ij<d)
	f += log(d/r_ij*sin(K*(r_ij-a))/sin(K*(d-a)));
      else
	f += 0;
    else{
      f=log(0.);
      break;
    }
  }
  
  if(f==log(0.))
    return f;

  double g=0;
  // g(alpha,lambda,r_ij)
  for(int j=0; j!=rows; j++)
    if(j==2)
      g+=lambda*sqr(pos(moved_particle,j));
    else
      g+=sqr(pos(moved_particle,j));
  
    
  //g = exp(-alpha*g);
  g = (-alpha*g);

  return g+f;

}


//=================================================================
// g(alpha,lambda,r_i)
//=================================================================
double PollsGauss:: valuePt(QickArray &pos){

  int cols, rows, dim;
  dim = pos.get_dimension_info(&cols,&rows);

  double g=0;
  // g(alpha,lambda,r_ij)
  for(int j=0; j!=rows; j++)
    if(j==2)
      g+=lambda*sqr(pos(moved_particle,j));
    else
      g+=sqr(pos(moved_particle,j));
  
    
  g = -alpha*g;

  return g;

}



//=================================================================
// f(a,r_ij)
//=================================================================
double PollsJastrow:: valuePt(QickArray &pos){

  int cols, rows, dim;
  dim = pos.get_dimension_info(&cols,&rows);

  if(!tmp_pos.is_initialized())
    tmp_pos.redim(rows);

  double f = 0;
  // f(a,r_ij)
  for(int moved_particle=0; moved_particle!=cols-1; moved_particle++)
    for(int i=moved_particle+1; i!=cols; i++){
      for(int j=0; j!=rows; j++)
	tmp_pos(j) = pos(moved_particle,j)-pos(i,j);
      double r_ij = 0;
      for(int j=0; j!=rows; j++){
	r_ij += sqr(tmp_pos(j));
      }
      if(r_ij > sqr(a) && r_ij < sqr(d))
	f += log(d/r_ij*sin(K*(r_ij-a))/sin(K*(d-a)));
      if(r_ij <= a){
	f = log(0.);
	break;
      }
    }
  
  return f;

}


//=================================================================
// g(alpha,lambda,r_i)f(a,r_ij)
//=================================================================
double PollsWaveAll:: valuePt(QickArray &pos){

  int cols, rows, dim;
  dim = pos.get_dimension_info(&cols,&rows);

  if(!tmp_pos.is_initialized())
    tmp_pos.redim(rows);

  double f = 0;


  // f(a,r_ij)
  for(int moved_particle=0; moved_particle!=cols-1; moved_particle++)
    for(int i=moved_particle+1; i!=cols; i++){
      for(int j=0; j!=rows; j++)
	tmp_pos(j) = pos(moved_particle,j)-pos(i,j);
      double r_ij = 0;
      for(int j=0; j!=rows; j++){
	r_ij += sqr(tmp_pos(j));
      }
      r_ij/=sqr(a);
      if(r_ij < sqr(d))
	f += log(sin(K*(sqrt(r_ij)-a))/sqrt(r_ij)*d/sin(K*(d-a)));
      if(r_ij <= a){
	f = log(0.);
	break;
      }
	//else
    }

  if(f==log(0.))
    return f;

  double g=0;

  // g(alpha,lambda,r_ij)
  for(int moved_particle=0; moved_particle!=cols; moved_particle++)
  for(int j=0; j!=rows; j++)
    if(j==2)
      g+=lambda*sqr(pos(moved_particle,j));
    else
      g+=sqr(pos(moved_particle,j));

  g = -alpha*g;
  return f+g;

}

//=================================================================
// 2*nabla psi / psi 
//=================================================================
QickArray& PollsQForce:: valuePt(QickArray &pos, bool dummy){

  int cols, rows, dim;
  dim = pos.get_dimension_info(&cols,&rows);

  if(!pos_ij.is_initialized())
    pos_ij.redim(rows);

  if(!dF.is_initialized()){
    dF.redim(cols,rows);
    dG.redim(cols,rows);
  }

  dF *= 0;

  double r_ij;

  for(int i=0; i!=cols; i++){

    for(int x=0; x!=rows; x++)
      if(x==2)
	dG(i,x) = lambda*pos(i,x);
      else
	dG(i,x) = pos(i,x);

    for(int j=0; j!=cols; j++) if(i!=j){

      r_ij=0;
      for(int x=0; x!=rows; x++){
	pos_ij(x) = pos(i,x)-pos(j,x);
	r_ij += sqr(pos_ij(x));
      }
      r_ij = sqrt(r_ij);

      if(r_ij<=a){
	dF*=0;
	return dF;
      }
      
      pos_ij *= du(r_ij)/r_ij;

      for(int x=0; x!=rows; x++)
	dF(i,x) += pos_ij(x);

    }

  }

  dG *= -2*alpha;

  dF += dG;

  dF *= 2;

  return dF;

}

double PollsLocalEnergy:: valuePt(QickArray &pos){

  int cols, rows, dim;
  dim = pos.get_dimension_info(&cols,&rows);

  if(!pos_ij.is_initialized()){
    pos_ij.redim(rows);
    pos_i.redim(rows);
    pos_ik.redim(rows);
  }

  double kinetic_e = 2*alpha*cols*(2+lambda);

  double e_lg2=0, e_lf1=0, e_lf2=0, e_lfg=0, r_ij, r_ik, du_ij;

  for(int i=0; i!=cols; i++){

    for(int x=0; x!=rows; x++)
      if(x==2){
	pos_i(x) = lambda*pos(i,x);
	e_lg2   += sqr(pos_i(x));
      }
      else{
	pos_i(x) = pos(i,x);
	e_lg2   += sqr(pos_i(x));
      }

    if(d) for(int j=0; j!=cols; j++) if(j!=i){

      r_ij=0;
      for(int x=0; x!=rows; x++){
	pos_ij(x) = pos(i,x)-pos(j,x);
	r_ij     += sqr(pos_ij(x));
      }
      r_ij = sqrt(r_ij);
      
      if(r_ij<d){
	du_ij   = du(r_ij)/r_ij;
	pos_ij *= du_ij;

	for(int x=0; x!=rows; x++)
	  e_lfg += pos_ij(x)*pos_i(x);
      
	e_lf1 += ddu(r_ij) + 2*du_ij;
      
	for(int k=0; k!=cols; k++) if(k!=i){
	  
	  r_ik=0;
	  for(int x=0; x!=rows; x++){
	    pos_ik(x) = pos(i,x)-pos(k,x);
	    r_ik     += sqr(pos_ik(x));
	  }
	  r_ik = sqrt(r_ik);
	  if(r_ik<d){
	    pos_ik *= du(r_ik)/r_ik;
	    
	    for(int x=0; x!=rows; x++)
	      e_lf2 += pos_ij(x)*pos_ik(x);
	  }
	}
      }
    }
  }

  kinetic_e += -sqr(2*alpha)*e_lg2 - e_lf1 - e_lf2 + 4*alpha*e_lfg;

  double potential_e = 0;

  for(int i=0; i!=cols; i++)
    for(int j=0; j!=rows; j++){
      if(j==2)
	potential_e += sqr(lambda*pos(i,j));
      else
	potential_e += sqr(pos(i,j));
    }

  return (kinetic_e+potential_e);

}

double PollsLocalEnergy:: valuePt(QickArray &pos, int dummy){

  int cols, rows, dim;
  dim = pos.get_dimension_info(&cols,&rows);

  if(!tmp_pos.is_initialized())
    tmp_pos.redim(rows);

  double r_ij;

  double kinetic_e=0;
  double du_ij;
  for(int i=0; i!=cols-1; i++)
    for(int j=i+1; j!=cols; j++){
      r_ij = 0;
      for(int x=0; x!=rows; x++)
	r_ij += sqr(pos(i,x)-pos(j,x));
      r_ij = sqrt(r_ij);
      if(r_ij<d){
	du_ij = du(r_ij);
	//kinetic_e += -2*du_ij/r_ij+sqr(du_ij)-ddu(r_ij);
	kinetic_e += sqr(du_ij);
      }
    }
      
  //kinetic_e *= .5;
  
  kinetic_e += alpha*cols*(2+lambda);

  double potential_e = 0;

  for(int i=0; i!=cols; i++)
    for(int j=0; j!=rows; j++){
      if(j==2)
	potential_e += sqr(lambda*pos(i,j));
      else
	potential_e += sqr(pos(i,j));
    }

  return kinetic_e+potential_e;

}

//=================================================================
// g(alpha,lambda,r_i)f(a,r_ij)
//=================================================================
double DuBoisWave:: valuePt(QickArray &pos){

  int cols, rows, dim;
  dim = pos.get_dimension_info(&cols,&rows);

  if(!tmp_pos.is_initialized())
    tmp_pos.redim(rows);
  
  double f = 0;
  
  // f(a,r_ij)
  if(a!=0) for(int i=moved_particle+1; i!=cols; i++){
    for(int j=0; j!=rows; j++)
      tmp_pos(j) = pos(moved_particle,j)-pos(i,j);
    double r_ij = 0;
    for(int j=0; j!=rows; j++){
      r_ij += sqr(tmp_pos(j));
    }
    r_ij = sqrt(r_ij);
    if(r_ij>a)
      //      if(r_ij<d)
	f += log(1-a/r_ij);
    else{
      f=log(0.);
      break;
    }
  }
  if(a!=0) for(int i=0; i!=moved_particle; i++){
    for(int j=0; j!=rows; j++)
      tmp_pos(j) = pos(moved_particle,j)-pos(i,j);
    double r_ij = 0;
    for(int j=0; j!=rows; j++){
      r_ij += sqr(tmp_pos(j));
    }
    r_ij = sqrt(r_ij);
    if(r_ij>a)
      //      if(r_ij<d)
	f += log(1-a/r_ij);
    else{
      f=log(0.);
      break;
    }
  }
  
  if(f==log(0.))
    return f;

  double g=0;
  // g(alpha,lambda,r_ij)
  for(int j=0; j!=rows; j++)
    if(j==2)
      g+=lambda*sqr(pos(moved_particle,j));
    else
      g+=sqr(pos(moved_particle,j));
  
    
  //g = exp(-alpha*g);
  g = (-alpha*g);

  return g+f;

}


//=================================================================
// g(alpha,lambda,r_i)
//=================================================================
double DuBoisGauss:: valuePt(QickArray &pos){

  int cols, rows, dim;
  dim = pos.get_dimension_info(&cols,&rows);
  
  double g=0;
  // g(alpha,lambda,r_ij)
  for(int j=0; j!=rows; j++)
    if(j==2)
      g+=lambda*sqr(pos(moved_particle,j));
    else
      g+=sqr(pos(moved_particle,j));
  
    
  g = (-alpha*g);

  return g;

}



//=================================================================
// f(a,r_ij)
//=================================================================
double DuBoisJastrow:: valuePt(QickArray &pos){

  int cols, rows, dim;
  dim = pos.get_dimension_info(&cols,&rows);

  if(!tmp_pos.is_initialized())
    tmp_pos.redim(rows);

  double f = 1;
  // f(a,r_ij)
  for(int moved_particle=0; moved_particle!=cols-1; moved_particle++)
    for(int i=moved_particle+1; i!=cols; i++){
      for(int j=0; j!=rows; j++)
	tmp_pos(j) = pos(moved_particle,j)-pos(i,j);
      double r_ij = 0;
      for(int j=0; j!=rows; j++){
	r_ij += sqr(tmp_pos(j));
      }
      r_ij = sqrt(r_ij);
      if(r_ij > a)
	//	if(r_ij < d)
	  f *= 1-a/r_ij;
      else{
	f = 0;
	break;
      }
    }
  
  return f;

}


//=================================================================
// g(alpha,lambda,r_i)f(a,r_ij)
//=================================================================
double DuBoisWaveAll:: valuePt(QickArray &pos){

  int cols, rows, dim;
  dim = pos.get_dimension_info(&cols,&rows);

  if(!tmp_pos.is_initialized())
    tmp_pos.redim(rows);

  double f = 0;

  // f(a,r_ij)
  if(a!=0) for(int moved_particle=0; moved_particle!=cols; moved_particle++)
    for(int i=moved_particle+1; i!=cols; i++){
      for(int j=0; j!=rows; j++)
	tmp_pos(j) = pos(moved_particle,j)-pos(i,j);
      double r_ij = 0;
      for(int j=0; j!=rows; j++){
	r_ij += sqr(tmp_pos(j));
      }
      r_ij = sqrt(r_ij);
      if(r_ij > a)
	//	if(r_ij < d)
	  f += log(1-a/r_ij);
      else{
	cerr << "klynk" << endl;
	f=log(0.);
	break;
      }
    }

  if(f==log(0.))
    return f;

  double g=0;

  // g(alpha,lambda,r_ij)
  for(int moved_particle=0; moved_particle!=cols; moved_particle++)
  for(int j=0; j!=rows; j++)
    if(j==2)
      g+=lambda*sqr(pos(moved_particle,j));
    else
      g+=sqr(pos(moved_particle,j));

  g = -alpha*g;

  return (f+g);

}

QickArray& DuBoisQForce:: valuePt(QickArray &pos, bool dummy){

  int cols, rows, dim;
  dim = pos.get_dimension_info(&cols,&rows);

  if(dim != 2) cerr << dim << " ";

  if(!pos_ij.is_initialized())
    pos_ij.redim(rows);

  if(!dF.is_initialized()){
    dF.redim(cols,rows);
    dG.redim(cols,rows);
  }

  dF *= 0;

  double r_ij;

  for(int i=0; i!=cols; i++){

    for(int x=0; x!=rows; x++)
      if(x==2)
	dG(i,x) = lambda*pos(i,x);
      else
	dG(i,x) = pos(i,x);

    if(a!=0) for(int j=0; j!=cols; j++) if(i!=j){

      r_ij=0;
      for(int x=0; x!=rows; x++){
	pos_ij(x) = pos(i,x)-pos(j,x);
	r_ij += sqr(pos_ij(x));
      }
      r_ij = sqrt(r_ij);

      if(r_ij<=a){
	dF*=0;
	return dF;
      }
      
      if(r_ij < d){
	pos_ij *= du(r_ij)/r_ij;
	
	for(int x=0; x!=rows; x++)
	  dF(i,x) += pos_ij(x);
      }
    }

  }

  dG *= -2*alpha;

  dF += dG;

  dF *= 2;

  return dF;

}

double DuBoisLocalEnergy:: valuePt(QickArray &pos){

  int cols, rows, dim;
  dim = pos.get_dimension_info(&cols,&rows);

  if(!pos_ij.is_initialized()){
    pos_ij.redim(rows);
    pos_i.redim(rows);
    pos_ik.redim(rows);
  }

  double kinetic_e = 2*alpha*cols*(2+lambda);

  double e_lg2=0, e_lf1=0, e_lf2=0, e_lfg=0, r_ij, r_ik, du_ij;

  for(int i=0; i!=cols; i++){

    for(int x=0; x!=rows; x++)
      if(x==2){
	pos_i(x) = lambda*pos(i,x);
	e_lg2   += sqr(pos_i(x));
      }
      else{
	pos_i(x) = pos(i,x);
	e_lg2   += sqr(pos_i(x));
      }

    if(a!=0) for(int j=0; j!=cols; j++) if(j!=i){
      
      r_ij=0;
      for(int x=0; x!=rows; x++){
	pos_ij(x) = pos(i,x)-pos(j,x);
	r_ij     += sqr(pos_ij(x));
      }
      r_ij = sqrt(r_ij);
      if(r_ij<=a){
	return 0.;
      }
      if(r_ij < d){
	du_ij   = du(r_ij)/r_ij;
	pos_ij *= du_ij;
	
	for(int x=0; x!=rows; x++)
	  e_lfg += pos_ij(x)*pos_i(x);
	
	e_lf1 += ddu(r_ij) + 2*du_ij;
	
	for(int k=0; k!=cols; k++) if(k!=i){
	  
	  r_ik=0;
	  for(int x=0; x!=rows; x++){
	    pos_ik(x) = pos(i,x)-pos(k,x);
	    r_ik     += sqr(pos_ik(x));
	  }
	  r_ik = sqrt(r_ik);
	  if(r_ik<=a){
	    return 0.;
	  }
	  if(r_ik < d){
	    pos_ik *= du(r_ik)/r_ik;
	    
	    for(int x=0; x!=rows; x++)
	      e_lf2 += pos_ij(x)*pos_ik(x);
	  }
	}
      }
    }
  }
  
  kinetic_e += -sqr(2*alpha)*e_lg2 - e_lf1 - e_lf2 + 4*alpha*e_lfg;
  
  double potential_e = 0;
  
  for(int i=0; i!=cols; i++)
    for(int j=0; j!=rows; j++){
      if(j==2)
	potential_e += sqr(lambda*pos(i,j));
      else
	potential_e += sqr(pos(i,j));
    }
  
  return (kinetic_e+potential_e);

}

double DuBoisLocalEnergy:: valuePt(QickArray &pos, int dummy){

  int cols, rows, dim;
  dim = pos.get_dimension_info(&cols,&rows);

  if(!tmp_pos.is_initialized())
    tmp_pos.redim(rows);

  double r_ij;

  double kinetic_e=0;
  double du_ij;
  if(a!=0) for(int i=0; i!=cols-1; i++)
    for(int j=i+1; j!=cols; j++){
      r_ij = 0;
      for(int x=0; x!=rows; x++)
	r_ij += sqr(pos(i,x)-pos(j,x));
      r_ij = sqrt(r_ij);
      if(r_ij < d){
	du_ij = du(r_ij);
	//kinetic_e += -2*du_ij/r_ij+sqr(du_ij)-ddu(r_ij);
	kinetic_e += sqr(du_ij);
      }
    }
      
  //kinetic_e *= .5;
  kinetic_e += alpha*cols*(2+lambda);

  double potential_e = 0;

  for(int i=0; i!=cols; i++)
    for(int j=0; j!=rows; j++){
      if(j==2)
	potential_e += sqr(lambda*pos(i,j));
      else
	potential_e += sqr(pos(i,j));
    }

  return kinetic_e+potential_e;

}

double DuBoisPotentialEnergy:: valuePt(QickArray &pos){

  int cols, rows, dim;
  dim = pos.get_dimension_info(&cols,&rows);

  double potential_e = 0;

  for(int i=0; i!=cols; i++)
    for(int j=0; j!=rows; j++){
      if(j==2)
	potential_e += sqr(lambda*pos(i,j));
      else
	potential_e += sqr(pos(i,j));
    }

  return potential_e;

}

//=================================================================
// g(alpha,lambda,r_i)f(a,r_ij)
//=================================================================
double DuBoisVortexWave:: valuePt(QickArray &pos){

  int cols, rows, dim;
  dim = pos.get_dimension_info(&cols,&rows);

  if(!tmp_pos.is_initialized())
    tmp_pos.redim(rows);
  
  double f = 0;
  
  // f(a,r_ij)
  for(int i=moved_particle+1; i!=cols; i++){
    for(int j=0; j!=rows; j++)
      tmp_pos(j) = pos(moved_particle,j)-pos(i,j);
    double r_ij = 0;
    for(int j=0; j!=rows; j++){
      r_ij += sqr(tmp_pos(j));
    }
    r_ij = sqrt(r_ij);
    if(r_ij>a)
      if(r_ij < d)
	f += log(1-a/r_ij);
    else{
      f=log(0.);
      break;
    }
  }
  for(int i=0; i!=moved_particle; i++){
    for(int j=0; j!=rows; j++)
      tmp_pos(j) = pos(moved_particle,j)-pos(i,j);
    double r_ij = 0;
    for(int j=0; j!=rows; j++){
      r_ij += sqr(tmp_pos(j));
    }
    r_ij = sqrt(r_ij);
    if(r_ij>a)
      if(r_ij<d)
	f += log(1-a/r_ij);
    else{
      f=log(0.);
      break;
    }
  }
  
  if(f==log(0.))
    return f;

  double g=0, rho_i=0;
  // g(alpha,lambda,r_ij)
  for(int j=0; j!=rows; j++)
    if(j==2)
      g+=lambda*sqr(pos(moved_particle,j));
    else{
      rho_i+=sqr(pos(moved_particle,j));
    }

  g += rho_i;
  g  = log(sqrt(rho_i))-alpha*g;

  return g+f;

}

//=================================================================
// g(alpha,lambda,r_i)f(a,r_ij)
//=================================================================
double DuBoisVortexWaveAll:: valuePt(QickArray &pos){

  int cols, rows, dim;
  dim = pos.get_dimension_info(&cols,&rows);

  if(!tmp_pos.is_initialized())
    tmp_pos.redim(rows);

  double F = 0;

  // F(a,r_ij)
  for(int moved_particle=0; moved_particle!=cols-1; moved_particle++)
    for(int i=moved_particle+1; i!=cols; i++){
      for(int j=0; j!=rows; j++)
	tmp_pos(j) = pos(moved_particle,j)-pos(i,j);
      double r_ij = 0;
      for(int j=0; j!=rows; j++){
	r_ij += sqr(tmp_pos(j));
      }
      r_ij = sqrt(r_ij);
      if(r_ij > a)
	if(r_ij < d)
	  F += log(1-a/r_ij);
      else{
	F = log(0.);
	break;
      }
    }

  if(F==log(0.))
    return F;

  double g=0, G=0, rho_i=0;
  // G(alpha,lambda,r_ij)
  for(int moved_particle=0; moved_particle!=cols; moved_particle++){
    g = rho_i = 0;
    for(int j=0; j!=rows; j++)
      if(j==2)
	g+=lambda*sqr(pos(moved_particle,j));
      else{
	rho_i+=sqr(pos(moved_particle,j));
	
      }	
    g += rho_i;
    G += log(sqrt(rho_i))-alpha*g;
  }

  return F+G;

}

//=================================================================
// 2*nabla psi / psi !?!?!?!?!?!?!?!?!?!?!?!?!?!?!?!?!?!?!?!?!?!?!?!?
//=================================================================
QickArray& DuBoisVortexQForce:: valuePt(QickArray &pos, bool dummy){

  int cols, rows, dim;
  dim = pos.get_dimension_info(&cols,&rows);

  if(!pos_ij.is_initialized())
    pos_ij.redim(rows);

  if(!dF.is_initialized()){
    dF.redim(cols,rows);
    dG.redim(cols,rows);
  }

  dF *= 0;

  double r_ij, xy_fac, rho_i;

  for(int i=0; i!=cols; i++){

    rho_i=0;
    for(int x=0; x!=2; x++)
      rho_i += sqr(pos(i,x));

    xy_fac = 1-1/(2*alpha*rho_i);

    for(int x=0; x!=rows; x++)
      if(x==2)
	dG(i,x) = lambda*pos(i,x);
      else
	dG(i,x) = xy_fac*pos(i,x);

    for(int j=0; j!=cols; j++) if(i!=j){

      r_ij=0;
      for(int x=0; x!=rows; x++){
	pos_ij(x) = pos(i,x)-pos(j,x);
	r_ij += sqr(pos_ij(x));
      }
      r_ij = sqrt(r_ij);
      
      if(r_ij<d){
	pos_ij *= du(r_ij)/r_ij;
	
	for(int x=0; x!=rows; x++)
	  dF(i,x) += pos_ij(x);
      }
    }

  }

  dG *= -2*alpha;

  dF += dG;

  dF *= 2;

  return dF;

}

double DuBoisVortexLocalEnergy:: valuePt(QickArray &pos, int dummy){


  int cols, rows, dim;
  dim = pos.get_dimension_info(&cols,&rows);

  if(!tmp_pos.is_initialized()){
    tmp_pos.redim(rows);
  }

  if(!pos_i.is_initialized()){
    pos_i.redim(rows);
  }

  double r_ij;

  double e_lg3=0;
  double kinetic_e=0;
  double du_ij;
  double rho_i;
  for(int i=0; i!=cols-1; i++){
    rho_i=0;
    for(int x=0; x!=rows; x++)
      if(x==2){
	pos_i(x) = lambda*pos(i,x);
      }
      else{
	pos_i(x) = pos(i,x);
	rho_i   += sqr(pos_i(x));
      }

    e_lg3 += 1/rho_i;
    for(int j=i+1; j!=cols; j++){
      r_ij = 0;
      for(int x=0; x!=rows; x++)
	r_ij += sqr(pos(i,x)-pos(j,x));
      r_ij = sqrt(r_ij);
      if(r_ij<d){
	du_ij = du(r_ij);
	kinetic_e += sqr(du_ij);
      }
    }
  }

  kinetic_e += alpha*cols*(2+lambda) - e_lg3;

  double potential_e = 0, tmp;

  for(int i=0; i!=cols; i++){
    rho_i=0;
    for(int j=0; j!=rows; j++){
      tmp = pos(i,j);
      if(j==2){
	potential_e += sqr(lambda*tmp);
      }
      else{
	rho_i += sqr(tmp);
      }
    }
    potential_e += rho_i + sqr(kappa)/rho_i;
  }

  return kinetic_e+potential_e;

}

double DuBoisVortexLocalEnergy:: valuePt(QickArray &pos){

  int cols, rows, dim;
  dim = pos.get_dimension_info(&cols,&rows);

  if(!pos_ij.is_initialized()){
    pos_ij.redim(rows);
    pos_i.redim(rows);
    pos_ik.redim(rows);
  }

  double kinetic_e = 2*alpha*cols*(4+lambda);

  double e_lg2=0, e_lg3=0, e_lf1=0, e_lf2=0, e_lfg1=0, e_lfg2=0;
  double r_ij, r_ik, rho_i, du_ij;

  for(int i=0; i!=cols; i++){

    rho_i=0;
    for(int x=0; x!=rows; x++)
      if(x==2){
	pos_i(x) = lambda*pos(i,x);
	e_lg2   += sqr(pos_i(x));
      }
      else{
	pos_i(x) = pos(i,x);
	rho_i   += sqr(pos_i(x));
      }

    e_lg2 += rho_i;

    e_lg3 += 1/rho_i;

    if(a) for(int j=0; j!=cols; j++) if(j!=i){

      r_ij=0;
      for(int x=0; x!=rows; x++){
	pos_ij(x) = pos(i,x)-pos(j,x);
	r_ij     += sqr(pos_ij(x));
      }
      r_ij = sqrt(r_ij);
      if(r_ij<d){
	du_ij   = du(r_ij)/r_ij;
	pos_ij *= du_ij;
	
	for(int x=0; x!=rows; x++){
	  e_lfg1 += pos_ij(x)*pos_i(x);
	  if(x!=2)
	    e_lfg2 += pos_ij(x)*pos_i(x)/rho_i;
	}
	e_lf1 += ddu(r_ij) + 2*du_ij;
	
	for(int k=0; k!=cols; k++) if(k!=i){
	  
	  r_ik=0;
	  for(int x=0; x!=rows; x++){
	    pos_ik(x) = pos(i,x)-pos(k,x);
	    r_ik     += sqr(pos_ik(x));
	  }
	  r_ik = sqrt(r_ik);
	  
	  if(r_ik < d){
	    pos_ik *= du(r_ik)/r_ik;
	    
	    for(int x=0; x!=rows; x++)
	      e_lf2 += pos_ij(x)*pos_ik(x);
	  }
	}
      }
    }

  }

  kinetic_e += -sqr(2*alpha)*e_lg2 - e_lg3 - e_lf1 - e_lf2 
    + 4*alpha*e_lfg1 - 2*e_lfg2;

  double potential_e = 0, tmp;

  for(int i=0; i!=cols; i++){
    rho_i=0;
    for(int j=0; j!=rows; j++){
      tmp = pos(i,j);
      if(j==2){
	potential_e += sqr(lambda*tmp);
      }
      else{
	rho_i += sqr(tmp);
      }
    }
    potential_e += rho_i + sqr(kappa)/rho_i;
  }

  return (kinetic_e+potential_e);

}

double DuBoisVortexPotentialEnergy:: valuePt(QickArray &pos){

  int cols, rows, dim;
  dim = pos.get_dimension_info(&cols,&rows);

  double potential_e = 0, tmp;

  for(int i=0; i!=cols; i++){
    double rho_i=0;
    for(int j=0; j!=rows; j++){
      tmp = pos(i,j);
      if(j!=2){
	rho_i += sqr(tmp);
      }
    }
    potential_e += sqr(kappa)/rho_i;
  }

  return potential_e;

}

double DuBoisVortexExcitationEnergy:: valuePt(QickArray &pos){

  int cols, rows, dim;
  dim = pos.get_dimension_info(&cols,&rows);

  double rho_i, excitation_e;

  excitation_e = 0;
  for(int i=0; i!=cols; i++){

    rho_i=0;
    for(int x=0; x!=rows; x++)
      if(x!=2)
	rho_i   += sqr(pos(i,x));
    
    excitation_e += 1/rho_i;
  }  

  return excitation_e*(1+kappa2);

}


//=================================================================
// g(alpha,lambda,r_i)f(a,r_ij)
//=================================================================
double ReattoVortexWave:: valuePt(QickArray &pos){

  int cols, rows, dim;
  dim = pos.get_dimension_info(&cols,&rows);

  if(!tmp_pos.is_initialized())
    tmp_pos.redim(rows);
  
  double f = 0;
  
  // f(a,r_ij)
  for(int i=moved_particle+1; i!=cols; i++){
    for(int j=0; j!=rows; j++)
      tmp_pos(j) = pos(moved_particle,j)-pos(i,j);
    double r_ij = 0;
    for(int j=0; j!=rows; j++){
      r_ij += sqr(tmp_pos(j));
    }
    r_ij = sqrt(r_ij);
    if(r_ij>a)
      if(r_ij < d)
	f += log(1-a/r_ij);
    else{
      f=log(0.);
      break;
    }
  }
  for(int i=0; i!=moved_particle; i++){
    for(int j=0; j!=rows; j++)
      tmp_pos(j) = pos(moved_particle,j)-pos(i,j);
    double r_ij = 0;
    for(int j=0; j!=rows; j++){
      r_ij += sqr(tmp_pos(j));
    }
    r_ij = sqrt(r_ij);
    if(r_ij>a)
      if(r_ij<d)
	f += log(1-a/r_ij);
    else{
      f=log(0.);
      break;
    }
  }
  
  if(f==log(0.))
    return f;

  double g=0, rho_i=0;
  // g(alpha,lambda,r_ij)
  for(int j=0; j!=rows; j++)
    if(j==2)
      g+=lambda*sqr(pos(moved_particle,j));
    else{
      rho_i+=sqr(pos(moved_particle,j));
    }

  g += rho_i;
  g  = u_vor(rho_i)-alpha*g;

  return g+f;

}

//=================================================================
// g(alpha,lambda,r_i)f(a,r_ij)
//=================================================================
double ReattoVortexWaveAll:: valuePt(QickArray &pos){

  int cols, rows, dim;
  dim = pos.get_dimension_info(&cols,&rows);

  if(!tmp_pos.is_initialized())
    tmp_pos.redim(rows);

  double F = 0;

  // F(a,r_ij)
  for(int moved_particle=0; moved_particle!=cols-1; moved_particle++)
    for(int i=moved_particle+1; i!=cols; i++){
      for(int j=0; j!=rows; j++)
	tmp_pos(j) = pos(moved_particle,j)-pos(i,j);
      double r_ij = 0;
      for(int j=0; j!=rows; j++){
	r_ij += sqr(tmp_pos(j));
      }
      r_ij = sqrt(r_ij);
      if(r_ij > a)
	if(r_ij < d)
	  F += log(1-a/r_ij);
      else{
	F = log(0.);
	break;
      }
    }

  if(F==log(0.))
    return F;

  double g=0, G=0, rho_i=0;
  // G(alpha,lambda,r_ij)
  for(int moved_particle=0; moved_particle!=cols; moved_particle++){
    g = rho_i = 0;
    for(int j=0; j!=rows; j++)
      if(j==2)
	g+=lambda*sqr(pos(moved_particle,j));
      else{
	rho_i+=sqr(pos(moved_particle,j));
	
      }	
    g += rho_i;
    G += u_vor(rho_i)-alpha*g;
  }

  return F+G;

}

//=================================================================
// 2*nabla psi / psi !?!?!?!?!?!?!?!?!?!?!?!?!?!?!?!?!?!?!?!?!?!?!?
//=================================================================
QickArray& ReattoVortexQForce:: valuePt(QickArray &pos, bool dummy){

  int cols, rows, dim;
  dim = pos.get_dimension_info(&cols,&rows);

  if(!pos_ij.is_initialized())
    pos_ij.redim(rows);

  if(!dF.is_initialized()){
    dF.redim(cols,rows);
    dG.redim(cols,rows);
  }

  dF *= 0;

  double r_ij, xy_fac, rho_i;

  for(int i=0; i!=cols; i++){

    rho_i=0;
    for(int x=0; x!=2; x++)
      rho_i += sqr(pos(i,x));

    xy_fac = 1-du_vor(rho_i)/(2*alpha*sqrt(rho_i));

    for(int x=0; x!=rows; x++)
      if(x==2)
	dG(i,x) = lambda*pos(i,x);
      else
	dG(i,x) = xy_fac*pos(i,x);

    for(int j=0; j!=cols; j++) if(i!=j){

      r_ij=0;
      for(int x=0; x!=rows; x++){
	pos_ij(x) = pos(i,x)-pos(j,x);
	r_ij += sqr(pos_ij(x));
      }
      r_ij = sqrt(r_ij);
      
      if(r_ij<d){
	pos_ij *= du(r_ij)/r_ij;
	
	for(int x=0; x!=rows; x++)
	  dF(i,x) += pos_ij(x);
      }
    }

  }

  dG *= -2*alpha;

  dF += dG;

  dF *= 2;

  return dF;

}

double ReattoVortexLocalEnergy:: valuePt(QickArray &pos, int dummy){


  int cols, rows, dim;
  dim = pos.get_dimension_info(&cols,&rows);

  if(!tmp_pos.is_initialized()){
    tmp_pos.redim(rows);
  }

  if(!pos_i.is_initialized()){
    pos_i.redim(rows);
  }

  double r_ij;

  double e_lg2=0, e_lg3=0;
  double kinetic_e=0;
  double du_ij;
  double rho_i;
  for(int i=0; i!=cols-1; i++){
    rho_i=0;
    for(int x=0; x!=rows; x++)
      if(x==2){
	pos_i(x) = lambda*pos(i,x);
	e_lg2   += sqr(pos_i(x));
      }
      else{
	pos_i(x) = pos(i,x);
	rho_i   += sqr(pos_i(x));
      }

    e_lg3 += du_vor(rho_i)/sqrt(rho_i);
    for(int j=i+1; j!=cols; j++){
      r_ij = 0;
      for(int x=0; x!=rows; x++)
	r_ij += sqr(pos(i,x)-pos(j,x));
      r_ij = sqrt(r_ij);
      if(r_ij<d){
	du_ij = du(r_ij);
	kinetic_e += sqr(du_ij);
      }
    }
  }

  kinetic_e += alpha*cols*(2+lambda) - e_lg3;

  double potential_e = 0, tmp;

  for(int i=0; i!=cols; i++){
    rho_i=0;
    for(int j=0; j!=rows; j++){
      tmp = pos(i,j);
      if(j==2){
	potential_e += sqr(lambda*tmp);
      }
      else{
	rho_i += sqr(tmp);
      }
    }
    potential_e += rho_i + sqr(kappa)/rho_i;
  }

  return kinetic_e+potential_e;

}

double ReattoVortexLocalEnergy:: valuePt(QickArray &pos){

  int cols, rows, dim;
  dim = pos.get_dimension_info(&cols,&rows);

  if(!pos_ij.is_initialized()){
    pos_ij.redim(rows);
    pos_i.redim(rows);
    pos_ik.redim(rows);
  }

  double kinetic_e = 2*alpha*cols*(2+lambda);

  double e_lg2=0, e_lg3=0, e_lf1=0, e_lf2=0, e_lfg1=0, e_lfg2=0;
  double r_ij, r_ik, rho_i, du_ij, du_vor_i;

  for(int i=0; i!=cols; i++){

    rho_i=0;
    for(int x=0; x!=rows; x++)
      if(x==2){
	pos_i(x) = lambda*pos(i,x);
	e_lg2   += sqr(pos_i(x));
      }
      else{
	pos_i(x) = pos(i,x);
	rho_i   += sqr(pos_i(x));
      }

    e_lg2 += rho_i;

    du_vor_i = du_vor(rho_i);
    rho_i = sqrt(rho_i);
    //e_lg3 += du_vor_i*(1/rho_i+rho_i/beta2-8*alpha);
    e_lg3 += du_vor_i*(2/rho_i-2*rho_i/beta2-4*alpha*rho_i);

    if(a) for(int j=0; j!=cols; j++) if(j!=i){

      r_ij=0;
      for(int x=0; x!=rows; x++){
	pos_ij(x) = pos(i,x)-pos(j,x);
	r_ij     += sqr(pos_ij(x));
      }
      r_ij = sqrt(r_ij);
      if(r_ij<d){
	du_ij   = du(r_ij)/r_ij;
	pos_ij *= du_ij;
	
	for(int x=0; x!=rows; x++){
	  e_lfg1 += pos_ij(x)*pos_i(x);
	  if(x!=2)
	    e_lfg2 += pos_ij(x)*pos_i(x)*du_vor_i/rho_i;
	}
	e_lf1 += ddu(r_ij) + 2*du_ij;
	
	for(int k=0; k!=cols; k++) if(k!=i){
	  
	  r_ik=0;
	  for(int x=0; x!=rows; x++){
	    pos_ik(x) = pos(i,x)-pos(k,x);
	    r_ik     += sqr(pos_ik(x));
	  }
	  r_ik = sqrt(r_ik);
	  
	  if(r_ik < d){
	    pos_ik *= du(r_ik)/r_ik;
	    
	    for(int x=0; x!=rows; x++)
	      e_lf2 += pos_ij(x)*pos_ik(x);
	  }
	}
      }
    }

  }

  kinetic_e += -sqr(2*alpha)*e_lg2 - e_lg3 - e_lf1 - e_lf2 
    + 4*alpha*e_lfg1 - 2*e_lfg2;

  double potential_e = 0, tmp;

  for(int i=0; i!=cols; i++){
    rho_i=0;
    for(int j=0; j!=rows; j++){
      tmp = pos(i,j);
      if(j==2){
	potential_e += sqr(lambda*tmp);
      }
      else{
	rho_i += sqr(tmp);
      }
    }
    potential_e += rho_i + sqr(kappa)/rho_i;
  }

  return (kinetic_e+potential_e);

}

double ReattoVortexPotentialEnergy:: valuePt(QickArray &pos){

  int cols, rows, dim;
  dim = pos.get_dimension_info(&cols,&rows);

  double potential_e = 0, tmp;

  for(int i=0; i!=cols; i++){
    double rho_i=0;
    for(int j=0; j!=rows; j++){
      tmp = pos(i,j);
      if(j!=2){
	rho_i += sqr(tmp);
      }
    }
    potential_e += sqr(kappa)/rho_i;
  }

  return potential_e;

}

double ReattoVortexExcitationEnergy:: valuePt(QickArray &pos){

  int cols, rows, dim;
  dim = pos.get_dimension_info(&cols,&rows);

  double rho_i, excitation_e, f_m;

  excitation_e = 0;
  for(int i=0; i!=cols; i++){

    rho_i=0;
    for(int x=0; x!=rows; x++)
      if(x!=2)
	rho_i   += sqr(pos(i,x));
    
    f_m = 1/(1-exp(-rho_i/beta2));

    excitation_e += kappa2/rho_i + (f_m*(f_m-2)+1)*4*rho_i/sqr(beta2);
  }  

  return excitation_e;

}

//=================================================================
// g(alpha,lambda,r_i)f(a,r_ij)
//=================================================================
double Reatto2VortexWave:: valuePt(QickArray &pos){

  int cols, rows, dim;
  dim = pos.get_dimension_info(&cols,&rows);

  if(!tmp_pos.is_initialized())
    tmp_pos.redim(rows);
  
  double f = 0;
  
  // f(a,r_ij)
  for(int i=moved_particle+1; i!=cols; i++){
    for(int j=0; j!=rows; j++)
      tmp_pos(j) = pos(moved_particle,j)-pos(i,j);
    double r_ij = 0;
    for(int j=0; j!=rows; j++){
      r_ij += sqr(tmp_pos(j));
    }
    r_ij = sqrt(r_ij);
    if(r_ij>a)
      if(r_ij < d)
	f += log(1-a/r_ij);
    else{
      f=log(0.);
      break;
    }
  }
  for(int i=0; i!=moved_particle; i++){
    for(int j=0; j!=rows; j++)
      tmp_pos(j) = pos(moved_particle,j)-pos(i,j);
    double r_ij = 0;
    for(int j=0; j!=rows; j++){
      r_ij += sqr(tmp_pos(j));
    }
    r_ij = sqrt(r_ij);
    if(r_ij>a)
      if(r_ij<d)
	f += log(1-a/r_ij);
    else{
      f=log(0.);
      break;
    }
  }
  
  if(f==log(0.))
    return f;

  double g=0, rho_i=0;
  // g(alpha,lambda,r_ij)
  for(int j=0; j!=rows; j++)
    if(j==2)
      g+=lambda*sqr(pos(moved_particle,j));
    else{
      rho_i+=sqr(pos(moved_particle,j));
    }

  g += rho_i;
  g  = u_vor(sqrt(rho_i))-alpha*g;

  return g+f;

}

//=================================================================
// g(alpha,lambda,r_i)f(a,r_ij)
//=================================================================
double Reatto2VortexWaveAll:: valuePt(QickArray &pos){

  int cols, rows, dim;
  dim = pos.get_dimension_info(&cols,&rows);

  if(!tmp_pos.is_initialized())
    tmp_pos.redim(rows);

  double F = 0;

  // F(a,r_ij)
  for(int moved_particle=0; moved_particle!=cols-1; moved_particle++)
    for(int i=moved_particle+1; i!=cols; i++){
      for(int j=0; j!=rows; j++)
	tmp_pos(j) = pos(moved_particle,j)-pos(i,j);
      double r_ij = 0;
      for(int j=0; j!=rows; j++){
	r_ij += sqr(tmp_pos(j));
      }
      r_ij = sqrt(r_ij);
      if(r_ij > a)
	if(r_ij < d)
	  F += log(1-a/r_ij);
      else{
	F = log(0.);
	break;
      }
    }

  if(F==log(0.))
    return F;

  double g=0, G=0, rho_i=0;
  // G(alpha,lambda,r_ij)
  for(int moved_particle=0; moved_particle!=cols; moved_particle++){
    g = rho_i = 0;
    for(int j=0; j!=rows; j++)
      if(j==2)
	g+=lambda*sqr(pos(moved_particle,j));
      else{
	rho_i+=sqr(pos(moved_particle,j));
	
      }	
    g += rho_i;
    G += u_vor(sqrt(rho_i))-alpha*g;
  }

  return F+G;

}

//=================================================================
// 2*nabla psi / psi !?!?!?!?!?!?!?!?!?!?!?!?!?!?!?!?!?!?!?!?!?!?!?
//=================================================================
QickArray& Reatto2VortexQForce:: valuePt(QickArray &pos, bool dummy){

  int cols, rows, dim;
  dim = pos.get_dimension_info(&cols,&rows);

  if(!pos_ij.is_initialized())
    pos_ij.redim(rows);

  if(!dF.is_initialized()){
    dF.redim(cols,rows);
    dG.redim(cols,rows);
  }

  dF *= 0;

  double r_ij, xy_fac, rho_i;

  for(int i=0; i!=cols; i++){

    rho_i=0;
    for(int x=0; x!=2; x++)
      rho_i += sqr(pos(i,x));

    rho_i = sqrt(rho_i);
    xy_fac = 1-du_vor(rho_i)/(2*alpha*rho_i);

    for(int x=0; x!=rows; x++)
      if(x==2)
	dG(i,x) = lambda*pos(i,x);
      else
	dG(i,x) = xy_fac*pos(i,x);

    for(int j=0; j!=cols; j++) if(i!=j){

      r_ij=0;
      for(int x=0; x!=rows; x++){
	pos_ij(x) = pos(i,x)-pos(j,x);
	r_ij += sqr(pos_ij(x));
      }
      r_ij = sqrt(r_ij);
      
      if(r_ij<d){
	pos_ij *= du(r_ij)/r_ij;
	
	for(int x=0; x!=rows; x++)
	  dF(i,x) += pos_ij(x);
      }
    }

  }

  dG *= -2*alpha;

  dF += dG;

  dF *= 2;

  return dF;

}

double Reatto2VortexLocalEnergy:: valuePt(QickArray &pos, int dummy){


  int cols, rows, dim;
  dim = pos.get_dimension_info(&cols,&rows);

  if(!tmp_pos.is_initialized()){
    tmp_pos.redim(rows);
  }

  if(!pos_i.is_initialized()){
    pos_i.redim(rows);
  }

  double r_ij;

  double e_lg2=0, e_lg3=0;
  double kinetic_e=0;
  double du_ij;
  double rho_i;
  for(int i=0; i!=cols-1; i++){
    rho_i=0;
    for(int x=0; x!=rows; x++)
      if(x==2){
	pos_i(x) = lambda*pos(i,x);
	e_lg2   += sqr(pos_i(x));
      }
      else{
	pos_i(x) = pos(i,x);
	rho_i   += sqr(pos_i(x));
      }

    rho_i = sqrt(rho_i);
    e_lg3 += du_vor(rho_i)/rho_i;
    for(int j=i+1; j!=cols; j++){
      r_ij = 0;
      for(int x=0; x!=rows; x++)
	r_ij += sqr(pos(i,x)-pos(j,x));
      r_ij = sqrt(r_ij);
      if(r_ij<d){
	du_ij = du(r_ij);
	kinetic_e += sqr(du_ij);
      }
    }
  }

  kinetic_e += alpha*cols*(2+lambda) - e_lg3;

  double potential_e = 0, tmp;

  for(int i=0; i!=cols; i++){
    rho_i=0;
    for(int j=0; j!=rows; j++){
      tmp = pos(i,j);
      if(j==2){
	potential_e += sqr(lambda*tmp);
      }
      else{
	rho_i += sqr(tmp);
      }
    }
    potential_e += rho_i + sqr(kappa)/rho_i;
  }

  return kinetic_e+potential_e;

}

double Reatto2VortexLocalEnergy:: valuePt(QickArray &pos){

  int cols, rows, dim;
  dim = pos.get_dimension_info(&cols,&rows);

  if(!pos_ij.is_initialized()){
    pos_ij.redim(rows);
    pos_i.redim(rows);
    pos_ik.redim(rows);
  }

  double kinetic_e = 2*alpha*cols*(2+lambda);

  double e_lg2=0, e_lg3=0, e_lf1=0, e_lf2=0, e_lfg1=0, e_lfg2=0;
  double r_ij, r_ik, rho_i, du_ij, du_vor_i;

  for(int i=0; i!=cols; i++){

    rho_i=0;
    for(int x=0; x!=rows; x++)
      if(x==2){
	pos_i(x) = lambda*pos(i,x);
	e_lg2   += sqr(pos_i(x));
      }
      else{
	pos_i(x) = pos(i,x);
	rho_i   += sqr(pos_i(x));
      }

    e_lg2 += rho_i;

    rho_i = sqrt(rho_i);

    du_vor_i = du_vor(rho_i);
    e_lg3 += du_vor_i*(1/rho_i-1/beta-4*alpha*rho_i);

    if(a) for(int j=0; j!=cols; j++) if(j!=i){

      r_ij=0;
      for(int x=0; x!=rows; x++){
	pos_ij(x) = pos(i,x)-pos(j,x);
	r_ij     += sqr(pos_ij(x));
      }
      r_ij = sqrt(r_ij);
      if(r_ij<d){
	du_ij   = du(r_ij)/r_ij;
	pos_ij *= du_ij;
	
	for(int x=0; x!=rows; x++){
	  e_lfg1 += pos_ij(x)*pos_i(x);
	  if(x!=2)
	    e_lfg2 += pos_ij(x)*pos_i(x)*du_vor_i/rho_i;
	}
	e_lf1 += ddu(r_ij) + 2*du_ij;
	
	for(int k=0; k!=cols; k++) if(k!=i){
	  
	  r_ik=0;
	  for(int x=0; x!=rows; x++){
	    pos_ik(x) = pos(i,x)-pos(k,x);
	    r_ik     += sqr(pos_ik(x));
	  }
	  r_ik = sqrt(r_ik);
	  
	  if(r_ik < d){
	    pos_ik *= du(r_ik)/r_ik;
	    
	    for(int x=0; x!=rows; x++)
	      e_lf2 += pos_ij(x)*pos_ik(x);
	  }
	}
      }
    }

  }

  kinetic_e += -sqr(2*alpha)*e_lg2 - e_lg3 - e_lf1 - e_lf2 
    + 4*alpha*e_lfg1 - 2*e_lfg2;

  double potential_e = 0, tmp;

  for(int i=0; i!=cols; i++){
    rho_i=0;
    for(int j=0; j!=rows; j++){
      tmp = pos(i,j);
      if(j==2){
	potential_e += sqr(lambda*tmp);
      }
      else{
	rho_i += sqr(tmp);
      }
    }
    potential_e += rho_i + sqr(kappa)/rho_i;
  }

  return (kinetic_e+potential_e);

}

double Reatto2VortexPotentialEnergy:: valuePt(QickArray &pos){

  int cols, rows, dim;
  dim = pos.get_dimension_info(&cols,&rows);

  double potential_e = 0, tmp;

  for(int i=0; i!=cols; i++){
    double rho_i=0;
    for(int j=0; j!=rows; j++){
      tmp = pos(i,j);
      if(j!=2){
	rho_i += sqr(tmp);
      }
    }
    potential_e += sqr(kappa)/rho_i;
  }

  return potential_e;

}

double Reatto2VortexExcitationEnergy:: valuePt(QickArray &pos){

  int cols, rows, dim;
  dim = pos.get_dimension_info(&cols,&rows);

  double rho_i, excitation_e, f_m;

  excitation_e = 0;
  for(int i=0; i!=cols; i++){

    rho_i=0;
    for(int x=0; x!=rows; x++)
      if(x!=2)
	rho_i   += sqr(pos(i,x));

    rho_i = sqrt(rho_i);

    f_m = 1/(1-exp(-rho_i/beta));

    excitation_e += kappa2/sqr(rho_i) + (f_m*(f_m-2)+1)/sqr(beta);
  }  

  return excitation_e;

}

//=================================================================
// exp(-sum_(i<j)(corr^5/r_ij^5))
//=================================================================
double LiqHe4Wave:: valuePt(QickArray &pos){

  int cols, rows, dim;
  dim = pos.get_dimension_info(&cols,&rows);
  QickArray tmp_pos(rows);
  double lh2 = sqr(l/2);
  double sum = 0;
  for(int i=0; i!=cols; i++)
    if(i!=moved_particle){
      for(int j=0; j!=rows; j++)
	tmp_pos(j) = fabs(pos(moved_particle,j)-pos(i,j));
      tmp_pos = distance(tmp_pos, l);
      double r_ij = 0;
      for(int j=0; j!=rows; j++){
	r_ij += sqr(tmp_pos(j));
      }
      if(r_ij > lh2)
	r_ij = lh2;

      sum += 1./(sqr(r_ij)*sqrt(r_ij));
    }

  return -.5*correlation5*sum;

}



//=================================================================
// exp(-sum_(i<j)(corr^5/r_ij^5))
//=================================================================
double LiqHe4WaveAll:: valuePt(QickArray &pos){

  int cols, rows, dim;
  dim = pos.get_dimension_info(&cols,&rows);
  QickArray tmp_pos(rows);
  double lh2 = sqr(l/2);
  double sum = 0;
  
  for(moved_particle=0; moved_particle!=cols; moved_particle++)
    for(int i=0; i!=cols; i++)
      if(i!=moved_particle){
	for(int j=0; j!=rows; j++)
	  tmp_pos(j) = fabs(pos(moved_particle,j)-pos(i,j));
	tmp_pos = distance(tmp_pos, l);
	double r_ij = 0;
	for(int j=0; j!=rows; j++){
	  r_ij += sqr(tmp_pos(j));
	}
	if(r_ij > lh2)
	  r_ij = lh2;
	
	sum += 1./(sqr(r_ij)*sqrt(r_ij));
      }
  
  std::cout << "sum: " << sum << std::endl;
  return -.5*correlation5*sum;

}


//=================================================================
// nabla psi / psi
//=================================================================
QickArray& LiqHe4QForce:: valuePt(QickArray &pos, bool dummy){

  int cols, rows, dim;
  dim = pos.get_dimension_info(&cols,&rows);
  QickArray tmp_pos(rows);
  double lh2 = sqr(l/2);
  if(!q_force.is_initialized())
    q_force.redim(cols,rows);
  q_force *= 0;
  for(int i=0; i!=cols; i++)
    for(int j=0; j!=cols; j++) if(i!=j){
      for(int k=0; k!=rows; k++)
	tmp_pos(k) = (pos(i,k)-pos(j,k));
      tmp_pos = distance(tmp_pos, l);
      double r_ij = 0;
      for(int k=0; k!=rows; k++){
	r_ij += sqr(tmp_pos(k));
      }
      if(r_ij > lh2)
	r_ij = lh2;
      tmp_pos /= sqr(r_ij)*r_ij*sqrt(r_ij);
      for(int k=0; k!=rows; k++)
	q_force(j,k) += tmp_pos(k);
    }

  q_force *= 5*correlation5;

  return q_force;

}

//=================================================================
// sum_(i<j) 1/r_ij^5
//=================================================================
double LiqHe4WfExp:: valuePt(QickArray& pos){

  int cols, rows, dim;
  dim = pos.get_dimension_info(&cols,&rows);
  QickArray tmp_pos(rows);
  double lh2 = sqr(l/2);
  
  double sum = 0;

  for(int i=0; i!=cols-1; i++)
    for(int j=i+1; j!=cols; j++){
      for(int k=0; k!=rows; k++)
	tmp_pos(k) = fabs(pos(i,k)-pos(j,k));
      tmp_pos = distance(tmp_pos, l);
      double r_ij = 0;
      for(int k=0; k!=rows; k++)
	r_ij += sqr(tmp_pos(k));
      if(r_ij > lh2)
	r_ij = lh2;
      sum += 1/(sqr(r_ij)*sqrt(r_ij));
    }

  return sum;
  
}

//=================================================================
// nabla psi / psi
//=================================================================
QickArray& LiqHe4QForceNum:: valuePt(QickArray &pos, bool dummy){

  int cols, rows, dim;
  dim = pos.get_dimension_info(&cols,&rows);
  if(!q_force.is_initialized())
    q_force.redim(cols,rows);
  Nabla nabla;
  LiqHe4WfExp wf_exp;

  QickArray params(2);

  params(0) = l;
  params(1) = correlation5;

  wf_exp.setParams(params);

  q_force = nabla(pos,wf_exp);

  q_force *= -correlation5;
  
  return q_force;

}

double LiqHe4LocalEnergy:: valuePt(QickArray &pos){

  int cols, rows, dum;
  dum = pos.get_dimension_info(&cols,&rows);
  QickArray tmp_pos(rows);

  const double kin_fac = 9.27647; //5*(hbar^2/m)/sigma^2
  const double four_epsilon = 40.88;

  double local_energy, kinetic_e, potential_e, z1, z2, z3;
  double dist, dist6, lh;
  local_energy = dist = 0;
  lh = sqr(l/2);
  kinetic_e = potential_e = 0;
  for(int particle_i=0; 
      particle_i!=cols-1; particle_i++){
    for(int particle_j=particle_i+1;
	particle_j!=cols; particle_j++){
      for(int dim=0; dim!=rows; dim++){
	tmp_pos(dim) = fabs(pos(particle_j,dim)-pos(particle_i,dim));
      }
      tmp_pos = distance(tmp_pos, l);
      dist = 0;
      for(int j=0; j!=rows; j++)
	dist += sqr(tmp_pos(j));
      if(dist<=lh) // cutoff correlation
	kinetic_e += kin_fac/(sqr(dist)*dist*sqrt(dist)); //kin/(dist^7)
      // printf("le %12.5E ", local_energy);
      dist = 0;

      // check image of x, y and z axis
      for(int i=-1; i!=1; i++){
	z1 = sqr(tmp_pos(0)+i*l);
	if(rows>1){ 
	  for(int j=-1; j!=1; j++){ 
	    z2 = sqr(tmp_pos(1)+j*l);
	    if(rows>2){ // 3 dimensions
	      for(int k=-1; k!=1; k++){
		z3 = sqr(tmp_pos(2)+k*l);
		dist = z1 + z2 + z3;
		//cerr << dist << endl;
		dist6=sqr(dist)*dist; //dist^6
		if(dist<=sqr(l)) // cutoff potential at L
		  potential_e += (1./dist6-1)/dist6;
	      }
	    }
	    else{ // 2 dimensions
	      dist = z1 + z2;
	      dist6=sqr(dist)*dist;
	      if(dist<=sqr(l)) // cutoff potential at L
		potential_e += (1./dist6-1)/dist6;
	    }
	  }
	}
	else{ // 1 dimension
	  dist6=sqr(z1)*z1;
	  if(dist<=sqr(l)) // cutoff potential at L
	    potential_e += (1./dist6-1)/dist6;
	}
      } // end for(i)
    } // end for(particle_j)
  } // end for(particle_i)

  kinetic_e *= correlation5;
  potential_e *= four_epsilon;
  //cerr << kinetic_e << " " << potential_e << endl;
  local_energy = kinetic_e+potential_e;
  return local_energy;  

}

//=================================================================
// V = -sum(number_particles/r_i) + sum(1/r_ij)
//=================================================================
double Potential:: valuePt(QickArray &pos){

  int cols, rows, dim;
  dim = pos.get_dimension_info(&cols,&rows);

  double potential_e=0;
  /* 
  ** adding potential energy; 
  ** V = -sum(number_particles/r_i) + sum(1/r_ij)
  */
  double r_single_particle;
  /* contribution from electron-proton potential  */
  for (int i = 0; i != cols; i++){
    r_single_particle = 0;
    for (int j = 0; j != rows; j++) {
      r_single_particle += sqr(pos(i,j));
    }
    potential_e -= cols/sqrt(r_single_particle);
  }
  
  /* contribution from electron-electron potential  */
  double r_ij, r_old_ik, r_old_jk;
  for (int i = 0; i != cols-1; i++) { 
    for (int j = i+1; j != cols; j++) {
      r_ij = 0;  
      for (int k = 0; k != rows; k++) { 
	r_old_ik = pos(i,k);
	r_old_jk = pos(j,k);
	r_ij += sqr(r_old_ik-r_old_jk);
      }

      potential_e += 1/sqrt(r_ij);
    }
  }

  return 2*potential_e;

}

//=================================================================
// g(alpha,lambda,r_i)f(a,r_ij)
//=================================================================
double MateuszWave:: valuePt(QickArray &pos){

  int cols, rows, dim;
  dim = pos.get_dimension_info(&cols,&rows);
  QickArray tmp_pos(rows);
  
  double f = 0;

  // f(a,r_ij)
  for(int i=0; i!=moved_particle; i++){
    for(int j=0; j!=rows; j++)
      tmp_pos(j) = pos(moved_particle,j)-pos(i,j);
    double r_ij = 0;
    for(int j=0; j!=rows; j++){
      r_ij += sqr(tmp_pos(j));
    }
    r_ij = sqrt(r_ij);
    f += log(1+a_1*exp(-sqr(r_ij/b_1))+a_2*exp(-sqr(r_ij/b_2)));
  }
  for(int i=moved_particle+1; i!=cols; i++){
    for(int j=0; j!=rows; j++)
      tmp_pos(j) = pos(moved_particle,j)-pos(i,j);
    double r_ij = 0;
    for(int j=0; j!=rows; j++){
      r_ij += sqr(tmp_pos(j));
    }
    r_ij = sqrt(r_ij);
    f += log(1+a_1*exp(-sqr(r_ij/b_1))+a_2*exp(-sqr(r_ij/b_2)));
  }
  
  double g=0;
  // g(alpha,lambda,r_ij)
  for(int j=0; j!=rows; j++)
    g+=sqr(pos(moved_particle,j));
  
    
  g = -alpha*g;

  return f+g;

}

//=================================================================
// g(alpha,lambda,r_i)f(a,r_ij)
//=================================================================
double MateuszWaveAll:: valuePt(QickArray &pos){

  int cols, rows, dim;
  dim = pos.get_dimension_info(&cols,&rows);
  QickArray tmp_pos(rows);

  double F = 0;

  // F(a,r_ij)
  for(int moved_particle=0; moved_particle!=cols-1; moved_particle++)
    for(int i=moved_particle+1; i!=cols; i++){
      for(int j=0; j!=rows; j++)
	tmp_pos(j) = pos(moved_particle,j)-pos(i,j);
      double r_ij = 0;
      for(int j=0; j!=rows; j++){
	r_ij += sqr(tmp_pos(j));
      }
      r_ij = sqrt(r_ij);
      F += log(1+a_1*exp(-sqr(r_ij/b_1))+a_2*exp(-sqr(r_ij/b_2)));
    }

  double g=0, G;
  for(int moved_particle=0; moved_particle!=cols; moved_particle++)
    // G(alpha,lambda,r_ij)
    for(int j=0; j!=rows; j++)
      g+=sqr(pos(moved_particle,j));
  
    
  G = -alpha*g;

  return F+G;

}


//=================================================================
// 2*nabla psi / psi !?!?!?!?!?!?!?!?!?!?!?!?!?!?!?!?!?!?!?!?!?!?!?!?
//=================================================================
QickArray& MateuszQForce:: valuePt(QickArray &pos, bool dummy){

  int cols, rows, dim;
  dim = pos.get_dimension_info(&cols,&rows);

  if(!pos_ij.is_initialized())
    pos_ij.redim(rows);

  if(!dF.is_initialized()){
    dF.redim(cols,rows);
    dG.redim(cols,rows);
  }

  dF *= 0;

  for(int i=0; i!=cols; i++){

    for(int x=0; x!=rows; x++)
      if(x==2)
	dG(i,x) = beta*pos(i,x);
      else
	dG(i,x) = pos(i,x);

  }

  dG *= -2*alpha;

  dF += dG;

  dF *= 2;

  return dF;

}

double MateuszLocalEnergy:: valuePt(QickArray &pos){

  int cols, rows, dim;
  dim = pos.get_dimension_info(&cols,&rows);

  QickArray pos_ij(rows);
  QickArray pos_i(rows);
  QickArray pos_ik(rows);

  double kinetic_e = 6*alpha*cols;

  double e_lg2=0, r_ij;

  for(int i=0; i!=cols; i++){

    for(int x=0; x!=rows; x++)
      if(x==2){
	pos_i(x) = lambda*pos(i,x);
	e_lg2   += sqr(pos_i(x));
      }
      else{
	pos_i(x) = pos(i,x);
	e_lg2   += sqr(pos_i(x));
      }

  }

  kinetic_e += -sqr(2*alpha)*e_lg2;

  double potential_e = 0;
  
  for(int i=0; i!=cols; i++)
    for(int j=0; j!=cols; j++) if(i!=j){
      r_ij = 0;
      for(int k=0; k!=rows; k++)
	r_ij += sqr(pos(i,k)-pos(j,k));
      r_ij = sqrt(r_ij);
      potential_e += (V_a*exp(-mu_a*r_ij)+V_r*exp(-mu_r*r_ij))/r_ij;

    }

  return (kinetic_e+potential_e);

}

double MateuszPotentialEnergy:: valuePt(QickArray &pos){

  int cols, rows, dim;
  dim = pos.get_dimension_info(&cols,&rows);

  double potential_e = 0, r_ij;

  for(int i=0; i!=cols; i++)
    for(int j=0; j!=cols; j++) if(i!=j){
      r_ij = 0;
      for(int k=0; k!=rows; k++)
	r_ij += sqr(pos(i,k)-pos(j,k));
      r_ij = sqrt(r_ij);
      potential_e += (V_a*exp(-mu_a*r_ij)+V_r*exp(-mu_r*r_ij))/r_ij;

    }


  return potential_e;

}

QickArray& Nabla:: valuePt(QickArray &pos, Func &func){

  const double h=1e-7;

  int cols,rows,dim;

  dim = pos.get_dimension_info(&cols,&rows);

  if(!dummy.is_initialized())
    dummy.redim(cols,rows);

  dummy *= 0;
  // using centered difference
  for(int i=0; i!=cols; i++)
    for(int j=0; j!=rows; j++){
      // f(x-h)
      pos(i,j) -= h;
      dummy(i,j) -= func(pos);
      // f(x+h)
      pos(i,j) += 2*h;
      dummy(i,j) += func(pos);
      // f(x)
      pos(i,j) -= h;
    }

  dummy /= 2*h;

  return dummy;

}


double Nabla2:: valuePt(Func &func, QickArray &pos){

  const double h=1e-4;

  int cols,rows,dim;

  dim = pos.get_dimension_info(&cols,&rows);

  double ans=0;

  for(int i=0; i!=cols; i++)
    for(int j=0; j!=rows; j++){
      // f(x-h)
      pos(i,j)-=h;
      ans += func(pos);
      // f(x+h)
      pos(i,j)+=2*h;
      ans += func(pos);
      // f(x)
      pos(i,j)-=h;
      ans -= 2*func(pos);
    }
  
  return ans/h/h;

}
