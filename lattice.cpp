
/*
/////////////////////////////////////////////////////////////////////////////

U(1) Lattice Gauge Theory Simulator | James Graham

The procedure here is much the same as it is in my first attempt at the problem
using python. The lattice is structured differently this time, which will make
my life harder but the computer's life easier.

/////////////////////////////////////////////////////////////////////////////
*/

#include <iostream>
#include <random>
#include <stdlib.h>
#include <typeinfo>

// declare update
int update(float* lattice, int Lx, int Ly, int Lt);

// declare field operators
float avg_plaquette(float* lattice, int Lx, int Ly, int Lt);
float jpc_plus(float* lattice, int Lx, int Ly, int t);
float jpc_minus(float* lattice, int Lx, int Ly, int t);
float flux(float* lattice, int Lx, int Ly, int t);

int main(){

  // Parameters for the size of the array

  const int Lx = 10;
  const int Ly = 10;
  const int Lt = 10;

  // Parameters for the simulation

  const float beta = 2.2;

  const int N_equilibration_configs = 20;
  const int N_configs_per_sample = 50;
  const int N_samples = 10;

  const int N_configs = N_configs_per_sample * N_samples;

  const int N_links = Lx * Ly * Lt * 3;

  // Parameters for selecting a new link

  const int N_V = 200;
  const float mu = 0;
  const float sigma = 0.5;

  // set up the RNG
  // seed??

  std::default_random_engine generator;
  std::normal_distribution<float> distribution(mu,sigma);

  // set up our V array

  float V [N_V];

  for (size_t i = 0; i < N_V; i += 2){
    float number = distribution(generator);
    V[i] = number;
    V[i+1] = -number;
  }

  // set up the lattice

  /*
  The coarsest graining will be in time
  Then in y
  Then in x
  And finally we put all the links next to each other for each x.
  */

  float lattice[N_links] = {0};
  
}

/*
This function is kinda gnarly and dirty because I will use it to actually change
the field configuration, rather than calling it and having it return a new link
that we put into the field configuration in lattice.cpp. I mean, if we did that,
we'd need to pass in coordinates and stuff and it would be more complicated and
proper and so on and so on

*sniffle*
*/

int update(float* lattice, int Lx, int Ly, int Lt){
  return 0;
}

/*
Here I have all the field operators
*/

float avg_plaquette(float* lattice, int Lx, int Ly, int Lt){
  return 0;
}

float jpc_plus(float* lattice, int Lx, int Ly, int t){
  return 0;
}

float jpc_minus(float* lattice, int Lx, int Ly, int t){
  return 0;
}

float flux(float* lattice, int Lx, int Ly, int t){
  return 0;
}






/*
///////////////////////////////////////////////////////////////////////////////

fin.

///////////////////////////////////////////////////////////////////////////////
*/
