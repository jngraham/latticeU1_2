
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

  // The coarsest graining will be in time
  // Then in y
  // Then in x
  // And finally we put all the links next to each other for each x.

  float lattice [N_links] = {0};

  return 0;
}









/*
///////////////////////////////////////////////////////////////////////////////

fin.

///////////////////////////////////////////////////////////////////////////////
*/
