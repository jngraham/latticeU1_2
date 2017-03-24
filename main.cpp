
/*
/////////////////////////////////////////////////////////////////////////////

U(1) Lattice Gauge Theory Simulator | James Graham

The procedure here is much the same as it is in my first attempt at the problem
using python. The lattice is structured differently this time, which will make
my life harder but the computer's life easier.

/////////////////////////////////////////////////////////////////////////////
*/

#include <iostream>
#include <fstream>

#include <random>

#include <typeinfo>
#include <cmath>

#include <stdlib.h>
#include <stdio.h>
#include <time.h>

#include "link_update.h"
#include "operators.h"

clock_t t1 = clock();

// Parameters for the size of the array

const int Lx = 10;
const int Ly = 10;
const int Lt = 10;

// Parameters for the simulation

const double beta = 2.2;

const int N_equilibration_configs = 2000;
const int N_configs_per_sample = 500;
const int N_samples = 10;

const int N_configs = N_configs_per_sample * N_samples;

const int N_links = Lx * Ly * Lt * 3;

int main(){

  // Parameters for selecting a new link

  const int N_V = 200;
  const double mu = 0;
  const double sigma = 1.2;

  // set up the RNG
  // seed??

  // std::mt19937 generator(time(0));
  std::default_random_engine generator;
  std::normal_distribution<double> gaussian_distribution(mu,sigma);

  // set up our V array

  double V [N_V];

  for (size_t i = 0; i < N_V; i += 2){
    double number = gaussian_distribution(generator);
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

  double lattice [N_links] = {0};

  // set up our data arrays

  double avg_plaquette_data [N_configs] = {0};
  double jpc_plus_data [Lt*N_configs] = {0};
  double jpc_minus_data [Lt*N_configs] = {0};
  double flux_data [Lt*N_configs] = {0};

  // update the field configuration so we "forget" the initial configuration

  for (size_t i = 0; i < N_equilibration_configs; i++){
    update(lattice, V, beta, Lx, Ly, Lt);
  }

  // for (int x = 0; x < Lx; x++){
  //   for (int y = 0; y < Ly; y++){
  //     for (int t = 0; t < Lt; t++){
  //
  //       std::cout << "x: " << x << ", y:" << y << ", t:" << t << ", x link value:" << lattice[3*x + 3*Lx*y + 3*Lx*Ly*t + 0] << ", y link value:" << lattice[3*x + 3*Lx*y + 3*Lx*Ly*t + 1] << ", t link value:" << lattice[3*x + 3*Lx*y + 3*Lx*Ly*t + 2] << "\n";
  //     }
  //   }
  // }

  // do the "science run"

  std::ofstream output;
  output.open("outputs.txt");

  for(size_t i = 0; i < N_configs; i++){
    update(lattice, V, beta, Lx, Ly, Lt);

    // calculate data

    avg_plaquette_data[i] = avg_plaquette(lattice, Lx, Ly, Lt);

  }

  clock_t t2 = clock();

  // return or otherwise output data

  double sum = 0;

  for(size_t i = 0; i < N_configs; i++){
    sum += avg_plaquette_data[i];
  }

  clock_t t3 = clock();

  output << "<cos U_p> = " << sum / N_configs << std::endl;
  // std::cout << "<cos U_p> = " << sum / N_configs << std::endl;

  output << "time 1: " << (double(t2) - double(t1)) / CLOCKS_PER_SEC << std::endl;
  // std::cout << "time 1: " << (double(t2) - double(t1)) / CLOCKS_PER_SEC << "\n";
  output << "time 2: " << (double(t3) - double(t2)) / CLOCKS_PER_SEC << std::endl;
  // std::cout << "time 2: " << (double(t3) - double(t2)) / CLOCKS_PER_SEC << "\n";

  output.close();
}

/*
///////////////////////////////////////////////////////////////////////////////

fin.

///////////////////////////////////////////////////////////////////////////////
*/
