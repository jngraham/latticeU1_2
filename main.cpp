
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
#include "globals.h"

clock_t t1 = clock();

int main(){

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

  double avg_plaquette_data [N_samples] = {0};

  // coarsest graining will be by sample, then by time

  double jpc_plus_data [Lt*N_samples] = {0};
  double jpc_minus_data [Lt*N_samples] = {0};
  double flux_data [Lt*N_samples] = {0};

  // declare the zero-time operators

  double this_avg_plaquette;
  double jpc_plus_zero;
  double jpc_minus_zero;
  double flux_zero;

  // update the field configuration so we "forget" the initial configuration

  std::cout << "initial average plaquette (should be 1): " << avg_plaquette(lattice) << "\n";

  for (size_t i = 0; i < N_equilibration_configs; i++){
    update(lattice, V);
  }

  // do the "science run"

  for (size_t i = 0; i < N_samples; i++){
    for (size_t j = 0; j < N_configs_per_sample; j++){

      update(lattice, V);

      // add the average plaquette for this configuration

      this_avg_plaquette = avg_plaquette(lattice);

      avg_plaquette_data[i] += this_avg_plaquette;

      // find the field operators for nt = 0

      jpc_plus_zero = jpc_plus(lattice, this_avg_plaquette, 0);
      jpc_minus_zero = jpc_minus(lattice, 0);
      flux_zero = flux(lattice, 0);

      // add the phi(0)*phi(0) to the data since we already know what that is

      jpc_plus_data[N_samples*i] += pow(jpc_plus_zero,2);
      jpc_minus_data[N_samples*i] += pow(jpc_minus_zero,2);
      flux_data[N_samples*i] += pow(flux_zero,2);

      // find the phi(t)*phi(0) for all other t in each configuration

      for (size_t t = 1; t < Lt; t++){
        jpc_plus_data[N_samples*i + t] += jpc_plus(lattice, this_avg_plaquette, t)*jpc_plus_zero;
        jpc_minus_data[N_samples*i + t] += jpc_minus(lattice, t)*jpc_minus_zero;
        flux_data[N_samples*i + t] += flux(lattice, t)*flux_zero;
      }

    }

    avg_plaquette_data[i] = avg_plaquette_data[i]/N_configs_per_sample;

    for (size_t t = 0; t < Lt; t++){
      jpc_plus_data[N_samples*i + t] += jpc_plus_data[N_samples*i + t]/N_configs_per_sample;
      jpc_minus_data[N_samples*i + t] += jpc_minus_data[N_samples*i + t]/N_configs_per_sample;
      flux_data[N_samples*i + t] += flux_data[N_samples*i + t]/N_configs_per_sample;
    }

    std::printf ("%i of %i samples complete\n", int(i+1), N_samples);

  }

  clock_t t2 = clock();

  // Output data to .csv so I can use excel or mathematica or python or whatever

  std::ofstream plaquette_output;
  std::ofstream mplus_output;
  std::ofstream mminus_output;
  std::ofstream flux_output;

  std::string params = "_" + std::to_string(Lx) + "x" + std::to_string(Ly) + "x" + std::to_string(Lt) + "_beta" + std::to_string(int(10*(beta-2)));
  plaquette_output.open("plaquette"+params+".csv");
  mplus_output.open("mplus"+params+".csv");
  mminus_output.open("mminus"+params+".csv");
  flux_output.open("flux"+params+".csv");

  plaquette_output << "Sample #,value,\n";

  mplus_output << "Sample #,";
  mminus_output << "Sample #,";
  flux_output << "Sample #,";

  for (size_t t = 0; t < Lt; t++){
    mplus_output << "time " << t << ",";
    mminus_output << "time " << t << ",";
    flux_output << "time " << t << ",";
  }

  mplus_output << "\n";
  mminus_output << "\n";
  flux_output << "\n";

  for (size_t i = 0; i < N_samples; i++){

    plaquette_output << i << "," << avg_plaquette_data[i] << ",\n";

    mplus_output << i << ",";
    mminus_output << i << ",";
    flux_output << i << ",";

    for (size_t t = 0; t < Lt; t++){
      mplus_output << jpc_plus_data[N_samples*i + t] << ",";
      mminus_output << jpc_minus_data[N_samples*i + t] << ",";
      flux_output << flux_data[N_samples*i + t] << ",";
    }

    mplus_output << "\n";
    mminus_output << "\n";
    flux_output << "\n";
  }

  plaquette_output.close();
  mplus_output.close();
  mminus_output.close();
  flux_output.close();

  clock_t t3 = clock();

  // output some data for my writing-code purposes

  double sum = 0;

  for(size_t i = 0; i < N_samples; i++){
    sum += avg_plaquette_data[i];
  }

  std::ofstream output;
  output.open("outputs.txt");

  output << "<cos U_p> = " << sum / N_samples << std::endl;
  std::cout << "<cos U_p> = " << sum / N_configs << std::endl;

  output << "time 1: " << (double(t2) - double(t1)) / CLOCKS_PER_SEC << std::endl;
  std::cout << "time 1: " << (double(t2) - double(t1)) / CLOCKS_PER_SEC << "\n";
  output << "time 2: " << (double(t3) - double(t2)) / CLOCKS_PER_SEC << std::endl;
  std::cout << "time 2: " << (double(t3) - double(t2)) / CLOCKS_PER_SEC << "\n";

  output.close();
}

/*
///////////////////////////////////////////////////////////////////////////////

fin.

///////////////////////////////////////////////////////////////////////////////
*/
