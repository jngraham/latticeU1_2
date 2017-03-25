
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

#include "globals.h"

#include "link_update.h"
#include "operators.h"
#include "write.h"

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
  // and in the case of the flux tube, by real and imaginary part

  double jpc_plus_data [Lt*N_samples] = {0};
  double jpc_minus_data [Lt*N_samples] = {0};
  double flux_data [Lt*N_samples*2] = {0};

  // declare the zero-time operators

  double this_avg_plaquette;
  double jpc_plus_zero;
  double jpc_minus_zero;
  double flux_zero [2];
  double flux_pair [2];

  // update the field configuration so we "forget" the initial configuration

  std::cout << "initial average plaquette (should be 1): " << avg_plaquette(lattice) << "\n";

  for (size_t i = 0; i < N_equilibration_configs; i++){
    update(lattice, V);
  }

  // do the "science run"

  for (size_t i = 0; i < N_samples; i++){

    // update the lattice then compute <\phi^\dagger(nt)\phi(0)> for all
    // our operators

    for (size_t j = 0; j < N_configs_per_sample; j++){

      update(lattice, V);

      // add the average plaquette for this configuration

      this_avg_plaquette = avg_plaquette(lattice);

      avg_plaquette_data[i] += this_avg_plaquette;

      // find the field operators for nt = 0
      // the jpc operators are real-valued

      jpc_plus_zero = jpc_plus(lattice, this_avg_plaquette, 0);
      jpc_minus_zero = jpc_minus(lattice, 0);
      flux(lattice, 0, &flux_zero[0], &flux_zero[1]);

      // add the phi(0)*phi(0) to the data since we already know what that is

      jpc_plus_data[Lt*i] += pow(jpc_plus_zero,2);
      jpc_minus_data[Lt*i] += pow(jpc_minus_zero,2);
      flux_data[2*Lt*i] += 1;

      // this line is not strictly necessary bc we initialize it to zero.
      // flux_data[2*Lt*i + 1] += 0;

      // find the phi(t)*phi(0) for all other t in each configuration

      for (size_t t = 1; t < Lt; t++){
        jpc_plus_data[Lt*i + t] += jpc_plus(lattice, this_avg_plaquette, t)*jpc_plus_zero;
        jpc_minus_data[Lt*i + t] += jpc_minus(lattice, t)*jpc_minus_zero;

        // the flux operator returns a complex exponential so we need to take
        // care of the real and imaginary pieces

        flux(lattice, t, &flux_pair[0], &flux_pair[1]);

        // real part
        flux_data[2*Lt*i + 2*t] += flux_pair[0]*flux_zero[0] + flux_pair[1]*flux_pair[1];

        // imaginary part
        flux_data[2*Lt*i + 2*t + 1] += flux_pair[0]*flux_zero[1] - flux_pair[1]*flux_zero[0];
      }

    }

    // thus far we have sums but
    avg_plaquette_data[i] = avg_plaquette_data[i]/N_configs_per_sample;

    for (size_t t = 0; t < Lt; t++){
      jpc_plus_data[Lt*i + t] = jpc_plus_data[Lt*i + t]/N_configs_per_sample;
      jpc_minus_data[Lt*i + t] = jpc_minus_data[Lt*i + t]/N_configs_per_sample;

      flux_data[2*Lt*i + 2*t] = flux_data[2*Lt*i + 2*t]/N_configs_per_sample;
      flux_data[2*Lt*i + 2*t + 1] = flux_data[2*Lt*i + 2*t + 1]/N_configs_per_sample;
    }

    std::printf ("%i of %i samples complete\n", int(i+1), N_samples);

  }

  clock_t t2 = clock();

  write(avg_plaquette_data, jpc_plus_data, jpc_minus_data, flux_data);

  clock_t t3 = clock();

  // output some data for my writing-code purposes

  double sum = 0;

  for(size_t i = 0; i < N_samples; i++){
    sum += avg_plaquette_data[i];
  }

  // std::ofstream output;
  // output.open("outputs.txt");

  // output << "<cos U_p> = " << sum / N_samples << std::endl;
  std::cout << "<cos U_p> = " << sum / N_samples << std::endl;

  // output << "time 1: " << (double(t2) - double(t1)) / CLOCKS_PER_SEC << std::endl;
  std::cout << "time 1: " << (double(t2) - double(t1)) / CLOCKS_PER_SEC << "\n";
  // output << "time 2: " << (double(t3) - double(t2)) / CLOCKS_PER_SEC << std::endl;
  std::cout << "time 2: " << (double(t3) - double(t2)) / CLOCKS_PER_SEC << "\n";

  // output.close();
}

/*
///////////////////////////////////////////////////////////////////////////////

fin.

///////////////////////////////////////////////////////////////////////////////
*/
