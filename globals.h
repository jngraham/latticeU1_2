
/*
/////////////////////////////////////////////////////////////////////////////

U(1) Lattice Gauge Theory Simulator | James Graham

Global parameters and the V set

/////////////////////////////////////////////////////////////////////////////
*/

#include <random>
// #include <stdlib.h>
// #include <stdio.h>
// #include <time.h>

// Parameters for the size of the array

const int Lx = 10;
const int Ly = 10;
const int Lt = 10;

// Parameters for the number of samples, field configurations, etc.

const double beta = 2.2;

const int N_equilibration_configs = 2000;
const int N_configs_per_sample = 500;
const int N_samples = 10;

const int N_configs = N_configs_per_sample * N_samples;

const int N_links = Lx * Ly * Lt * 3;

// Parameters for V

const int N_V = 200;
const double mu = 0;
const double sigma = 1.2;

// Set up the RNG

// std::mt19937 generator(time(0));
extern std::default_random_engine generator;
