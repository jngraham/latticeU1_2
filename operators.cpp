
/*
/////////////////////////////////////////////////////////////////////////////

U(1) Lattice Gauge Theory Simulator | James Graham

Here I write the functions for the field operators on the lattice

/////////////////////////////////////////////////////////////////////////////
*/

#include <cmath>

#include "operators.h"
#include "globals.h"

double avg_plaquette(double* array){

  double plaquette_sum = 0;

  int N_plaquettes = 3*Lx*Ly*Lt;

  double p_xy = 0;
  double p_xt = 0;
  double p_yt = 0;

  for (int x = 0; x < Lx; x++){
    for (int y = 0; y < Ly; y++){
      for (int t = 0; t < Lt; t++){

        p_xy = array[3*x + 3*Lx*y + 3*Lx*Ly*t + 0] + array[3*((x+1)%Lx) + 3*Lx*y + 3*Lx*Ly*t + 1] - array[3*x + 3*Lx*((y+1)%Ly) + 3*Lx*Ly*t + 0] - array[3*x + 3*Lx*y + 3*Lx*Ly*t + 1];
        p_xt = array[3*x + 3*Lx*y + 3*Lx*Ly*t + 0] + array[3*((x+1)%Lx) + 3*Lx*y + 3*Lx*Ly*t + 2] - array[3*x + 3*Lx*y + 3*Lx*Ly*((t+1)%Lt) + 0] - array[3*x + 3*Lx*y + 3*Lx*Ly*t + 2];
        p_yt = array[3*x + 3*Lx*y + 3*Lx*Ly*t + 1] + array[3*x + 3*Lx*((y+1)%Ly) + 3*Lx*Ly*t + 2] - array[3*x + 3*Lx*y + 3*Lx*Ly*((t+1)%Lt) + 1] - array[3*x + 3*Lx*y + 3*Lx*Ly*t + 2];

        plaquette_sum = plaquette_sum + cos(p_xy) + cos(p_xt) + cos(p_yt);
      }
    }
  }

  return plaquette_sum / N_plaquettes;
}

double jpc_plus(double* array, int t){
  return 0;
}

double jpc_minus(double* array, int t){
  return 0;
}

double flux(double* array, int t){
  return 0;
}
