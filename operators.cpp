
/*
/////////////////////////////////////////////////////////////////////////////

U(1) Lattice Gauge Theory Simulator | James Graham

Here I write the functions for the field operators on the lattice

/////////////////////////////////////////////////////////////////////////////
*/

#include <cmath>

#include "operators.h"

float avg_plaquette(float* array, int Lx, int Ly, int Lt){

  float plaquette_sum = 0;

  int N_plaquettes = 3*Lx*Ly*Lt;

  float p_xy = 0;
  float p_xt = 0;
  float p_yt = 0;

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

float jpc_plus(float* array, int Lx, int Ly, int t){
  return 0;
}

float jpc_minus(float* array, int Lx, int Ly, int t){
  return 0;
}

float flux(float* array, int Lx, int Ly, int t){
  return 0;
}
