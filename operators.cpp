
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

  int xnext;
  int ynext;
  int tnext;

  double plaquette_sum = 0;

  int N_plaquettes = 3*Lx*Ly*Lt;

  double p_xy = 0;
  double p_xt = 0;
  double p_yt = 0;

  for (int x = 0; x < Lx; x++){

    xnext = (x+1)%Lx;

    for (int y = 0; y < Ly; y++){

      ynext = (y+1)%Ly;

      for (int t = 0; t < Lt; t++){

        tnext = (t+1)%Lt;

        p_xy = array[3*x + 3*Lx*y + 3*Lx*Ly*t + 0] + array[3*xnext + 3*Lx*y + 3*Lx*Ly*t + 1] - array[3*x + 3*Lx*ynext + 3*Lx*Ly*t + 0] - array[3*x + 3*Lx*y + 3*Lx*Ly*t + 1];
        p_xt = array[3*x + 3*Lx*y + 3*Lx*Ly*t + 0] + array[3*xnext + 3*Lx*y + 3*Lx*Ly*t + 2] - array[3*x + 3*Lx*y + 3*Lx*Ly*tnext + 0] - array[3*x + 3*Lx*y + 3*Lx*Ly*t + 2];
        p_yt = array[3*x + 3*Lx*y + 3*Lx*Ly*t + 1] + array[3*x + 3*Lx*ynext + 3*Lx*Ly*t + 2] - array[3*x + 3*Lx*y + 3*Lx*Ly*tnext + 1] - array[3*x + 3*Lx*y + 3*Lx*Ly*t + 2];

        plaquette_sum = plaquette_sum + cos(p_xy) + cos(p_xt) + cos(p_yt);
      }
    }
  }

  return plaquette_sum / N_plaquettes;
}

double jpc_plus(double* array, double VEV, int t){
  return 0;
}

double jpc_minus(double* array, int t){

  int xnext;
  int ynext;

  double sum = 0;

  double p_xy = 0;

  for (int x = 0; x < Lx; x++){

    xnext = (x+1)%Lx;

    for(int y = 0; y < Ly; y++){

      ynext = (y+1)%Ly;

      p_xy = array[3*x + 3*Lx*y + 3*Lx*Ly*t + 0] + array[3*xnext + 3*Lx*y + 3*Lx*Ly*t + 1] - array[3*x + 3*Lx*ynext + 3*Lx*Ly*t + 0] - array[3*x + 3*Lx*y + 3*Lx*Ly*t + 1];

      sum += sin(p_xy);
    }
  }

  return sum;
}

double flux(double* array, int t, double* re_ptr, double* im_ptr){

  double phase_sum = 0;
  double re_sum = 0;
  double im_sum = 0;

  for (int y = 0; y < Ly; y++){
    for (int x = 0; x < Lx; x++){
      phase_sum += array[3*x + 3*Lx*y + 3*Lx*Ly*t + 0];
    }

    re_sum += cos(phase_sum);
    im_sum += sin(phase_sum);

    phase_sum = 0;
  }

  *re_ptr = re_sum;
  *im_ptr = im_sum;

  return 0;
}
