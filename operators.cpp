
/*
/////////////////////////////////////////////////////////////////////////////

U(1) Lattice Gauge Theory Simulator | James Graham

Here I write the functions for the field operators on the lattice

/////////////////////////////////////////////////////////////////////////////
*/

#include <iostream>
#include <cmath>

#include "operators.h"
#include "globals.h"

double avg_plaquette(double* array){

  int xthis;
  int xnext;

  int ythis;
  int ynext;

  int tthis;
  int tnext;

  double plaquette_sum = 0;

  int N_plaquettes = 3*Lx*Ly*Lt;

  double p_xy = 0;
  double p_xt = 0;
  double p_yt = 0;

  for (int x = 0; x < Lx; x++){

    xthis = 3*x;
    xnext = 3*((x+1)%Lx);

    for (int y = 0; y < Ly; y++){

      ythis = 3*Lx*y;
      ynext = 3*Lx*((y+1)%Ly);

      for (int t = 0; t < Lt; t++){

        tthis = 3*Lx*Ly*t;
        tnext = 3*Lx*Ly*((t+1)%Lt);

        p_xy = array[xthis + ythis + tthis] + array[xnext + ythis + tthis + 1] - array[xthis + ynext + tthis] - array[xthis + ythis + tthis + 1];
        p_xt = array[xthis + ythis + tthis] + array[xnext + ythis + tthis + 2] - array[xthis + ythis + tnext] - array[xthis + ythis + tthis + 2];
        p_yt = array[xthis + ythis + tthis + 1] + array[xthis + ynext + tthis + 2] - array[xthis + ythis + tnext + 1] - array[xthis + ythis + tthis + 2];

        plaquette_sum = plaquette_sum + cos(p_xy) + cos(p_xt) + cos(p_yt);
      }
    }
  }

  return plaquette_sum / N_plaquettes;
}

double jpc_plus(double* array, double VEV, int t){

  int xthis;
  int xnext;

  int ythis;
  int ynext;

  int tthis = 3*Lx*Ly*t;

  double sum = 0;

  double p_xy = 0;

  for (int x = 0; x < Lx; x++){

    xthis = 3*x;
    xnext = 3*((x+1)%Lx);

    for(int y = 0; y < Ly; y++){

      ythis = 3*y;
      ynext = 3*((y+1)%Ly);

      p_xy = array[xthis + ythis + tthis] + array[xnext + ythis + tthis + 1] - array[xthis + ynext + tthis] - array[xthis + ythis + tthis + 1];

      sum += cos(p_xy) - VEV;
    }
  }

  return sum;
}

double jpc_minus(double* array, int t){

  int xthis;
  int xnext;

  int ythis;
  int ynext;

  int tthis = 3*Lx*Ly*t;

  double sum = 0;

  double p_xy = 0;

  for (int x = 0; x < Lx; x++){

    xthis = 3*x;
    xnext = 3*((x+1)%Lx);

    for(int y = 0; y < Ly; y++){

      ythis = 3*y;
      ynext = 3*((y+1)%Ly);

      p_xy = array[xthis + ythis + tthis] + array[xnext + ythis + tthis + 1] - array[xthis + ynext + tthis] - array[xthis + ythis + tthis + 1];

      sum += sin(p_xy);
    }
  }

  return sum;
}

double flux(double* array, int t, double* re_ptr, double* im_ptr){

  double phase_sum = 0;
  double re_sum = 0;
  double im_sum = 0;

  int ythis;
  int tthis = 3*Lx*Ly*t;

  for (int y = 0; y < Ly; y++){

    ythis = 3*Lx*y;

    for (int x = 0; x < Lx; x++){

      // less efficient to declare an int xthis and then define it, since we only use it once per loop.
      phase_sum += array[3*x + ythis + tthis];
    }

    re_sum += cos(phase_sum);
    im_sum += sin(phase_sum);

    phase_sum = 0;
  }

  *re_ptr = re_sum;
  *im_ptr = im_sum;

  return 0;
}
