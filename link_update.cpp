
/*
/////////////////////////////////////////////////////////////////////////////

U(1) Lattice Gauge Theory Simulator | James Graham

Here I write the function(s) to update the lattice

/////////////////////////////////////////////////////////////////////////////
*/

#include <iostream>
#include <random>
#include <cmath>
#include <stdlib.h>

#include "link_update.h"

/*
This function is kinda gnarly and dirty because I will use it to actually change
the field configuration, rather than calling it and having it return a new link
that we put into the field configuration in lattice.cpp. I mean, if we did that,
we'd need to pass in coordinates and stuff and it would be more complicated and
proper and so on and so on

*sniffle*
*/

std::default_random_engine update_generator;
std::uniform_int_distribution<int> int_distribution(0,199);
std::uniform_real_distribution<double> real_distribution(0.0,1.0);

int xnext;
int xprev;
int ynext;
int yprev;
int tnext;
int tprev;

double old_link;
double new_link;

double staple1;
double staple2;
double staple3;
double staple4;

double old_action;
double new_action;

double C;
double z;

int update_corners(double* array, double* gauss, double beta, int Lx, int Ly, int Lt);
int update_edges(double* array, double* gauss, double beta, int Lx, int Ly, int Lt);
int update_bulk(double* array, double* gauss, double beta, int Lx, int Ly, int Lt);

int update(double* array, double* gauss, double beta, int Lx, int Ly, int Lt){

  // for the corners of the lattice (8 of them)

  // update_corners(array, gauss, beta, Lx, Ly, Lt);

  // for the edges of the lattice (12 of them)

  // update_edges(array, gauss, beta, Lx, Ly, Lt);

  // for the bulk of the lattice

  update_bulk(array, gauss, beta, Lx, Ly, Lt);

  return 0;
}

int update_corners(double* array, double* gauss, double beta, int Lx, int Ly, int Lt){
  return 0;
}

int update_edges(double* array, double* gauss, double beta, int Lx, int Ly, int Lt){
  return 0;
}

int update_bulk(double* array, double* gauss, double beta, int Lx, int Ly, int Lt){

  int N_links = 3*Lx*Ly*Lt;
  // int N_acceptances = 0;

  for (int x = 0; x < Lx; x++){

    // get x neighbour coordinates

    xnext = (x+1)%Lx;
    xprev = (x+Lx-1)%Lx;

    for (int y = 0; y < Ly; y++){

      // get y neighbour coordinates

      ynext = (y+1)%Ly;
      yprev = (y+Ly-1)%Ly;

      for (int t = 0; t < Lt; t++){

        // get t neighbour coordinates

        tnext = (t+1)%Lt;
        tprev = (t+Lt-1)%Lt;

        // update the x link

        old_link = array[3*x + 3*Lx*y + 3*Lx*Ly*t + 0];
        new_link = old_link + gauss[int_distribution(update_generator)];
        // new_link = old_link + gauss[int_distribution(generator)];

        staple1 = array[3*xnext + 3*Lx*y + 3*Lx*Ly*t + 1] - array[3*x + 3*Lx*ynext + 3*Lx*Ly*t + 0] - array[3*x + 3*Lx*y + 3*Lx*Ly*t + 1];
        staple2 = - array[3*xnext + 3*Lx*yprev + 3*Lx*Ly*t + 1] - array[3*x + 3*Lx*yprev + 3*Lx*Ly*t + 0] + array[3*x + 3*Lx*yprev + 3*Lx*Ly*t + 1];
        staple3 = array[3*xnext + 3*Lx*y + 3*Lx*Ly*t + 2] - array[3*x + 3*Lx*y + 3*Lx*Ly*tnext + 0] - array[3*x + 3*Lx*y + 3*Lx*Ly*t + 2];
        staple4 = - array[3*xnext + 3*Lx*y + 3*Lx*Ly*tprev + 2] - array[3*x + 3*Lx*y + 3*Lx*Ly*tprev + 0] + array[3*x + 3*Lx*y + 3*Lx*Ly*tprev + 2];

        old_action = cos(old_link + staple1) + cos(old_link + staple2) + cos(old_link + staple3) + cos(old_link + staple4);
        new_action = cos(new_link + staple1) + cos(new_link + staple2) + cos(new_link + staple3) + cos(new_link + staple4);

        C = exp(- beta * new_action)/exp(- beta * old_action);

        z = real_distribution(update_generator);
        // z = real_distribution(generator);

        if (z < C){
          // N_acceptances++;
          array[3*x + 3*Lx*y + 3*Lx*Ly*t + 0] = new_link;
        }

        // update the y link

        old_link = array[3*x + 3*Lx*y + 3*Lx*Ly*t + 1];
        new_link = old_link + gauss[int_distribution(update_generator)];
        // new_link = old_link + gauss[int_distribution(generator)];

        staple1 = array[3*x + 3*Lx*ynext + 3*Lx*Ly*t + 0] - array[3*xnext + 3*Lx*y + 3*Lx*Ly*t + 1] - array[3*x + 3*Lx*y + 3*Lx*Ly*t + 0];
        staple2 = - array[3*xprev + 3*Lx*ynext + 3*Lx*Ly*t + 0] - array[3*xprev + 3*Lx*y + 3*Lx*Ly*t + 1] + array[3*xprev + 3*Lx*y + 3*Lx*Ly*t + 0];
        staple3 = array[3*x + 3*Lx*ynext + 3*Lx*Ly*t + 2] - array[3*x + 3*Lx*y + 3*Lx*Ly*tnext + 1] - array[3*x + 3*Lx*y + 3*Lx*Ly*t + 2];
        staple4 = - array[3*x + 3*Lx*ynext + 3*Lx*Ly*tprev + 2] - array[3*x + 3*Lx*y + 3*Lx*Ly*tprev + 0] + array[3*x + 3*Lx*y + 3*Lx*Ly*t + 2];

        old_action = cos(old_link + staple1) + cos(old_link + staple2) + cos(old_link + staple3) + cos(old_link + staple4);
        new_action = cos(new_link + staple1) + cos(new_link + staple2) + cos(new_link + staple3) + cos(new_link + staple4);

        C = exp(- beta * new_action)/exp(- beta * old_action);

        z = real_distribution(update_generator);
        // z = real_distribution(generator);

        if (z < C){
          // N_acceptances++;
          array[3*x + 3*Lx*y + 3*Lx*Ly*t + 1] = new_link;
        }

        // update the t link

        old_link = array[3*x + 3*Lx*y + 3*Lx*Ly*t + 2];
        new_link = old_link + gauss[int_distribution(update_generator)];
        // new_link = old_link + gauss[int_distribution(generator)];

        staple1 = array[3*x + 3*Lx*y + 3*Lx*Ly*tnext + 0] - array[3*xnext + 3*Lx*y + 3*Lx*Ly*t + 2] - array[3*x + 3*Lx*y + 3*Lx*Ly*t + 0];
        staple2 = - array[3*xprev + 3*Lx*y + 3*Lx*Ly*tnext + 0] - array[3*xprev + 3*Lx*y + 3*Lx*Ly*t + 2] + array[3*xprev + 3*Lx*y + 3*Lx*Ly*t + 0];
        staple3 = array[3*x + 3*Lx*y + 3*Lx*Ly*tnext + 1] - array[3*x + 3*Lx*ynext + 3*Lx*Ly*t + 2] - array[3*x + 3*Lx*y + 3*Lx*Ly*t + 1];
        staple4 = - array[3*x + 3*Lx*yprev + 3*Lx*Ly*tnext + 1] - array[3*x + 3*Lx*yprev + 3*Lx*Ly*t + 2] + array[3*x + 3*Lx*yprev + 3*Lx*Ly*t + 1];

        old_action = cos(old_link + staple1) + cos(old_link + staple2) + cos(old_link + staple3) + cos(old_link + staple4);
        new_action = cos(new_link + staple1) + cos(new_link + staple2) + cos(new_link + staple3) + cos(new_link + staple4);

        C = exp(- beta * new_action)/exp(- beta * old_action);

        z = real_distribution(update_generator);
        // z = real_distribution(generator);

        if (z < C){
          // N_acceptances++;
          array[3*x + 3*Lx*y + 3*Lx*Ly*t + 2] = new_link;
        }

      }
    }
  }

  // std::cout << "acceptance rate: " << double(N_acceptances) / double(N_links) << std::endl;

  return 0;
}
