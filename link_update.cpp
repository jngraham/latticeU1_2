
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
#include "globals.h"

/*
This function is kinda gnarly and dirty because I will use it to actually change
the field configuration, rather than calling it and having it return a new link
that we put into the field configuration in lattice.cpp. I mean, if we did that,
we'd need to pass in coordinates and stuff and it would be more complicated and
proper and so on and so on

*sniffle*
*/

std::default_random_engine update_generator;
std::uniform_int_distribution<int> int_distribution(0,N_V-1);
std::uniform_real_distribution<double> real_distribution(0.0,1.0);

int update(double* array, double* gauss){

  int xthis;
  int xnext;
  int xprev;

  int ythis;
  int ynext;
  int yprev;

  int tthis;
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

  // int N_acceptances = 0;

  for (int x = 0; x < Lx; x++){

    // get x neighbour coordinates

    xthis = 3*x;
    xnext = 3*((x+1)%Lx);
    xprev = 3*((x+Lx-1)%Lx);

    for (int y = 0; y < Ly; y++){

      // get y neighbour coordinates

      ythis = 3*Lx*y;
      ynext = 3*Lx*((y+1)%Ly);
      yprev = 3*Lx*((y+Ly-1)%Ly);

      for (int t = 0; t < Lt; t++){

        // get t neighbour coordinates

        tthis = 3*Lx*Ly*t;
        tnext = 3*Lx*Ly*((t+1)%Lt);
        tprev = 3*Lx*Ly*((t+Lt-1)%Lt);

        // update the x link

        old_link = array[xthis + ythis + tthis];
        new_link = old_link + gauss[int_distribution(update_generator)];
        // new_link = old_link + gauss[int_distribution(generator)];

        staple1 = array[xnext + ythis + tthis + 1] - array[xthis + ynext + tthis] - array[xthis + ythis + tthis + 1];
        staple2 = - array[xnext + yprev + tthis + 1] - array[xthis + yprev + tthis] + array[xthis + yprev + tthis + 1];
        staple3 = array[xnext + ythis + tthis + 2] - array[xthis + ythis + tnext] - array[xthis + ythis + tthis + 2];
        staple4 = - array[xnext + ythis + tprev + 2] - array[xthis + ythis + tprev] + array[xthis + ythis + tprev + 2];

        old_action = cos(old_link + staple1) + cos(old_link + staple2) + cos(old_link + staple3) + cos(old_link + staple4);
        new_action = cos(new_link + staple1) + cos(new_link + staple2) + cos(new_link + staple3) + cos(new_link + staple4);

        // note that if the new action is smaller than the old action, then
        // C > 1 so we accept the new link.
        C = exp(- beta * new_action)/exp(- beta * old_action);

        z = real_distribution(update_generator);
        // z = real_distribution(generator);

        if (z < C){
          // N_acceptances++;
          array[xthis + ythis + tthis] = new_link;
        }

        // update the y link

        old_link = array[xthis + ythis + tthis + 1];
        new_link = old_link + gauss[int_distribution(update_generator)];
        // new_link = old_link + gauss[int_distribution(generator)];

        staple1 = array[xthis + ynext + tthis] - array[xnext + ythis + tthis + 1] - array[xthis + ythis + tthis];
        staple2 = - array[xprev + ynext + tthis] - array[xprev + ythis + tthis + 1] + array[xprev + ythis + tthis];
        staple3 = array[xthis + ynext + tthis + 2] - array[xthis + ythis + tnext + 1] - array[xthis + ythis + tthis + 2];
        staple4 = - array[xthis + ynext + tprev + 2] - array[xthis + ythis + tprev + 1] + array[xthis + ythis + tprev + 2];

        old_action = cos(old_link + staple1) + cos(old_link + staple2) + cos(old_link + staple3) + cos(old_link + staple4);
        new_action = cos(new_link + staple1) + cos(new_link + staple2) + cos(new_link + staple3) + cos(new_link + staple4);

        C = exp(- beta * new_action)/exp(- beta * old_action);

        z = real_distribution(update_generator);
        // z = real_distribution(generator);

        if (z < C){
          // N_acceptances++;
          array[xthis + ythis + tthis + 1] = new_link;
        }

        // update the t link

        old_link = array[xthis + ythis + tthis + 2];
        new_link = old_link + gauss[int_distribution(update_generator)];
        // new_link = old_link + gauss[int_distribution(generator)];

        staple1 = array[xthis + ythis + tnext] - array[xnext + ythis + tthis + 2] - array[xthis + ythis + tthis];
        staple2 = - array[xprev + ythis + tnext] - array[xprev + ythis + tthis + 2] + array[xprev + ythis + tthis];
        staple3 = array[xthis + ythis + tnext + 1] - array[xthis + ynext + tthis + 2] - array[xthis + ythis + tthis + 1];
        staple4 = - array[xthis + yprev + tnext + 1] - array[xthis + yprev + tthis + 2] + array[xthis + yprev + tthis + 1];

        old_action = cos(old_link + staple1) + cos(old_link + staple2) + cos(old_link + staple3) + cos(old_link + staple4);
        new_action = cos(new_link + staple1) + cos(new_link + staple2) + cos(new_link + staple3) + cos(new_link + staple4);

        C = exp(- beta * new_action)/exp(- beta * old_action);

        z = real_distribution(update_generator);
        // z = real_distribution(generator);

        if (z < C){
          // N_acceptances++;
          array[xthis + ythis + tthis + 2] = new_link;
        }

      }
    }
  }

  // std::cout << "acceptance rate: " << double(N_acceptances) / double(N_links) << std::endl;

  return 0;
}
