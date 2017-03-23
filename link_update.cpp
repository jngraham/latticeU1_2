
/*
/////////////////////////////////////////////////////////////////////////////

U(1) Lattice Gauge Theory Simulator | James Graham

Here I write the function(s) to update the lattice

/////////////////////////////////////////////////////////////////////////////
*/

#include <random>
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

int update(float* lattice, int Lx, int Ly, int Lt){
  return 0;
}
