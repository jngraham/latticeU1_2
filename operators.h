
/*
/////////////////////////////////////////////////////////////////////////////

U(1) Lattice Gauge Theory Simulator | James Graham

Here I declare the functions for the field operators on the lattice

/////////////////////////////////////////////////////////////////////////////
*/

float avg_plaquette(float* lattice, int Lx, int Ly, int Lt);

/*
We have field oeprators that find whatever it is we want at the appropriate time.
We shall call these functions from lattice.cpp when we want to populate our data
*/

float jpc_plus(float* lattice, int Lx, int Ly, int t);

float jpc_minus(float* lattice, int Lx, int Ly, int t);

float flux(float* lattice, int Lx, int Ly, int t);
