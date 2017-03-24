
/*
/////////////////////////////////////////////////////////////////////////////

U(1) Lattice Gauge Theory Simulator | James Graham

Here I declare the functions for the field operators on the lattice

/////////////////////////////////////////////////////////////////////////////
*/

double avg_plaquette(double* lattice);

/*
We have field oeprators that find whatever it is we want at the appropriate time.
We shall call these functions from lattice.cpp when we want to populate our data
*/

double jpc_plus(double* lattice);

double jpc_minus(double* lattice);

double flux(double* lattice);
