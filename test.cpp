
#include <iostream>
#include <stdio.h>      /* printf, NULL */
#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */

int main ()
{
  const int Lx = 4;
  const int Ly = 4;
  const int Lt = 6;

  double lattice [3*Lx*Ly*Lt] = {};

  for (int x = 0; x < Lx; x++){
    for (int y = 0; y < Ly; y++){
      for (int t = 0; t < Lt; t++){
        lattice[3*x + 3*Lx*y + 3*Lx*Ly*t + 0] = x;
        lattice[3*x + 3*Lx*y + 3*Lx*Ly*t + 1] = y;
        lattice[3*x + 3*Lx*y + 3*Lx*Ly*t + 2] = t;
      }
    }
  }

  for (size_t i = 0; i < 3*Lx*Ly*Lt; i++){
    std::cout << lattice[i] << " ";

    if ((i+1)%(3*Lx)==0){
      std::cout << std::endl;
    }
  }

  std::cout << std::endl;

}
