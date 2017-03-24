
#include <iostream>

int bar(double* array, int index){

  array[index] = 2;

  return 0;
}

int main(){
  double foo [10] = {};

  bar(foo, 3);

  std::cout << foo[3] << std::endl;
}
