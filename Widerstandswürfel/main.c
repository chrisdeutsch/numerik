#include "matrix.h"
#include <stdio.h>

int main() {
	MATRIX m = matrix_alloc(2, 2);
  int permutation[2] = {1, 2};
  
  m.elem[0][0] = 4;
  m.elem[0][1] = 3;
  m.elem[1][0] = 6;
  m.elem[1][1] = 3;
  
  matrix_print(&m);
  printf("\n");
  
  LU_decomp(&m, permutation);
  
  matrix_print(&m);
  printf("\n");
  
  printf("%i, %i\n", permutation[0], permutation[1]);
  
  return 0;
}
