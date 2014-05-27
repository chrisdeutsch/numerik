#include "matrix.h"
#include <stdio.h>

int main() {
	MATRIX m = matrix_alloc(3, 1);
  MATRIX n = matrix_alloc(1, 3);
  MATRIX result;
  
  m.elem[0][0] = 1;
  m.elem[1][0] = 2;
  m.elem[2][0] = 3;
  
  n.elem[0][0] = 3;
  n.elem[0][1] = 2;
  n.elem[0][2] = 1;
  
  matrix_print(&m);
  printf("\n");
  matrix_print(&n);
  printf("\n");
  
  result = matrix_mult(&m, &n);
  
  matrix_print(&result);
  printf("\n");
  
  return 0;
}