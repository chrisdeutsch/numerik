#include "matrix.h"

int main() {
	MATRIX m = matrix_alloc(3, 1);
  MATRIX n = matrix_alloc(1, 3);
  MATRIX result;
  
  matrix_set(&m, 0, 0, 1);
  matrix_set(&m, 1, 0, 2);
  matrix_set(&m, 2, 0, 3);
  
  matrix_set(&n, 0, 0, 3);
  matrix_set(&n, 0, 1, 2);
  matrix_set(&n, 0, 2, 1);
  
  
  matrix_print(&m);
  matrix_print(&n);
  
  result = matrix_mult(&m, &n);
  
  matrix_print(&result);
  return 0;
}