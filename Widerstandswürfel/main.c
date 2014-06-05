#include "numerik_bespin_deutsch_linalg.h"
#include "gls.h"
#include <stdio.h>



int main() {
  MATRIX *m = matrix_alloc(6);
  VECTOR *b = vector_alloc(6);
  VECTOR *sol = vector_alloc(6);
  
  double R[12] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
  
  cube_diag(m, b, R);  
  
  matrix_print(m);
  printf("\n");
  
  
  linear_solve(m, b, sol);
  
  vector_print(sol);
  printf("\n");
  
  matrix_free(m);
  vector_free(b);
  vector_free(sol);
  
  return 0;
}
