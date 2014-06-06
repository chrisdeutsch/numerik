#include "numerik_bespin_deutsch_linalg.h"
#include "gls.h"
#include <stdio.h>

void table(void (*func)(MATRIX*, VECTOR*, double*), int dim, int resistor, double start, double stop, double step);

int main() {
  int i;
  for (i = 0; i < 12; i++) {
    table(oktahedron_edge, 8, i, 0.0, 100.0, 0.1);
  }
  
  
  return 0;
}

void table(void (*func)(MATRIX*, VECTOR*, double*), int dim, int resistor, double start, double stop, double step) {
  MATRIX *m = matrix_alloc(dim);
  VECTOR *b = vector_alloc(dim);
  VECTOR *sol = vector_alloc(dim);
  
  double R[12] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
  
  R[resistor] = start;
  
  printf("\nR%i\tR\n", resistor + 1);
  while (R[resistor] <= stop) {
    func(m, b, R);
    if ( linear_solve(m, b, sol) == 0 ) {
      printf("%f\t%f\n", R[resistor], 1.0 / sol->elem[dim - 1]);
    }
    
    R[resistor] += step;  
  }
  
  matrix_free(m);
  vector_free(b);
  vector_free(sol);
}
