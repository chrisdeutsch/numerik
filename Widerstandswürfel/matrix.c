#include "matrix.h"
#include <stdlib.h>
#include <stdio.h>

MATRIX matrix_alloc(int n, int m) {
  MATRIX ret;
  
  /* Fehlercheck fehlt */
  ret.data = malloc(n * m * sizeof(double));
  ret.n = n;
  ret.m = m;
  
  return ret;  
}

void matrix_free(MATRIX A) {
  /* Fehlercheck */
  free(A.data);
}

void matrix_print(MATRIX A) {
  int i, j;
  for(i = 0; i < n; i++) {
    for(j = 0; j < m; j++) {
      printf("test");
    }
  }
  
}

/* Access */
double matrix_get(MATRIX *A, int i, int j) {
  return A->data[i * A->m + j];
}

void matrix_set(MATRIX *A, int i, int j, double value) {
  A->data[i * A->m + j] = value;
}

/* Operationen */
MATRIX matrix_mult(MATRIX A, MATRIX B);