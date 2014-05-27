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

void matrix_init(MATRIX *A, double value) {
  int i, j;
  for (i = 0; i < A->n; i++) {
    for (j = 0; j < A->m; j++) {
      matrix_set(A, i, j, value);
    }
  }
}

void matrix_free(MATRIX *A) {
  /* Fehlercheck */
  free(A->data);
}

void matrix_print(MATRIX *A) {
  int i, j;
  for (i = 0; i < A->n; i++) {
    for (j = 0; j < A->m; j++) {
      printf("%f ", matrix_get(A, i, j));
    }
    printf("\n");
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
MATRIX matrix_mult(MATRIX *A, MATRIX *B) {
  /* (L x M) * (M x N) = (L x n) */
  MATRIX ret;
  int i, j, k;
  double sum;
  
  /* Kompatibilitaet */
  if (A->m != B->n) {
    printf("Matrizen nicht kompatibel\n");
  }
  
  ret = matrix_alloc(A->n, B->m);
  
  for (i = 0; i < A->n; i++) {
    for (j = 0; j < B->m; j++) {
      sum = 0;
      for (k = 0; k < A-> m; k++) {
        sum += matrix_get(A, i, k) * matrix_get(B, k, j);
      }
      matrix_set(&ret, i, j, sum);
    }
  }
  
  return ret;
}
