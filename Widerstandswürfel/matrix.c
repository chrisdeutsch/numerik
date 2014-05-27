#include "matrix.h"
#include <stdlib.h>
#include <stdio.h>

MATRIX matrix_alloc(int n, int m) {
  MATRIX ret;
  int i;
  
  /* Fehlercheck fehlt */
  ret.data = malloc(n * m * sizeof(double));
  ret.n = n;
  ret.m = m;
  
  /* GET Pointerstruktur */
  ret.elem = malloc(n * sizeof(double *));
  for (i = 0; i < n; i++) {
    ret.elem[i] = ret.data + i * m;
  }
  
  return ret;  
}

void matrix_init(MATRIX *A, double value) {
  int i, j;
  for (i = 0; i < A->n; i++) {
    for (j = 0; j < A->m; j++) {
      A->elem[i][j] = value;
    }
  }
}

void matrix_free(MATRIX *A) {
  /* Fehlercheck */
  free(A->elem);
  free(A->data);
}

void matrix_print(MATRIX *A) {
  int i, j;
  for (i = 0; i < A->n; i++) {
    for (j = 0; j < A->m; j++) {
      printf("%f ", A->elem[i][j]);
    }
    printf("\n");
  }
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
        sum += A->elem[i][k] * B->elem[k][j];
      }
      ret.elem[i][j] = sum;
    }
  }
  
  return ret;
}
