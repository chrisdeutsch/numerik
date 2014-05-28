#include "matrix.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

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

void matrix_swap_row(MATRIX *A, int i, int j) {
  double *temp;
  
  temp = A->elem[i];
  A->elem[i] = A->elem[j];
  A->elem[j] = temp;
}

void LU_decomp(MATRIX *A, int *permutation) {
  int i, j, k;
  int piv;
  int temp;
  
  for (i = 0; i < A->n; i++) {            
    for (j = i; j < A->n; j++) {
      for (k = 0; k < i; k++) {
        A->elem[i][j] -= A->elem[i][k] * A->elem[k][j];
      }
    }
    for (j = i + 1; j < A->n; j++) {
      for (k = 0; k < i; k++) {
        A->elem[j][i] -= A->elem[j][k] * A->elem[k][i];
      }
      A->elem[j][i] /= A->elem[i][i];
    }
  }
}

int pivot(MATRIX *A, int k) {
  int i, j;
  int piv = k;
  double max = 0;
  double temp, sum;
  
  if (k < A->n - 1) {
    for (i = k; i < A->n; i++) {
      sum = 0;
      for (j = k; j < A->n; j++) {
        sum += fabs(A->elem[i][j]);
      }
      temp = A->elem[i][k] / sum;
    
      if (temp > max) {
        max = temp;
        piv = i;
      }
    }
  }
  
  return piv;
}

VECTOR vector_alloc(int n) {
  VECTOR ret;
  
  /* Fehlercheck */
  ret.elem = malloc(n * sizeof(double));
  
  return ret;
}

void vector_init(VECTOR *v, double value) {
  int i;
  for (i = 0; i < v->n; i++) {
    v->elem[i] = value;
  }
}

void vector_free(VECTOR *v) {
  free(v->elem);
}

void vector_print(VECTOR *v) {
  int i;
  for (i = 0; i < v->n; i++) {
    printf("%f\n", v->elem[i]);
  }
}

VECTOR matrix_vector_mult(MATRIX *A, VECTOR *v) {
  VECTOR ret;
  int i, j;
  double sum;
  
  if (A->m != v->n) {
    printf("Nicht kompatibel\n");
  }
  
  ret = vector_alloc(A->n);
  
  
  for (i = 0; i < A->n; i++) {
    sum = 0;
    for (j = 0; j < A->m; j++) {
      break; /*WEITER*/
    }
  }
  
    return ret;
}
