#include "matrix.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

MATRIX *matrix_alloc(int n) {
  int i;
  MATRIX *ret;
  
  /* Allokiert ein Matrix-Struct */
  if ( NULL == (ret = malloc(sizeof(MATRIX))) ) {
    return NULL;
  }
  
  /* Speicher fuer die Matrixelemente */
  if ( NULL == (ret->data = malloc(n * n * sizeof(double))) ) {
    free(ret);
    return NULL;
  }
  ret->n = n;
  
  /* Darstellung als 2D-Array */
  if ( NULL == (ret->elem = malloc(n * sizeof(double*))) ) {
    free(ret);
    free(ret->data);
    return NULL;
  }
  for (i = 0; i < n; i++) {
    ret->elem[i] = ret->data + i * n;
  }
  
  return ret;
}

void matrix_init(MATRIX *A, double value) {
  int i;
  int n = A->n;
  
  for (i = 0; i < n * n; i++) {
    A->data[i] = value;
  }
}

void matrix_free(MATRIX *A) {
  free(A->data);
  free(A->elem);
  free(A);
}

void matrix_print(MATRIX *A) {
  int i, j;
  int n = A->n;
  
  for (i = 0; i < n; i++) {
    for (j = 0; j < n; j++) {
      printf("%f ", A->elem[i][j]);
    }
    printf("\n");
  }
}

void matrix_swap_row(MATRIX *A, int i, int j) {
  double *temp;
  int n = A->n;
  
  if ( i < n && j < n ) {
    temp = A->elem[i];
    A->elem[i] = A->elem[j];
    A->elem[j] = temp;
  }
}

VECTOR *vector_alloc(int n) {
  VECTOR *ret;
  
  if ( NULL == (ret = malloc(sizeof(VECTOR))) ) {
    return NULL;
  }
  
  if ( NULL == (ret->elem = malloc(n * sizeof(double))) ) {
    free(ret);
    return NULL;
  }
  ret->n = n;
  
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
  free(v);
}

void vector_print(VECTOR *v) {
  int i;
  for (i = 0; i < v->n; i++) {
    printf("%f\n", v->elem[i]);
  }
}

int LU_decomp(MATRIX *A, int *permutation) {
  int i, j, k;
  int n = A->n;
  int piv, temp;
  
  /* Spalte der Gauß-Elimination */
  for (i = 0; i < n; i++) {
    /* Berechnung des Pivot-Elements */
    piv = pivot_row(A, i);
    if ( piv != i ) {
      /* Tauschen der Zeilen */
      matrix_swap_row(A, i, piv);
      /* Merken der Vertauschung */
      temp = permutation[i];
      permutation[i] = permutation[piv];
      permutation[piv] = temp;
    }
    
    /* Berechnung der oberen Dreiecksmatrix U */
    for (j = i; j < n; j++) {
      for (k = 0; k < i; k++) {
        A->elem[i][j] -= A->elem[i][k] * A->elem[k][j];
      }
    }
    
    /* Ueberpruefung ob das Diagonalelement ausreichend von 0 verschieden ist,
     * also die Matrix nicht singulaer bzw. fast singulaer ist */
    if ( fabs(A->elem[i][i]) < 1E-10 ) {
      return -1;
    }
    
    /* Berechnung der unteren Dreiecksmatrix L */
    for (j = i + 1; j < n; j++) {
      for (k = 0; k < i; k++) {
        A->elem[j][i] -= A->elem[j][k] * A->elem[k][i];
      }
      A->elem[j][i] /= A->elem[i][i];
    }
  }
  
  return 0;
}

int pivot_row(MATRIX *A, int k) {
  int i, j;
  int n = A->n;
  int piv;
  double max = 0;
  double temp, sum;
  
  /* Im letzten Schritt der Gauß-Elimination ist kein Zeilentausch moeglich */
  if ( k == n - 1 ) return k;
  
  /* Iteration ueber alle moeglichen Pivot-Zeilen */
  for (i = k; i < n; i++) {
    sum = 0;
    /* Berechnung des Terms, der maximiert werden soll */
    for (j = k; j < n; j++) {
      sum += fabs(A->elem[i][j]);
    }
    temp = A->elem[i][k] / sum;
    
    if ( temp > max ) {
      max = temp;
      piv = i;
    }
  }
  
  return piv;
}

int LU_solve(MATRIX *LU, VECTOR *Pb, VECTOR *sol) {
  /* Es ist LUx = Pb zu loesen. */
  VECTOR *y = vector_alloc(Pb->n);
  
  /* Loese zuerst Ly = Pb mit y = Ux durch Vorwaertssubstitution */
  LU_forward_sub(LU, Pb, y);
  /* Loese Ux = y durch Rueckwaertssubstitution */
  LU_back_sub(LU, y, sol);
  
  vector_free(y);
  
  return 0;
}

int LU_forward_sub(MATRIX *LU, VECTOR *b, VECTOR *sol) {
  int i, j;
  int n = LU->n;
  
  if ( n != b->n ) return -1;
  
  for (i = 0; i < n; i++) {
    sol->elem[i] = b->elem[i];
    for (j = 0; j < i; j++) {
      sol->elem[i] -= LU->elem[i][j] * sol->elem[j];
    }
    /* Es wird hier nicht durch L[i][i] geteilt, da die Diagonalelemente der
     * unteren Dreickecksmatrix 1 sind. */
  }
  
  return 0;
}

int LU_back_sub(MATRIX *LU, VECTOR *b, VECTOR *sol) {
  int i, j;
  int n = LU->n;
  
  if ( n != b->n ) return -1;
  
  for (i = n - 1; i >= 0; i--) {
    sol->elem[i] = b->elem[i];
    for (j = i + 1; j < n; j++) {
      sol->elem[i] -= LU->elem[i][j] * sol->elem[j];
    }
    sol->elem[i] /= LU->elem[i][i];
  }
  
  return 0;
}

int linear_solve(MATRIX *A, VECTOR *b, VECTOR *sol) {
  int i;
  int n = A->n;
  int *permutation = malloc(n * sizeof(int));
  VECTOR *Pb = vector_alloc(n);
  
  if (A->n != b->n) return -1;
  
  /* Fuellt das Permutationsarray mit {0, 1, ..., n-1} */
  for (i = 0; i < n; i++) {
    permutation[i] = i;
  }
  
  /* LU-Zerlegung der Matrix */
  LU_decomp(A, permutation);
  
  /* !!!DEBUG MESSAGE!!! */
  printf("{");
  for (i = 0; i < n - 1; i++) {
    printf("%i, ", permutation[i]);
  }
  printf("%i}\n", permutation[n-1]);
  
  /* Permutation von b */
  for (i = 0; i < n; i++) {
    Pb->elem[i] = b->elem[permutation[i]];
  }
  
  /* Loesung des identischen Gleichungssystems LUx = Pb */
  LU_solve(A, Pb, sol);
  
  return 0;
}
