#ifndef _MATRIX_H_
#define _MATRIX_H_

typedef struct {
  double *data;
  int n, m;
} MATRIX;

/* Infrastruktur */
MATRIX matrix_alloc(int n, int m);

void matrix_free(MATRIX A);

void matrix_print(MATRIX A);

/* Access */
double matrix_get(MATRIX *A, int i, int j);

double matrix_set(MATRIX *A, int i, int j, double value);

/* Operationen */
MATRIX matrix_mult(MATRIX A, MATRIX B);


#endif
