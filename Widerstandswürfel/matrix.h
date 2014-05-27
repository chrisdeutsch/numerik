#ifndef _MATRIX_H_
#define _MATRIX_H_

typedef struct {
  double *data;
  double **elem;
  int n, m;
} MATRIX;

/* Infrastruktur */
MATRIX matrix_alloc(int n, int m);

void matrix_init(MATRIX *A, double value);

void matrix_free(MATRIX *A);

void matrix_print(MATRIX *A);

/* Operationen */
MATRIX matrix_mult(MATRIX *A, MATRIX *B);


#endif
