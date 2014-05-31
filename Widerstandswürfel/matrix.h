#ifndef _MATRIX_H_
#define _MATRIX_H_

typedef struct {
  double *data;
  double **elem;
  int n, m;
} MATRIX;

typedef struct {
  double *elem;
  int n;
} VECTOR;

/* Infrastruktur */
MATRIX matrix_alloc(int n, int m);

void matrix_init(MATRIX *A, double value);

void matrix_free(MATRIX *A);

void matrix_print(MATRIX *A);


/* Operationen */
MATRIX matrix_mult(MATRIX *A, MATRIX *B);

void matrix_swap_row(MATRIX *A, int i, int j);

int LU_decomp(MATRIX *A, int *permutation);

int pivot(MATRIX *A, int k);

/* Löst ein LU aufgespaltenes System  L*U*x = P*b dabei ist Pb bereits permutiert*/
int LU_solve(MATRIX *LU, VECTOR *Pb, VECTOR *sol);

int LU_forward_sub(MATRIX *LU, VECTOR *b, VECTOR *sol);

int LU_back_sub(MATRIX *LU, VECTOR *b, VECTOR *sol);

int linear_solve(MATRIX *A, VECTOR *b, VECTOR *sol);

/* Vektorkram */
VECTOR vector_alloc(int n);

void vector_init(VECTOR *v, double value);

void vector_free(VECTOR *v);

void vector_print(VECTOR *v);

#endif
