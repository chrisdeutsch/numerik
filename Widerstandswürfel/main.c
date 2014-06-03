#include "matrix.h"
#include <stdio.h>

void cube_diag(MATRIX *m, VECTOR *b, double *R);
void oktahedron(MATRIX *m, VECTOR *b, double *R);
void oktahedron_edge(MATRIX *m, VECTOR *b, double *R);

int main() {
  MATRIX *m = matrix_alloc(8);
  VECTOR *b = vector_alloc(8);
  VECTOR *sol = vector_alloc(8);
  
  double R[12] = {5, 5, 2, 5, 2, 2, 2, 2, 5, 2, 2, 5};
  
  oktahedron_edge(m, b, R);
  
  
  matrix_print(m);
  printf("\n");
  
  
  linear_solve(m, b, sol);
  
  vector_print(sol);
  printf("\n");
  
  return 0;
}

void cube_diag(MATRIX *m, VECTOR *b, double *R) {
  /* Inhomogenitaet */
  b->elem[0] = 0;
  b->elem[1] = 0;
  b->elem[2] = 0;
  b->elem[3] = 0;
  b->elem[4] = 0;
  b->elem[5] = 1;

  /* erste Zeile */
  m->elem[0][0] = R[0] + R[1] + R[2] + R[3];
  m->elem[0][1] = -R[3];
  m->elem[0][2] = 0;
  m->elem[0][3] = -R[2];
  m->elem[0][4] = -R[1];
  m->elem[0][5] = 0;

  /* zweite Zeile */
  m->elem[1][0] = -R[3];
  m->elem[1][1] = R[3] + R[5] + R[6] + R[9];
  m->elem[1][2] = -R[9];
  m->elem[1][3] = -R[6];
  m->elem[1][4] = -R[5];
  m->elem[1][5] = -R[5];

  /* dritte Zeile */
  m->elem[2][0] = 0;
  m->elem[2][1] = -R[9];
  m->elem[2][2] = R[8] + R[9] + R[10] + R[11];
  m->elem[2][3] = -R[10];
  m->elem[2][4] = -R[8];
  m->elem[2][5] = -R[8] - R[11];

  /* vierte Zeile */
  m->elem[3][0] = -R[2];
  m->elem[3][1] = -R[6];
  m->elem[3][2] = -R[10];
  m->elem[3][3] = R[2] + R[6] + R[7] + R[10];
  m->elem[3][4] = 0;
  m->elem[3][5] = 0;

  /* fünfte Zeile */
  m->elem[4][0] = -R[1];
  m->elem[4][1] = -R[5];
  m->elem[4][2] = -R[8];
  m->elem[4][3] = 0;
  m->elem[4][4] = R[1] + R[4] + R[5] + R[8];
  m->elem[4][5] = R[5] + R[8];

  /* sechste Zeile */
  m->elem[5][0] = 0;
  m->elem[5][1] = -R[5];
  m->elem[5][2] = -R[8] - R[11];
  m->elem[5][3] = 0;
  m->elem[5][4] = R[5] + R[8];
  m->elem[5][5] = R[5] + R[8] + R[11];
}

void oktahedron(MATRIX *m, VECTOR *b, double *R) {
  vector_init(b, 0.0);
  
  /* Inhomogenitaet */
  b->elem[0] = 1;
  
  /* erste Zeile */
  m->elem[0][0] = 0;
  m->elem[0][1] = 0;
  m->elem[0][2] = -R[3];
  m->elem[0][3] = 0;
  m->elem[0][4] = 0;
  m->elem[0][5] = -R[11];
  m->elem[0][6] = 0;
  m->elem[0][7] = R[3] + R[11];
  
  /* zweite Zeile */
  m->elem[1][0] = R[0] + R[1] + R[4];
  m->elem[1][1] = -R[1];
  m->elem[1][2] = 0;
  m->elem[1][3] = -R[4];
  m->elem[1][4] = 0;
  m->elem[1][5] = 0;
  m->elem[1][6] = -R[5];
  m->elem[1][7] = 0;
  
  /* dritte Zeile */
  m->elem[2][0] = -R[1];
  m->elem[2][1] = R[1] + R[2] + R[5];
  m->elem[2][2] = -R[2];
  m->elem[2][3] = 0;
  m->elem[2][4] = -R[5];
  m->elem[2][5] = 0;
  m->elem[2][6] = -R[5];
  m->elem[2][7] = 0;
  
  /* vierte Zeile */
  m->elem[3][0] = 0;
  m->elem[3][1] = -R[2];
  m->elem[3][2] = R[2] + R[3] + R[6];
  m->elem[3][3] = 0;
  m->elem[3][4] = 0;
  m->elem[3][5] = -R[6];
  m->elem[3][6] = -R[6];
  m->elem[3][7] = -R[3];
  
  /* fuenfte Zeile */
  m->elem[4][0] = -R[4];
  m->elem[4][1] = 0;
  m->elem[4][2] = 0;
  m->elem[4][3] = R[9] + R[4] + R[8];
  m->elem[4][4] = -R[9];
  m->elem[4][5] = 0;
  m->elem[4][6] = R[4];
  m->elem[4][7] = 0;
  
  /* sechste Zeile */
  m->elem[5][0] = 0;
  m->elem[5][1] = -R[5];
  m->elem[5][2] = 0;
  m->elem[5][3] = -R[9];
  m->elem[5][4] = R[5] + R[9] + R[10];
  m->elem[5][5] = -R[10];
  m->elem[5][6] = R[5];
  m->elem[5][7] = 0;
  
  /* siebte Zeile */
  m->elem[6][0] = 0;
  m->elem[6][1] = 0;
  m->elem[6][2] = -R[6];
  m->elem[6][3] = 0;
  m->elem[6][4] = -R[10];
  m->elem[6][5] = R[6] + R[10] + R[11];
  m->elem[6][6] = R[6];
  m->elem[6][7] = -R[11];
  
  /* achte Zeile */
  m->elem[7][0] = -R[4];
  m->elem[7][1] = -R[5];
  m->elem[7][2] = -R[6];
  m->elem[7][3] = R[4];
  m->elem[7][4] = R[5];
  m->elem[7][5] = R[6];
  m->elem[7][6] = R[4] + R[5] + R[6] + R[7];
  m->elem[7][7] = 0;
}

void oktahedron_edge(MATRIX *m, VECTOR *b, double *R) {
  vector_init(b, 0.0);
  
  /* Inhomogenitaet */
  b->elem[0] = 1;
  
  /* erste Zeile */
  m->elem[0][0] = 0;
  m->elem[0][1] = 0;
  m->elem[0][2] = -R[3];
  m->elem[0][3] = 0;
  m->elem[0][4] = 0;
  m->elem[0][5] = 0;
  m->elem[0][6] = 0;
  m->elem[0][7] = R[3];
  
  /* zweite Zeile */
  m->elem[1][0] = R[0] + R[1] + R[4];
  m->elem[1][1] = -R[1];
  m->elem[1][2] = 0;
  m->elem[1][3] = -R[4];
  m->elem[1][4] = 0;
  m->elem[1][5] = 0;
  m->elem[1][6] = -R[5];
  m->elem[1][7] = 0;
  
  /* dritte Zeile */
  m->elem[2][0] = -R[1];
  m->elem[2][1] = R[1] + R[2] + R[5];
  m->elem[2][2] = -R[2];
  m->elem[2][3] = 0;
  m->elem[2][4] = -R[5];
  m->elem[2][5] = 0;
  m->elem[2][6] = -R[5];
  m->elem[2][7] = 0;
  
  /* vierte Zeile */
  m->elem[3][0] = 0;
  m->elem[3][1] = -R[2];
  m->elem[3][2] = R[2] + R[3] + R[6];
  m->elem[3][3] = 0;
  m->elem[3][4] = 0;
  m->elem[3][5] = -R[6];
  m->elem[3][6] = -R[6];
  m->elem[3][7] = -R[3];
  
  /* fuenfte Zeile */
  m->elem[4][0] = -R[4];
  m->elem[4][1] = 0;
  m->elem[4][2] = 0;
  m->elem[4][3] = R[9] + R[4] + R[8];
  m->elem[4][4] = -R[9];
  m->elem[4][5] = 0;
  m->elem[4][6] = R[4];
  m->elem[4][7] = 0;
  
  /* sechste Zeile */
  m->elem[5][0] = 0;
  m->elem[5][1] = -R[5];
  m->elem[5][2] = 0;
  m->elem[5][3] = -R[9];
  m->elem[5][4] = R[5] + R[9] + R[10];
  m->elem[5][5] = -R[10];
  m->elem[5][6] = R[5];
  m->elem[5][7] = 0;
  
  /* siebte Zeile */
  m->elem[6][0] = 0;
  m->elem[6][1] = 0;
  m->elem[6][2] = -R[6];
  m->elem[6][3] = 0;
  m->elem[6][4] = -R[10];
  m->elem[6][5] = R[6] + R[10] + R[11];
  m->elem[6][6] = R[6];
  m->elem[6][7] = 0;
  
  /* achte Zeile */
  m->elem[7][0] = -R[4];
  m->elem[7][1] = -R[5];
  m->elem[7][2] = -R[6];
  m->elem[7][3] = R[4];
  m->elem[7][4] = R[5];
  m->elem[7][5] = R[6];
  m->elem[7][6] = R[4] + R[5] + R[6] + R[7];
  m->elem[7][7] = 0;
}