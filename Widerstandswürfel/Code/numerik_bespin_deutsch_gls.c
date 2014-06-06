#include "numerik_bespin_deutsch_gls.h"

void cube_diag(MATRIX *A, VECTOR *b, double *R) {
  /* Inhomogenitaet */
  b->elem[0] = 0;
  b->elem[1] = 0;
  b->elem[2] = 0;
  b->elem[3] = 0;
  b->elem[4] = 0;
  b->elem[5] = 1;
  
  /* Koeffizientenmatrix */
  /* erste Zeile */
  A->elem[0][0] =  R[0] + R[1] + R[2] + R[3];
  A->elem[0][1] = -R[3];
  A->elem[0][2] =  0;
  A->elem[0][3] = -R[2];
  A->elem[0][4] = -R[1];
  A->elem[0][5] =  0;

  /* zweite Zeile */
  A->elem[1][0] = -R[3];
  A->elem[1][1] =  R[3] + R[5] + R[6] + R[9];
  A->elem[1][2] = -R[9];
  A->elem[1][3] = -R[6];
  A->elem[1][4] = -R[5];
  A->elem[1][5] = -R[5];

  /* dritte Zeile */
  A->elem[2][0] =  0;
  A->elem[2][1] = -R[9];
  A->elem[2][2] =  R[8] + R[9] + R[10] + R[11];
  A->elem[2][3] = -R[10];
  A->elem[2][4] = -R[8];
  A->elem[2][5] = -R[8] - R[11];

  /* vierte Zeile */
  A->elem[3][0] = -R[2];
  A->elem[3][1] = -R[6];
  A->elem[3][2] = -R[10];
  A->elem[3][3] =  R[2] + R[6] + R[7] + R[10];
  A->elem[3][4] =  0;
  A->elem[3][5] =  0;

  /* fünfte Zeile */
  A->elem[4][0] = -R[1];
  A->elem[4][1] = -R[5];
  A->elem[4][2] = -R[8];
  A->elem[4][3] =  0;
  A->elem[4][4] =  R[1] + R[4] + R[5] + R[8];
  A->elem[4][5] =  R[5] + R[8];

  /* sechste Zeile */
  A->elem[5][0] =  0;
  A->elem[5][1] = -R[5];
  A->elem[5][2] = -R[8] - R[11];
  A->elem[5][3] =  0;
  A->elem[5][4] =  R[5] + R[8];
  A->elem[5][5] =  R[5] + R[8] + R[11];
}

void cube_facediag(MATRIX *A, VECTOR *b, double *R) {
  /* Inhomogenitaet */
  b->elem[0] = 0;
  b->elem[1] = 0;
  b->elem[2] = 0;
  b->elem[3] = 0;
  b->elem[4] = 0;
  b->elem[5] = 1;
  
  /* Koeffizientenmatrix */
  /* erste Zeile */
  A->elem[0][0] =  R[0] + R[1] + R[2] + R[3];
  A->elem[0][1] = -R[3];
  A->elem[0][2] =  0;
  A->elem[0][3] = -R[2];
  A->elem[0][4] = -R[1];
  A->elem[0][5] =  0;

  /* zweite Zeile */
  A->elem[1][0] = -R[3];
  A->elem[1][1] =  R[3] + R[5] + R[6] + R[9];
  A->elem[1][2] = -R[9];
  A->elem[1][3] = -R[6];
  A->elem[1][4] = -R[5];
  A->elem[1][5] = -R[5];

  /* dritte Zeile */
  A->elem[2][0] =  0;
  A->elem[2][1] = -R[9];
  A->elem[2][2] =  R[8] + R[9] + R[10] + R[11];
  A->elem[2][3] = -R[10];
  A->elem[2][4] = -R[8];
  A->elem[2][5] = -R[8];

  /* vierte Zeile */
  A->elem[3][0] = -R[2];
  A->elem[3][1] = -R[6];
  A->elem[3][2] = -R[10];
  A->elem[3][3] =  R[2] + R[6] + R[7] + R[10];
  A->elem[3][4] =  0;
  A->elem[3][5] =  0;

  /* fünfte Zeile */
  A->elem[4][0] = -R[1];
  A->elem[4][1] = -R[5];
  A->elem[4][2] = -R[8];
  A->elem[4][3] =  0;
  A->elem[4][4] =  R[1] + R[4] + R[5] + R[8];
  A->elem[4][5] =  R[5] + R[8];

  /* sechste Zeile */
  A->elem[5][0] =  0;
  A->elem[5][1] = -R[5];
  A->elem[5][2] = -R[8];
  A->elem[5][3] =  0;
  A->elem[5][4] =  R[5] + R[8];
  A->elem[5][5] =  R[5] + R[8];
}

void cube_edge(MATRIX *A, VECTOR *b, double *R) {
  /* Inhomogenitaet */
  b->elem[0] = 0;
  b->elem[1] = 0;
  b->elem[2] = 0;
  b->elem[3] = 0;
  b->elem[4] = 0;
  b->elem[5] = 1;
  
  /* Koeffizientenmatrix */
  /* erste Zeile */
  A->elem[0][0] =  R[0] + R[1] + R[2] + R[3];
  A->elem[0][1] = -R[3];
  A->elem[0][2] =  0;
  A->elem[0][3] = -R[2];
  A->elem[0][4] = -R[1];
  A->elem[0][5] =  0;

  /* zweite Zeile */
  A->elem[1][0] = -R[3];
  A->elem[1][1] =  R[3] + R[5] + R[6] + R[9];
  A->elem[1][2] = -R[9];
  A->elem[1][3] = -R[6];
  A->elem[1][4] = -R[5];
  A->elem[1][5] = 0;

  /* dritte Zeile */
  A->elem[2][0] =  0;
  A->elem[2][1] = -R[9];
  A->elem[2][2] =  R[8] + R[9] + R[10] + R[11];
  A->elem[2][3] = -R[10];
  A->elem[2][4] = -R[8];
  A->elem[2][5] = -R[11];

  /* vierte Zeile */
  A->elem[3][0] = -R[2];
  A->elem[3][1] = -R[6];
  A->elem[3][2] = -R[10];
  A->elem[3][3] =  R[2] + R[6] + R[7] + R[10];
  A->elem[3][4] =  0;
  A->elem[3][5] =  0;

  /* fünfte Zeile */
  A->elem[4][0] = -R[1];
  A->elem[4][1] = -R[5];
  A->elem[4][2] = -R[8];
  A->elem[4][3] =  0;
  A->elem[4][4] =  R[1] + R[4] + R[5] + R[8];
  A->elem[4][5] =  0;

  /* sechste Zeile */
  A->elem[5][0] =  0;
  A->elem[5][1] =  0;
  A->elem[5][2] = -R[11];
  A->elem[5][3] =  0;
  A->elem[5][4] =  0;
  A->elem[5][5] =  R[11];
}

void octahedron(MATRIX *A, VECTOR *b, double *R) {
  /* Inhomogenitaet */
  b->elem[0] = 1;
  b->elem[1] = 0;
  b->elem[2] = 0;
  b->elem[3] = 0;
  b->elem[4] = 0;
  b->elem[5] = 0;
  b->elem[6] = 0;
  b->elem[7] = 0;
  
  /* Koeffizientenmatrix */
  /* erste Zeile */
  A->elem[0][0] =  0;
  A->elem[0][1] =  0;
  A->elem[0][2] = -R[3];
  A->elem[0][3] =  0;
  A->elem[0][4] =  0;
  A->elem[0][5] = -R[11];
  A->elem[0][6] =  0;
  A->elem[0][7] =  R[3] + R[11];
  
  /* zweite Zeile */
  A->elem[1][0] =  R[0] + R[1] + R[4];
  A->elem[1][1] = -R[1];
  A->elem[1][2] =  0;
  A->elem[1][3] = -R[4];
  A->elem[1][4] =  0;
  A->elem[1][5] =  0;
  A->elem[1][6] = -R[4];
  A->elem[1][7] =  0;
  
  /* dritte Zeile */
  A->elem[2][0] = -R[1];
  A->elem[2][1] =  R[1] + R[2] + R[5];
  A->elem[2][2] = -R[2];
  A->elem[2][3] =  0;
  A->elem[2][4] = -R[5];
  A->elem[2][5] =  0;
  A->elem[2][6] = -R[5];
  A->elem[2][7] =  0;
  
  /* vierte Zeile */
  A->elem[3][0] =  0;
  A->elem[3][1] = -R[2];
  A->elem[3][2] =  R[2] + R[3] + R[6];
  A->elem[3][3] =  0;
  A->elem[3][4] =  0;
  A->elem[3][5] = -R[6];
  A->elem[3][6] = -R[6];
  A->elem[3][7] = -R[3];
  
  /* fuenfte Zeile */
  A->elem[4][0] = -R[4];
  A->elem[4][1] =  0;
  A->elem[4][2] =  0;
  A->elem[4][3] =  R[4] + R[8] + R[9];
  A->elem[4][4] = -R[9];
  A->elem[4][5] =  0;
  A->elem[4][6] =  R[4];
  A->elem[4][7] =  0;
  
  /* sechste Zeile */
  A->elem[5][0] =  0;
  A->elem[5][1] = -R[5];
  A->elem[5][2] =  0;
  A->elem[5][3] = -R[9];
  A->elem[5][4] =  R[5] + R[9] + R[10];
  A->elem[5][5] = -R[10];
  A->elem[5][6] =  R[5];
  A->elem[5][7] =  0;
  
  /* siebte Zeile */
  A->elem[6][0] =  0;
  A->elem[6][1] =  0;
  A->elem[6][2] = -R[6];
  A->elem[6][3] =  0;
  A->elem[6][4] = -R[10];
  A->elem[6][5] =  R[6] + R[10] + R[11];
  A->elem[6][6] =  R[6];
  A->elem[6][7] = -R[11];
  
  /* achte Zeile */
  A->elem[7][0] = -R[4];
  A->elem[7][1] = -R[5];
  A->elem[7][2] = -R[6];
  A->elem[7][3] =  R[4];
  A->elem[7][4] =  R[5];
  A->elem[7][5] =  R[6];
  A->elem[7][6] =  R[4] + R[5] + R[6] + R[7];
  A->elem[7][7] =  0;
}

void octahedron_edge(MATRIX *A, VECTOR *b, double *R) {
  /* Inhomogenitaet */
  b->elem[0] = 1;
  b->elem[1] = 0;
  b->elem[2] = 0;
  b->elem[3] = 0;
  b->elem[4] = 0;
  b->elem[5] = 0;
  b->elem[6] = 0;
  b->elem[7] = 0;
  
  /* Koeffizientenmatrix */
  /* erste Zeile */
  A->elem[0][0] =  0;
  A->elem[0][1] =  0;
  A->elem[0][2] = -R[3];
  A->elem[0][3] =  0;
  A->elem[0][4] =  0;
  A->elem[0][5] =  0;
  A->elem[0][6] =  0;
  A->elem[0][7] =  R[3];
  
  /* zweite Zeile */
  A->elem[1][0] =  R[0] + R[1] + R[4];
  A->elem[1][1] = -R[1];
  A->elem[1][2] =  0;
  A->elem[1][3] = -R[4];
  A->elem[1][4] =  0;
  A->elem[1][5] =  0;
  A->elem[1][6] = -R[4];
  A->elem[1][7] =  0;
  
  /* dritte Zeile */
  A->elem[2][0] = -R[1];
  A->elem[2][1] =  R[1] + R[2] + R[5];
  A->elem[2][2] = -R[2];
  A->elem[2][3] =  0;
  A->elem[2][4] = -R[5];
  A->elem[2][5] =  0;
  A->elem[2][6] = -R[5];
  A->elem[2][7] =  0;
  
  /* vierte Zeile */
  A->elem[3][0] =  0;
  A->elem[3][1] = -R[2];
  A->elem[3][2] =  R[2] + R[3] + R[6];
  A->elem[3][3] =  0;
  A->elem[3][4] =  0;
  A->elem[3][5] = -R[6];
  A->elem[3][6] = -R[6];
  A->elem[3][7] = -R[3];
  
  /* fuenfte Zeile */
  A->elem[4][0] = -R[4];
  A->elem[4][1] =  0;
  A->elem[4][2] =  0;
  A->elem[4][3] =  R[4] + R[8] + R[9];
  A->elem[4][4] = -R[9];
  A->elem[4][5] =  0;
  A->elem[4][6] =  R[4];
  A->elem[4][7] =  0;
  
  /* sechste Zeile */
  A->elem[5][0] =  0;
  A->elem[5][1] = -R[5];
  A->elem[5][2] =  0;
  A->elem[5][3] = -R[9];
  A->elem[5][4] =  R[5] + R[9] + R[10];
  A->elem[5][5] = -R[10];
  A->elem[5][6] =  R[5];
  A->elem[5][7] =  0;
  
  /* siebte Zeile */
  A->elem[6][0] =  0;
  A->elem[6][1] =  0;
  A->elem[6][2] = -R[6];
  A->elem[6][3] =  0;
  A->elem[6][4] = -R[10];
  A->elem[6][5] =  R[6] + R[10] + R[11];
  A->elem[6][6] =  R[6];
  A->elem[6][7] =  0;
  
  /* achte Zeile */
  A->elem[7][0] = -R[4];
  A->elem[7][1] = -R[5];
  A->elem[7][2] = -R[6];
  A->elem[7][3] =  R[4];
  A->elem[7][4] =  R[5];
  A->elem[7][5] =  R[6];
  A->elem[7][6] =  R[4] + R[5] + R[6] + R[7];
  A->elem[7][7] =  0;
}
