#include "matrix.h"
#include <stdio.h>

int main() {
  MATRIX m = matrix_alloc(6, 6);
  int permutation[6] = {1, 2, 3, 4, 5, 6};
  
  double R[12] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
  
  matrix_init(&m, 0.0);
  
  /*erste Zeile*/
  m.elem[0][0] =  R[0] + R[1] + R[2] + R[3];
  m.elem[0][1] = -R[3];
  m.elem[0][2] =  0;
  m.elem[0][3] = -R[2];
  m.elem[0][4] = -R[1];
  m.elem[0][5] =  0;

  /*zweite Zeile*/
  m.elem[1][0] = -R[3];
  m.elem[1][1] =  R[3] + R[5] + R[6] + R[9];
  m.elem[1][2] = -R[9];
  m.elem[1][3] = -R[6];
  m.elem[1][4] = -R[5];
  m.elem[1][5] = -R[5];

  /*dritte Zeile*/
  m.elem[2][0] =  0;
  m.elem[2][1] = -R[9];
  m.elem[2][2] =  R[8] + R[9] + R[10] + R[11];
  m.elem[2][3] = -R[10];
  m.elem[2][4] = -R[8];
  m.elem[2][5] = -(R[8] + R[11]);

  /*vierte Zeile*/
  m.elem[3][0] = -R[2];
  m.elem[3][1] = -R[6];
  m.elem[3][2] = -R[10];
  m.elem[3][3] =  R[2] + R[6] + R[7] + R[10];
  m.elem[3][4] =  0;
  m.elem[3][5] =  0;

  /*fünfte Zeile*/
  m.elem[4][0] = -R[1];
  m.elem[4][1] = -R[5];
  m.elem[4][2] = -R[8];
  m.elem[4][3] =  0;
  m.elem[4][4] =  R[1] + R[4] + R[5] + R[8];
  m.elem[4][5] =  R[5] + R[8];

  /*sechste Zeile*/
  m.elem[5][0] =  0;
  m.elem[5][1] = -R[5];
  m.elem[5][2] = -R[8] - R[11];
  m.elem[5][3] =  0;
  m.elem[5][4] =  R[5] + R[8];
  m.elem[5][5] =  R[5] + R[8] + R[11];

  matrix_print(&m);
  printf("\n");
  
  if (LU_decomp(&m, permutation) == -1) {
    printf("Matrix ist (fast) singulaer!\n");
  } else {
  
    matrix_print(&m);
    printf("\n");
    
    printf("{%i, %i, %i, %i, %i, %i}\n", permutation[0], permutation[1], permutation[2], permutation[3], permutation[4], permutation[5]);
  }
  return 0;
}
