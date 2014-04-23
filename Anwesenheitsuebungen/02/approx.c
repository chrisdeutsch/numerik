#include <stdio.h>
#include <stdlib.h>
#include <math.h>

typedef struct {
  double x;
  double y;
} tuple;

double *newton_coeff(tuple *stuetz, int n);

int main() {
  int count = 4;

  tuple *s = malloc(count * sizeof(tuple));
  for (int i = 0; i < count; i++) {
    s->x = i;
    s->y = sinh(i);
  }

  double *coeff = NULL;

  coeff = newton_coeff(s, count);

  for (int i = 0; i < count; i++) {
    printf("C%i: %f\n", i, coeff[i]);
  }

  free(coeff);
  return 0;
}

double *newton_coeff(tuple *stuetz, int n) {
  double *ret;
  ret = malloc(n * sizeof(tuple));

  // Matrix initialisieren
  double **mat;
  mat = malloc(n * sizeof(double *));      // Zeilen
  *mat = malloc((n + 1) * sizeof(double)); // Spalten

  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n + 1; j++) {
      mat[i][j] = 0;
    }
  }

  // Stuetzwerte einlesen
  for (int i = 0; i < n; i++) {
    mat[i][0] = stuetz->x;
    mat[i][1] = stuetz->y;
  }

  // Berechnen
  for (int j = 2; j < n + 1; j++) {
    for (int i = 0; i <= n - j; i++) {
      mat[i][j] =
          (mat[i + 1][j - 1] - mat[i][j - 1]) / (mat[i + 1][0] - mat[i][0]);
    }
  }

  /*
  mat[i][j] i zeilen j spalte
  */
  for (int i = 0; i < n; i++) {
    ret[i] = mat[0][i + 1];
  }

  free(*mat);
  free(mat);
  return ret;
}
