#include <stdio.h>
#include <stdlib.h>

typedef struct {
  double x;
  double y;
} tuple;

double *newton_coeff(tuple *stuetz, int n);

int main() {
  tuple *s = NULL;
  int count = 0;
  double *coeff = NULL;

  coeff = newton_coeff(s, count);

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
    for (int i = 0; i < n - j; i++) {
      break;
    }
  }

  /*
  mat[i][j] i zeilen j spalte
  */

  free(*mat);
  free(mat);
  return ret;
}
