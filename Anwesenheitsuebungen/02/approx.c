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

  return ret;
}
