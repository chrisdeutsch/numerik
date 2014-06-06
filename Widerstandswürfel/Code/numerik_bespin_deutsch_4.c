#include "numerik_bespin_deutsch_linalg.h"
#include "numerik_bespin_deutsch_gls.h"
#include <stdio.h>

void table(void (*func)(MATRIX*, VECTOR*, double*), int dim, int resistor, double start, double stop, double step);

int main(int argc, char **argv) {
  void (*func[5])(MATRIX*, VECTOR*, double*);
  char *label[5];
  
  int geometry, resistor;
  double start, stop, step;
  
  func[0] = cube_diag;
  label[0] = "Raumdiagonale eines Wuerfels";
  
  func[1] = cube_facediag;
  label[1] = "Flaechendiagonale eines Wuerfels";
  
  func[2] = cube_edge;
  label[2] = "Kante eines Wuerfels";
  
  func[3] = octahedron;
  label[3] = "gegenueberliegende Spitzen eines Oktaeders";
  
  func[4] = octahedron_edge;
  label[4] = "Kante eines Oktaeders";
  
  if ( argc != 6 ) {
    printf("Benutzung:\n");
    return -1;
  }
  
  
  if ( sscanf(argv[1], "%i", &geometry) != 1 ||
       sscanf(argv[2], "%i", &resistor) != 1 ||
       sscanf(argv[3], "%lf", &start) != 1 ||
       sscanf(argv[4], "%lf", &stop) != 1 ||
       sscanf(argv[5], "%lf", &step) != 1 ) {
    printf("Fehler in den uebergebenen Argumenten!\n");
    return -1;
  }
  
  printf("## %s ##\n", label[geometry]);
  printf("Geo: %i\n", geometry);
  printf("R: %i\n", resistor);
  printf("Start: %f\n", start);
  printf("Stop: %f\n", stop);
  printf("Step: %f\n", step);
  printf("%i\n", geometry > 2 ? 8 : 6);
  table(func[geometry], geometry > 2 ? 8 : 6, resistor, start, stop, step);
  
  
  return 0;
}

void table(void (*func)(MATRIX*, VECTOR*, double*), int dim, int resistor, double start, double stop, double step) {
  MATRIX *m = matrix_alloc(dim);
  VECTOR *b = vector_alloc(dim);
  VECTOR *sol = vector_alloc(dim);
  
  double R[12] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
  
  if ( m == NULL || b == NULL || sol == NULL ) {
    matrix_free(m);
    vector_free(b);
    vector_free(sol);
    printf("Probleme bei der Allokierung von Speicher!\n");
    return;
  }
  
  R[resistor] = start;
  
  printf("# R%i\t\tR_E\n", resistor + 1);
  while (R[resistor] <= stop) {
    func(m, b, R);
    if ( linear_solve(m, b, sol) == 0 ) {
      printf("%f\t%f\n", R[resistor], 1.0 / sol->elem[dim - 1]);
    }
    
    R[resistor] += step;  
  }
  
  matrix_free(m);
  vector_free(b);
  vector_free(sol);
}
