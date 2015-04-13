/* gcc -o numerik_5 -O2 numerik_bespin_deutsch_sparse_matrix.c numerik_bespin_deutsch_poisson.c numerik_bespin_deutsch_5.c -lm */
/* Christian Bespin, Christopher Deutsch */

/* Programmaufruf: Erklaerung bei Aufruf des Programms ohne Argumente */

#include <stdio.h>
#include "numerik_bespin_deutsch_sparse_matrix.h"
#include "numerik_bespin_deutsch_poisson.h"

int main(int argc, char **argv) {
  GRID *grid = NULL;
  SPARSE_MATRIX *A = NULL;
  VECTOR *b = NULL, *x = NULL;
  int i, ret = -1, a = -1, sym = -1;
  
  if (argc != 3) {
    printf("Benutzung: %s a sym\n"
           "a: Gitterabstand in Hundersteln\n"
           "sym: 0 keine Symmetrie\n"
           "     1 Symmetrie\n", argv[0]);
  }
  
  if (sscanf(argv[1], "%i", &a) != 1 ||
      sscanf(argv[2], "%i", &sym) != 1) {
    printf("Auslesen der Programmargumente fehlgeschlagen\n");
    return -1;
  }
  
  /* Erstelle die Geometrie des Problems in "grid" */
  printf("Diskretisierung der Geometrie...\n");
  if (sym == 0) {
    ret = geometry(&grid, a);
  } else if (sym == 1) {
    ret = geometry_sym(&grid, a);
  } else {
    printf("Symmetrieparameter ungueltig\n");
    return -1;
  }
  
  if (ret == -1) {
    printf("Fehler bei der Allokierung des Speichers fuer die Geometrie\n");
    return -1;
  } else if (ret == -2) {
    printf("Die Geometrie ist bei gegebenem Gitterabstand nicht darstellbar.\n"
           "Bitte ein a (in Hundersteln), welches 50 ohne Rest teilt.\n"
           "Ohne Symmetriebetrachtung muss a nur 100 ohne Rest teilen.\n");
    return -1;
  }
  
  /* Erstelle das aus der Diskretisierung folgende Gleichungsystem in der Koef-
   * fizientenmatrix "A" und der Inhomogenitaet "b" */
  printf("Aufstellen des Gleichungssystems...\n");
  ret = setup_gls(grid, &A, &b);
  if (ret == -1) {
    printf("Fehler bei der Allokierung des Speichers fuer das Gleichungssystem\n");
    return -1;
  } else if (ret == -2) {
    printf("Fehler bei der Allokierung des Speichers der duennen Matrix\n");
    return -1;
  }
  
  
  /* Naeherungsloesung fuer das Gauss-Seidel-Verfahren */
  printf("Loesen des Gleichungssystems...\n");
  x = vector_alloc(grid->eq_count);
  if (x == NULL) {
    printf("Fehler bei der Allokierung des Speichers des Loesungsvektors\n");
  }
  for (i = 0; i < grid->eq_count; i++) {
    x->elem[i] = 0.24;
  }
  
  /* Loese das Gleichungssystem mit dem Gauß-Seidel-Verfahren. "x" enthaelt
   * zunaechst den Startvektor des Iterationsverfahrens und nachher die Loesung
   * des Gleichungssystems A x = b */
  gauss_seidel(A, b, x, 1E-6);
  
  /* Ordnet die berechnete Loesung wieder in die Geometrie ein */
  enter_solution(grid, x);
  
  /* Output */
  if (sym == 0) {
    mathematica_output(grid, "bespin_deutsch_poisson_loesung.txt");
    printf("Ausgabe geschrieben in bespin_deutsch_poisson_loesung.txt\n");
  } else if (sym == 1) {
    mathematica_sym_output(grid, "bespin_deutsch_poisson_loesung.txt");
    printf("Ausgabe geschrieben in bespin_deutsch_poisson_loesung.txt\n");
  }
  
  /* Allokierten Speicher freigeben */
  grid_free(grid);
  matrix_free(A);
  vector_free(b);
  vector_free(x);
  
  return 0;
}
