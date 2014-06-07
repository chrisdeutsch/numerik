/* gcc -o numerik_4 -O2 numerik_bespin_deutsch_linalg.c numerik_bespin_deutsch_gls.c numerik_bespin_deutsch_4.c -lm */
/* Christian Bespin, Christopher Deutsch */

/* Programmaufruf: Erklaerung bei Aufruf des Programms ohne Argumente */

#include "numerik_bespin_deutsch_linalg.h"
#include "numerik_bespin_deutsch_gls.h"
#include <stdio.h>

/* Erstellt eine Tabelle vo */
void resistance_table(void (*GLS)(MATRIX*, VECTOR*, double*), int size,
                      int resistor, double start, double stop, double step);

int main(int argc, char **argv) {
  void (*GLS[5])(MATRIX *, VECTOR *, double *);
  char *label[5];
  
  int geometry, resistor;
  int i;
  double start, stop, step;
  
  /* Zugehoerige Gleichungssysteme
   * implementiert in numerik_bespin_deutsch_gls.c */
  GLS[0] = cube_diag;
  label[0] = "Raumdiagonale eines Wuerfels";
  
  GLS[1] = cube_facediag;
  label[1] = "Flaechendiagonale eines Wuerfels";
  
  GLS[2] = cube_edge;
  label[2] = "Kante eines Wuerfels";
  
  GLS[3] = octahedron_diag;
  label[3] = "Raumdiagonale eines Oktaeders";
  
  GLS[4] = octahedron_edge;
  label[4] = "Kante eines Oktaeders";
  
  /* Ueberprueft die Anzahl der Programmargumente und erklaert die Benutzung */
  if ( argc != 6 ) {
    printf("Benutzung:\n");
    printf("%s geometry resistor start stop step\n\n", argv[0]);
    printf("geometry: Geometrie des Widerstandsnetzwerks\n");
    for (i = 0; i < 5; i++) {
      printf("  %i: %s\n", i, label[i]);
    }
    printf("\nresistor: Der zu aendernde Widerstand. Dabei enspricht\n"
           "          resistor=i R(i+1) in den Skizzen der PDF\n\n");
    printf("start, stop, step: Berechnet den Ersatzwiderstand des Netzwerks\n"
           "                   fuer Aenderung des Widerstands \"resistor\"\n"
           "                   von \"start\" bis \"stop\" in Schritten von "
           "\"step\"\n"
           "                   Einheit: Ohm\n");
    
    return -1;
  }
  
  /* Liesst die Programmargumente ein und ueberprueft auf Fehler */
  if ( sscanf(argv[1], "%i", &geometry) != 1 ||
       sscanf(argv[2], "%i", &resistor) != 1 ||
       sscanf(argv[3], "%lf", &start) != 1 ||
       sscanf(argv[4], "%lf", &stop) != 1 ||
       sscanf(argv[5], "%lf", &step) != 1 ) {
    printf("Fehler in den uebergebenen Argumenten!\n");
    return -1;
  }
  
  /* Ueberprueft die Programmargumente auf Richtigkeit */
  if ( geometry < 0 || geometry > 4 ) {
    printf("Ungueltige Geometrie\n");
    return -1;
  }
  if ( resistor < 0 || resistor > 11 ) {
    printf("Ungueltiger Widerstand\n");
    return -1;
  }
  if ( start < 0 ) {
    printf("Der Startwert fuer den Widerstand muss groesser 0 sein\n");
    return -1;
  }
  if ( stop < start ) {
    printf("Der Endwert fuer den Widerstand muss"
           " groesser als der Startwert sein\n");
    return -1;
  }
  if ( step <= 0 ) {
    printf("Die Schrittweite muss groesser als 0 sein\n");
    return -1;
  }
  
  /* Tabellenkopf */
  printf("# %s\n", label[geometry]);
  printf("# Alle anderen Widerstaende R = 1 Ohm\n");
  
  /* Berechnung und Ausgabe der Ersatzwiderstaende */
  resistance_table(GLS[geometry], geometry > 2 ? 8 : 6,
                   resistor, start, stop, step);
  
  return 0;
}

void resistance_table(void (*GLS)(MATRIX*, VECTOR*, double*), int size,
                      int resistor, double start, double stop, double step) {
  /* Allokiert Speicher fuer Koeffizientenmatrix, Inhomogenitaet und Loesung */
  MATRIX *M = matrix_alloc(size);
  VECTOR *b = vector_alloc(size);
  VECTOR *sol = vector_alloc(size);
  
  /* Array der Widerstaende R[i] entspricht dabei R(i+1) in der PDF */
  double R[12] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
  
  /* Ueberprueft ob alle Speicherallokierungen erfolgreich waren */
  if ( M == NULL || b == NULL || sol == NULL ) {
    matrix_free(M);
    vector_free(b);
    vector_free(sol);
    printf("Probleme bei der Allokierung von Speicher!\n");
    return;
  }
  
  /* setzt den veraenderlichen Widerstand auf seinen Startwert */
  R[resistor] = start;
  
  /* Tabellenkopf */
  printf("# R%i / Ohm\tR_E / Ohm\n", resistor + 1);
  
  while (R[resistor] <= stop) {
    /* Setzt "M" und "b" auf das Gleichungssystem fuer die gegebenen Widerstands-
     * werte in dem Array "R" */
    GLS(M, b, R);
    
    /* Loest das Gleichungssystem und speichert die Loesung in "sol" */
    if ( linear_solve(M, b, sol) == 0 ) {
      /* Gleichungssystem erfolgreich geloest: Berechnung des Ersatzwiderstandes
       * durch R = U/I_ges, I_ges ist jeweis das letzte Element des Loesungs-
       * vektors und U wurde 1V gewaehlt */
      printf("%f\t%f\n", R[resistor], 1.0 / sol->elem[size - 1]);
    }
    
    R[resistor] += step;
  }
  
  matrix_free(M);
  vector_free(b);
  vector_free(sol);
}
