/* gcc -O2 -std=c99 -Wall -pedantic -o madelung madelung_konstante.c -lm */
/* Christian Bespin, Christopher Deutsch */

#include <stdio.h>
#include <math.h>

/* Literaturwert nach Y. Sakamoto */
#define MAD_CONST_3D 1.7475645946331822
/* Hier noch eine vernünftige Quelle */
#define MAD_CONST_2D 1.61554

/* Berechnet die Madelungkonstante für einen Würfel der Kantenlänge 2n wobei die
   Länge
   auf den Gitterabstand normalisiert ist */
double madelung_3d(int n);

/* Berechnet die Madelungkonstante für ein Quadrat der Seitenlänge 2n wobei die
   Länge
   auf den Gitterabstand normalisiert ist */
double madelung_2d(int n);

/* Berechnet den Betrag des Vektors (x,y,z) */
double dist(int x, int y, int z);

/* Vorzeichenfunktion */
int sign_z(int x, int y, int z);

/* Beschreibung was eigentlich in der main-Funktion berechnet wird */
int main() {
  printf("Output:\n");
  return 0;
}

double madelung_3d(int n) {
  double sum = 0; // wo ist der Unterschied zwischen sum und rest?
  double rest = 0;

  for (int m = 1; m <= n;
       m++) { // Madelungkonst von Würfel mit Kantenlänge 2 + Madelungkonst von
              // Würfel mit Kantenlänge 4 + ...?
    /* Den Rest vom letzten Durchgang aufaddieren und danach resetten */
    sum += rest;
    rest = 0;

    /* Fläche (welche Fläche wird berechnet?)
       sum: 6 (Sym.) * 1/2 (Evjens) = 3
       rest: 6 (Sym.) * 1/2 (Rest) = 3
    */
    for (int y = -m + 1; y < m; y++) {
      for (int x = -m + 1; x < m; x++) {
        sum += 3 * sign_z(x, y, m) / dist(x, y, m);
        rest += 3 * sign_z(x, y, m) / dist(x, y, m);
      }
    }

    /* Kanten (welche Kante wird berechnet?)
       sum: 12 (Sym.) * 1/4 (Evjens) = 3
       rest: 12 (Sym.) * 3/4 (Rest) = 9
    */
    for (int x = -m + 1; x < m; x++) {
      sum += 3 * sign_z(x, m, m) / dist(x, m, m);
      rest += 9 * sign_z(x, m, m) / dist(x, m, m);
    }

    /* Ecken (welche Ecke wird berechnet?)
       sum: 8 (Sym.) * 1/8 (Evjens) = 1
       rest: 8 (Sym.) * 7/8 (Rest) = 7
    */
    sum += sign_z(m, m, m) / dist(m, m, m);
    rest += 7 * sign_z(m, m, m) / dist(m, m, m);
  }

  return sum;
}

double madelung_2d(int n) {
  double sum = 0;
  double rest = 0;

  for (int m = 1; m <= n; m++) {
    /* Den Rest vom letzten Durchgang aufaddieren und danach resetten */
    sum += rest;
    rest = 0;

    /* Kante (welche Kante wird berechnet?)
       sum: 4 (Sym.) * 1/2 (Evjens) = 2
       rest: 4 (Sym.) * 1/2 (Rest) = 2
    */
    for (int x = -m + 1; x < m; x++) { // War das nicht Kanten 1/2 gewichtet?
      sum += 2 * sign_z(x, m, 0) / dist(x, m, 0);
      rest += 2 * sign_z(x, m, 0) / dist(x, m, 0);
    }

    /* Ecken * 4 wegen Symmetrie (welche Ecke wird berechnet?)
       sum: 4 (Sym.) * 1/4 (Evjens) = 1
       rest: 4 (Sym.) * 3/4 (Rest) = 3
    */
    sum += sign_z(m, m, 0) / dist(m, m, 0);
    rest += 3 * sign_z(m, m, 0) / dist(m, m, 0);
  }

  return sum;
}

double dist(int x, int y, int z) { return sqrt(x * x + y * y + z * z); }

/* Es gilt pow(-1, x + y + z) == (x + y + z) % 2 ? -1 : 1
   In Worten: Wenn der Exponent nicht gerade ist -> -1
              sonst -> 1 */
int sign_z(int x, int y, int z) { return -(x + y + z) % 2 ? -1 : 1; }
