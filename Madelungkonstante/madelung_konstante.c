/* gcc -O2 -std=c99 -Wall -pedantic -o madelung madelung_konstante.c -lm */
/* Christian Bespin, Christopher Deutsch */

#include <stdio.h>
#include <math.h>

/* Literaturwert nach Y. Sakamoto */
#define MAD_CONST_3D 1.7475645946331822
/* Hier noch eine vernünftige Quelle */
#define MAD_CONST_2D 1.6155426267128247

/* Berechnet die Madelungkonstante für einen Würfel; bricht ab wenn die Änderung
 * vom vorigen Durchlauf kleiner als epsilon ist */
double madelung_3d(double epsilon);

/* Berechnet die Madelungkonstante für einen Würfel; bricht ab wenn die Änderung
 * vom vorigen Durchlauf kleiner als epsilon ist !!!!(Vermutlich keine hinreichende Bedingung)!!!! */
double madelung_2d(double epsilon);

/* Berechnet den Betrag des Vektors (x,y,z) */
double dist(int x, int y, int z);

/* Vorzeichenfunktion */
int sign_z(int x, int y, int z);

/* Beschreibung was eigentlich in der main-Funktion berechnet wird */
int main() {
  printf("Berechnung der Madelung-Konstante\n");
  printf("###################################\n");
  printf("3D: \n");
  printf("alpha\t\tdelta\n");
  double mad = madelung_3d(1E-5);
  printf("%.12f\t%.2E\n", mad, fabs(mad - MAD_CONST_3D)/MAD_CONST_3D);
  printf("###################################\n");
  printf("2D: \n");
  printf("alpha\t\tdelta\n");
  mad = madelung_2d(1E-5);
  printf("%.12f\t%.2E\n", mad, fabs(mad - MAD_CONST_2D)/MAD_CONST_2D);
  
  return 0;
}

double madelung_3d(double epsilon) {
  double sum = 0; // wo ist der Unterschied zwischen sum und residual?
  double residual = 0;
  double prev = 1E100; // Vorheriger Schleifendurchgang: -1, da die for Schleife
                       // betreten werden soll

  for (int m = 1; fabs(sum - prev) > epsilon; m++) {
    prev = sum;
    /* Den Rest vom letzten Durchgang aufaddieren und danach resetten */
    sum += residual;
    residual = 0;

    /* Fläche (welche Fläche wird berechnet?)
       sum: 6 (Sym.) * 1/2 (Evjens) = 3
       residual: 6 (Sym.) * 1/2 (Rest) = 3
    */
    for (int y = -m + 1; y < m; y++) {
      for (int x = -m + 1; x < m; x++) {
        sum += 3 * sign_z(x, y, m) / dist(x, y, m);
        residual += 3 * sign_z(x, y, m) / dist(x, y, m);
      }
    }

    /* Kanten (welche Kante wird berechnet?)
       sum: 12 (Sym.) * 1/4 (Evjens) = 3
       residual: 12 (Sym.) * 3/4 (Rest) = 9
    */
    for (int x = -m + 1; x < m; x++) {
      sum += 3 * sign_z(x, m, m) / dist(x, m, m);
      residual += 9 * sign_z(x, m, m) / dist(x, m, m);
    }

    /* Ecken (welche Ecke wird berechnet?)
       sum: 8 (Sym.) * 1/8 (Evjens) = 1
       residual: 8 (Sym.) * 7/8 (Rest) = 7
    */
    sum += sign_z(m, m, m) / dist(m, m, m);
    residual += 7 * sign_z(m, m, m) / dist(m, m, m);
  }

  return sum;
}

double madelung_2d(double epsilon) {
  double sum = 0;
  double residual = 0;

  for (int m = 1; fabs(sum - MAD_CONST_2D) > epsilon; m++) {
    /* Den Rest vom letzten Durchgang aufaddieren und danach resetten */
    sum += residual;
    residual = 0;

    /* Kante (welche Kante wird berechnet?)
       sum: 4 (Sym.) * 1/2 (Evjens) = 2
       residual: 4 (Sym.) * 1/2 (Rest) = 2
    */
    for (int x = -m + 1; x < m; x++) {
      sum += 2 * sign_z(x, m, 0) / dist(x, m, 0);
      residual += 2 * sign_z(x, m, 0) / dist(x, m, 0);
    }

    /* Ecken * 4 wegen Symmetrie (welche Ecke wird berechnet?)
       sum: 4 (Sym.) * 1/4 (Evjens) = 1
       residual: 4 (Sym.) * 3/4 (Rest) = 3
    */
    sum += sign_z(m, m, 0) / dist(m, m, 0);
    residual += 3 * sign_z(m, m, 0) / dist(m, m, 0);
  }

  return sum;
}

double dist(int x, int y, int z) { return sqrt(x * x + y * y + z * z); }

/* Es gilt pow(-1, x + y + z) == (x + y + z) % 2 ? -1 : 1
   In Worten: Wenn der Exponent nicht gerade ist -> -1
              sonst -> 1 */
int sign_z(int x, int y, int z) { return -((x + y + z) % 2 ? -1 : 1); }
