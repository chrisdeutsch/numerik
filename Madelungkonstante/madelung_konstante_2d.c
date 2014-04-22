/*
HISTORY:
*/

/* gcc -O2 -std=c99 -Wall -pedantic -o madelung madelung_konstante_2d.c -lm */
/* Christian Bespin, Christopher Deutsch */

#include <stdio.h>
#include <math.h>

/* Berechnet die Madelungkonstante für einen Quadrat der Seitenlänge 2n wobei
   die Länge
   auf den Gitterabstand normalisiert ist */
double evjens_method(int n);

/* Beschreibung was eigentlich in der main-Funktion berechnet wird */
int main() {
  printf("Madelungkonstante: %f\n", evjens_method(100));
  return 0;
}

double evjens_method(int n) {
  double sum = 0;
  double rest = 0;

  /* effizient der pow-Funktion diskutieren (könnte durch eine Modulo Operation
     eliminiert werden
     macht aber den Code unklarer */
  for (int i = 1; i <= n; i++) {
    /* Den Rest vom letzten Durchgang aufaddieren und danach resetten */
    sum += rest;
    rest = 0;

    /*
    magic numbers müssen noch entfernt bzw. erklärt werden
    */

    /* Kante  * 4 wegen Symmetrie (welche Kante wird berechnet?) */
    for (int x = -i + 1; x < i; x++) {
      sum += 2 * pow(-1, x + i) / sqrt(x * x + i * i);
      rest += 2 * pow(-1, x + i) / sqrt(x * x + i * i);
    }

    /* Ecken * 4 wegen Symmetrie (welche Ecke wird berechnet?)*/
    sum += pow(-1, i + i) / sqrt(i * i + i * i);
    rest += 3 * pow(-1, i + i) / sqrt(i * i + i * i);
  }

  return sum;
}
