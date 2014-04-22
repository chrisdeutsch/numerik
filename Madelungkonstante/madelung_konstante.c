/* gcc -O2 -std=c99 -Wall -pedantic -o madelung madelung_konstante.c -lm */
/* Christian Bespin, Christopher Deutsch */

#include <stdio.h>
#include <math.h>

/* Berechnet die Madelungkonstante für einen Würfel der Kantenlänge 2n wobei die Länge 
   auf den Gitterabstand normalisiert ist */
double evjens_method(int n);

/* Beschreibung was eigentlich in der main-Funktion berechnet wird */
int main() {
  /* Literaturwert nach Y. Sakamoto */
  const double lit = 1.7475645946331822;
  double mad;

  for (int i = 1; i < 10; i++) {
    mad = -evjens_method(i);
    printf("Madelung-Konstante: %.12f \t n: %i \t Delta_rel: %.2E\n", mad, i,
           fabs(mad - lit) / lit);
  }
  
  for (int i = 1; i < 10; i++) {
    mad = -evjens_method(10 * i);
    printf("Madelung-Konstante: %.12f \t n: %i \t Delta_rel: %.2E\n", mad, 10 * i,
           fabs(mad - lit) / lit);
  }
  return 0;
}

double evjens_method(int n) {
  double sum = 0;
  double rest = 0;
  
  /* effizient der pow-Funktion diskutieren (könnte durch eine Modulo Operation eliminiert werden
     macht aber den Code unklarer */
  for (int i = 1; i <= n; i++) {
    /* Den Rest vom letzten Durchgang aufaddieren und danach resetten */
    sum += rest;
    rest = 0;

    /*
    magic numbers müssen noch entfernt bzw. erklärt werden
    */

    /* Fläche  * 6 wegen Symmetrie (welche Fläche wird berechnet?) */
    for (int y = -i + 1; y < i; y++) {
      for (int x = -i + 1; x < i; x++) {
        sum += 3 * pow(-1, x + y + i) / sqrt(x * x + y * y + i * i);
        rest += 3 * pow(-1, x + y + i) / sqrt(x * x + y * y + i * i);
      }
    }

    /* Kanten * 12 wegen Symmetrie (welche Kante wird berechnet?) */
    for (int x = -i + 1; x < i; x++) {
      sum += 3 * pow(-1, x + i + i) / sqrt(x * x + i * i + i * i);
      rest += 9 * pow(-1, x + i + i) / sqrt(x * x + i * i + i * i);
    }

    /* Ecken * 8 wegen Symmetrie (welche Ecke wird berechnet?)*/
    sum += pow(-1, i + i + i) / sqrt(i * i + i * i + i * i);
    rest += 7 * pow(-1, i + i + i) / sqrt(i * i + i * i + i * i);
  }

  return sum;
}
