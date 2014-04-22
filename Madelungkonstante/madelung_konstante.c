/* gcc -O2 -std=c99 -Wall -pedantic -o madelung madelung_konstante.c -lm */
/* Christian Bespin, Christopher Deutsch */

#include <stdio.h>
#include <math.h>

/* Berechnet die Madelungkonstante für einen Würfel der Kantenlänge 2n wobei die Länge 
   auf den Gitterabstand normalisiert ist */
double madelung_3d(int n);

/* Berechnet die Madelungkonstante für ein Quadrat der Seitenlänge 2n wobei die Länge
   auf den Gitterabstand normalisiert ist */
double madelung_2d(int n);


/* Beschreibung was eigentlich in der main-Funktion berechnet wird */
int main() {
  /* Literaturwert nach Y. Sakamoto */
  const double lit = 1.7475645946331822;
  double mad;

  for (int i = 1; i < 10; i++) {
    mad = -madelung_3d(i);
    printf("Madelung-Konstante: %.12f \t n: %i \t Delta_rel: %.2E\n", mad, i,
           fabs(mad - lit) / lit);
  }
  
  for (int i = 1; i < 10; i++) {
    mad = -madelung_3d(10 * i);
    printf("Madelung-Konstante: %.12f \t n: %i \t Delta_rel: %.2E\n", mad, 10 * i,
           fabs(mad - lit) / lit);
  }
  
  printf("Madelungkonstante: %f\n", -madelung_2d(100));
  
  return 0;
}

double madelung_3d(int n) {
  double sum = 0;
  double rest = 0;
  
  /* effizient der pow-Funktion diskutieren (könnte durch eine Modulo Operation eliminiert werden
     macht aber den Code unklarer */
  for (int m = 1; m <= n; m++) {
    /* Den Rest vom letzten Durchgang aufaddieren und danach resetten */
    sum += rest;
    rest = 0;

    /*
    magic numbers müssen noch entfernt bzw. erklärt werden
    */

    /* Fläche  * 6 wegen Symmetrie (welche Fläche wird berechnet?) */
    for (int y = -m + 1; y < m; y++) {
      for (int x = -m + 1; x < m; x++) {
        sum += 3 * pow(-1, x + y + m) / sqrt(x * x + y * y + m * m);
        rest += 3 * pow(-1, x + y + m) / sqrt(x * x + y * y + m * m);
      }
    }

    /* Kanten * 12 wegen Symmetrie (welche Kante wird berechnet?) */
    for (int x = -m + 1; x < m; x++) {
      sum += 3 * pow(-1, x + m + m) / sqrt(x * x + m * m + m * m);
      rest += 9 * pow(-1, x + m + m) / sqrt(x * x + m * m + m * m);
    }

    /* Ecken * 8 wegen Symmetrie (welche Ecke wird berechnet?)*/
    sum += pow(-1, m + m + m) / sqrt(m * m + m * m + m * m);
    rest += 7 * pow(-1, m + m + m) / sqrt(m * m + m * m + m * m);
  }

  return sum;
}


double madelung_2d(int n) {
  double sum = 0;
  double rest = 0;

  /* effizient der pow-Funktion diskutieren (könnte durch eine Modulo Operation
     eliminiert werden
     macht aber den Code unklarer */
  for (int m = 1; m <= n; m++) {
    /* Den Rest vom letzten Durchgang aufaddieren und danach resetten */
    sum += rest;
    rest = 0;

    /*
    magic numbers müssen noch entfernt bzw. erklärt werden
    */

    /* Kante  * 4 wegen Symmetrie (welche Kante wird berechnet?) */
    for (int x = -m + 1; x < m; x++) {
      sum += 2 * pow(-1, x + m) / sqrt(x * x + m * m);
      rest += 2 * pow(-1, x + m) / sqrt(x * x + m * m);
    }

    /* Ecken * 4 wegen Symmetrie (welche Ecke wird berechnet?)*/
    sum += pow(-1, m + m) / sqrt(m * m + m * m);
    rest += 3 * pow(-1, m + m) / sqrt(m * m + m * m);
  }

  return sum;
}
