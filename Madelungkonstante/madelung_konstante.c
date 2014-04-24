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

/* Berechnet den Betrag des Vektors (x,y,z) */
double dist(int x, int y, int z);

/* Vorzeichenfunktion */
int sign_z(int x, int y, int z);


/* Beschreibung was eigentlich in der main-Funktion berechnet wird */
int main() {
  /* Literaturwert nach Y. Sakamoto */
  const double lit = 1.7475645946331822;
  double mad;

  for (int i = 1; i < 10; i++) {
    mad = madelung_3d(i);
    printf("Madelung-Konstante: %.12f \t n: %i \t Delta_rel: %.2E\n", mad, i,
           fabs(mad - lit) / lit);
  }
  
  for (int i = 1; i < 10; i++) {
    mad = madelung_3d(10 * i);
    printf("Madelung-Konstante: %.12f \t n: %i \t Delta_rel: %.2E\n", mad, 10 * i,
           fabs(mad - lit) / lit);
  }
  
  printf("Madelungkonstante: %f\n", madelung_2d(100));
  
  return 0;
}

double madelung_3d(int n) {
  double sum = 0; //wo ist der Unterschied zwischen sum und rest?
  double rest = 0;
  
  /* effizient der pow-Funktion diskutieren (könnte durch eine Modulo Operation eliminiert werden
     macht aber den Code unklarer */
  for (int m = 1; m <= n; m++) { //Madelungkonst von Würfel mit Kantenlänge 2 + Madelungkonst von Würfel mit Kantenlänge 4 + ...?
    /* Den Rest vom letzten Durchgang aufaddieren und danach resetten */
    sum += rest;
    rest = 0;

    /*
    magic numbers müssen noch entfernt bzw. erklärt werden
    */

    /* Fläche  * 6 wegen Symmetrie (welche Fläche wird berechnet?) */
    for (int y = -m + 1; y < m; y++) {
      for (int x = -m + 1; x < m; x++) {
        sum += 3 * sign_z(x, y, m) / dist(x, y, m);
        rest += 3 * sign_z(x, y, m) / dist(x, y, m);
      }
    }

    /* Kanten * 12 wegen Symmetrie (welche Kante wird berechnet?) */
    for (int x = -m + 1; x < m; x++) {
      sum += 3 * sign_z(x, m, m) / dist(x, m, m);
      rest += 9 * sign_z(x, m, m) / dist(x, m, m);
    }

    /* Ecken * 8 wegen Symmetrie (welche Ecke wird berechnet?)*/
    sum += sign_z(m, m, m) / dist(m, m, m);
    rest += 7 * sign_z(m, m, m) / dist(m, m, m);
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
    for (int x = -m + 1; x < m; x++) { //War das nicht Kanten 1/2 gewichtet?
      sum += 2 * sign_z(x, m, 0) / dist(x, m, 0);
      rest += 2 * sign_z(x, m, 0) / dist(x, m, 0);
    }

    /* Ecken * 4 wegen Symmetrie (welche Ecke wird berechnet?)*/
    sum += sign_z(m, m, 0) / dist(m, m, 0);
    rest += 3 * sign_z(m, m, 0) / dist(m, m, 0);
  }

  return sum;
}

double dist(int x, int y, int z) {
  return sqrt(x*x + y*y + z*z);
}

int sign_z(int x, int y, int z) {
  return -pow(-1, x + y + z);
}
