/* gcc -O2 -std=c99 -Wall -pedantic -o madelung madelung_konstante.c -lm */
/* Christian Bespin, Christopher Deutsch */

#include <stdio.h>
#include <math.h>

/* Literaturwert nach Y. Sakamoto */
#define MAD_CONST_3D 1.7475645946331822
/* Literaturwert von OEIS */
#define MAD_CONST_2D 1.6155426267128247

/* Berechnet die Madelungkonstante fuer einen Wuerfel; bricht ab wenn die
 * Aenderung vom vorigen Durchlauf kleiner als epsilon ist */
double madelung_3d(double epsilon);

/* Berechnet die Madelungkonstante für ein Quadrat; bricht ab wenn die
 * Abweichung vom Literaturwert kleiner als epsilon ist */
double madelung_2d(double epsilon);

/* Abstand der Gitterstelle (x,y,z) vom Ursprung (Längen sind mit dem Gitter-
 * abstand normalisiert) */
double dist(int x, int y, int z);

/* Berechnet das Vorzeichen des Ions an der Gitterstelle (x,y,z) (Es wird
 * angenommen, dass wir die Madelungkonstante fuer ein Chlorion berechnen,
 * da wir so eine positive Madelung-Konstante erhalten) */
int sign_z(int x, int y, int z);

/* Beschreibung was eigentlich in der main-Funktion berechnet wird */
int main() {
  double mad;
  
  printf("#############################################\n");
  printf("Berechnung der Madelung-Konstante\n");
  printf("#############################################\n");
  printf("Im dreidimensionalen Fall: (Genauigkeit: 10E-5)\n");
  printf("Madelung-Konstante\trelativer Fehler\n");
  
  mad = madelung_3d(1E-5);
  
  printf("%.12f\t\t%.2E\n", mad, fabs(mad - MAD_CONST_3D)/MAD_CONST_3D);
  printf("#############################################\n");
  printf("Im zweidimensionalen Fall: (Genauigkeit: 10E-5) \n");
  printf("Madelung-Konstante\trelativer Fehler\n");
  
  mad = madelung_2d(1E-5);
  
  printf("%.12f\t\t%.2E\n", mad, fabs(mad - MAD_CONST_2D)/MAD_CONST_2D);
  
  return 0;
}

double madelung_3d(double epsilon) {
  /* In dieser Variable wird der Wert der Madelung-Konstante gespeichert */
  double mconst = 0;
  /* Diese Variable speichert den Rest der Gewichtung */
  double residual = 0;
  
  /* Speichert den Wert der Madelung-Konstante des vorigen Schleifendurchlaufs
   * damit die Aenderung seit dem Durchlauf berechnet werden kann
   * Der Startwert wurde so gewählt, dass die for-Schleife wenigstens einmal 
   * betreten wird */
  double prev = 1E100; 
  
  /* Diese Schleife zaehlt die aktuell hinzugefuegte Schale um das Zentralion.
   * Waehrend jedem Durchgang entspricht dies einem Wuerfel mit Kantenlaenge 2m
   * und Mittelpunkt im Ursprung (0,0,0). Die for-Schleife wird beendet wenn die
   * Differenz von vorigem Durchlauf zu aktuellen Durchlauf >= epsilon ist. */
  for (int m = 1; fabs(mconst - prev) > epsilon; m++) {
    /* Hier wird der die Madelung-Konstante des letzten Schleifendurchlaufs
     * zwischengespeichert */
    prev = mconst;
    /* Der Rest des letzten Durchlaufs wird aufaddiert, damit der Innenraum des
     * Wuerfels voll gewichtet wird. Danach wird der Rest für den aktuellen
     * Durchlauf wieder auf 0 gesetzt */
    mconst += residual;
    residual = 0;
    
    /* Ab hier wird unter Symmetriebetrachtung der Beitrag der einzelnen Stuecke
     * der neuen Schale auf die Madelung-Konstante aufsummiert. Dabei wird stets
     * der Rest der Gewichtung in der Variable residual aufsummiert */

    /* Flaechen (es wird die Flaeche (x,y,m) mit x,y = -m + 1, ... , m - 1
     * berechnet)
     * mconst: 6 (Sym.) * 1/2 (Evjens) = 3
     * residual: 6 (Sym.) * 1/2 (Rest) = 3 */
    for (int y = -m + 1; y < m; y++) {
      for (int x = -m + 1; x < m; x++) {
        mconst += 3 * sign_z(x, y, m) / dist(x, y, m);
        residual += 3 * sign_z(x, y, m) / dist(x, y, m);
      }
    }

    /* Kanten (es wird die Kante (x,m,m) mit x = -m + 1, ... , m - 1 berechnet)
     * mconst: 12 (Sym.) * 1/4 (Evjens) = 3
     * residual: 12 (Sym.) * 3/4 (Rest) = 9 */
    for (int x = -m + 1; x < m; x++) {
      mconst += 3 * sign_z(x, m, m) / dist(x, m, m);
      residual += 9 * sign_z(x, m, m) / dist(x, m, m);
    }

    /* Ecken (es wird die Ecke (m,m,m) berechnet)
     * mconst: 8 (Sym.) * 1/8 (Evjens) = 1
     * residual: 8 (Sym.) * 7/8 (Rest) = 7 */
    mconst += sign_z(m, m, m) / dist(m, m, m);
    residual += 7 * sign_z(m, m, m) / dist(m, m, m);
  }

  return mconst;
}

double madelung_2d(double epsilon) {
  /* In dieser Variable wird der Wert der Madelung-Konstante gespeichert */
  double mconst = 0;
  /* Diese Variable speichert den Rest der Gewichtung */
  double residual = 0;
  
  /* Diese Schleife vergroeßert die Seitenlaenge des Quadrat 
  	  bei jedem Durchgang um 2, bis mconst genau genug ist
  */
  for (int m = 1; fabs(mconst - MAD_CONST_2D) > epsilon; m++) {
    /* Den Rest vom letzten Durchgang aufaddieren und danach resetten */
    mconst += residual;
    residual = 0;

    /* Kante (hier Kante in x-Richtung mit y=m, z=0))
       mconst: 4 (Sym.) * 1/2 (Evjens) = 2
       residual: 4 (Sym.) * 1/2 (Rest) = 2
    */
    for (int x = -m + 1; x < m; x++) {
      mconst += 2 * sign_z(x, m, 0) / dist(x, m, 0);
      residual += 2 * sign_z(x, m, 0) / dist(x, m, 0);
    }

    /* Ecken * 4 wegen Symmetrie (hier Ecke (m, m, 0))
       mconst: 4 (Sym.) * 1/4 (Evjens) = 1
       residual: 4 (Sym.) * 3/4 (Rest) = 3
    */
    mconst += sign_z(m, m, 0) / dist(m, m, 0);
    residual += 3 * sign_z(m, m, 0) / dist(m, m, 0);
  }

  return mconst;
}

double dist(int x, int y, int z) {
  return sqrt(x * x + y * y + z * z);
}

/* Umsetzung der Gleichung (4) in der beiliegenden pdf-Datei
 * aequivalent zu -pow(-1, x + y + z) aber effizienter 
 * gerade Exponenten: -1
 * ungerade Exponenten: +1 */
int sign_z(int x, int y, int z) {
  if((x + y + z) % 2 == 0)
    return -1;
  else
    return 1;
}
