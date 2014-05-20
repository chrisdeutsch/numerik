#include "numerik_bespin_deutsch_mathfunctions.h"
#include <math.h>
#include <stdio.h>

int find_root(function f, double x1, double x2, double epsilon, double *root) {
  /* Es wird davon ausgegangen, dass eine Nullstelle mit Vorzeichenwechsel
   * im Intervall [x1,x2] vorliegt und die Funktion stetig ist. Es werden die
   * Funktionswerte am Rand und in der Mitte des Intervalls berechnet */
  double f1, f2;
  double mid, fmid;
  
  f1 = f.func(x1, f.args);
  f2 = f.func(x2, f.args);
  
  /* Ueberpruefung ob das Intervall sinnvoll gewaehlt wurde */
  if (f1 * f2 > 0) {
    return 1;
  }
  
  /* Bisektion solange bis die Nullstelle mit einem absoluten Fehler <= epsilon
   * bestimmt wurde. */
  while (fabs(f1 - f2) > epsilon) {
    mid = 0.5 * (x1 + x2);
    fmid = f.func(mid, f.args);
    
    /* Bisektion: das Intervall wurde halbiert und der Funktionswert der Mitte
     * berechnet. Einige Faelle (Rest analog):
     * Fall (+,+,-): Die Nullstelle liegt im rechten Teilintervall 
     * Fall (+,-,-): Die Nullstelle liegt im linken Teilintervall */
    if (f1 * fmid > 0) {
      x1 = mid;
      f1 = fmid;
    } else {
      x2 = mid;
      f2 = fmid;
    }
  }
  
  *root = 0.5 * (x1 + x2);
  return 0;
}

double integrate(function f, double a, double b, double epsilon, int rdepth) {
  /* Halbiert das Intervall [a,b] zum ersten Mal und berechnet Funktionswerte 
   * an den Stuetzstellen */
  double h = 0.5 * (b - a);
  double mid = a + h;
  
  double fa = f.func(a, f.args);
  double fmid = f.func(mid, f.args);
  double fb = f.func(b, f.args);
  
  /* Berechnet die erste Naeherung nach Simpson */
  double simp = (fa + 4*fmid + fb) * h / 3;
  
  /* Leitet die Rekursion ein. Dabei werden alle schon berechneten Funktions-
   * werte weitergegeben, damit die Anzahl der Funktionsaufrufe minimiert wird.
   * Es wird auch der aktuelle Wert des Integrals ueberreicht, damit verglichen
   * werden kann, welchen Einfluss die Verfeinerung auf das Ergebnis hat */
  return recursive_simpson(f, simp, a, mid, b, fa, fmid, fb, epsilon, rdepth);
}

double recursive_simpson(function f, double prev_simp,
                         double x0, double x1, double x2,
                         double f0, double f1, double f2,
                         double epsilon, int rdepth) {
  rdepth--;
  /* Wir haben zwei Intervalle [x0,x1] und [x1,x2] mit den Funktionswerten f0, 
   * f1 und f2. Es wird die naechste Verfeinerung nach Simpson berechnet. */
  double h = 0.5 * (x1 - x0);
  /* Berechne den Beitrag des Intervalls [x0,x1] zum Integral, dazu muss die
   * Stuetzstelle in der Haelfte des Intervalls berechnet werden: */
  double x01 = h + x0;
  double f01 = f.func(x01, f.args);
  
  double simp01 = (f0 + 4*f01 + f1) * h / 3;
  
  /* Berechne den Beitrag des Intervalls [x1,x2] zum Integral: */
  double x12 = h + x1;
  double f12 = f.func(x12, f.args);
  
  double simp12 = (f1 + 4*f12 + f2) * h / 3;
  
  /* Neue Abschaetzung des Integrals: */
  double simp = simp01 + simp12;
  
  /* Vergleich des neuen Schaetzwertes mit dem Alten. Sollte die Differenz
   * groesser als epsilon sein, werden die Intervalle halbiert und der Prozess
   * erneut auf die halbierten Intervalle durchgefuehrt. */
  if (rdepth > 0 && fabs(prev_simp - simp) > epsilon) {
    return recursive_simpson(f, simp01, x0, x01, x1, f0, f01, f1, 0.5*epsilon, rdepth)
           + recursive_simpson(f, simp12, x1, x12, x2, f1, f12, f2, 0.5*epsilon, rdepth);
  } else {
    return simp;
  }
}
