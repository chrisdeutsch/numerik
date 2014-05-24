#ifndef _MATHFUNCTIONS_H_
#define _MATHFUNCTIONS_H_

/* Funktionen mit beliebigen Argumenten */
typedef struct {
  double (*func)(double x, void *args);
  void *args;
} function;

/* Berechnet die Nullstelle einer stetigen Funktion f auf dem geklammerten
 * Intervall [x1,x2] mit absolutem Fehler <= epsilon und speichert diese in 
 * root */
int find_root(function f, double x1, double x2, double epsilon, double *root);

/* Integriert die Funktion f rekursiv mit der adaptiven Simpsonmethode
 * von a bis b.
 * rdepth: maximale Rekursionstiefe 
 * epsilon: Abbruchbedingung (Abweichung von der letzten Naeherung) */
double integrate(function f, double a, double b, double epsilon, int rdepth);

/* Rekursionsfunktion zur Integration */
double recursive_simpson(function f, double prev_simp,
                         double x0, double x1, double x2,
                         double f0, double f1, double f2,
                         double epsilon, int rdepth);

#endif
