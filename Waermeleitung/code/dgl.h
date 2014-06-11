#ifndef _DGL_H
#define _DGL_H

typedef struct {
  double (*func)(double x, double y, void *args);
  void *args;
} function;

/* Quellen & Senkenfunktion f, Randbedingung Phi, Gitterabstand h und Lösungsarray sol */
int Poisson(function f, function phi, double h, double *sol);

#endif
