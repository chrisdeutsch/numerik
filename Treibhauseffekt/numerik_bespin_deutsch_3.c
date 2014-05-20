/* gcc -o numerik_3 -O2 numerik.c main.c */
/* Namen */

#include <stdio.h>
#include <math.h>
#include "numerik_bespin_deutsch_mathfunctions.h"

/* physikalische Konstanten */
const double c_h = 6.62606957E-34;
const double c_k = 1.3806488E-23;
const double c_c = 299792458;
const double c_pi = 3.1415926535897932384626433832795;

const double c_S2 = 2.45E-11;
const double c_S3 = 2.74E-10;
const double c_gamma = 3.00E+10;
const double c_nu2 = 2.00E+13;
const double c_nu3 = 7.04E+13;

const double c_n0 = 7.3E+25;
const double c_TS = 5750;
const double c_epsilon_s = 5.5E-6;


double sigma(double nu);
double rho(double nu, double T);
/* analytisch */
double Z(double T);

/*  */
double f(double nu, double n_schl);
double diff(double nu, void *args);
double epsilon(double T, double n_schl);
double eqn(double TE, void *args);

int main() {
  double root;

  double args[2];
  
  function F;
  F.func = eqn;
  
  int i;
  
  for (i = 1; i <= 1000; i++) {
    args[0] = i / 10.0;
    args[1] = c_epsilon_s * (2 - epsilon(c_TS, args[0])) * pow(c_TS, 4);
    
    F.args = args;
    find_root(F, 250, 350, 0.1, &root);
    printf("%f\t%f\n", args[0], root);
  }
  
  return 0;
}

double rho(double nu, double T) {
  return 8 * c_pi * c_h / pow(c_c, 3) * pow(nu, 3) / ( exp(c_h / c_k * nu / T ) -1 );
}

double sigma(double nu) {
  return c_S2 / c_pi * c_gamma / ((nu - c_nu2) * (nu - c_nu2) + c_gamma * c_gamma)
         + c_S3 / c_pi * c_gamma / ((nu - c_nu3) * (nu - c_nu3) + c_gamma * c_gamma);
}

double Z(double T) {
  return 8 * c_pi * c_h / pow(c_c, 3) * pow(c_pi, 4) / (15 * pow(c_h / c_k / T, 4));
}

double f(double nu, double n_schl) {
  return exp(-n_schl * c_n0 * sigma(nu));
}

double diff(double nu, void *args) {
  double T = *(double*)args;
  double n_schl = *((double*)args + 1);
  return rho(nu, T) * (1 - f(nu, n_schl));
}

double epsilon(double T, double n_schl) {
  double args[2] = {T, n_schl};
  
  function F;
  F.func = diff;
  F.args = args;
  
  /* Die gewaehlten Grenzen sind gut fuer T = 5750 sowie T im 300er Bereich */
  return integrate(F, 1E12, 1E15, 1E-11, 20) / Z(T);
}

double eqn(double TE, void *args) {
  double n_schl = *(double*)args;
  double rhs = *((double*)args + 1);
  return (2 - epsilon(TE, n_schl))*pow(TE, 4) - rhs;
}
