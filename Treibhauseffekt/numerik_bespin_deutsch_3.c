/* gcc -o numerik_3 -O2 numerik_bespin_deutsch_mathfunctions.c numerik_bespin_deutsch_3.c */
/* Christian Bespin, Christopher Deutsch */

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

/* Wirkungsquerschnitt */
double sigma(double nu);
/* Irgendwas-Energiedichte */
double rho(double nu, double T);
/* analytisch integrierte Energiedichte */
double Z(double T);

/* Integrand der Emissivitaet: *args = {dbl T, dbl n} */
double epsilon_integrand(double nu, void *args);
/* Berechnet die Emissivitaet (das Integral ist in den Grenzen a, b um an die
 * Berechnung an den jeweiligen Parametersatz T,n anpassen zu koennen) */
double epsilon(double T, double n, double a, double b);

/* Die eigentlich zu loesende nichtlineare Gleichung: *args = {dbl n, dbl rhs} */
double eqn(double TE, void *args);

int main() {
  double root;
  int i;
  
  /* Verpackt die zu loesende Gleichung */
  function F;
  double args[2];
  F.func = eqn;
  
  for (i = 1; i <= 1000; i++) {
    args[0] = i / 10.0;
    args[1] = c_epsilon_s * (2 - epsilon(c_TS, args[0], 1E12, 1E15)) * pow(c_TS, 4);
    
    F.args = args;
    find_root(F, 250, 350, 0.1, &root);
    printf("%f\t%f\n", args[0], root);
  }
  
  return 0;
}

double rho(double nu, double T) {
  return 8 * c_pi * c_h * pow(nu / c_c, 3) / ( exp(c_h / c_k * nu / T ) -1 );
}

double sigma(double nu) {
  double gamma_sqred = c_gamma * c_gamma;
  double nu2_diff = nu - c_nu2;
  double nu3_diff = nu - c_nu3;
  
  return c_gamma / c_pi * (c_S2 / (nu2_diff * nu2_diff + gamma_sqred)
                           + c_S3 / (nu3_diff * nu3_diff + gamma_sqred));

}

double Z(double T) {
  /* Form der analytischen Gleichung mit nur einem pow(x, 4)-Aufruf */
  return 8.0 / 15.0 * c_pi * c_h * c_c * pow(c_pi * c_k * T / c_h / c_c, 4);
}

double epsilon_integrand(double nu, void *args) {
  /* Extrahiert die Argumente */
  double T = *(double*)args;
  double n = *((double*)args + 1);
  
  /*  */
  double f = exp(-n * c_n0 * sigma(nu));
  
  return rho(nu, T) * (1 - f);
}

double epsilon(double T, double n, double a, double b) {
  function F;
  double args[2];
  /* in den Integranden zu substituierende Groessen */
  args[0] = T;
  args[1] = n;
  
  /* Wrapper fuer Funktion + Argumente */
  F.func = epsilon_integrand;
  F.args = args;
  
  /* Die gewaehlten Grenzen sind gut fuer T = 300 K */
  return integrate(F, a, b, 1E-11, 20) / Z(T);
}

double eqn(double TE, void *args) {
  double n = *(double*)args;
  double rhs = *((double*)args + 1);
  
  return (2 - epsilon(TE, n, 1E12, 1E14)) * pow(TE, 4) - rhs;
}
