/* gcc -o numerik_3 -O2 numerik_bespin_deutsch_mathfunctions.c numerik_bespin_deutsch_3.c */
/* Christian Bespin, Christopher Deutsch */

#include <stdio.h>
#include <math.h>
#include "numerik_bespin_deutsch_mathfunctions.h"

/* physikalische Konstanten */
const double k_h = 6.62606957E-34;
const double k_kB = 1.3806488E-23;
const double k_c0 = 299792458;
const double k_pi = 3.1415926535897932384626433832795;

/* Konstanten der Breit-Wigner Formel */
const double k_S2 = 2.45E-11;
const double k_S3 = 2.74E-10;
const double k_gamma = 3.00E+10;
const double k_nu2 = 2.00E+13;
const double k_nu3 = 7.04E+13;

/* CO2 Gehalt der Atmosphaere */
const double k_n0 = 7.3E+25;

/* Sonnenparameter */
const double k_TS = 5750;
const double k_epsilon_s = 5.5E-6;

/* spektrale Energiedichte */
double rho(double nu, double T);

/* Wirkungsquerschnitt (Breit-Wigner) */
double sigma(double nu);

/* analytisch integrierte Energiedichte */
double Z(double T);

/* Integrand der Emissivitaet: *args = {double T, double n} */
double epsilon_integrand(double nu, void *args);

/* Berechnet die Emissivitaet (das Integral ist in den Grenzen a, b um an die
 * Berechnung an den jeweiligen Parametersatz T,n anpassen zu koennen) */
double epsilon(double T, double n, double a, double b);

/* Die zu loesende nichtlineare Gleichung fuer das thermische Gleichgewicht.
 * Argumente: *args = {double n, double rhs} */
double equilibrium_eqn(double TE, void *args);


int main(int argc, char **argv) {
  double equi_temp;
  int i;
  
  /* Verpackt die zu loesende Gleichung */
  function F;
  double args[2];
  F.func = equilibrium_eqn;
  F.args = args;
  
  for (i = 1; i <= 1000; i++) {
    args[0] = i / 10.0;
    args[1] = k_epsilon_s * (2 - epsilon(k_TS, args[0], 1E12, 1E15)) * pow(k_TS, 4);
    
    find_root(F, 250, 350, 0.01, &equi_temp);
    printf("%f\t%f\n", args[0], equi_temp);
  }
  
  return 0;
}

double rho(double nu, double T) {
  return 8 * k_pi * k_h * pow(nu / k_c0, 3) / ( exp(k_h / k_kB * nu / T ) -1 );
}

double sigma(double nu) {
  double gamma_sq= k_gamma * k_gamma;
  double nu2_diff = nu - k_nu2;
  double nu3_diff = nu - k_nu3;
  
  return k_gamma / k_pi * (k_S2 / (nu2_diff * nu2_diff + gamma_sq)
                           + k_S3 / (nu3_diff * nu3_diff + gamma_sq));

}

double Z(double T) {
  /* Form der analytischen Gleichung mit nur einem pow(x, 4)-Aufruf */
  return 8.0 / 15.0 * k_pi * k_h * k_c0 * pow(k_pi * k_kB * T / k_h / k_c0, 4);
}

double epsilon_integrand(double nu, void *args) {
  /* Extrahiert die Argumente */
  double T = *(double*)args;
  double n = *((double*)args + 1);
  
  /*  */
  double f = exp(-n * k_n0 * sigma(nu));
  
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

double equilibrium_eqn(double TE, void *args) {
  double n = *(double*)args;
  double rhs = *((double*)args + 1);
  
  return (2 - epsilon(TE, n, 1E12, 1E14)) * pow(TE, 4) - rhs;
}
