/* gcc -o numerik_3 -O2 numerik_bespin_deutsch_mathfunctions.c numerik_bespin_deutsch_3.c -lm */
/* Christian Bespin, Christopher Deutsch */

/* Programmaufruf: ./numerik_3 n_start, n_stop, n_step */
/* Das Programm gibt T_E(n) von n = n_start bis n = n_stop in Schritten von n_step aus */
/* Alternativ: ./numerik_3 n fuer den Wert von T_E(n) */

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

/* CO2 Konzentration der Atmosphaere */
const double k_N0 = 7.3E+25;

/* Strahlungseigenschaften der Sonne */
const double k_T_sun = 5750;
const double k_emissivity_sun = 5.5E-6;


/* spektrale Energiedichte */
double rho(double nu, double T);

/* Wirkungsquerschnitt (Breit-Wigner) */
double sigma(double nu);

/* analytisch integrierte Energiedichte */
double Z(double T);

/* Integrand der Emissivitaet: *args = {double T, double n} */
double emissivity_integrand(double nu, void *args);

/* Berechnet die Emissivitaet (das Integral ist in den Grenzen a, b um an die
 * Berechnung an den jeweiligen Parametersatz T,n anpassen zu koennen) */
double emissivity(double T, double n, double a, double b);

/* Die zu loesende nichtlineare Gleichung fuer das thermische Gleichgewicht.
 * Argumente: *args = {double n, double rhs}
 * Dabei ist rhs die rechte Seite der Gleichung (2) (PDF) 
 * Die Funktion liefert den Wert von lhs - rhs, sodass nur die Nullstelle
 * gefunden werden muss */
double equilibrium_eqn(double TE, void *args);


int main(int argc, char **argv) {
  double n_start, n_stop, n_step;
  double equi_temp;
  
  /* struct der zu loesenden Gleichung */
  function F;
  double args[2];
  
  
  /* ---Auswertung der Programm-Parameter--- */
  if (argc != 4 && argc != 2) {
    printf("\nBenutzung:\n");
    printf("Berechnung der Gleichgewichtstemperatur fuer ein n:\n"
           "%s n\n\n", argv[0]);
    printf("Tabelle von n = start bis stop mit step Schrittbreite:\n"
           "%s start stop step\n", argv[0]);
    return 1;
  }
  
  sscanf(argv[1], "%lf", &n_start);
  if (n_start > 100) {
    printf("n_start sollte kleiner als 100 sein\n");
    return 1;
  }
  
  if (argc == 2) {
    n_stop = n_start;
    n_step = 1;
  } else {
    sscanf(argv[2], "%lf", &n_stop);
    sscanf(argv[3], "%lf", &n_step);
    if (n_start < 0) {
      printf("Ungueltiger Startwert\n");
      return 1;
    }
    if (n_step <= 0) {
      printf("Ungueltige Schrittweite\n");
      return 1;
    }
    if (n_stop < n_start) {
      printf("n_stop sollte groesser sein als n_start\n");
      return 1;
    }
    if (n_stop > 100) {
      printf("# Warnung: n_stop > 100; setze n_stop = 100\n");
      n_stop = 100;
    }
  }
  /* ---Auswertung der Programm-Parameter--- */
  
  
  F.func = equilibrium_eqn;
  F.args = args;
  
  printf("# n\tT_E [K]\n");  
  while (n_start <= n_stop) {
    args[0] = n_start;
    /* rechte Seite der Gleichung (2) (PDF); unabhaengig von T_E, daher
     * Berechnung ausserhalb der Bisektion */
    args[1] = k_emissivity_sun * (2 - emissivity(k_T_sun, args[0], 1E13, 1E15))
                 * pow(k_T_sun, 4);
    
    find_root(F, 250, 350, 0.01, &equi_temp);
    printf("%.3f\t%.3f\n", args[0], equi_temp);
    
    n_start += n_step;
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

double emissivity_integrand(double nu, void *args) {
  /* Extrahiert die Argumente */
  double T = *(double*)args;
  double n = *((double*)args + 1);
  
  /* Wahrscheinlichkeit fuer das Entweichen eines Photons */
  double f = exp(-n * k_N0 * sigma(nu));
  
  return rho(nu, T) * (1 - f);
}

double emissivity(double T, double n, double a, double b) {
  function F;
  double args[2];
  
  /* in den Integranden zu substituierende Groessen */
  args[0] = T;
  args[1] = n;
  
  F.func = emissivity_integrand;
  F.args = args;
  
  return integrate(F, a, b, 1E-11, 20) / Z(T);
}

double equilibrium_eqn(double TE, void *args) {
  double n = *(double*)args;
  double rhs = *((double*)args + 1);
  
  return (2 - emissivity(TE, n, 1E12, 1E14)) * pow(TE, 4) - rhs;
}
