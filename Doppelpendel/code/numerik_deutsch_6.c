/* Christopher Deutsch */
/* gcc -o numerik_6 -O2 numerik_deutsch_ode_solver.c numerik_deutsch_fft.c numerik_deutsch_6.c -lm */

/* Verwendung: Ausfuehrliche Erklaerung, wenn das Programm ohne Argumente aufgerufen wird */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include "numerik_deutsch_ode_solver.h"
#include "numerik_deutsch_fft.h"

/* Funktionsargumente:
 * theta_1 = y[0], omega_1 = y[1]
 * theta_2 = y[2], omega_2 = y[3]
 * Parameter:
 * g = params[0], mu = params[1], L_1 = params[2], L_2 = params[3] */
double d_theta_1(double t, double *y, double *params);
double d_omega_1(double t, double *y, double *params);
double d_theta_2(double t, double *y, double *params);
double d_omega_2(double t, double *y, double *params);

/* Speichert die Loesung in der Datei "filename" und berechnet die verallgemei-
 * nerten Impulse der Koordinaten sowie die kartesischen Koordinaten fuer die
 * Punktmassen */
void save_ode_solution(ODE_SOLUTION *solution, double m1, double m2,
                       double L1, double L2, char *filename);

/* Berechnet das Leistungsspektrum der DFT (Betraege der Fourierkoeffizienten)
 * und die zu jedem Koeffizienten korrespondierende Kreisfrequenz */
void save_power_spectrum(int n, double delta,
                         double complex *f1, double complex *f2,
                         char *filename);

int main(int argc, char **argv) {
  int i;
  
  /* Zaehlt die eingelesenen Argumente */
  int input_cnt;
  
  /* physikalische Parameter (Default-Werte) */
  double g = 9.81;
  double m1 = 1.0;
  double m2 = 1.0;
  double L1 = 1.0;
  double L2 = 1.0;
  
  /* Laenge der Zeitentwicklung und Zeitschritt */
  double t;
  double h;
  
  /* FFT-Variablen */
  int n, r;
  double complex *f1, *f2;
  
  /* Anfangswert des DGL-Sys. */
  double y0[4];
  
  /* DGL-Sys. definieren: */
  ODE_SOLUTION *solution;
  ODE_SYSTEM system;
  ODE eqns[4];
  double params[4];
  
  eqns[0].dydt = d_theta_1;
  eqns[1].dydt = d_omega_1;
  eqns[2].dydt = d_theta_2;
  eqns[3].dydt = d_omega_2;
  
  system.eqns = eqns;
  system.dimension = 4;
  
  /* Parameter in den DGLn setzen */
  for (i = 0; i < 4; i++) {
    eqns[i].params = params;
  }
  
  /* Programmargumente einlesen und verarbeiten */
  input_cnt = 0;
  if ((argc != 7) && (argc != 12)) {
    printf("Benutzung:\n"
           "%s t h theta1 omega1 theta2 omega2 (g m1 m2 L1 L2)\n\n"
           "Die Klammer enthaelt optionale Argumente\n"
           "t: Laenge der Zeitentwicklung\n"
           "h: Zeitschritt\n"
           "theta/omega: Anfangsbedingung des Pendels (Bezeichnung PDF)\n"
           "g: Gravitationsbeschleunigung\n"
           "m: Massen\n"
           "L1: Laengen der Pendelstange\n\n"
           "Beispielaufruf:\n"
           "%s 10 0.001 3.14 0 3.14 0\nEntwicklung ueber 10 s in Schritten von 0.001 s mit den beiden\n"
           "Winkeln gleich Pi und ohne anfaengliche Winkelgeschwindigkeit)\n\n",
           argv[0], argv[0]);
    return -1;
  }
  input_cnt += sscanf(argv[1], "%lf", &t);
  input_cnt += sscanf(argv[2], "%lf", &h);
  input_cnt += sscanf(argv[3], "%lf", y0);
  input_cnt += sscanf(argv[4], "%lf", y0+1);
  input_cnt += sscanf(argv[5], "%lf", y0+2);
  input_cnt += sscanf(argv[6], "%lf", y0+3);
  
  if (input_cnt != 6) {
    printf("Es konnten nicht alle Argumente eingelesen werden\n");
    return -1;
  }
  
  /* Optionale Argumente einlesen */
  if (argc == 12) {
    input_cnt += sscanf(argv[7], "%lf", &g);
    input_cnt += sscanf(argv[8], "%lf", &m1);
    input_cnt += sscanf(argv[9], "%lf", &m2);
    input_cnt += sscanf(argv[10], "%lf", &L1);
    input_cnt += sscanf(argv[11], "%lf", &L2);
    
    if (input_cnt != 11) {
      printf("Es konnten nicht alle Argumente eingelesen werden\n");
      return -1;
    }
  }
  
  /* Ausgabe welche Parameter zum Berechnen verwendet werden */
  printf("Verwendeter Parametersatz:\n"
         "Laenge der Zeitentwicklung und Zeitschritt:\n"
         "t = %.3f s;  h = %f s\n\n"
         "Anfangswerte:\n"
         "theta1 = %.4f rad;  omega1 = %.4f rad/s\n"
         "theta2 = %.4f rad;  omega2 = %.4f rad/s\n\n"
         "Optionale Parameter:\n"
         "g = %.3f m s^-2;  m1 = %.3f kg;  m2 = %.3f kg\n"
         "L1 = %.3f m;  L2 = %.3f m\n\n",
         t, h, y0[0], y0[1], y0[2], y0[3], g, m1, m2, L1, L2);
  
  /* Die eingelesenen Parameter werden in die Parameterliste der Differential-
   * gleichung eingetragen */
  params[0] = g;
  params[1] = m2 / (m1 + m2);
  params[2] = L1;
  params[3] = L2;
  
  /* Loesen des DGL-Sys. */
  printf("Loesen des Differentialgleichungssystems...\n");
  solution = rk4_solve(&system, y0, 0, t, h);
  
  /* Berechnete Loesung speichern (hier werden noch einige Berechnungen durch-
   * gefuehrt, wie der verallgemeinerte Impuls, Trajektorie in kartesischen
   * koordinaten etc.) */
  printf("Speichern der Loesung in \"numerik_deutsch_ode_solution.txt\"...\n");
  save_ode_solution(solution, m1, m2, L1, L2, "numerik_deutsch_ode_solution.txt");
  
  
  /* #### FFT #### */
  
  /* Berechnet die maximale Anzahl n = 2^r an Datenpunkten fuer eine
   * Radix-2-FFT */
  r = 0, n = 1;
  while ((n <<= 1) <= solution->t_count) r++;
  n = 1 << r;
  
  /* Arrays komplexer Zahlen fuer die Fourierkoeffizienten */
  f1 = malloc(n * sizeof(double complex));
  f2 = malloc(n * sizeof(double complex));
  if (f1 == NULL || f2 == NULL) {
    printf("Speicher fuer die Fourierkoeffizienten konnte nicht allokiert werden\n");
    return -1;
  }
  
  /* Schreibe die Loesung fuer theta1 und theta2 in f1 / f2 */
  for (i = 0; i < n; i++) {
    f1[i] = solution->y[0][i];
    f2[i] = solution->y[2][i];
  }
  
  /* Berechnet die DFT der ersten 2^r Datenpunkte fuer theta1 und theta2 */
  printf("Berechnen der FFT...\n");
  if (fft(r, f1) == FFT_ALLOC_ERROR ||
      fft(r, f2) == FFT_ALLOC_ERROR) {
    printf("Speicher zur FFT konnte nicht allokiert werden\n");
    return -1;
  }
  
  /* Speichert das Leistungsspektrum der DFT */
  printf("Speichern der DFT in \"numerik_deutsch_power_spectrum.txt\"...\n");
  save_power_spectrum(n, h, f1, f2, "numerik_deutsch_power_spectrum.txt");
  
  /* Speicher wieder freigeben */
  ode_solution_free(solution);
  free(f1);
  free(f2); 
  
  return 0;
}

double d_theta_1(double t, double *y, double *params) {
  double omega_1 = y[1];
  return omega_1;
}

double d_omega_1(double t, double *y, double *params) {
  double theta_1 = y[0];
  double omega_1 = y[1];
  double theta_2 = y[2];
  double omega_2 = y[3];
  
  double g = params[0];
  double mu = params[1];
  double L_1 = params[2];
  double L_2 = params[3];
  
  double numerator, denominator;
  
  double cosdiff = cos(theta_1 - theta_2);
  
  
  numerator = (0.5 * mu - 1) * g * sin(theta_1)
              - 0.5 * mu * g * sin(theta_1 - 2 * theta_2)
              - mu * sin(theta_1 - theta_2) * (L_1 * omega_1 * omega_1 * cosdiff
                                               + L_2 * omega_2 * omega_2);
  
  denominator = L_1 * (1 - mu * cosdiff * cosdiff);
  
  return numerator / denominator;
}

double d_theta_2(double t, double *y, double *params) {
  double omega_2 = y[3];
  return omega_2;
}

double d_omega_2(double t, double *y, double *params) {
  double theta_1 = y[0];
  double omega_1 = y[1];
  double theta_2 = y[2];
  double omega_2 = y[3];
  
  double g = params[0];
  double mu = params[1];
  double L_1 = params[2];
  double L_2 = params[3];
  
  double numerator, denominator;
  
  double cosdiff = cos(theta_1 - theta_2);
  
  numerator = sin(theta_1 - theta_2) * (g * cos(theta_1)
                                        + L_1 * omega_1 * omega_1
                                        + mu * L_2 * omega_2 * omega_2 * cosdiff);
  
  denominator = L_2 * (1 - mu * cosdiff * cosdiff);
  
  return numerator / denominator;
}

void save_ode_solution(ODE_SOLUTION *solution, double m1, double m2,
                       double L1, double L2, char *filename) {
  int i;
  
  /* Anzahl der Datenpunkte */
  int t_count = solution->t_count;
  
  double *t = solution->t;
  double **y = solution->y;
  
  /* verallgemeinerte Impulse */
  double p1, p2;
  double cosdiff;
  
  /* kartesische Koordinaten */
  double x1, x2, y1, y2;
  /* Sinus/Cosinus der Winkel */
  double s1, c1, s2, c2;
  
  FILE *file = fopen(filename, "w");
  if (file == NULL) {
    printf("Konnte die Datei %s nicht erstellen\n", filename);
    return;
  }
  
  /* Tabellenkopf */
  fprintf(file, "t\ttheta1[rad]\tomega1[rad/s]\ttheta2[rad]\tomega2[rad/s]\tp1[kg m^2 s^-1]\tp2[kg m^2 s^-1]\tx1[m]\ty1[m]\tx2[m]\ty2[m]\n");
  
  /* Loesung */
  for (i = 0; i < t_count; i++) {
    /* Berechne die verallgemeinerten Impulse */
    cosdiff = cos(y[0][i] - y[2][i]);
    p1 = (m1 + m2) * L1 * L1 * y[1][i] + m2 * L1 * L2 * y[3][i] * cosdiff;
    p2 = m2 * L1 * L2 * y[1][i] * cosdiff + m2 * L2 * L2 * y[3][i];
    
    /* Berechne die Koordinaten im kartesischen KS */
    s1 = sin(y[0][i]), c1 = cos(y[0][i]);
    s2 = sin(y[2][i]), c2 = cos(y[2][i]);
    
    x1 = L1 * s1;
    y1 = -L1 * c1;
    x2 = L1 * s1 + L2 * s2;
    y2 = -L1 * c1 - L2 * c2;
    
    /* Schreiben in die Datei */
    fprintf(file, "%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n",
            t[i], y[0][i], y[1][i], y[2][i], y[3][i], p1, p2, x1, y1, x2, y2);
  }
  
  fclose(file);
}

void save_power_spectrum(int n, double delta,
                         double complex *f1, double complex *f2,
                         char *filename) {
  int i;
  FILE *file;
  
  file = fopen(filename, "w");
  if (file == NULL) {
    printf("Konnte die Datei %s nicht erstellen\n", filename);
    return;
  }
  
  /* Tabellenkopf */
  fprintf(file, "k\tomega[rad/s]\t|g_k(theta_1)|[rad]\t|g_k(theta_2)|[rad]\n");
  
  for (i = 0; i < n; i++) {
    fprintf(file, "%i\t%f\t%f\t%f\n",
            i, 2 * fft_pi * i/(n * delta), cabs(f1[i]), cabs(f2[i]));
  }
  
  fclose(file);
}