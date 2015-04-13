#ifndef _ODE_SOLVER_H
#define _ODE_SOLVER_H

/* NOTATION:
 * y_i'[t] = f_i(t, y_0, y_1, ..., y_n)  fuer i = 0, ..., n */

/* Eine Differentialgleichung eines Systems von n-Differentialgleichungen. 
 * Die Funktion dydt implementiert die rechte Seite von y_i'[t] (vgl. Notation).
 * Die y_0, y_1, ..., y_n-Abhaengigkeiten werden ueber das Array durch y[0],
 * y[1], ..., y[n] implementiert. Ausserdem koennen zusaetzliche Parameter in
 * "params" uebergeben werden. */
typedef struct {
  double (*dydt)(double t, double *y, double *params);
  double *params;
} ODE;

/* Repraesentation eines Differentialgleichungssystem aus "dimension"-Glei-
 * chungen. "eqns" ist das Array von "ODE" Objekten. */
typedef struct {
  ODE *eqns;
  int dimension;
} ODE_SYSTEM;

/* Temporaere Arrays fuer das Runge-Kutta-Verfahren 4. Ordnung, um staendige
 * Allokierungen bei der Entwicklung des DGL-Sys. zu vermeiden. */
typedef struct {
  int dimension;
  
  double *k_1;
  double *k_2;
  double *k_3;
  double *k_4;
  double *y;
} RK4_WORKSPACE;

/* Struktur zum Speichern der Loesung:
 * dimension: Anzahl der Gleichungen und Funktionen des Differentialgleichungs-
 *            systems
 * t_count: Anzahl der berechneten Loesungspunkte
 * t: Werte der Veraenderlichen
 * y: Funktionswerte der Funktionen zu dem jeweiligen Wert der Veraenderlichen
 *    in "t" (y[i][j] ist der Wert der i-ten Funktion zum "Zeitpunkt" t[j]) */
typedef struct {
  int dimension;
  int t_count;
  
  /* Werte der Veraenderlichen */
  double *t;
  /* Werte der Funktionen */
  double **y;
} ODE_SOLUTION;


/* Allokiert den Workspace zum Entwickeln der DGL */
RK4_WORKSPACE *rk4_workspace_alloc(int dimension);

/* Gibt den Speicher des Workspaces wieder frei */
void rk4_workspace_free(RK4_WORKSPACE *workspace);

/* Entwickelt das System fuer die Anfangsbedingung (y0, t0) auf (y1, t0 + h) */
void rk4_evolve(ODE_SYSTEM *system, RK4_WORKSPACE *space,
               double *y0, double t0, double h, double *y1);

/* Allokiert Speicher fuer die Loesung eines Differentialgleichungssystem mit
 * "dimension"-Gleichungen und "t_count" Loesungspunkten */
ODE_SOLUTION *ode_solution_alloc(int dimension, int t_count);

/* Gibt den Speicher der Loesung wieder frei */
void ode_solution_free(ODE_SOLUTION *sol);

/* Loest das Differentialgleichungssystem von "t0" bis "t1" in Schritten von "h"
 * mit der Anfangsbedingung "y0" und gibt die Loesung zurueck. */
ODE_SOLUTION *rk4_solve(ODE_SYSTEM *system, double *y0,
                        double t0, double t1, double h);

#endif
