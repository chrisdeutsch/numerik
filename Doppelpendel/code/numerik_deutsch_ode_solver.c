#include "numerik_deutsch_ode_solver.h"
#include <stdlib.h>
#include <assert.h>

RK4_WORKSPACE *rk4_workspace_alloc(int dimension) {
  RK4_WORKSPACE *ret = malloc(sizeof(RK4_WORKSPACE));
  if (ret == NULL) return NULL;
  
  ret->dimension = dimension;
  
  /* Es wird nur ein Array allokiert und dementsprechend die Pointer gesetzt, um
   * 5 Arrays gleicher Groesse zu erhalten */
  ret->k_1 = malloc(5 * dimension * sizeof(double));
  if (ret->k_1 == NULL) {
    free(ret);
    return NULL;
  }
  
  ret->k_2 = ret->k_1 + dimension;
  ret->k_3 = ret->k_2 + dimension;
  ret->k_4 = ret->k_3 + dimension;
  ret->y = ret->k_4 + dimension;
  
  return ret;
}

void rk4_workspace_free(RK4_WORKSPACE *workspace) {
  free(workspace->k_1);
  free(workspace);
}

void rk4_evolve(ODE_SYSTEM *system, RK4_WORKSPACE *workspace,
               double *y0, double t0, double h, double *y1) {
  int i, j;
  
  double (*f_i)(double t, double *y, double *params);
  double *params;
  int dimension = system->dimension;
  
  /* Benoetigte temporaere Arrays */
  double *y_temp = workspace->y;
  double *k_1 = workspace->k_1;
  double *k_2 = workspace->k_2;
  double *k_3 = workspace->k_3;
  double *k_4 = workspace->k_4;
  
  assert(system->dimension == workspace->dimension);
  
  /* Berechnung von k_1 */
  for (i = 0; i < dimension; i++) {
    f_i = (system->eqns[i]).dydt;
    params = (system->eqns[i]).params;
    
    k_1[i] = h * f_i(t0, y0, params);
  }
  
  /* Berechnung von k_2 */
  for (i = 0; i < dimension; i++) {
    f_i = (system->eqns[i]).dydt;
    params = (system->eqns[i]).params;
    
    /* Berechne "y_temp = y0 + 0.5 * k_1" fuer alle Gleichungen */
    for (j = 0; j < dimension; j++) {
      y_temp[j] = y0[j] + 0.5 * k_1[j];
    }
    
    k_2[i] = h * f_i(t0 + 0.5 * h, y_temp, params);
  }
  
  /* Berechnung von k_3 */
  for (i = 0; i < dimension; i++) {
    f_i = (system->eqns[i]).dydt;
    params = (system->eqns[i]).params;
    
    /* Berechne y_temp = y0 + 0.5 * k_2 fuer alle Gleichungen */
    for (j = 0; j < dimension; j++) {
      y_temp[j] = y0[j] + 0.5 * k_2[j];
    }
    
    k_3[i] = h * f_i(t0 + 0.5 * h, y_temp, params);
  }
  
  /* Berechnung von k_4 */
  for (i = 0; i < dimension; i++) {
    f_i = (system->eqns[i]).dydt;
    params = (system->eqns[i]).params;
    
    /* Berechne y_temp = y0 + k_3 fuer alle Gleichungen */
    for (j = 0; j < dimension; j++) {
      y_temp[j] = y0[j] + k_3[j];
    }
    
    k_4[i] = h * f_i(t0 + 0.5 * h, y_temp, params);
  }
  
  /* Berechnet den Loesungsvektor "y1" nach der Zeitentwicklung um "h" durch:
   * y_(n+1) = y_n + h/6 * (k1 + 2*k2 + 2*k3 + k4) fuer jede Gleichung */
  for (i = 0; i < dimension; i++) {
    y1[i] = y0[i] + (k_1[i] + 2 * k_2[i] + 2 * k_3[i] + k_4[i]) / 6.0;
  }
}

ODE_SOLUTION *ode_solution_alloc(int dimension, int t_count) {
  int i;
  ODE_SOLUTION *ret = malloc(sizeof(ODE_SOLUTION));
  if (ret == NULL) {
    return NULL;
  }
  
  ret->t = malloc(t_count * sizeof(double));
  if (ret->t == NULL) {
    free(ret);
    return NULL;
  }
  
  ret->y = malloc(dimension * sizeof(double*));
  if (ret->y == NULL) {
    free(ret->t);
    free(ret);
    return NULL;
  }
  
  for (i = 0; i < dimension; i++) {
    ret->y[i] = malloc(t_count * sizeof(double));
    
    if (ret->y[i] == NULL) {
      for (; i >= 0; i--) {
        free(ret->y[i]);
      }
      free(ret->t);
      free(ret);
      return NULL;
    }
  }
  
  /* Datenfelder setzen: */
  ret->dimension = dimension;
  ret->t_count = t_count;
  
  return ret;
}

void ode_solution_free(ODE_SOLUTION *sol) {
  int i;
  
  free(sol->t);
  for (i = 0; i < sol->dimension; i++) {
    free(sol->y[i]);
  }
  free(sol->y);
  free(sol);
}

ODE_SOLUTION *rk4_solve(ODE_SYSTEM *system, double *y0, double t0, double t1, double h) {
  int i, j;
  int dimension = system->dimension;
  int t_steps = (t1 - t0) / h;
  
  ODE_SOLUTION *sol;
  RK4_WORKSPACE *workspace;
  double *y, *y_next, *swap_temp;
  
  /* Allokiert Speicher fuer die Loesung (t_steps + 1, da die Anfangsbedingung
   * mitgespeichert wird) */
  sol = ode_solution_alloc(dimension, t_steps + 1);
  
  /* Workspace */
  workspace = rk4_workspace_alloc(dimension);
  if (workspace == NULL) {
    ode_solution_free(sol);
    return NULL;
  }
  
  /* Temporaere Loesung */
  y = malloc(dimension * sizeof(double));
  y_next = malloc(dimension * sizeof(double));
  if (y == NULL || y_next == NULL) {
    free(y);
    free(y_next);
    ode_solution_free(sol);
    rk4_workspace_free(workspace);
    return NULL;
  }
  
  /* Anfangsbedingung eintragen */
  sol->t[0] = t0;
  for (i = 0; i < dimension; i++) {
    y[i] = sol->y[i][0] = y0[i];
  }
  
  for (i = 1; i <= t_steps; i++) {
    rk4_evolve(system, workspace, y, t0, h, y_next);
    t0 += h;
    
    /* Loesung eintragen */
    sol->t[i] = t0;
    for (j = 0; j < dimension; j++) {
      sol->y[j][i] = y_next[j];
    }
    
    /* Tausche Pointer um kopieren der einzelnen Werte zu vermeiden */
    swap_temp = y;
    y = y_next;
    y_next = swap_temp;
  }
  
  free(y);
  free(y_next);
  rk4_workspace_free(workspace);
  return sol;
}
