#ifndef _POISSON_H
#define _POISSON_H

#include "numerik_bespin_deutsch_sparse_matrix.h"

/* Klassifizierung der Punkte der Diskretisierung:
 * REGULAR: Innerer Punkt mit unbekanntem u
 * DIRICHLET: auesserer Punkte mit bekanntem u
 * NEUMANN: auesserer Punkt mit bekannter Ableitung in X/Y-Richtung
 * NONE: aeusserer Punkt */
typedef enum {
  REGULAR,
  DIRICHLET,
  NEUMANN_Y,
  NEUMANN_X,
  NEUMANN_XY,
  NONE
} PointType;

typedef struct {
  /* LSG. diskrete Poisson Gleichung */
  double u;
  /* RHS Poisson Gleichung */
  double f;
  /* Klassifizierung des Punktes */
  PointType type;
  /* Position im Gleichungssystem */
  int position;
  } GridPoint;

typedef struct {
  int m;
  int n;
  /* Gitterabstand */
  double h;
  
  /* Anzahl der Gleichungen im GLS */
  int eq_count;
  
  /* Array von GridPoints */
  GridPoint *grid;
  /* Pointerstruktur zum einfachen Zugriff auf die Elemente (grid->elem[m][n]) */
  GridPoint **elem;
} GRID;

/* Allokierung eines (m x n)-Grids mit Gitterabstand h. Bei Fehlgeschlagener
 * Allokierung wird "NULL" zurueckgegeben. */
GRID *grid_alloc(int m, int n, double h);

/* Freigeben des Speichers des Grids */
void grid_free(GRID *grid);

/* Erstellt die Geometrie mit dem Gitterabstand "h_100" in Hundersteln */
int geometry(GRID **grid, int h_100);

/* Erstellt die Geometrie mit dem Gitterabstand "h_100" in Hundersteln unter 
 * Beachtung der Symmetrie */
int geometry_sym(GRID **grid, int h_100);

/* Erstellt das Gleichungssystem zur Diskretisierung */
int setup_gls(GRID *grid, SPARSE_MATRIX **A, VECTOR **b);

/* Fuegt den Nachbarpunkt gemaess seines Typs in das Gleichungssystem ein */
int handle_neighbor(GridPoint *neighbor, SPARSE_MATRIX *A, VECTOR *b, int k, int times);

/* Traegt die Loesung in das Grid ein */
void enter_solution(GRID *grid, VECTOR *sol);

void mathematica_output(GRID *grid, char *filename);

void mathematica_sym_output(GRID *grid, char *filename);

#endif
