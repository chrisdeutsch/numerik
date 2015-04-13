#include "numerik_bespin_deutsch_poisson.h"
#include <stdio.h>
#include <stdlib.h>

GRID *grid_alloc(int m, int n, double h) {
  int i;
  GRID *ret = malloc(sizeof(GRID));
  
  if (ret == NULL) return NULL;
  
  ret->m = m;
  ret->n = n;
  ret->h = h;
  
  /* 1-dim. Array von GridPoints */
  ret->grid = malloc(m * n * sizeof(GridPoint));
  if (ret->grid == NULL) {
    free(ret);
    return NULL;
  }
  
  /* Pointerstruktur zum 2-dim. Zugriff auf die Elemente */
  ret->elem = malloc(m * sizeof(GridPoint*));
  if (ret->elem == NULL) {
    free(ret->grid);
    free(ret);
    return NULL;
  }
  
  for (i = 0; i < m; i++) {
    ret->elem[i] = ret->grid + i * n;
  }
  
  return ret;
}

void grid_free(GRID *grid) {
  free(grid->grid);
  free(grid->elem);
  free(grid);
}

int geometry(GRID **grid, int h_100) {
  int i, j;
  
  double h = (double)h_100 / 100;
  
  /* Dimension des Grids */
  int m = 400 / h_100 + 1;
  int n = 500 / h_100 + 1;
  
  /* Division von b, c durch a */
  int b_a = 200 / h_100;
  int c_a = 100 / h_100;
  
  GridPoint **elem;
  
  /* Ueberpruefen ob die Geometrie bei gegebenem Gitterabstand darstellbar ist */
  if (100 % h_100 != 0) {
    return -2;
  }
  
  *grid = grid_alloc(m, n, h);
  
  /* Ueberpruefen ob die Allokierung erfolgreich war */
  if (grid == NULL) {
    return -1;
  }
  
  elem = (*grid)->elem;
  
  /* Bilde zuerst ein Rechteck und nehme dann die Furchen aus der Geometrie */
  /* Kante links & rechts */
  for (i = 0; i < m; i++) {
    elem[i][0].u = 0.22;
    elem[i][0].f = 0;
    elem[i][0].type = DIRICHLET;
    
    elem[i][n-1].u = 0.22;
    elem[i][n-1].f = 0;
    elem[i][n-1].type = DIRICHLET;
  }
  
  /* Unterkante und Oberkante */
  for (i = 1; i < n - 1; i++) {
    elem[0][i].f = 0.95;
    elem[0][i].type = NEUMANN_Y;

    elem[m-1][i].u = 0.22;
    elem[m-1][i].f = 0;
    elem[m-1][i].type = DIRICHLET;
  }
  
  /* Mitte */
  for (i = 1; i < m - 1; i++) {
    for (j = 1; j < n - 1; j++) {
      elem[i][j].f = 0;
      elem[i][j].type = REGULAR;
    }
  }
  
  /* Furchen */
  for (i = b_a + 1; i < m; i++) {
    for (j = c_a + 1; j < 2 * c_a; j++) {
      elem[i][j].u = 0.0;
      elem[i][j].f = 0;
      elem[i][j].type = NONE;
      
      elem[i][j+2*c_a].u = 0.0;
      elem[i][j+2*c_a].f = 0;
      elem[i][j+2*c_a].type = NONE;
    }
  }
  /* Raender der Furchen */
  for (i = 0; i < 4; i++) {
    for (j = b_a; j < m; j++) {
      elem[j][(i+1)*c_a].u = 0.22;
      elem[j][(i+1)*c_a].f = 0;
      elem[j][(i+1)*c_a].type = DIRICHLET;
    }
  }
  for (i = 1; i < c_a; i++) {
    elem[b_a][c_a+i].u = 0.22;
    elem[b_a][c_a+i].f = 0;
    elem[b_a][c_a+i].type = DIRICHLET;
    
    elem[b_a][i+3*c_a].u = 0.22;
    elem[b_a][i+3*c_a].f = 0;
    elem[b_a][i+3*c_a].type = DIRICHLET;
  }
  
  return 0;
}

int geometry_sym(GRID **grid, int h_100) {
  int i, j;
  
  double h = (double)h_100 / 100;
  
  /* Dimension des Grids */
  int m = 400 / h_100 + 1;
  int n = 250 / h_100 + 1;
  
  /* Division von b, c durch a */
  int b_a = 200 / h_100;
  int c_a = 100 / h_100;
  
  GridPoint **elem;
  
  /* Ueberpruefen ob die Geometrie bei gegebenem Gitterabstand darstellbar ist */
  if (50 % h_100 != 0) {
    return -2;
  }
  
  *grid = grid_alloc(m, n, h);
  
  /* Ueberpruefen ob die Allokierung erfolgreich war */
  if (grid == NULL) {
    return -1;
  }
  
  elem = (*grid)->elem;
  
  /* Kante links & rechts */
  for (i = 0; i < m; i++) {
    elem[i][0].u = 0.22;
    elem[i][0].f = 0;
    elem[i][0].type = DIRICHLET;
    
    elem[i][n-1].f = 0;
    elem[i][n-1].type = NEUMANN_X;
  }
  /* Die untere Ecke ist Neumannbedingung in X und Y Richtung */
  elem[0][n-1].f = 0;
  elem[0][n-1].type = NEUMANN_XY;
  
  /* Die obere Ecke ist Dirichletbedingung */
  elem[m-1][n-1].u = 0.22;
  elem[m-1][n-1].f = 0;
  elem[m-1][n-1].type = DIRICHLET;
  
  /* Unterkante und Oberkante */
  for (i = 1; i < n - 1; i++) {
    elem[0][i].f = 0.95;
    elem[0][i].type = NEUMANN_Y;

    elem[m-1][i].u = 0.22;
    elem[m-1][i].f = 0;
    elem[m-1][i].type = DIRICHLET;
  }
  
  /* Mitte */
  for (i = 1; i < m - 1; i++) {
    for (j = 1; j < n - 1; j++) {
      elem[i][j].f = 0;
      elem[i][j].type = REGULAR;
    }
  }
  
  /* Furche */
  for (i = b_a + 1; i < m; i++) {
    for (j = c_a + 1; j < 2 * c_a; j++) {
      elem[i][j].u = 0.0;
      elem[i][j].f = 0;
      elem[i][j].type = NONE;
    }
  }
  /* Raender der Furche */
  for (i = 0; i < 2; i++) {
    for (j = b_a; j < m; j++) {
      elem[j][(i+1)*c_a].u = 0.22;
      elem[j][(i+1)*c_a].f = 0;
      elem[j][(i+1)*c_a].type = DIRICHLET;
    }
  }
  for (i = 1; i < c_a; i++) {
    elem[b_a][c_a+i].u = 0.22;
    elem[b_a][c_a+i].f = 0;
    elem[b_a][c_a+i].type = DIRICHLET;
  }
  
  return 0;
}

int setup_gls(GRID *grid, SPARSE_MATRIX **A, VECTOR **b) {
  int i, j, k;
  double h = grid->h;
  int m = grid->m, n = grid->n;
  int eq_count;
  
  SPARSE_MATRIX *A_temp = NULL;
  VECTOR *b_temp = NULL;
  GridPoint *current = NULL;
  
  /* Zaehlt die Anzahl der Gleichungen im Gleichungssystem (jeder innere Punkt
   * und jeder Punkt mit Neumann-Randbedingung liefert eine Gleichung). Gleich-
   * zeitig wird die Position der Unbekannten im Gleichungssystem gespeichert, 
   * um die Loesung des Systems der Geometrie zuordnen zu koennen. */
  eq_count = 0;
  for (i = 0; i < grid->m; i++) {
    for (j = 0; j < grid->n; j++) {
      current = &(grid->elem[i][j]);
      if (current->type == REGULAR ||
          current->type == NEUMANN_Y ||
          current->type == NEUMANN_X ||
          current->type == NEUMANN_XY) {
        current->position = eq_count;
        eq_count++;
      }
    }
  }
  
  grid->eq_count = eq_count;
  
  /* Koeffizientenmatrix und Inhomogenitaet des Gleichungssystems.
   * Temporaere Variablen, damit A bzw. b nicht staendig dereferenziert werden
   * muss. */
  A_temp = matrix_alloc(eq_count);
  b_temp = vector_alloc(eq_count);
  
  *A = A_temp;
  *b = b_temp;
  
  if (A_temp == NULL || b_temp == NULL) return -1;

  /* Initialisiere null */
  for (i = 0; i < eq_count; i++) {
    b_temp->elem[i] = 0;
  }
  
  /* Stelle fuer jeden Punkt der Geometrie die Gleichung auf.
   * 0 <= k < eq_count zaehlt die aktuell betrachtete Zeile im GLS */
  k = 0;
  for (i = 0; i < m; i++) {
    for (j = 0; j < n; j++) {
      /* Innere Punkte in das Gleichungssystem eintragen */
      if (grid->elem[i][j].type == REGULAR) {
        /* Betrachteter Punkt der Geometrie (vlg. diskretisierter 2-dim.
         * Laplace-Operator) */
        if (matrix_set(A_temp, k, k, 4) != 0) return -2;
        b_temp->elem[k] += h * h * grid->elem[i][j].f;
        
        /* Nachbarn des Punktes */
        /* Oben */
        current = &(grid->elem[i+1][j]);
        if (handle_neighbor(current, A_temp, b_temp, k, 1) != 0) return -2;
        
        /* Unten */
        current = &(grid->elem[i-1][j]);
        if (handle_neighbor(current, A_temp, b_temp, k, 1) != 0) return -2;
        
        /* Links */
        current = &(grid->elem[i][j-1]);
        if (handle_neighbor(current, A_temp, b_temp, k, 1) != 0) return -2;
        
        /* Rechts */
        current = &(grid->elem[i][j+1]);
        if (handle_neighbor(current, A_temp, b_temp, k, 1) != 0) return -2;
        
        k++;
      } 
      /* Punkte mit Neumann-Randbedingung in das Gleichungssystem eintragen */
      else if (grid->elem[i][j].type == NEUMANN_Y) {
        /* Betrachteter Punkt der Geometrie (vgl. diskretisierter 2-dim.
         * Laplace-Operator) */
        if (matrix_set(A_temp, k, k, 4) != 0) return -2;
        b_temp->elem[k] += h * h * grid->elem[i][j].f;
        
        /* Nachbarn des Punktes */
        /* Oben
         * dieser Punkt wird 2-fach gezaehlt, da das Gitter gemaess der
         * Neumann-Bedingung um eine Reihe nach unten fortgesetzt werden kann
         * (vgl. PDF) */
        current = &(grid->elem[i+1][j]);
        if (handle_neighbor(current, A_temp, b_temp, k, 2) != 0) return -2;
        
        /* Links */
        current = &(grid->elem[i][j-1]);
        if (handle_neighbor(current, A_temp, b_temp, k, 1) != 0) return -2;
        
        /* Rechts */
        current = &(grid->elem[i][j+1]);
        if (handle_neighbor(current, A_temp, b_temp, k, 1) != 0) return -2;
        k++;
      }
      /* Symmetrieachse */
      else if (grid->elem[i][j].type == NEUMANN_X) {
        /* Betrachteter Punkt der Geometrie (vgl. diskretisierter 2-dim.
         * Laplace-Operator) */
        if (matrix_set(A_temp, k, k, 4) != 0) return -2;
        b_temp->elem[k] += h * h * grid->elem[i][j].f;
        
        /* Nachbarn des Punktes */
        /* Oben */
        current = &(grid->elem[i+1][j]);
        if (handle_neighbor(current, A_temp, b_temp, k, 1) != 0) return -2;
        
        /* Unten */
        current = &(grid->elem[i-1][j]);
        if (handle_neighbor(current, A_temp, b_temp, k, 1) != 0) return -2;
        
        /* Links
         * dieser Punkt wird 2-fach gezaehlt, da die Symmetrie effektiv eine
         * Neumann-Bedingung mit der Ableitung in x-Richtung gleich 0 ist. */
        current = &(grid->elem[i][j-1]);
        if (handle_neighbor(current, A_temp, b_temp, k, 2) != 0) return -2;
        k++;
      }
      /* Punkt mit Neumann-Randbedingung in X und Y Richtung
       * (tritt an Symmetrieachse auf) */
      else if (grid->elem[i][j].type == NEUMANN_XY) {
        /* Betrachteter Punkt der Geometrie (vgl. diskretisierter 2-dim.
         * Laplace-Operator) */
        if (matrix_set(A_temp, k, k, 4) != 0) return -2;
        b_temp->elem[k] += h * h * grid->elem[i][j].f;
        
        /* Nachbarn des Punktes */
        /* Oben 
         * 2-fach, da dies eine Neumann-Bedingung in y-Richtung ist */
        current = &(grid->elem[i+1][j]);
        if (handle_neighbor(current, A_temp, b_temp, k, 2));
        
        /* Links
         * 2-fach, da dies eine Neumann-Bedingung in x-Richtung ist (Sym.) */
        current = &(grid->elem[i][j-1]);
        if (handle_neighbor(current, A_temp, b_temp, k, 2));
        k++;
      }
    }
  }
  
  return 0;
}

int handle_neighbor(GridPoint *neighbor, SPARSE_MATRIX *A, VECTOR *b, int k, int factor) {
  /* Wenn der Nachbar ein innerer Punkt oder eine Neumann-Randbedingung ist,
   * ist der Wert fuer u unbekannt und muss in das GLS aufgenommen werden */
  if (neighbor->type == REGULAR ||
      neighbor->type == NEUMANN_Y ||
      neighbor->type == NEUMANN_X ||
      neighbor->type == NEUMANN_XY) {
    if (matrix_set(A, k, neighbor->position, -factor) != 0) return -1;
  }
  /* Hat der Nachbar eine Dirichlet-Randbedingung, ist der Wert fuer u bekannt
   * und kann auf die rechte Seite der Gleichung gebracht werden. */
  else if (neighbor->type == DIRICHLET) {
    b->elem[k] += factor * neighbor->u;
  }
  return 0;
}

void enter_solution(GRID *grid, VECTOR *sol) {
  int i, j;
  int m = grid->m, n = grid->n;
  
  for (i = 0; i < m; i++) {
    for (j = 0; j < n; j++) {
      if (grid->elem[i][j].type == REGULAR ||
          grid->elem[i][j].type == NEUMANN_Y ||
          grid->elem[i][j].type == NEUMANN_X ||
          grid->elem[i][j].type == NEUMANN_XY) {
        grid->elem[i][j].u = sol->elem[grid->elem[i][j].position];
      }
    }
  }
}

void mathematica_output(GRID *grid, char *filename) {
  int i, j;
  FILE *file = fopen(filename, "w");
  int m = grid->m, n = grid->n;
  
  fprintf(file, "{");
  for (i = m - 1; i >= 0; i--) {
    fprintf(file, "{");
    for (j = 0; j < n; j++) {
      fprintf(file, "%f", grid->elem[i][j].u);
      if (j < n - 1) fprintf(file, ", ");
    }
    if (i != 0) fprintf(file, "},\n");
    else fprintf(file, "}");
  }
  fprintf(file,"}");
  fclose(file);
}

void mathematica_sym_output(GRID *grid, char *filename) {
  int i, j;
  FILE *file = fopen(filename, "w");
  int m = grid->m, n = grid->n;
  
  fprintf(file, "{");
  for (i = m - 1; i >= 0; i--) {
    fprintf(file, "{");
    for (j = 0; j < n; j++) {
      fprintf(file, "%f", grid->elem[i][j].u);
      fprintf(file, ", ");
    }
    for (j-=2; j >= 0; j--) {
      fprintf(file, "%f", grid->elem[i][j].u);
      if (j > 0) fprintf(file, ", ");
    }
    if (i != 0) fprintf(file, "},\n");
    else fprintf(file, "}");
  }
  fprintf(file,"}");
  fclose(file);
}
