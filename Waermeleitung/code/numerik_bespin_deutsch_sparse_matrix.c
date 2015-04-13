#include "numerik_bespin_deutsch_sparse_matrix.h"
#include <stdlib.h>
#include <math.h>

SPARSE_MATRIX *matrix_alloc(int n) {
  int i;
  SPARSE_MATRIX *ret;
  
  /* Allokiert ein SPARSE_MATRIX-Struct */
  ret = malloc(sizeof(SPARSE_MATRIX));
  if (ret == NULL) {
    return NULL;
  }
  
  ret->n = n;
  
  /* Allokiert ein Array fuer die verketteten Listen */
  ret->row = malloc(n * sizeof(NODE*));
  if (ret->row == NULL) {
    free(ret);
    return NULL;
  }
  
  /* Initialisiere die Listen als "leer" */
  for (i = 0; i < n; i++) {
    ret->row[i] = NULL;
  }
  
  return ret;
}

double matrix_get(SPARSE_MATRIX *M, int m, int n) {
  /* Liste der n-ten Zeile */
  NODE *current = M->row[m];
  /* Geht durch die Liste und sucht nach dem Spaltenelement */
  while (current != NULL) {
    if (current->column == n) return current->value;
    current = current->next;
  }
  /* Wurde das Element nicht in der Liste gefunden war es null */
  return 0;
}

int matrix_set(SPARSE_MATRIX *M, int m, int n, double value) {
  /* neues Listenelement allokieren */
  NODE *new = malloc(sizeof(NODE*));
  
  if (NULL == new) return -1;
  
  /* fuegt neues Element an den Anfang der Liste ein */
  new->next = M->row[m];
  new->column = n;
  new->value = value;
  
  M->row[m] = new;
  
  return 0;
}

void matrix_free(SPARSE_MATRIX *M) {
  int i;
  NODE *current, *temp;
  
  /* Zeilen der Matrix */
  for (i = 0; i < M->n; i++) {
    current = M->row[i];
    /* Geht durch die Liste und gibt jeden Knoten frei */
    while (current != NULL) {
      temp = current->next;
      free(current);
      current = temp;
    }
  }
  free(M->row);
  free(M);
}

int gauss_seidel(SPARSE_MATRIX *M, VECTOR *b, VECTOR *sol, double epsilon) {
  /* Maximum der Abweichung einer Variable zwischen zwei aufeinander folgenden
   * Iterationsschritten */
  double delta;
  /* Wert des Loesungselements bei der letzten Iteration (zur Berechnung von
   * "delta") */
  double prev;
  /* Wert des Diagonalelements der jeweiligen Zeile */
  double diag_elem = 0;
  /* Dimension der Matrix */
  int n = M->n;
  
  NODE *current;
  int k;
  
  /* Dimensionskonflikt zwischen Matrix und den Vektoren */
  if (n != b->n || n != sol->n) return -1;
  
  /* Iteration des Gauss-Seidel-Verfahrens bis die maximale Abweichung "delta"
   * kleiner als "epsilon" */
  do {
    delta = 0;
    
    /* Iteration der k-ten Variable des GLS */
    for (k = 0; k < n; k++) {
      current = M->row[k];
      prev = sol->elem[k];
      
      /* Es wird vermieden temporaere Variablen anzulegen und stattdessen direkt
       * im Loesungsvektor gerechnet */
      sol->elem[k] = b->elem[k];
      
      /* Geht durch die nicht verschwindenden Elemente in der k-ten Zeile der 
       * Matrix. Dabei wird das Diagonalelement der Zeile in "diag_elem" ge-
       * speichert und alle anderen werden gemaess dem Gauss-Seidel-Verfahren
       * verarbeitet. */
      while (current != NULL) {
        if (current->column == k) {
          diag_elem = current->value;
        } else {
          /* Hier wird zum Teil schon mit Variablen im naechsten Iterations-
           * schritt gerechnet (sol->elem[l] mit l < k), was zu einer
           * schnelleren Konvergenz fuehrt */
          sol->elem[k] -= current->value * sol->elem[current->column];
        }
        current = current->next;
      }
      
      /* es wird davon ausgegangen, dass "diag_elem" ungleich null */
      sol->elem[k] /= diag_elem;
      
      /* Maximum der Abweichung */
      if (delta < fabs(prev - sol->elem[k])) {
        delta = fabs(prev - sol->elem[k]);
      }
    }
  } while (delta > epsilon);
  
  return 0;
}


VECTOR *vector_alloc(int n) {
  VECTOR *ret;
  
  if ( NULL == (ret = malloc(sizeof(VECTOR))) ) {
    return NULL;
  }
  
  if ( NULL == (ret->elem = malloc(n * sizeof(double))) ) {
    free(ret);
    return NULL;
  }
  ret->n = n;
  
  return ret;
}

void vector_free(VECTOR *v) {
  free(v->elem);
  free(v);
}
