#ifndef _SPARSE_MATRIX_H
#define _SPARSE_MATRIX_H

/* Knoten einer einfach verketteten Liste. */
typedef struct NODE {
  /* Pointer auf den naechsten Knoten */
  struct NODE *next;
  /* Datenfelder */
  int column;
  double value;
} NODE;

/* Duenne (n x n)-Matrix:
 * Jede Zeile der Matrix wird als verkettete Liste der Nicht-Null-Elemente dar-
 * gestellt. Da jede Zeile der Matrix mindestens ein Element hat, wird beim
 * allokieren der Matrix sofort ein n-dimenensionales Array von verketteten
 * Listen (NODE **row) angelegt. */
typedef struct {
  /* Dimension der Matrix */
  int n;
  /* Verkettete Listen der Zeilen der Matrix */
  NODE **row;
} SPARSE_MATRIX;

/* Vector-Struct fuer n-dimensionale Vektoren
 * elem: Speicherblock fuer Vektorelemente
 *       v->elem[i] fuer das i-te Element */
typedef struct {
  int n;
  double *elem;
} VECTOR;


/* Allokiert eine duenne (n x n)-Matrix
 * Rueckgabewert:
 * NULL: Allokierung fehlgeschlagen */
SPARSE_MATRIX *matrix_alloc(int n);

/* Gibt das Element M[m][n] aus */
double matrix_get(SPARSE_MATRIX *M, int m, int n);

/* Setzt das Element M[m][n] auf den Wert "value"
 * Sollte nicht zum mehrfachen Setzen desselben Elements verwendet werden.
 * Rueckgabewert:
 * 0: Erfolg
 * -1: Allokierung fehlgeschlagen */
int matrix_set(SPARSE_MATRIX *M, int m, int n, double value);

/* Gibt den Speicher der gesamten Matrix wieder frei */
void matrix_free(SPARSE_MATRIX *M);

/* Loest das Gleichungssystem M x = b und speichert den Loesungsvektor in "sol".
 * Der Startvektor fuer das iterative Gauss-Seidel-Verfahren wird mit "sol"
 * uebergeben. Es darf kein Element der Hauptdiagonalen verschwinden.
 * Rueckgabewert:
 * 0: Werfolg
 * -1: Dimensionskonflikt */
int gauss_seidel(SPARSE_MATRIX *M, VECTOR *b, VECTOR *sol, double epsilon);

/* Allokiert einen n-dimensionalen Vektor */
VECTOR *vector_alloc(int n);

/* Gibt den Speicher des Vektors wieder frei */
void vector_free(VECTOR *v);

#endif
