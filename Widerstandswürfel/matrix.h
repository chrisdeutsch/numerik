#ifndef _MATRIX_H_
#define _MATRIX_H_

/* Matrix-Struct fuer (n x n)-Matrizen 
 * data: Speicherblock fuer die Matrixelemente
 * elem: Pointer auf den Zeilenanfangspointer fuer einfachen Zugriff
 *       M->elem[Zeile][Spalte] */
typedef struct {
  double *data;
  double **elem;
  int n;
} MATRIX;

/* Vector-Struct fuer n-dimensionale Vektoren
 * elem: Speicherblock fuer Vektorelemente
 *       Zugriff: v->elem[i] fuer das i-te Element */
typedef struct {
  double *elem;
  int n;
} VECTOR;


/* Allokiert eine (n x n)-Matrix */
MATRIX *matrix_alloc(int n);

/* Initialisiert eine quad. Matrix mit dem Wert "value" */
void matrix_init(MATRIX *A, double value);

/* Gibt den Speicher der Matrix wieder frei */
void matrix_free(MATRIX *A);

/* Gibt die Matrix auf dem Schirm aus (DEBUG) */
void matrix_print(MATRIX *A);

/* Vertauschung der Zeilen i, j einer Matrix */
void matrix_swap_row(MATRIX *A, int i, int j);


/* Allokiert einen n-dimensionalen Vektor */
VECTOR *vector_alloc(int n);

/* Initialisiert einen Vektor mit dem Wert "value" */
void vector_init(VECTOR *v, double value);

/* Gibt den Speicher des Vektors wieder frei */
void vector_free(VECTOR *v);

/* Gibt den Vektor auf dem Schirm aus (DEBUG) */
void vector_print(VECTOR *v);


/* LU/LR-Zerlegung der Matrix A mit Pivotisierung
 * Die Permutierung der Zeilen wird analog auf dem Array "permutation" durchge-
 * fuehrt. Es wird fuer die LU-Zerlegung keine neue Matrix angelegt, sondern das
 * Ergebnis der Zerlegung direkt in A gespeichert. Dabei entspricht U dem oberen
 * Dreieck von A und L dem Unteren mit 1-en auf der Diagonalen.
 *
 * Rueckgabewert:
 * 0: Erfolgreiche Zerlegung
 * -1: Matrix ist (fast) singulaer */
int LU_decomp(MATRIX *A, int *permutation);

/* Liefert den Zeilenindex des besten Pivot-Elements im k-ten Schritt der Gauß-
 * Elimination */
int pivot_row(MATRIX *A, int k);

/* Loest das System L U x = P b, dabei ist "LU" die LU-Zerlegung nach
 * "LU_decomp" und "Pb" der nach der Pivotisierung permutierte Vektor b.
 * Der Loesungsvektor wird in "sol" gespeichert. */
int LU_solve(MATRIX *LU, VECTOR *Pb, VECTOR *sol);

/* Fuehrt Vorwaerts-/Rueckwaerts-Substitution auf die LU-Zerlegung "LU" mit der
 * Inhomogenitaet "b" durch und speichert das Ergebnis in "sol"
 *
 * Rueckgabewert:
 * 0: Erfolg
 * -1: Dimensionskonflikt zwischen Matrix und Vektor */
int LU_forward_sub(MATRIX *LU, VECTOR *b, VECTOR *sol);
int LU_back_sub(MATRIX *LU, VECTOR *b, VECTOR *sol);

/* Loest das lineare Gleichungssystem Ax = b und speichert den Loesungsvektor in
 * "sol"
 *
 * Rueckgabewert:
 * 0: Erfolg
 * -1: Dimensionskonflikt zwischen A und b
 * -2: Matrix ist (fast) singulaer */
int linear_solve(MATRIX *A, VECTOR *b, VECTOR *sol);

#endif
