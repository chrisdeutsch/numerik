#ifndef _GLS_H
#define _GLS_H

#include "numerik_bespin_deutsch_linalg.h"

/* Gleichungssysteme aus der Kirchhoff'schen Maschenregel
 *
 * Die Funktionen erhalten ein Array von Widerstandswerten (12 Stueck) und
 * weisen der Matrix A die Koeffizientenmatrix fuer die uebergebenen
 * Widerstaende zu. Außerdem wird die Inhomogenitaet fuer das jeweilige
 * Gleichungssystem gesetzt. 
 *
 * Fuer die Gleichungssysteme und die Geometrien sei auf den Anhang der PDF
 * verwiesen. */

/* ueber die Raumdiagonale des Wuerfels */
void cube_diag(MATRIX *A, VECTOR *b, double *R);

/* ueber die Flaechendiagonale des Wuerfels */
void cube_facediag(MATRIX *A, VECTOR *b, double *R);

/* ueber eine Kante des Wuerfels */
void cube_edge(MATRIX *A, VECTOR *b, double *R);

/* ueber zwei gegenueberliegende Spitzen (Raumdiagonale) des Oktaeder */
void octahedron_diag(MATRIX *A, VECTOR *b, double *R);

/* ueber eine Kante des Oktaeders */
void octahedron_edge(MATRIX *A, VECTOR *b, double *R);

#endif
