#ifndef _GLS_H
#define _GLS_H

#include "numerik_bespin_deutsch_linalg.h"

/* ueber R6, R9 und R12 */
void cube_diag(MATRIX *m, VECTOR *b, double *R);
/* ueber R6 und R9 */
void cube_facediag(MATRIX *m, VECTOR *b, double *R);
/* ueber R12 */
void cube_edge(MATRIX *m, VECTOR *b, double *R);
/* siehe Skizze */
void oktahedron(MATRIX *m, VECTOR *b, double *R);
/* ueber R4 */
void oktahedron_edge(MATRIX *m, VECTOR *b, double *R);

#endif
