#ifndef _FFT_H
#define _FFT_H

#include <complex.h>

extern const double fft_pi;

typedef enum {
  FFT_SUCCESS,
  FFT_ALLOC_ERROR
} FFT_ERR;

/* Radix-2 FFT: Berechnet die diskrete Fouriertransformation der n = 2^r -
 * Funktionswerte in "f". Dabei werden die Koeffizienten in richtiger Ordnung in
 * f gespeichert.
 * Verwendete Konvention fuer die Fouriertransformation:
 *   g_k = 1/Sqrt[n] * Sum[f_j * Exp[-I * 2*pi/n * k * j], {j, 0, n-1}]
 *   fuer k = 0, 1, ..., n - 1 */
FFT_ERR fft(int r, double complex *f);

/* Bringt die 2^r Koeffizienten in "f" in die vorgesehene Reihenfolge */
FFT_ERR fft_rearrange(int r, double complex *f);

#endif
