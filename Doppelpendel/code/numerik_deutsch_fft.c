#include "numerik_deutsch_fft.h"
#include <stdlib.h>
#include <math.h>

const double fft_pi = 3.1415926535897932384626433832795;

FFT_ERR fft(int r, double complex *f) {
  int i, j, k;
  int m, K, n;
  int a, b;
  double complex *w;
  
  /* Berechne n = 2^r durch r-maligen Bitshift */
  n = 1 << r;
  
  /* Lookup-Table fuer die mehrfach auftretenden Exponentialfaktoren */
  w = malloc(n * sizeof(double complex));
  if (w == NULL) {
    return FFT_ALLOC_ERROR;
  }
  
  /* Exponentialfaktoren */
  w[0] = 1;
  w[1] = cexp(-2. * fft_pi / n * I);
  for (i = 2; i < n; i++) {
    w[i] = w[i-1] * w[1];
  }
  
  /* Berechnung der Fourierkoeffizienten */
  m = n / 2;
  K = 1;
  
  for (i = 0; i < r; i++) {
    for (k = 0; k < K; k++) {
      for (j = 0; j < m; j++) {
        a = 2 * k * m + j;
        b = a + m;
        
        f[a] = f[a] + f[b];
        f[b] = w[K*j] * (f[a] - 2. * f[b]);
      }
    }
    m = m / 2;
    K = 2 * K;
  }
  
  /* Normierung */
  for (i = 0; i < n; i++) {
    f[i] = f[i] / sqrt(n);
  }
  
  free(w);
  /* Reihenfolge wiederherstellen */
  return fft_rearrange(r, f);
}

FFT_ERR fft_rearrange(int r, double complex *f) {
  int i;
  int l, m, n;
  
  int *index_table;
  double complex *g;
  
  /* 2^r mit Bitshift */
  n = 1 << r;
  
  /* Lookup-Table fuer die korrekte Reihenfolge */
  index_table = malloc(n * sizeof(int));
  /* Zwischenspeicher zum Vertauschen der Reihenfolge */
  g = malloc(n * sizeof(double complex));
  
  if (index_table == NULL || g == NULL) {
    free(index_table);
    free(g);
    return FFT_ALLOC_ERROR;
  }
  
  /* Lookup-Table anlegen */
  l = 0;
  for (i = 0; i < n - 1; i++) {
    index_table[i] = l;
    m = n / 2;
    while (m <= l) {
      l = l - m;
      m = m / 2;
    }
    l = l + m;
  }
  index_table[n-1] = n-1;
  
  /* f in g zwischenspeichern */
  for (i = 0; i < n; i++) {
    g[i] = f[i];
  }
  
  /* von g in f mit richtiger Reihenfolge speichern */
  for (i = 0; i < n; i++) {
    f[i] = g[index_table[i]];
  }
  
  free(index_table);
  free(g);
  
  return FFT_SUCCESS;
}
