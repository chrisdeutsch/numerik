/* Kompilieranweisungen */
/* Christian Bespin, Christopher Deutsch */
/* TODO: Fakultätsfaktoren ausmultiplizieren */

#include <stdio.h>
#include <math.h>

/* Kelvin-Funktionen und ihre Ableitungen */
double ber(double x);
double bei(double x);
double d_ber(double x);
double d_bei(double x);

/* Die folgenden Funktionen werden für die asymptotische Näherung benötigt 
   vgl. englisches Wikipedia und Abramowitz, Stegun (der ker/kei Anteil wurde
   vernachlässigt) */
double f0(double x);
double g0(double x);
double d_f0(double x);
double d_g0(double x);

/* Berechnet zunächst Werte für Kupfer Radius im cgs-System: Einheit cm */
void table(double I_0, double sigma, double mu,
           double omega, double radius, int N);

/* vergleicht die Ergebnisse der Funktion dbl func(dbl) mit den idealwerten in
   der Datei hinter filename (Format: x-Wert f(x)-Wert). Sollte der relative
   Fehler größer als epsilon sein meldet sich die Funktion */
void test_func(double (*func)(double), char *filename, double epsilon);

/* Konvergenzeinstellung */
const double kThreshold = 10;
const double kEpsilon = 1E-6;

/* Vorberechnete Konstanten (pi wegen ansi konformität M_PI nicht immer vorhanden) */
const double kSqrt2 = 1.4142135623730950488016887242097;
const double kPi = 3.1415926535897932384626433832795;

int main() {
  double I_0 = 1; /* cgs-System Einheit: Fr/s */
  double sigma = 5.356E+17; /* cgs-System Einheit: 1/s */
  double omega = 1E+6; /* Einheit: 1/s */
  double mu = 1; /* Einheit: keine */
  table(I_0, sigma, mu, omega, 0.1, 49);
  
  return 0;
}

double ber(double x) {
  if (fabs(x) < kThreshold) {
    /* Reihendarstellung */
    double sum = 1;
    double summand = 1;
    const double kFactor = pow(x, 4) / 16;
    
    int k = 1;
    
    while (fabs(summand) > kEpsilon * fabs(sum)) {
      summand *= -kFactor / ((2*k - 1) * (2*k - 1) * 2*k * 2*k);
      sum += summand;
      k++;
    }
    
    return sum;
  } else {
    /* asymptotische Näherung */
    double alpha = x / kSqrt2 - kPi / 8;
    double factor = exp(x / kSqrt2) / sqrt(2 * kPi * x);
    
    return factor * (f0(x) * cos(alpha) + g0(x) * sin(alpha));
  }
}

double bei(double x) {
  if (fabs(x) < kThreshold) {
    /* Reihendarstellung */
    double sum = x*x / 4;
    double summand = sum;
    double factor = pow(x, 4) / 16;
    
    int k = 1;
    
    while (fabs(summand) > kEpsilon * fabs(sum)) {
      summand *= -factor / (2*k * 2*k * (2*k + 1) * (2*k + 1));
      sum += summand;
      k++;
    }
    
    return sum;
  } else {
    /* asymptotische Näherung */
    double alpha = x / kSqrt2 - kPi / 8;
    double factor = exp(x / kSqrt2) / sqrt(2 * kPi * x);
    
    return factor * (f0(x) * sin(alpha) - g0(x) * cos(alpha));
  }
}

double d_ber(double x) {
  if (fabs(x) < kThreshold) {
    /* Reihendarstellung */
    double sum = -pow(x, 3) / 16;
    double summand = sum;
    double factor = pow(x, 4) / 16;
    
    int k = 2;
    
    while (fabs(summand) > kEpsilon * fabs(sum)) {
      summand *= -factor / ((2*k - 2) * (2*k - 1) * (2*k - 1) * 2*k);
      sum += summand;
      k++;
    }
    
    return sum;
  } else {
    /* asymptotische Näherung */
    double alpha = x / kSqrt2 - kPi / 8;
    double factor = exp(x / kSqrt2) / sqrt(2 * kPi * x);
    /* factor abgeleitet nach x */
    double d_factor = factor * (1 / kSqrt2 - 1 / (2 * x));
    
    /* Sinus/Cosinus-Berechnung */
    double sin_a = sin(alpha);
    double cos_a = cos(alpha);
    
    /* im wesentlichen Produktregel */
    return d_factor * (f0(x) * cos_a + g0(x) * sin_a) +
           factor * (d_f0(x) * cos_a + d_g0(x) * sin_a -
                     f0(x) * sin_a / kSqrt2 + g0(x) * cos_a / kSqrt2);
  }
}

double d_bei(double x) {
  if (fabs(x) < kThreshold) {
    /* Reihendarstellung */
    double sum = x / 2;
    double summand = sum;
    const double kFactor = pow(x, 4) / 16;
    
    int k = 1;
    
    while (fabs(summand) > kEpsilon * fabs(sum)) {
      summand *= -kFactor / ((2*k - 1) * 2*k * 2*k * (2*k + 1));
      sum += summand;
      k++;
    }
    
    return sum;
  } else {
    /* asymptotische Näherung */
    double alpha = x / kSqrt2 - kPi / 8;
    double factor = exp(x / kSqrt2) / sqrt(2 * kPi * x);
    /* factor abgeleitet nach x */
    double d_factor = factor * (1 / kSqrt2 - 1 / (2 * x));
    
    /* Sinus/Cosinus-Berechnung */
    double sin_a = sin(alpha);
    double cos_a = cos(alpha);
    
    /* im wesentlichen Produktregel */
    return d_factor * (f0(x) * sin_a - g0(x) * cos_a) +
           factor * (d_f0(x) * sin_a - d_g0(x) * cos_a +
                     f0(x) * cos_a / kSqrt2 + g0(x) * sin_a / kSqrt2);
  }
}

double f0(double x) {
  double sum = 1;
  double factor = 1;
  
  int k = 1;
  
  /* factor Vergleich, da Sinus null wird */
  while (fabs(factor) > kEpsilon * fabs(sum)) {
    factor *= (2*k - 1) * (2*k - 1) / (8 * k * x);
    sum += cos(k * kPi / 4) * factor;
    k++;
  }
  
  return sum;
}

double g0(double x) {
  double sum = 0;
  double factor = 1;
  
  int k = 1;
  
  /* factor Vergleich, da Sinus null wird */
  while (fabs(factor) > kEpsilon * fabs(sum)) {
    factor *= (2*k - 1) * (2*k - 1) / (8 * k * x);
    sum += sin(k * kPi / 4) * factor;
    k++;
  }
  
  return sum;
}

double d_f0(double x) {
  double sum = -1/(8 * x * x * kSqrt2);
  double factor = -1/(8 * x * x);
  
  int k = 2;
  
  while (fabs(factor) > kEpsilon * fabs(sum)) {
    factor *= (2*k - 1) * (2*k - 1) /(8 * (k - 1) * x);
    sum += cos(k * kPi / 4) * factor;
    k++;
  }
  
  return sum;
}

double d_g0(double x) {
  double sum = -1/(8 * x * x * kSqrt2);
  double factor = -1/(8 * x * x);
  
  int k = 2;
  
  while (fabs(factor) > kEpsilon * fabs(sum)) {
    factor *= (2*k - 1) * (2*k - 1) /(8 * (k - 1) * x);
    sum += sin(k * kPi / 4) * factor;
    k++;
  }
  
  return sum;
}

void table(double I_0, double sigma, double mu,
           double omega, double radius, int N) {
  /* Einheit 1/cm */
  double kappa = 2 * sqrt(kPi * sigma * mu * omega) / 2.99792458E10;
  
  /* Vorberechnete Werte die in in jedem Schleifendurchlauf gleich sind */
  double factor = I_0 * kappa / (2 * kPi * radius); /* cgs-System Einheit: Fr/s/cm^2 */
  double denominator = d_ber(kappa*radius) * d_ber(kappa*radius)
                       + d_bei(kappa*radius) * d_bei(kappa*radius);
  
  int i;
  double step = radius / (N - 1);
  
  printf("#rho\tamplitude\tphase\n");
  for(i = 0; i < N; i++) {
    double rho = i*step;
    
    double real = ( ber(kappa*rho) * d_bei(kappa*radius)
                    - bei(kappa*rho)* d_ber(kappa*radius) ) / denominator;
    double imag = ( ber(kappa*rho) * d_ber(kappa*radius)
                    + bei(kappa*rho) * d_bei(kappa*radius) ) / denominator;
    
    double amplitude = factor * sqrt(real*real + imag*imag);
    double phase = atan2(imag, real);
    
    printf("%.4f\t%f\t%f\n", i*step, amplitude, phase);
  }
}

void test_func(double (*func)(double), char *filename, double epsilon) {
  FILE *file = fopen(filename, "r");
  
  double x, chkf;
  double del_f;
  int count = 0;
  int failed = 0;
  
  if (file == NULL) {
    printf("Konnte die Datei %s nicht öffnen.", filename);
    return;
  }
  
  printf("Fehlercheck(%s):\nx\t\tdel\n", filename);
  while(fscanf(file, "%lf %lE", &x, &chkf) != EOF) {
    count++;
    del_f = fabs((func(x) - chkf)/chkf);
    
    if (del_f > epsilon) {
      failed++;
      printf("%f\t%E\n", x, del_f);
    }
  }
  printf("%i von %i Tests nicht bestanden (epsilon = %.1E)\n\n\n",
         failed, count, epsilon);
  
  fclose(file);
}
