/* Kompilieranweisungen */
/* Christian Bespin, Christopher Deutsch */
/* TODO: mehrfache Sinus/Cosinusberechnung in d_ber/d_bei entfernen */

#include <stdio.h>
#include <math.h>

/* Kelvin-Funktionen und ihre Ableitungen */
double ber(double x);
double bei(double x);
double d_ber(double x);
double d_bei(double x);

/* vergleicht die Ergebnisse der Funktion dbl func(dbl) mit den idealwerten in
   der Datei hinter filename (Format: x-Wert f(x)-Wert). Sollte der relative
   Fehler größer als epsilon sein meldet sich die Funktion */
void test_func(double (*func)(double), char *filename, double epsilon);

/* Die folgenden Funktionen werden für die asymptotische Näherung benötigt 
   vgl. englisches Wikipedia und Abramowitz, Stegun (der ker/kei Anteil wurde
   vernachlässigt) */
double f0(double x);
double g0(double x);
double d_f0(double x);
double d_g0(double x);

/* Konvergenzeinstellung */
const double kThreshold = 10;
const double kEpsilon = 1E-6;

/* Vorberechnete Konstanten (pi wegen ansi konformität M_PI nicht immer vorhanden) */
const double kSqrt2 = 1.4142135623730950488016887242097;
const double kPi = 3.1415926535897932384626433832795;

int main() {
  test_func(ber, "literatur/ber.tsv", 1E-7);
  test_func(bei, "literatur/bei.tsv", 1E-7);
  test_func(d_ber, "literatur/dber.tsv", 1E-7);
  test_func(d_bei, "literatur/dbei.tsv", 1E-7);
  return 0;
}

double ber(double x) {
  if (fabs(x) < kThreshold) {
    /* Reihendarstellung */
    double sum = 1;
    double summand = 1;
    const double kFactor = pow(x, 4) / 16;
    
    int i = 2;
    
    while (fabs(summand) > kEpsilon * fabs(sum)) {
      summand *= -kFactor / ((i - 1) * (i - 1) * i * i);
      sum += summand;
      i += 2;
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
    const double kFactor = pow(x, 4) / 16;
    
    int i = 3;
    
    while (fabs(summand) > kEpsilon * fabs(sum)) {
      summand *= -kFactor / ((i - 1) * (i - 1) * i * i);
      sum += summand;
      i += 2;
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
    const double kFactor = pow(x, 4) / 16;
    
    int i = 4;
    
    while (fabs(summand) > kEpsilon * fabs(sum)) {
      summand *= -kFactor / ((i - 2) * (i - 1) * (i - 1) * i);
      sum += summand;
      i += 2;
    }
    
    return sum;
  } else {
    /* asymptotische Näherung */
    double alpha = x / kSqrt2 - kPi / 8;
    double factor = exp(x / kSqrt2) / sqrt(2 * kPi * x);
    /* factor abgeleitet nach x */
    double d_factor = factor * (1 / kSqrt2 - 1 / (2 * x));
    
    /* im wesentlichen Produktregel */
    return d_factor * (f0(x) * cos(alpha) + g0(x) * sin(alpha)) +
           factor * (d_f0(x) * cos(alpha) + d_g0(x) * sin(alpha) -
                     f0(x) * sin(alpha) / kSqrt2 + g0(x) * cos(alpha) / kSqrt2);
  }
}

double d_bei(double x) {
  if (fabs(x) < kThreshold) {
    /* Reihendarstellung */
    double sum = x / 2;
    double summand = sum;
    const double kFactor = pow(x, 4) / 16;
    
    int i = 2;
    
    while (fabs(summand) > kEpsilon * fabs(sum)) {
      summand *= -kFactor / ((i - 1) * i * i * (i + 1));
      sum += summand;
      i += 2;
    }
    
    return sum;
  } else {
    /* asymptotische Näherung */
    double alpha = x / kSqrt2 - kPi / 8;
    double factor = exp(x / kSqrt2) / sqrt(2 * kPi * x);
    /* factor abgeleitet nach x */
    double d_factor = factor * (1 / kSqrt2 - 1 / (2 * x));
    
    /* im wesentlichen Produktregel */
    return d_factor * (f0(x) * sin(alpha) - g0(x) * cos(alpha)) +
           factor * (d_f0(x) * sin(alpha) - d_g0(x) * cos(alpha) +
                     f0(x) * cos(alpha) / kSqrt2 + g0(x) * sin(alpha) / kSqrt2);
  }
}

double f0(double x) {
  double sum = 1;
  double factor = 1;
  
  int i = 1;
  
  /* factor Vergleich, da Sinus null wird */
  while (fabs(factor) > kEpsilon * fabs(sum)) {
    factor *= (2*i - 1) * (2*i - 1) / (8 * i * x);
    sum += cos(i * kPi / 4) * factor;
    i++;
  }
  
  return sum;
}

double g0(double x) {
  double sum = 0;
  double factor = 1;
  
  int i = 1;
  
  /* factor Vergleich, da Sinus null wird */
  while (fabs(factor) > kEpsilon * fabs(sum) || (i - 1)%4 == 0) {
    factor *= (2*i - 1) * (2*i - 1) / (8 * i * x);
    sum += sin(i * kPi / 4) * factor;
    i++;
  }
  
  return sum;
}

double d_f0(double x) {
  double sum = -1/(8 * x * x * kSqrt2);
  double factor = -1/(8 * x * x);
  
  int i = 2;
  
  while (fabs(factor) > kEpsilon * fabs(sum)) {
    factor *= (2*i - 1) * (2*i - 1) /(8 * (i - 1) * x);
    sum += cos(i * kPi / 4) * factor;
    i++;
  }
  
  return sum;
}

double d_g0(double x) {
  double sum = -1/(8 * x * x * kSqrt2);
  double factor = -1/(8 * x * x);
  
  int i = 2;
  
  while (fabs(factor) > kEpsilon * fabs(sum)) {
    factor *= (2*i - 1) * (2*i - 1) /(8 * (i - 1) * x);
    sum += sin(i * kPi / 4) * factor;
    i++;
  }
  
  return sum;
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
