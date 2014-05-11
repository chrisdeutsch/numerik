/* gcc -O2 numerik_bespin_deutsch_2.c -o numerik_2 -lm */
/* Christian Bespin, Christopher Deutsch */

#include <stdio.h>
#include <math.h>

/* Kelvin-Funktionen und ihre Ableitungen */
double ber(double x);
double bei(double x);
double d_ber(double x);
double d_bei(double x);

/* Die folgenden Funktionen werden für die asymptotische Näherung benötigt 
   vgl. Abramowitz, Stegun (der ker/kei Anteil wurde vernachlaessigt) */
double f0(double x);
double g0(double x);
double d_f0(double x);
double d_g0(double x);

/* Berechnet eine Wertetabelle der Stromverteilung mit N gleichverteilten Werten
 * auf das Intervall 0 bis rho_0; Einheiten der Parameter im cgs-System */
void table(double I_0, double sigma, double mu,
           double omega, double rho_0, int N);

/* vergleicht die Ergebnisse der Funktion dbl func(dbl) mit den Idealwerten in
 * der Datei hinter filename (Format: x-Wert f(x)-Wert). Sollte der relative
 * Fehler groesser als epsilon sein meldet sich die Funktion 
 * Benutzung: test_func(ber, "mathematica_vglswerte/ber.tsv", 1E-6); */
void test_func(double (*func)(double), char *filename, double epsilon);

/* Konvergenzeinstellung der Reihenberechnung für die Kelvin-Funktionen 
 * kThreshold: für Argumente groesser als 10 wird die asymptotische Näherung
 *             verwendet 
 * kEpsilon: ist der letzte Summand der Reihe kleiner als kEpsilon * sum, wobei
 *           sum die aktuelle Partialsumme ist, wird die Reihenberechnung ab-
 *           gebrochen */
const double kThreshold = 10;
const double kEpsilon = 1E-6;

/* Vorberechnete Konstanten */
const double kSqrt2 = 1.4142135623730950488016887242097;
const double kPi = 3.1415926535897932384626433832795;

int main() {
  /* Standardwerte fuer Kupfer */
  double I_0 = 1;           /* Strom [Fr/s] */
  double sigma = 5.356E+17; /* Leitfaehigkeit [1/s] */
  double omega = 1E+6;      /* Kreisfrequenz des Wechselstroms [1/s] */
  double mu = 0.999994;     /* Permeabilitaet des Leiters [ ] */
  double rho_0 = 0.1;       /* Radius des Leiters [cm] */
  int N;                    /* Anzahl der berechneten Werte */
  
  /* Ein/Ausgabe: die Berechnung findet in der table-Funktion statt */
  int choice;
  
  test_func(ber, "mathematica_vglswerte/ber.tsv", 1E-6);
  test_func(bei, "mathematica_vglswerte/bei.tsv", 1E-6);
  test_func(d_ber, "mathematica_vglswerte/dber.tsv", 1E-6);
  test_func(d_bei, "mathematica_vglswerte/dbei.tsv", 1E-6);
  
  printf("# Berechnung der Stromdichteverteilung in einem zylindischen Leiter\n");
  printf("# [1]: Standard-Werte fuer Kupferdraht vom Radius 1mm und Strom I_0 = 1 Fr/s\n");
  printf("# [2]: benutzerdefinierte Parameter\n");
  printf("# Auswahl: ");
  if (scanf("%i", &choice) != 1) return 1;
  
  if (choice == 1) {
    printf("\n# Benutze Standard-Werte fuer Kupfer\n");
    
    printf("# Kreisfrequenz des Wechselstroms in Hz: ");
    if (scanf("%lf", &omega) != 1) return 1;
    omega = fabs(omega);
    
    printf("# Anzahl der berechneten Werte N: ");
    if (scanf("%i", &N) != 1) return 1;
    N = N < 2 ? 2 : N; /* N muss mind. 2 sein */
    
    table(I_0, sigma, mu, omega, rho_0, N);
  } else if (choice == 2) {
    printf("\n# benutzerdefinierter Parametersatz\n");
    
    printf("# Strom I_0 in Fr/s: ");
    if (scanf("%lf", &I_0) != 1) return 1;
    I_0 = fabs(I_0);
    
    printf("# Leitfaehigkeit sigma in 1/s: ");
    if (scanf("%lf", &sigma) != 1) return 1;
    sigma = fabs(sigma);
    
    printf("# Permeabilitaet mu: ");
    if (scanf("%lf", &mu) != 1) return 1;
    
    printf("# Leiterradius rho_0 in cm: ");
    if (scanf("%lf", &rho_0) != 1) return 1;
    rho_0 = fabs(rho_0);
    
    printf("# Kreisfrequenz des Wechselstroms in Hz: ");
    if (scanf("%lf", &omega) != 1) return 1;
    omega = fabs(omega);
    
    printf("# Anzahl der berechneten Werte N: ");
    if (scanf("%i", &N) != 1) return 1;
    N = N < 2 ? 2 : N; /* N muss mind. 2 sein */
    
    table(I_0, sigma, mu, omega, rho_0, N);
  } else {
    printf("Eingabe ungültig\n");
    return 1;
  }
  
  return 0;
}

double ber(double x) {
  if (fabs(x) < kThreshold) {
    /* Reihendarstellung */
    double sum = 1;
    double summand = 1;
    double factor = pow(x, 4) / 16;
    
    int k = 1;
    
    while (fabs(summand) > kEpsilon * fabs(sum)) {
      /* berechnet den naechsten Summanden aus dem Vorigen */
      summand *= -factor / ((2*k - 1) * (2*k - 1) * 2*k * 2*k);
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
    /* sum enthält bereits den ersten Summanden daher Start mit k = 1 */
    double sum = x*x / 4;
    double summand = sum;
    double factor = pow(x, 4) / 16;
    
    int k = 1;
    
    while (fabs(summand) > kEpsilon * fabs(sum)) {
      /* berechnet den naechsten Summanden aus dem Vorigen */
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
    /* sum enthält bereits den ersten Summanden daher Start mit k = 2 */
    double sum = -pow(x, 3) / 16;
    double summand = sum;
    double factor = pow(x, 4) / 16;
    
    int k = 2;
    
    while (fabs(summand) > kEpsilon * fabs(sum)) {
      /* berechnet den naechsten Summanden aus dem Vorigen */
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
    /* sum enthält bereits den ersten Summanden daher Start mit k = 1 */
    double sum = x / 2;
    double summand = sum;
    double factor = pow(x, 4) / 16;
    
    int k = 1;
    
    while (fabs(summand) > kEpsilon * fabs(sum)) {
      /* berechnet den naechsten Summanden aus dem Vorigen */
      summand *= -factor / ((2*k - 1) * 2*k * 2*k * (2*k + 1));
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
  
  /* vgl. Abschnitt 2.3.3 in der PDF */
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
  
  /* vgl. Abschnitt 2.3.3 in der PDF */
  while (fabs(factor) > kEpsilon * fabs(sum)) {
    factor *= (2*k - 1) * (2*k - 1) / (8 * k * x);
    sum += sin(k * kPi / 4) * factor;
    k++;
  }
  
  return sum;
}

double d_f0(double x) {
  /* sum enthält bereits den ersten Summanden daher Start mit k = 2 */
  double sum = -1/(8 * x * x * kSqrt2);
  double factor = -1/(8 * x * x);
  
  int k = 2;
  
  /* vgl. Abschnitt 2.3.3 in der PDF */
  while (fabs(factor) > kEpsilon * fabs(sum)) {
    factor *= (2*k - 1) * (2*k - 1) /(8 * (k - 1) * x);
    sum += cos(k * kPi / 4) * factor;
    k++;
  }
  
  return sum;
}

double d_g0(double x) {
  /* sum enthält bereits den ersten Summanden daher Start mit k = 2 */
  double sum = -1/(8 * x * x * kSqrt2);
  double factor = -1/(8 * x * x);
  
  int k = 2;
  
  /* vgl. Abschnitt 2.3.3 in der PDF */
  while (fabs(factor) > kEpsilon * fabs(sum)) {
    factor *= (2*k - 1) * (2*k - 1) /(8 * (k - 1) * x);
    sum += sin(k * kPi / 4) * factor;
    k++;
  }
  
  return sum;
}

void table(double I_0, double sigma, double mu,
           double omega, double rho_0, int N) { 
  /* Vorberechnete Werte die in in jedem Schleifendurchlauf gleich sind */
  double kappa = 2 * sqrt(kPi * sigma * mu * omega) / 2.99792458E10;
  double factor = I_0 * kappa / (2 * kPi * rho_0);

  /* Die Ableitungen der Kelvin-Funktion haben immer dasselbe Argument: */
  double d_ber_rho0 = d_ber(kappa * rho_0);
  double d_bei_rho0 = d_bei(kappa * rho_0);
  
  double denominator = d_ber_rho0 * d_ber_rho0 + d_bei_rho0 * d_bei_rho0;
  
  /* Laufvariable und Schrittgroesse der rho-Werte in der Wertetabelle */
  int i;
  double step = rho_0 / (N - 1);
  
  if (fabs(kappa * rho_0) > 1000) {
    printf("\n\n# WARNUNG: Die implementierten Kelvin-Funktionen wurden nur "
           "fuer Argumente |x| < 1000 getestet. Es kann zu Ueberlaeufen kommen\n");
  }
  
  printf("\n# Verwendete Parameter:\n");
  printf("# Strom: I_0 = %E Fr/s\n", I_0);
  printf("# Leitfaehigkeit: sigma = %E 1/s\n", sigma);
  printf("# Permeabilitaet: mu = %E\n", mu);
  printf("# Wechselstromkreisfrequenz: omega = %E 1/s\n", omega);
  printf("# Leiterradius: rho_0 = %f cm\n\n", rho_0);
 
  printf("#rho[cm]\t\t|j|[Fr/s/cm^2]\t\tphi[rad]\n");
  
  /* Berechnung nach den Formeln 18 - 20 in der PDF */
  for(i = 0; i < N; i++) {
    double rho = i * step;
    
    /* Mehrfachberechnung von ber/bei vermeiden: */
    double ber_rho = ber(kappa * rho);
    double bei_rho = bei(kappa * rho);
    
    double real = (ber_rho * d_bei_rho0 - bei_rho * d_ber_rho0) / denominator;
    
    double imag = (ber_rho * d_ber_rho0 + bei_rho * d_bei_rho0) / denominator;
    
    double amplitude = factor * sqrt(real*real + imag*imag);
    double phase = atan2(imag, real);
    
    printf("%f\t\t%f\t\t%f\n", rho, amplitude, phase);
  }
}

void test_func(double (*func)(double), char *filename, double epsilon) {
  FILE *file = fopen(filename, "r");
  
  double x, chkf;
  double del_f;
  int count = 0;
  int failed = 0;
  
  if (file == NULL) {
    printf("Konnte die Datei %s nicht öffnen.\n", filename);
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
