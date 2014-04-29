#include <stdio.h>
#include <math.h>

double ber(double x);
double bei(double x);

int main() {
  printf("Hallo Welt\n");
  printf("%f \n", ber(2));

  return 0;
}

double ber(double x) {
  double epsilon = 1E-6;
  double summe = 1;
  double a = -pow(x, 4) / 4;

  for (int k = 2; fabs(a) > epsilon; k++) {
    summe += a;
    a *= -pow(x, 4) / 16 * 1 / (pow((2 * k + 2), 2) * pow((2 * k + 1), 2));
  }
  
  return summe;
}
