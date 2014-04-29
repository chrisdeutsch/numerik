#include <stdio.h>
#include <math.h>

double ber(double x);
double bei(double x);

int main() {
  printf("Hallo Welt\n");

  return 0;
}

double bei(double x) {
  const double xhalf = x / 2;
  
  double sum = xhalf * xhalf;
  double summand = xhalf * xhalf;

  double epsilon = 1E-6;
  int i = 3;

  while (fabs(summand) > epsilon) {
    summand *= -xhalf * xhalf * xhalf * xhalf / ((i - 1) * (i - 1) * i * i);
    i += 2;
    sum += summand;
  }

  return sum;
}