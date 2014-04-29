#include <stdio.h>
#include <math.h>

double ber(double x);
double bei(double x);

int main() {
  printf("ber(100) = %.10E\n", ber(100));
  printf("bei(100) = %.10E\n", bei(100));
  return 0;
}

double ber(double x) {
   const double factor = pow(x,4)/16;

   double sum = 1;
   double summand = 1;

   double epsilon = 1E-6;
   int i = 2;

   while (fabs(summand) > epsilon) {
     summand *= -factor / ((i - 1) * (i - 1) * i * i);
     i += 2;
     sum += summand;
   }

   return sum;
}

double bei(double x) {
  const double factor = pow(x,4)/16;

  double sum = x*x/4;
  double summand = sum;

  double epsilon = 1E-6;
  
  int i = 3;
  while (fabs(summand) > epsilon) {
    summand *= -factor / ((i - 1) * (i - 1) * i * i);
    i += 2;
    sum += summand;
  }

  return sum;
}
