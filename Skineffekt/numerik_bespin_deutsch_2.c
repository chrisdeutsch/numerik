#include <stdio.h>
#include <math.h>

#define PI 3.14

double ber(double x);
double bei(double x);

double f1(double x);
double g1(double x);

int main() {
  printf("ber(100) = %.10E\n", ber(100));
  printf("bei(100) = %.10E\n", bei(100));
  
  printf("f(100) = %.5E\n", f1(1000));
  printf("g(100) = %.5E\n", g1(1000));
  
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

double f1(double x) {
  double sum = 1;
  double factor = 1;
  double epsilon = 1E-4;
  
  int k = 1;
  while(fabs(factor) > epsilon) {
    factor *= (2 * k - 1) * (2 * k - 1) / (8 * k * x);
    
    sum += cos(k * PI / 4) * factor;
    k++;
  }
  
  return sum;
}

double g1(double x) {
  double sum = 0;
  double factor = 1;
  double epsilon = 1E-4;
  
  int k = 1;
  while(fabs(factor) > epsilon) {
    factor *= (2 * k - 1) * (2 * k - 1) / (8 * k * x);
    
    sum += sin(k * PI / 4) * factor;
    k++;
  }
  
  return sum;  
}
