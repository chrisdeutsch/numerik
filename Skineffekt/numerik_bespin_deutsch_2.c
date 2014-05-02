#include <stdio.h>
#include <math.h>

double ber(double x);
double bei(double x);

double f1(double x);
double g1(double x);

void test(double epsilon);

int main() {
	double epsilon = 1E-5;
	
  printf("ber(100) = %.10E\n", ber(100));
  printf("bei(100) = %.10E\n", bei(100));
  test(epsilon);
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
  
}

double g1(double x) {
  
}

void test(double epsilon) {
	FILE *file = fopen("berbei_min.tsv", "r");
	double x, ckbei, ckber;
	
	while(fscanf(file, "%lf \t %lE \t %lE", &x, &ckber, &ckbei) != EOF) {
		double delta = fabs((ber(x)-ckber)/ckber);
		if (delta<epsilon) {
			printf("%f %lE %lE\n", delta, ckber, ber(x));
		}
	}
	fclose(file);
}
