#include <stdio.h>
#include <math.h>

#define PI 3.14

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


  printf("f(100) = %.5E\n", f1(1000));
  printf("g(100) = %.5E\n", g1(1000));

  return 0;
}

double ber(double x) {
  if (fabs(x) <= 100) {
    const double factor = pow(x, 4) / 16;

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
  } else {
    return exp(x / sqrt(2)) / sqrt(2 * PI * x) *
           (f1(x) * cos(x / sqrt(2) - PI / 8) +
            g1(x) * sin(x / sqrt(2) - PI / 8));
  }
}

double bei(double x) {
  const double factor = pow(x, 4) / 16;

  double sum = x * x / 4;
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
  while (fabs(factor) > epsilon) {
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
  while (fabs(factor) > epsilon) {
    factor *= (2 * k - 1) * (2 * k - 1) / (8 * k * x);

    sum += sin(k * PI / 4) * factor;
    k++;
  }

  return sum;
}

void test(double epsilon) {
	FILE *file = fopen("berbei_min.tsv", "r");
	double x, ckbei, ckber;
	
	printf("x \t \t ber(x) \t ber Lit\t delta ber \t bei(x) \t Lit bei \t delta bei\n");
	while(fscanf(file, "%lf \t %lE \t %lE", &x, &ckber, &ckbei) != EOF) {
		double deltaber = fabs((ber(x)-ckber)/ckber);
		double deltabei = fabs((bei(x)-ckbei)/ckbei);
		if (deltaber<epsilon || deltabei<epsilon) {
			printf("%f \t %lE \t %lE \t %lE \t %lE \t %f \t %f\n", x, ber(x), ckber, deltaber, bei(x), ckbei, deltabei);
		}
	}
	fclose(file);
}
