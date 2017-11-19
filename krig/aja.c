#include <math.h>
#include <stdio.h>
#include <string.h>
#define SQUARE(x) ((x)*(x))


/* -------------------- */
/* covariance functions */
/* -------------------- */

double C_covScalingFactor(const char *type) {
	if (strcmp(type, "gauss") == 0) return(sqrt(2.)/2.);
	else if (strcmp(type, "matern3_2") == 0) return(sqrt(3.)); 
	else if (strcmp(type, "matern5_2") == 0) return(sqrt(5.));
	else return(1.);
}


double C_covWhiteNoise(const double *x1, const int *n1, const double *x2, const int *n2, const int *d, const int *i1, const int *i2, const double *param, const double *scaling_factor, const double *var) {
  double s = 0.;
  for (int k = 0; k < *d; k++) {
    s += fabs((x1[*i1 + *n1 * k] - x2[*i2 + *n2 * k]));
  }
  if (s < 0.000000000000001) return(*var);
  else return(0.);
}

double C_covGauss(const double *x1, const int *n1, const double *x2, const int *n2, const int *d, const int *i1, const int *i2, const double *param, const double *scaling_factor, const double *var) {
  double s = 0.;
  for (int k = 0; k < *d; k++) {
    s += SQUARE((x1[*i1 + *n1 * k] - x2[*i2 + *n2 * k])/ (param[k] / *scaling_factor));
  }
  return(exp(-s) * *var);
}

