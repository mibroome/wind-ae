#include "nr_util.h"

void mmid(double y[], double dydx[], int nvar, double xs, double htot,
          int nstep, double yout[],
          void (*derivs)(double, double[], double[])) {
  int n, i;
  double x, swap, h2, h, *ym, *yn;

  ym = vector(1, nvar);
  yn = vector(1, nvar);
  h  = htot/nstep;

  for (i = 1; i <= nvar; i++) {
    ym[i] = y[i];
    yn[i] = y[i]+h*dydx[i];
  }

  x  = xs+h;
  (*derivs)(x, yn, yout);
  h2 = 2.0*h;

  for (n = 2; n <= nstep; n++) {
    for (i = 1; i <= nvar; i++) {
      swap  = ym[i]+h2*yout[i];
      ym[i] = yn[i];
      yn[i] = swap;
    }
    x += h;
    (*derivs)(x, yn, yout);
  }

  for (i = 1; i <= nvar; i++) {
    yout[i] = 0.5*(ym[i]+yn[i]+h*yout[i]);
  }

  free_vector(yn, 1, nvar);
  free_vector(ym, 1, nvar);

  return;
}
