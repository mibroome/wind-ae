#include <math.h>
#include <stdio.h>
#include "nr_util.h"

#define MAXSTP 1E6
#define TINY 1.0e-30

extern int kmax, kount;
extern double *xp, **yp, dxsav;

void odeint(double ystart[], int nvar, double x1, double x2, double eps,
            double h1, double hmin, int *nok, int *nbad,
            void (*derivs)(double, double [], double []),
            void (*rkqs)(double [], double [], int, double *, double, double,
                         double [], double *, double *,
                         void (*)(double, double [], double []))) {
  long nstp;
  int i;
  double xsav, x, hnext, hdid, h;
  double *yscal, *y, *dydx;

  yscal = vector(1, nvar);
  y = vector(1, nvar);
  dydx  = vector(1, nvar);
  x = x1;
  h = SIGN(h1, x2-x1);
  *nok = (*nbad) = kount = 0;
  for (i = 1; i <= nvar; i++) {
    y[i] = ystart[i];
  }
  if (kmax > 0) {
    xsav = x-dxsav*2.0;
  }
  for (nstp = 1; nstp <= MAXSTP; nstp++) {
    (*derivs)(x, y, dydx);
    for (i = 1; i <= nvar; i++) {
      yscal[i] = fabs(y[i])+fabs(dydx[i]*h)+TINY;
    }
    if (kmax > 0 && kount < kmax-1 && fabs(x-xsav) > fabs(dxsav)) {
      xp[++kount] = x;
      for (i = 1; i <= nvar; i++) {
        yp[i][kount] = y[i];
      }
      xsav = x;
    }
    if ((x+h-x2)*(x+h-x1) > 0.0) {
      h = x2-x;
    }
    (*rkqs)(y, dydx, nvar, &x, h, eps, yscal, &hdid, &hnext, derivs);
    if (hdid == h) {
      ++(*nok);
    }
    else {
      ++(*nbad);
    }
    if ((x-x2)*(x2-x1) >= 0.0) {
      for (i = 1; i <= nvar; i++) {
        ystart[i] = y[i];
      }
      if (kmax) {
        xp[++kount] = x;
        for (i = 1; i <= nvar; i++) {
          yp[i][kount] = y[i];
        }
      }
      free_vector(dydx, 1, nvar);
      free_vector(y, 1, nvar);
      free_vector(yscal, 1, nvar);
      return;
    }
    if (fabs(hnext) <= hmin) {
      nrerror("Step size too small in odeint");
    }
    h = hnext;
  }
  nrerror("Too many steps in routine odeint");

  return;
}
#undef MAXSTP
#undef TINY
