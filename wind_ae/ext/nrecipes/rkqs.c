#include <math.h>
#include <stdio.h>
#include "nr_prototypes.h"
#include "nr_util.h"
#include "../../src/defs.h"

#define SAFETY 0.9
#define PGROW -0.2
#define PSHRNK -0.25
#define ERRCON 1.89e-4
double **d, *x;

void rkqs(double y[], double dydx[], int n, double *x, double htry, double eps,
    double yscal[], double *hdid, double *hnext,
    void(*derivs)(double, double [], double []))
{
    void rkck(double y[], double dydx[], int n, double x, double h,
        double yout[], double yerr[], void(*derivs)(double, double [], double []));
//     void rk4(double y[], double dydx[], int n, double x, double h,
//         double yout[], void(*derivs)(double, double [], double []));
    int i;
    double errmax, h, xnew, *yerr, *ytemp;

    yerr = vector(1, n);
    ytemp = vector(1, n);
    h = htry;
    for (;;){
        rkck(y, dydx, n, *x, h, ytemp, yerr, derivs);
//         rk4(y, dydx, n, *x, h, ytemp, derivs);
        errmax = 0.0;
        for (i = 1; i <= n; i++) errmax = FMAX(errmax, fabs(yerr[i] / yscal[i]));
        errmax /= eps;
        if (errmax > 1.0) {
            h = SAFETY*h*pow(errmax, PSHRNK);
            if (h < 0.1*h) h *= 0.1;
            xnew = (*x) + h;
            if (xnew == *x) nrerror("stepsize underflow in rkqs");
            continue;
        }
        else {
            if (errmax > ERRCON) *hnext = SAFETY*h*pow(errmax, PGROW);
            else *hnext = 5.0*h;
            *x += (*hdid = h);
            for (i = 1; i <= n; i++) y[i] = ytemp[i];
            break;
        }
    }
        free_vector(ytemp, 1, n);
        free_vector(yerr, 1, n);
}
