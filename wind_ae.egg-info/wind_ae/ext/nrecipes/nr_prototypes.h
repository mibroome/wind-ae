#ifndef _NR_H_
#define _NR_H_

/* bksub.c */
void bksub(int ne, int nb, int jf, int k1, int k2, double ***c);

/* bsstep.c */
void bsstep(double y[], double dydx[], int nv, double *xx, double htry,
            double eps, double yscal[], double *hdid, double *hnext,
            void (*derivs)(double, double [], double []));

/* rk4.c (added 9/20/23)*/
void rk4(double y[], double dydx[], int n, double x, double h, double yout[],
         void(*derivs)(double, double[], double[]));


/* rkck.c (added 9/12/23)*/
void rkck(double y[], double dydx[], int n, double x, double h, double yout[], 
          double yerr[], void (*derivs)(double, double [], double []));

/* rkqs.c (added 9/12/23)*/
void rkqs(double y[], double dydx[], int n, double *x, double htry, double eps,
          double yscal[], double *hdid, double *hnext,
          void(*derivs)(double, double [], double []));

/* mmid.c */
void mmid(double y[], double dydx[], int nvar, double xs, double htot,
          int nstep, double yout[], void (*derivs)(double, double[], double[]));

/* rkdumb.c (added 9/20/23)*/
void rkdumb(double vstart[], int nvar, double x1, double x2, long nstep,
            void(*derivs)(double, double[], double[]));

/* odeint.c */
void odeint(double ystart[], int nvar, double x1, double x2,
            double eps, double h1, double hmin, int *nok, int *nbad,
            void (*derivs)(double, double [], double []),
            void (*rkqs)(double [], double [], int, double *, double, double,
                         double [], double *, double *,
                         void (*)(double, double [], double [])));

/* pinvs.c */
void pinvs(int ie1, int ie2, int je1, int jsf, int jc1, int k,
           double ***c, double **s);

/* polint.c */
void polint(double xa[], double ya[], int n, double x, double *y, double *dy);

/* pzextr.c */
void pzextr(int iest, double xest, double yest[], double yz[], double dy[],
            int nv);

/* red.c */
void red(int iz1, int iz2, int jz1, int jz2, int jm1, int jm2, int jmf,
         int ic1, int jc1, int jcf, int kc, double ***c, double **s);

/* solvde.c */
void solvde(int itmax, double conv, double slowc, double scalv[],
            int indexv[], int ne, int nb, int m, double **y, double ***c,
            double **s);

#endif /* _NR_H_ */
