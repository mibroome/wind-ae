/*============================================================================*
 *! \file relax.c                                                             *
 *  \brief Computes the solution in the relaxation domain for the two-point   *
 *         boundary valued problem                                            *
 *============================================================================*/

/* Standard C */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
/* My Headers */
#include "defs.h"
#include "wind.h"
#include "globals.h"
#include "prototypes.h"

/*----------------------------------------------------------------------------*
 *======================== PRIVATE FUNCTION PROTOTYPES =======================*
 *----------------------------------------------------------------------------*/
static void init_relax(double ***y_p, double ***s_p, double ****c_p,
                       int *itmax_p, double *conv_p, double*slowc_p);
static void free_relax(double **y, double **s, double ***c);
static void setup_indices(int *indexv);
static void setup_scales(double *scalv);
static void set_initial_guess(double *x, double *y[NE+1],
                              EQNVARS *equationvars_p);
static void set_relax_soln(double *x, double **y, EQNVARS *equationvars_p);

/*============================================================================*
 *=========================== PUBLIC FUNCTIONS ===============================*
 *============================================================================*
 * relax - Handles the relaxation routine and updates solution array          *
 *----------------------------------------------------------------------------*/

void relax(EQNVARS *equationvars_p) {
  /* variables required for solvde */
  double **y, **s, ***c;
  int itmax, indexv[NE+1];
  double conv, slowc, scalv[NE+1];

  /* Setup */
  /* Initialize variables used by solvede */
  init_relax(&y, &s, &c, &itmax, &conv, &slowc);
//       printf("INIT RELAX done\n");
  /* Set indices to be compatible with Numerical Recipe's boundary condition
     requirements */
  setup_indices(indexv);
//       printf("SETUP INDICES done\n");
  /* Set scales of the equation variables for use in computing the convergence
     criterion */
  setup_scales(scalv);
//       printf("SETUP SCALES done\n");
  /* the initial guess is the input to the function relax */
  set_initial_guess(x, y, equationvars_p);
//     printf("SET INITIAL GUESS done\n");

  /* Call the relaxation routine */
  printf("Starting relaxation routine...\n");
  solvde(itmax, conv, slowc, scalv, indexv, NE, NB, M, y, c, s);
  printf("Finished relaxation routine!\n");

  /* Set the arrays input to the function relax to the solution */
  set_relax_soln(x, y, equationvars_p);

  /* Free memory */
  free_relax(y, s, c);

  return;
}

/*============================================================================*
 *=========================== PRIVATE FUNCTIONS ==============================*
 *============================================================================*
 * init_relax        - inits variables used for relaxation                    *
 * free_relax        - frees variables used for relaxation                    *
 * setup_indices     - remaps the indices to conform to solvde boundary setup *
 * setup_scales      - sets scales for each variable used for convergence     *
 * set_initial_guess - sets the initial guess for relaxation                  *
 * set_relax_soln    - remaps the relaxation solution to our solution array   *
 *----------------------------------------------------------------------------*/

/*============================================================================*
 *! \fn static void init_relax(double ***y_p, double ***s_p, double ****c_p,  *
 *                             int *itmax_p, double *conv_p, double*slowc_p)  *
 *                              int jsf)                                      *
 *  \brief Initialize the variables used by the relaxation routine solvede    *
 *----------------------------------------------------------------------------*/

static void init_relax(double ***y_p, double ***s_p, double ****c_p,
                       int *itmax_p, double *conv_p, double*slowc_p) {
  int i, j;
  double **y, **s, ***c;

  y = calloc_2d_array_gross(1, NE, 1, M);
  s = calloc_2d_array_gross(1, NE, 1, (2*NE+1));
  c = calloc_3d_array_gross(1, NE, 1, (NE-NB+1), 1, (M+1));
  *itmax_p = g_itmax;
  *conv_p  = CONV;
  *slowc_p = SLOWC;

  *y_p = y;
  *s_p = s;
  *c_p = c;

  /* set to zero to start for cleanliness */
  for (i = 1; i <= NE; i++) {
    for (j = 1; j <= (2*NE+1); j++) {
      s[i][j] = 0.0;
    }
  }

  return;
}

/*============================================================================*
 *! \fn static void free_relax(double **y, double **s, double ***c)           *
 *  \brief Free memory allocated for the variables used by solvede            *
 *----------------------------------------------------------------------------*/

static void free_relax(double **y, double **s, double ***c) {
  free_2d_array_gross(y, 1, 1);
  free_2d_array_gross(s, 1, 1);
  free_3d_array_gross(c, 1, 1, 1);

  return;
}

/*============================================================================*
 *! \fn static void setup_indices(int *indexv)                                *
 *  \brief Set the indices of the equation variables such that the interior   *
 *         boundary conditions variables come first, and the outer boundary   *
 *         condition variables go last (see Numerical Recipes for discussion) *
 *----------------------------------------------------------------------------*/

static void setup_indices(int *indexv) {
  /* NR solvde convention: y[j] stores the variable for which indexv[j]
   * is the s-matrix column.  y-slot layout:
   *   y[1]=v, y[2]=z, y[3]=rho, y[4]=T, y[5]=F,
   *   y[6+j]=Ys[j], y[6+N+j]=Ncol[j]
   *
   * s-column assignments:
   *   Cols 1..NB   = lower BC vars: rho(1), Ys[j](2..N+1), T(N+2), F(N+3=NB)
   *   Cols NB+1..NE= upper BC vars: Ncol[j](N+4..2N+3), v(2N+4), z(2N+5)
   */
  int i, j;

  for (i = 1; i <= NE; i++) indexv[i] = i;   /* default identity */

  indexv[1] = 4+2*NSPECIES;   /* v   -> col 4+2N (upper BC) */
  indexv[2] = 5+2*NSPECIES;   /* z   -> col 5+2N (upper BC) */
  indexv[3] = 1;               /* rho -> col 1    (lower BC) */
  indexv[4] = 2+NSPECIES;      /* T   -> col 2+N  (lower BC) */
  indexv[5] = 3+NSPECIES;      /* F   -> col 3+N  (lower BC, = NB) */
  for (j = 0; j < NSPECIES; j++) {
    indexv[6+j]          = 2+j;           /* Ys[j]   -> col 2+j  (lower BC) */
    indexv[6+NSPECIES+j] = 4+NSPECIES+j;  /* Ncol[j] -> col 4+N+j (upper BC) */
  }

  return;
}

/*============================================================================*
 *! \fn static void setup_scales(double *scalv)                               *
 *  \brief Scales for variables to determine the convergence criterion.       *
 *         Scales are defined in relax.h.                                     *
 *----------------------------------------------------------------------------*/
static void setup_scales(double *scalv) {
  int j;
  scalv[1] = VSCALE;
  scalv[2] = ZSCALE;
  scalv[3] = (double)g_rhoscale;
  scalv[4] = TEMPSCALE;
  scalv[5] = FSCALE;              /* F = kappa*dT/dr */
  for (j = 0; j < NSPECIES; j++) {
      scalv[6+j]          = YSSCALE;
      scalv[6+NSPECIES+j] = NCOLSCALE;
  }

  return;
}

/*============================================================================*
 *! \fn static void set_initial_guess(double *x, double *y[NE+1],             *
 *                                    EQNVARS *equationvars_p)                *
 *  \brief Set the solvde variables to have the inital guess provided         *
 *----------------------------------------------------------------------------*/

static void set_initial_guess(double *x, double *y[NE+1],
                              EQNVARS *equationvars_p) {
  int j,i;

  for (j = 1; j <= M; j++) {
    x[j]    = equationvars_p->q[INPTS+j-1];
    y[1][j] = equationvars_p->v[INPTS+j-1];
    y[2][j] = equationvars_p->z[INPTS+j-1];
    y[3][j] = equationvars_p->rho[INPTS+j-1];
    y[4][j] = equationvars_p->T[INPTS+j-1];
    y[5][j] = 0.0;  /* F = 0: conduction will be ramped in by the wrapper */
    for (i = 0; i < NSPECIES; i++) {
      y[6+i][j]          = equationvars_p->Ys[INPTS+j-1][i];
      y[6+NSPECIES+i][j] = equationvars_p->Ncol[INPTS+j-1][i];
    }

  }
    
  return;
}

/*============================================================================*
 *! \fn static void set_relax_soln(double *x, double **y,                     *
 *                                 EQNVARS *equationvars_p)                   *
 *  \brief Take the relaxed solution and map it to our solution array.        *
 *----------------------------------------------------------------------------*/

static void set_relax_soln(double *x, double **y, EQNVARS *equationvars_p) {
  int k, i;

  for (k = 0; k < M; k++) {
    equationvars_p->z[INPTS+k]   = y[2][k+1];
    equationvars_p->v[INPTS+k]   = y[1][k+1];
    equationvars_p->rho[INPTS+k] = y[3][k+1];
    equationvars_p->T[INPTS+k]   = y[4][k+1];
    equationvars_p->Fp[INPTS+k]  = y[5][k+1];
    for (i = 0; i < NSPECIES; i++) {
      equationvars_p->Ys[INPTS+k][i]   = y[6+i][k+1];
      equationvars_p->Ncol[INPTS+k][i] = y[6+NSPECIES+i][k+1];
    }
    equationvars_p->r[INPTS+k] = x[k+1]*y[2][k+1] + parameters.Rmin;
    equationvars_p->q[INPTS+k] = x[k+1];

    /* dT/dr: F/kappa_norm when conduction is on; FD of T when off */
    if (parameters.conduction) {
      I_EQNVARS gv;
      gv.q   = x[k+1];
      gv.v   = y[1][k+1];
      gv.z   = y[2][k+1];
      gv.rho = y[3][k+1];
      gv.T   = y[4][k+1];
      for (i = 0; i < NSPECIES; i++) {
        gv.Ys[i]   = y[6+i][k+1];
        gv.Ncol[i] = y[6+NSPECIES+i][k+1];
      }
      double kn = get_kappa_norm(&gv);
      equationvars_p->dTdr[INPTS+k] = y[5][k+1] / kn;
    } else {
      if (k > 0) {
        double z_avg   = 0.5*(y[2][k+1] + y[2][k]);
        double delta_r = z_avg * (x[k+1] - x[k]);
        equationvars_p->dTdr[INPTS+k] = (delta_r > 0.0)
                                        ? (y[4][k+1] - y[4][k]) / delta_r
                                        : 0.0;
      } else {
        equationvars_p->dTdr[INPTS+k] = 0.0;
      }
    }
  }

  return;
}
