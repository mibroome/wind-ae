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
  *itmax_p = ITMAX;
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
    //set such that it is easy to list and index by the indices (e.g., Ys1 and Ys2 are now adjacent in indices [5+j] and the values
    //of these indices = 2+j as is required by Numerical Recipes so that the interior boundary conditions come first.
  int i,j,k,m;

  for (i = 1; i <= NE; i++) {
    indexv[i] = i;
  }
  indexv[1] = 3+2*NSPECIES;
  indexv[2] = 4+2*NSPECIES;
  indexv[3] = 1;
  indexv[4] = 2+NSPECIES;
  for (j=0; j<NSPECIES; j++){
      k = 5+j;
      m = 5+j+NSPECIES;
      indexv[k] = 2+j; 
      indexv[m] = 3+NSPECIES+j;
  }

  return;
}

/*============================================================================*
 *! \fn static void setup_scales(double *scalv)                               *
 *  \brief Scales for variables to determine the convergence criterion.       *
 *         Scales are defined in relax.h.                                     *
 *----------------------------------------------------------------------------*/
static void setup_scales(double *scalv) {
  int j,k,m;
  /* Note: not using the indexv's */
  scalv[1] = VSCALE;
  scalv[2] = ZSCALE;
  scalv[3] = RHOSCALE;
  scalv[4] = TEMPSCALE;
  for (j=0; j<NSPECIES; j++) {
      k = 5+j;
      m = 5+NSPECIES+j;
      scalv[k] = FPSCALE;
      scalv[m] = NCOLSCALE;      
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
  int j,i,k,m;

  for (j = 1; j <= M; j++) {
    x[j] = equationvars_p->q[INPTS+j-1];
    y[1][j] = equationvars_p->v[INPTS+j-1];
    y[2][j] = equationvars_p->z[INPTS+j-1];
    y[3][j] = equationvars_p->rho[INPTS+j-1];
    y[4][j] = equationvars_p->T[INPTS+j-1];
    for (i=0; i<NSPECIES; i++){
        k=i+5;
        m=i+5+NSPECIES;
        y[k][j] = equationvars_p->Ys[INPTS+j-1][i];
        y[m][j] = equationvars_p->Ncol[INPTS+j-1][i];
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
  int k,j,l,m;

  for (k = 0; k < M; k++) {
    equationvars_p->z[INPTS+k]    = y[2][k+1];
    equationvars_p->v[INPTS+k]    = y[1][k+1]; //FIX HERE
    equationvars_p->rho[INPTS+k]  = y[3][k+1];
    equationvars_p->T[INPTS+k]    = y[4][k+1];
    for (j=0; j<NSPECIES; j++) {
        l=j+5;
        m=j+5+NSPECIES;
        equationvars_p->Ys[INPTS+k][j]   = y[l][k+1];
        equationvars_p->Ncol[INPTS+k][j] = y[m][k+1];
    }
    equationvars_p->r[INPTS+k]    = x[k+1]*y[2][k+1]+parameters.Rmin;
    equationvars_p->q[INPTS+k]    = x[k+1];
  }

  return;
}
