/*============================================================================*
 *! \file intode.c (integrate ordinary differential equations)                *
 *  \brief Integrates the solution outside the relaxation domain, which spans *
 *         from Rmin to the critical point. Goes outwards to Rmax, and goes   *
 *         inwards to RIN.                                                    *
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

/* global variables for odeint */
int kmax, kount;
double *xp, **yp, dxsav;
/* global variables for rkdumb */
double *xx, **y;
/* static global variables */
static double thez;
static I_EQNVARS critvars;

/*----------------------------------------------------------------------------*
 *======================== PRIVATE FUNCTION PROTOTYPES =======================*
 *----------------------------------------------------------------------------*/
static void ode_system(double x, double y[], double dydx[]);
static void init_odeint(int indx, double **ystart, EQNVARS equationvars);
static void free_odeint(double *ystart);
static void set_odeint_soln(int i, double end, double *ystart,
                            EQNVARS *equationvars_p);

/*============================================================================*
 *=========================== PUBLIC FUNCTIONS ===============================*
 *============================================================================*
 * integrate_ode - integrate the system of odes outside the relaxation domain *
 *----------------------------------------------------------------------------*/

/*============================================================================*
 *! \fn int integrate_ode(EQNVARS *equationvars_p)                            *
 *  \brief With the relaxation solution in hand, we now have the critical     *
 *         point and need not worry about futher boundary conditions. Thus,   *
 *         we expand the the solution outside the relaxation domain, which is *
 *         no longer a two-point boundary value problem, and use standard ode *
 *         solves for our inital value problem (specifically Bulirsch-Stoer)  *
 *----------------------------------------------------------------------------*/

#define DEBUG_IN 0
#define DEBUG_OUT 0
#define NSTEPS 50000
void integrate_ode(EQNVARS *equationvars_p) {
  int i,j,k, nok, nbad;//, ne;
  double *ystart;
  double rsp, start, end, delta;

  /* Integrate outwards */
  if (parameters.integrate_outward) {
    printf("Integrating outwards...\n");
    /* Initialize variables used by odeint */
    init_odeint(INPTS+M-1, &ystart, *equationvars_p);
    rsp = equationvars_p->r[INPTS+M-1];
    if (parameters.Rmax-rsp < 0.1*rsp) {
      delta = 0.1*rsp/ADDPTS;
    }
    else {
      delta = (parameters.Rmax-rsp)/ADDPTS;
    }
    if ((parameters.Rmax-rsp)/ADDPTS > 0.1) {
      printf("Outward stepping probably too large.\n"
             "  Suggest running smaller Rmax or larger ADDPTS.\n");
    }
    end = rsp;
    for (i = 0; i < ADDPTS; i++) {
      start = end;
      end   = start+delta;
        /* Not adapted to array-based multispecies approach, so commented out*/
// #if DEBUG_OUT
//       printf("outward step %d [%e, %e] {%d, %e, %e, %e}\n",
//              i, start, end, NUM_EQNS, ACCURACY, STEP_GUESS, STEP_MIN);
//       printf("  rho:%e, v:%e, Ys:%e, Ys2:%e, Ncol:%e, Ncol2:%e, T:%e\n",
//              ystart[1], ystart[2], ystart[4], ystart[5], ystart[6], ystart[7], ystart[4]);
// #endif
        
      /* Call the integration routine */
      /* Numerical Recipes' Bulirsch-Stoer Method (adaptive steps)*/
      odeint(ystart, NUM_EQNS, start, end, ACCURACY, STEP_GUESS, STEP_MIN, &nok,
             &nbad, ode_system, bsstep);
      /* Numerical Recipes' Runge-Kutta (not 4) (adaptive steps)*/
//       odeint(ystart, NUM_EQNS, start, end, ACCURACY, STEP_GUESS, STEP_MIN, &nok,
//              &nbad, ode_system, rkqs);
     
       /* Numerical Recipes' Runge-Kutta 4 (Defined number of steps, NSTEPS)*/        
//       rkdumb(ystart, NUM_EQNS, start, end,NSTEPS,ode_system);
//       for (ne=1; ne<NUM_EQNS+1; ne++){
//            ystart[ne] = y[ne][NSTEPS+1];
//       } //neccessary for rkdumb to reset ystart to the output for each step
        

      /* Fix unphysical results */
      for (j=0; j<NSPECIES; j++){
          k=j+4;
          ystart[k] = ystart[k] <= 1.0 ? ystart[k] : 1.0;
          ystart[k] = ystart[k] >= 0.0 ? ystart[k] : 0.0;
      }
        
      /* Set the arrays input to the function integrate_past_sp to the solution */
      set_odeint_soln(INPTS+M+i, end, ystart, equationvars_p);
//        For r[Mach[i]<last Mach] call linearize
    }
  }

  /* Free memory */
//   if (parameters.integrate_outward || parameters.integrate_inward) {
  if (parameters.integrate_outward) {
    free_odeint(ystart);
  }

  return;
}
#undef DEBUG_IN
#undef DEBUG_OUT

/*============================================================================*
 *=========================== PRIVATE FUNCTIONS ==============================*
 *============================================================================*
 * ode_system      - the system of ode to be integrated                       *
 * init_odeint     - initalize the ode solver variables                       *
 * free_odeint     - free all the ode solver arrays                           *
 * set_odeint_soln - map the ode solver solution to solution array            *
 *----------------------------------------------------------------------------*/

/*============================================================================*
 *! \fn static void ode_system(double x, double y[], double dydx[])           *
 *  \brief System of ode that odeint solves                                   *
 *  \note Must be of function form taken by odeint, which is currently        *
 *        derivs(double x, double y[], double dydx[])                         *
 *----------------------------------------------------------------------------*/

static void ode_system(double x, double y[], double dydx[]) {
  I_EQNVARS vars;
    double temp_dydx[NSPECIES];
    int j,k,m;

  vars.r    = x;
  vars.rho  = y[1];
  vars.v    = y[2];
  vars.T    = y[3];
  for (j=0; j<NSPECIES; j++){
      k = j+4;
      m = j+4+NSPECIES;
      vars.Ys[j]  = y[k];
      vars.Ncol[j] = y[m];
  }
  vars.q    = (x-parameters.Rmin)/thez;
  vars.z    = thez;

  /* d(Ncol)/dr */
  get_dNcoldr(temp_dydx, vars);
  for (j=0; j<NSPECIES; j++){
       m = j+4+NSPECIES;
        dydx[m] = temp_dydx[j]; 
   }
    
      /* d(Ys)/dr */
  get_dYsdr(temp_dydx, vars,1e4,1); // Must come second.
  for (j=0; j<NSPECIES; j++){
       k = j+4;
       dydx[k] = temp_dydx[j]; 
   }

  /* dv/dr */
  get_dvdr(&dydx[2], vars);

  /* d(rho)/dr */
  get_drhodr(&dydx[1], vars, dydx[2]);

  /* dT/dr */
  get_dTdr(&dydx[3], vars, dydx[1], temp_dydx,1e4); //temp must = dYsdr array

  return;
}

/*============================================================================*
 *! \fn static void init_odeint(int indx, double **ystart,                    *
 *                              EQNVARS equationvars)                         *
 *  \brief Initalizes the values used for our ode solver routine              *
 *----------------------------------------------------------------------------*/

static void init_odeint(int indx, double **ystart, EQNVARS equationvars) {
  /* Number of variables we hardwired assign by hand */
  /* Take care to also modify ode_system() if altering */
  const int assigned = 3+2*NSPECIES; /*CHANGED for multispecies*/
  int j,k,m;
    

  /* set up variables for odeint */
  kmax  = MAX_STORED_STEPS;
  dxsav = 0.1;

  /* we'll be indexing starting at 1 instead of 0 (gross) */
  *ystart = calloc_1d_array_gross(1, NUM_EQNS);
  xp = calloc_1d_array_gross(1, MAX_STORED_STEPS);
  yp = calloc_2d_array_gross(1, NUM_EQNS, 1, MAX_STORED_STEPS);
  xx = calloc_1d_array_gross(1, NSTEPS+1);
  y = calloc_2d_array_gross(1, NUM_EQNS, 1, NSTEPS+1);

  if (assigned != NUM_EQNS) {
    fprintf(stderr, "ERROR: initializing odeint with (%d) variables, differs "
            "from NUM_EQNS:%d\n", assigned, NUM_EQNS);
    exit(301);
  }
  /* Assign indices 1 thru assigned */
  (*ystart)[1] = equationvars.rho[indx];
  (*ystart)[2] = equationvars.v[indx];
  (*ystart)[3] = equationvars.T[indx];
   for (j=0; j<NSPECIES; j++){
        k = j+4;
        m = j+4+NSPECIES;
        (*ystart)[k] = equationvars.Ys[indx][j];
//        printf("Ncol[%d][%d] = %.2e\n",indx,j,equationvars.Ncol[indx][j]);
        (*ystart)[m] = equationvars.Ncol[indx][j];
    }

  thez = equationvars.z[indx];

  critvars.r    = equationvars.r[INPTS+M-1];
  critvars.rho  = equationvars.rho[INPTS+M-1];
  critvars.v    = equationvars.v[INPTS+M-1];
  critvars.T    = equationvars.T[INPTS+M-1];
   for (j=0; j<NSPECIES; j++){
      critvars.Ys[j]   = equationvars.Ys[INPTS+M-1][j];
      critvars.Ncol[j] = equationvars.Ncol[INPTS+M-1][j];
    }
  critvars.q    = equationvars.q[INPTS+M-1];
  critvars.z    = equationvars.z[INPTS+M-1];

  return;
}

/*============================================================================*
 *! \fn static void free_odeint(double *ystart)                               *
 *  \brief Free memory allocated for the variables used by odeint             *
 *----------------------------------------------------------------------------*/

static void free_odeint(double *ystart) {
  free_1d_array_gross(xp, 1);
  free_1d_array_gross(xx, 1);
  free_1d_array_gross(ystart, 1);
  free_2d_array_gross(yp, 1, 1);
  free_2d_array_gross(y, 1, 1);

  return;
}

/*============================================================================*
 *! \fn static void set_odeint_soln(int i, double end, double *ystart,        *
 *                                  EQNVARS *equationvars_p)                  *
 *  \brief Sets the equationsvar array with the ode solver solution at grid i *
 *----------------------------------------------------------------------------*/

static void set_odeint_soln(int i, double end, double *ystart,
                            EQNVARS *equationvars_p) {
  int j,k,m;
    
  equationvars_p->rho[i]  = ystart[1];
  equationvars_p->v[i]    = ystart[2];
  equationvars_p->T[i]    = ystart[3];
  for (j=0; j<NSPECIES; j++){
      k=j+4;
      m=j+4+NSPECIES;
//       if (ystart[k]>1.0){
//           equationvars_p->Ys[i][j]   = 1.0;  //in xray version having problem with unphysical neutral frac (Ys > 1). HERE
//       }
//       else{
//           equationvars_p->Ys[i][j]   = ystart[k];
//       }
      equationvars_p->Ys[i][j]   = ystart[k];
      equationvars_p->Ncol[i][j] = ystart[m];  
  }
  equationvars_p->r[i]    = end;
  equationvars_p->z[i]    = equationvars_p->z[INPTS];
  equationvars_p->q[i]    = ((equationvars_p->r[i]-equationvars_p->r[INPTS])
                             /equationvars_p->z[i]);
  return;
}
