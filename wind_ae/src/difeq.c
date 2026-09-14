/*============================================================================*
 *! \file difeq.c                                                             *
 *  \brief Equations and derivatives to be solved for by via relaxation.      *
 *         The analytic finite difference Jacobian of the system of equations,*
 *         or s, is burdensome to calculate and code, so numerical results are*
 *         implemented in its place.                                          *
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
/*=========================== SETTING s FUNCTIONS ============================*/
static void set_bcs_base(double **s, double **y, int ne, int indexv[], int jsf);
static void set_bcs_sp(double **s, double **y, int ne, int indexv[], int jsf);
static void set_interior_eqns(double **s, double **y, int ne, int indexv[],
                              int jsf, int k);
/*======================= NUMERICAL FD FOR s FUNCTIONS =======================*/
static void get_bc_evals(VARLIST *varmod_p, double *eqn_eval_p,
                         VARLIST *eqn_pl_p, VARLIST *eqn_mi_p, int k, double *x,
                         double **y, double (*ymod)[2], int eqnnum);
/* Column-by-column helper: perturb one variable, evaluate ALL interior eqns */
static void get_col_evals(int var_idx, int side, double pertscale,
                          int k, double *x, double **y, double (*ymod)[2],
                          double col_deriv[]);
static void zero_ymod(double (*ymod)[2]);

/*============================================================================*
 *=========================== PUBLIC FUNCTIONS ===============================*
 *============================================================================*
 * difeq - Called by solvede, used to evaluate s matrix                       *
 *----------------------------------------------------------------------------*/

/*============================================================================*
 *! \fn void difeq(int k, int k1, int k2, int jsf, int is1, int isf,          *
 *                 int indexv[], int ne, double **s, double **y)              *
 *  \brief used by solvede to evaluate the finite difference Jacobian s       *
 *----------------------------------------------------------------------------*/
#define DIAGNOSTIC_MATRICES 0               //prints the matrices that numerical recipes uses for relaxation method
void difeq(int k, int k1, int k2, int jsf, int is1, int isf, int indexv[],
           int ne, double **s, double **y) {
    
    #if DIAGNOSTIC_MATRICES  //CURRENTLY ONLY SET UP FOR NSPECIES=2 
    char rownames[4+2*NSPECIES_MAX][20] = {             //each row of the matrix corresponds to one
                         "     v_eqn()",            //of the finite difference equations
                         "     rho_eqn()",          //ORDER MATTERS (see Numerical Recipes in C)
                         "     z_eqn()",
                         "     Ncol_eqn()",
                         "     Ncol2_eqn()",
                         "     Ys_eqn()",
                         "     Ys2_eqn()",
                         "     T_eqn()"
                     };
          char base_rownames[4+2*NSPECIES_MAX][15] = {
                     "",
                     "",
                     "",
                     "",
                     "     T_BC",
                     "     Ys_BC",
                     "     Ys2_BC",
                     "     rho_BC"
                 };
          char sp_rownames[4+2*NSPECIES_MAX][20] = {
                 "     Bernoulli_SP",
                 "     v_crit_SP",
                 "     Ncol1_SP",
                 "     Ncol2_sp",
                 "",
                 "",
                 "",
                 ""
             };
     char colnames[100] = "rho  Ys   Ys2   T  Ncol Ncol2   v   z  |   rho   Y   Ys2   T  Ncol Ncol2   v   z \n";
      // each column gives the result of the row's finite difference equation f(x) for different 'x' variables
    #endif
    
  if (k == k1) { /* boundary conditions at first point */
    set_bcs_base(s, y, ne, indexv, jsf);
    #if DIAGNOSTIC_MATRICES
    int rows, cols;
    printf("Base BCs \n");
      printf("%s",colnames);
    for(rows = 1; rows <= ne; rows++){
        for(cols = 1; cols <= 2*ne+1; cols++){
            printf("%.2e ", (s[rows][cols]));
            if (cols==ne){printf("|");}
            if (cols==2*ne){printf("|");}
            if (cols==2*ne+1){printf("%s",base_rownames[rows-1]);}
        }
        printf("\n");
    }
    #endif
    /* Good place to do once per relaxation iteration calculations */
    linearize_dvdr_crit(x, y);
  }
  else if (k > k2) { /* boundary conditions at last point */
    set_bcs_sp(s, y, ne, indexv, jsf);
    #if DIAGNOSTIC_MATRICES
    int rows, cols;
    printf("Sonic Point BCs \n");
    printf("%s",colnames);
    for(rows = 1; rows <= ne; rows++){
        for(cols = 1; cols <= 2*ne+1; cols++){
            printf("%.2e ", (s[rows][cols]));
            if (cols==ne){printf("|");}
            if (cols==2*ne){printf("|");}
            if (cols==2*ne+1){printf("%s",sp_rownames[rows-1]);}
        }
        printf("\n");
    }
    #endif
  }
  else { /* interior points */
    set_interior_eqns(s, y, ne, indexv, jsf, k);
//     printf("k = %d \n",k);
      
    #if DIAGNOSTIC_MATRICES
    int rows, cols;
    if(k==k1+1){
        printf("First Interior Point Matrix \n");
        printf("%s",colnames);
        for(rows = 1; rows <= ne; rows++){
            for(cols = 1; cols <= 2*ne+1; cols++){
                printf("%.2e ", (s[rows][cols]));
                if (cols==ne){printf("|");}
                if (cols==2*ne){printf("|");}
                if (cols==2*ne+1){printf("%s",rownames[rows-1]);}
            }
            printf("\n");
        }
     }
    if(k==k2){
        printf("Last Interior Point Matrix \n");
        printf("%s",colnames);
        for(rows = 1; rows <= ne; rows++){
            for(cols = 1; cols <= 2*ne+1; cols++){
                printf("%.2e ", (s[rows][cols]));
                if (cols==ne){printf("|");}
                if (cols==2*ne){printf("|");}
                if (cols==2*ne+1){printf("%s",rownames[rows-1]);}

            }
            printf("\n");
        }
     }
  #endif
  }

  return;
}

/*============================================================================*
 *=========================== PRIVATE FUNCTIONS ==============================*
 *============================================================================*
 * set_bcs_base      - sets s at the inner boundary                           *
 * set_bcs_sp        - sets s at the outer boundary                           *
 * set_interior_eqns - sets s at k-th interior point (column-by-column)      *
 * get_bc_evals      - calculates the boundary numerical jacobian             *
 * get_col_evals     - perturbs one var, returns all interior eqn derivs      *
 * zero_ymod         - zeroes the ymod matrix                                 *
 *----------------------------------------------------------------------------*/

/*============================================================================*
 *! \fn static int set_bcs_base(double **s, double **y, int ne, int indexv[], *
 *                              int jsf)                                      *
 *  \brief Sets s at the inner boundary  
 * Finite Difference Equations (e.g. E_1(v)) are numbered as follows:
 *     As a function of:
 *
 * Dependent Variables 
 * indexv[1] - v
 * indexv[2] - z     - see Murray-Clay et al. 2009 for more on "z"
 * indexv[3] - rho
 * indexv[4] - T
 * indexv[5] - Ys of species 1
 * indexv[6] - Ys of species 2
 ... etc.
 * indexv[5+NSPECIES] - Ncol of species 1
 * indexv[6+NSPECIES] - Ncol of species 2
 ... etc.
 
 Indexv[#] equals the actual index value of the variable in the y array.
 y is fed to the Numerical Recipes relaxation method ode solver in ext/Nrecipes.
 Per Numerical Recipes in C 2nd ed.:
 The values in y must be ordered such that the variables in the early slots 
 correspond to those for which there is an interior boundary condition 
 (e.g. y[1] = rho, because rho(r=r_min) is a BC; y[2] = Ncol, b/c Ncol(r_min) in BC, etc.).
 In the indev language used in the code the former would be y[indexv[3]]=rho.
 *----------------------------------------------------------------------------*/

static void set_bcs_base(double **s, double **y, int ne, int indexv[],
                         int jsf) {
  /* NE = 5+2*NSPECIES; last BC row = NE */
  int j, i, last = 2*NSPECIES+5;

  /* rho(Rmin) fixed — y[3]=rho */
  s[last][ne+indexv[1]] = 0.0;
  s[last][ne+indexv[2]] = 0.0;
  s[last][ne+indexv[3]] = 1.0;
  s[last][ne+indexv[4]] = 0.0;
  s[last][ne+indexv[5]] = 0.0;
  for (j=0; j<NSPECIES; j++){
      s[last][ne+indexv[j+6]]          = 0.0;
      s[last][ne+indexv[j+6+NSPECIES]] = 0.0;
  }
  s[last][jsf] = y[3][1] - parameters.rho_rmin;

  /* Ys[j](Rmin) fixed — y[6+j]=Ys[j] */
  for (j=0; j<NSPECIES; j++){
      int row = last-j-1;
      s[row][ne+indexv[1]] = 0.0;
      s[row][ne+indexv[2]] = 0.0;
      s[row][ne+indexv[3]] = 0.0;
      s[row][ne+indexv[4]] = 0.0;
      s[row][ne+indexv[5]] = 0.0;
      for (i=0; i<NSPECIES; i++){
          s[row][ne+indexv[i+6+NSPECIES]] = 0.0;
          s[row][ne+indexv[i+6]]          = (i==j) ? 1.0 : 0.0;
      }
      s[row][jsf] = y[6+j][1] - parameters.Ys_rmin[j];
  }

  /* T(Rmin) fixed — y[4]=T */
  s[last-NSPECIES-1][ne+indexv[1]] = 0.0;
  s[last-NSPECIES-1][ne+indexv[2]] = 0.0;
  s[last-NSPECIES-1][ne+indexv[3]] = 0.0;
  s[last-NSPECIES-1][ne+indexv[4]] = 1.0;
  s[last-NSPECIES-1][ne+indexv[5]] = 0.0;
  for (j=0; j<NSPECIES; j++){
      s[last-NSPECIES-1][ne+indexv[j+6]]          = 0.0;
      s[last-NSPECIES-1][ne+indexv[j+6+NSPECIES]] = 0.0;
  }
  s[last-NSPECIES-1][jsf] = y[4][1] - parameters.T_rmin;

  /* F(Rmin) = 0 — y[5]=F */
  {
    int row = last-NSPECIES-2;
    s[row][ne+indexv[1]] = 0.0;
    s[row][ne+indexv[2]] = 0.0;
    s[row][ne+indexv[3]] = 0.0;
    s[row][ne+indexv[4]] = 0.0;
    s[row][ne+indexv[5]] = 1.0;
    for (j=0; j<NSPECIES; j++){
        s[row][ne+indexv[j+6]]          = 0.0;
        s[row][ne+indexv[j+6+NSPECIES]] = 0.0;
    }
    s[row][jsf] = y[5][1];   /* residual: F[1] - 0 = 0 */
  }

  return;
}

/*============================================================================*
 *! \fn static void set_bcs_sp(double **s, double **y, int ne, int indexv[],  *
 *                             int jsf)                                       *
 *  \brief Sets s at the outer boundary                                       *
 *----------------------------------------------------------------------------*/

static void set_bcs_sp(double **s, double **y, int ne, int indexv[], int jsf) {
  double eqn;
  VARLIST eqn_pl, eqn_mi;
  VARLIST varmod;
  double ymod[NE_MAX+1][2];
  int j,i;

  /* Numerator (Bernoulli) condition at the sonic point */
  get_bc_evals(&varmod, &eqn, &eqn_pl, &eqn_mi, M, x, y, ymod, SPCRITEQN);

  s[1][ne+indexv[3]] = (eqn_pl.rho-eqn_mi.rho)/varmod.rho;
  s[1][ne+indexv[1]] = (eqn_pl.v-eqn_mi.v)/varmod.v;
  s[1][ne+indexv[2]] = (eqn_pl.z-eqn_mi.z)/varmod.z;
  s[1][ne+indexv[4]] = (eqn_pl.T-eqn_mi.T)/varmod.T;
  s[1][ne+indexv[5]] = 0.0;  /* F: spcrit_eqn has no dependence on F */
  for (j=0; j<NSPECIES; j++){
      s[1][ne+indexv[j+6]] = (eqn_pl.Ys[j]-eqn_mi.Ys[j])/varmod.Ys[j];
      s[1][ne+indexv[j+6+NSPECIES]] = (eqn_pl.Ncol[j]-eqn_mi.Ncol[j])/varmod.Ncol[j];
//         printf("Bernoulli condition: Ncol: (%.17e - %.17e) / %.17e     species=%d\n",eqn_pl.Ncol[j],eqn_mi.Ncol[j],varmod.Ncol[j],j);
  }
  s[1][jsf] = eqn;
//       printf("Bernoulli condtion: residual: %.17e\n", eqn);

  /* velocity condition at the sonic point */
  get_bc_evals(&varmod, &eqn, &eqn_pl, &eqn_mi, M, x, y, ymod, SPVEQN);

  s[2][ne+indexv[3]] = (eqn_pl.rho-eqn_mi.rho)/varmod.rho;
  s[2][ne+indexv[1]] = (eqn_pl.v-eqn_mi.v)/varmod.v;
  s[2][ne+indexv[2]] = (eqn_pl.z-eqn_mi.z)/varmod.z;
  s[2][ne+indexv[4]] = (eqn_pl.T-eqn_mi.T)/varmod.T;
  s[2][ne+indexv[5]] = 0.0;  /* F: spv_eqn has no dependence on F */
  for (j=0; j<NSPECIES; j++){
      s[2][ne+indexv[j+6]] = (eqn_pl.Ys[j]-eqn_mi.Ys[j])/varmod.Ys[j];
      s[2][ne+indexv[j+6+NSPECIES]] = (eqn_pl.Ncol[j]-eqn_mi.Ncol[j])/varmod.Ncol[j];
//       printf("Velocity condition: Ncol: (%.17e - %.17e) / %.17e     species=%d\n",eqn_pl.Ncol[j],eqn_mi.Ncol[j],varmod.Ncol[j],j);
  }
  s[2][jsf] = eqn;
//   printf("Velocity condtion: residual: %.17e\n", eqn);
  /* When BREEZEPARAM = 1, we have a transsonic solution;
     when between 0 and 1, a breeze */

  /* the column density is fixed at the sonic point for H */
  for (j=0; j<NSPECIES; j++){
      s[3+j][ne+indexv[1]] = 0.0;
      s[3+j][ne+indexv[2]] = 0.0;
      s[3+j][ne+indexv[3]] = 0.0;
      s[3+j][ne+indexv[4]] = 0.0;
      s[3+j][ne+indexv[5]] = 0.0;  /* F: Ncol BC has no F dependence */
      for (i=0; i<NSPECIES; i++){
          s[3+j][ne+indexv[i+6]]          = 0.0;
          s[3+j][ne+indexv[i+6+NSPECIES]] = (i==j) ? 1.0 : 0.0;
      }
     s[3+j][jsf] = y[j+6+NSPECIES][M]-parameters.Ncol_sp[j];

  }
//     printf("y[5][M] = %.4e \n",y[5][M]);

  return;
}

/*============================================================================*
 *! \fn static void set_interior_eqns(double **s, double **y, int ne,         *
 *                                    int indexv[], int jsf, int k)           *
 *  \brief Sets s at the k-th interior point                                  *
 *----------------------------------------------------------------------------*/

static void set_interior_eqns(double **s, double **y, int ne, int indexv[],
                              int jsf, int k) {
  /*------------------------------------------------------------------------*
   * Column-by-column Jacobian strategy                                      *
   *                                                                         *
   * Old code (row-by-row): for each of (3+2*N) equations, call             *
   *   get_eqn_evals → (4+8*N) eval_eqn calls, each potentially hitting    *
   *   calc_gql_rates.  Total: O(N^2) GLQ calls.                            *
   *                                                                         *
   * New code (column-by-column): for each of (4+4*N) variable perturbations*
   *   call get_col_evals → 2 eval_interior_eqns calls.                    *
   *   eval_interior_eqns evaluates ALL equations in one pass with one GLQ  *
   *   result (shared by all eqns via the calc_gql_rates cache).            *
   *   Total GLQ calls: O(N) — one per Ncol perturbation direction.         *
   *                                                                         *
   * The Ncol rows are filled analytically (no numerical perturbation):     *
   *   Ncol_eqn[j] = Ncol[j][k] - Ncol[j][k-1]                            *
   *               - dNdr[j](avgvars) * delta_r                             *
   *   dNdr[j] = -Ys[j] * n_abs_j(rho, Ys, chain) * Rp/NCOL0              *
   * This eliminates 4*N eval_interior_eqns calls and their GLQ misses.    *
   *------------------------------------------------------------------------*/

  double ymod[NE_MAX+1][2];
  /* col_deriv[row] = (f(+h) - f(-h)) from one column perturbation.
   * Rows: 0=v, 1=rho, 2..2+N-1=Ncol[j], 2+N..2+2N-1=ion[j], 2+2N=T, 3+2N=F */
  double col_pl[3+2*NSPECIES_MAX+1], col_mi[3+2*NSPECIES_MAX+1];
  double res0[3+2*NSPECIES_MAX+1]; /* unperturbed residuals */
  /* Analytic partial derivatives of dNdr[j] wrt each variable:             */
  double dNdr[NSPECIES_MAX];       /* dNdr at avg point (unperturbed)          */
  double dNdr_drho[NSPECIES_MAX];  /* d(dNdr[j])/d(rho_avg)                   */
  double dNdr_dYs[NSPECIES_MAX][NSPECIES_MAX]; /* d(dNdr[j])/d(Ys_avg[i])         */
  I_EQNVARS avgvars;
  double ravg, delta_r, delta_z, h;
  int i, j, last=5+2*NSPECIES, parent;  /* NE = 5+2*NSPECIES */

  /* Build average variables (same as before) */
  avgvars.q     = 0.5*(x[k]+x[k-1]);
  avgvars.v     = 0.5*(y[1][k]+y[1][k-1]);
  avgvars.z     = 0.5*(y[2][k]+y[2][k-1]);
  avgvars.rho   = 0.5*(y[3][k]+y[3][k-1]);
  avgvars.T     = 0.5*(y[4][k]+y[4][k-1]);
  for (j=0; j<NSPECIES; j++){
      avgvars.Ys[j]    = 0.5*(y[j+6][k]+y[j+6][k-1]);
      avgvars.Ncol[j]  = 0.5*(y[j+6+NSPECIES][k]+y[j+6+NSPECIES][k-1]);
  }
  get_rad(&ravg, &avgvars);
  avgvars.r    = ravg;

  /* delta_r = z_avg * (q[k] - q[k-1]) — used by analytic Ncol partials */
  delta_r = avgvars.z * (x[k] - x[k-1]);
  delta_z = x[k] - x[k-1];   /* for the z-column's effect on delta_r */

  /*----------------------------------------------------------------------*
   * Pre-compute unperturbed residuals (one GLQ call, sets the cache)     *
   *----------------------------------------------------------------------*/
  zero_ymod(ymod);
  eval_interior_eqns(k, x, y, ymod, res0);

  /*----------------------------------------------------------------------*
   * Pre-compute analytic dNdr and its partial derivatives                *
   *                                                                       *
   * dNdr[j] = -Ys_avg[j] * n_abs_j(rho_avg, Ys_avg) * Rp/NCOL0         *
   *                                                                       *
   * For an independent species j (chain_parent[j] < 0):                  *
   *   n_abs_j/rho = HX[j]/atomic_mass[j]   (constant)                   *
   *   d(dNdr[j])/d(rho_avg) = -Ys[j]*HX[j]/m[j]*Rp/NCOL0*(RHO0)        *
   *   d(dNdr[j])/d(Ys[i])   = 0  if i!=j                                *
   *   d(dNdr[j])/d(Ys[j])   = -HX[j]/m[j]*rho*RHO0*Rp/NCOL0            *
   *                                                                       *
   * For a child species j (chain_parent[j] = p):                         *
   *   n_abs_j/rho = Ys[j]*(1-Ys[p])*HX[p]/m[p]                         *
   *   d(dNdr[j])/d(rho_avg) = -Ys[j]*(1-Ys[p])*HX[p]/m[p]*RHO0*Rp/NCOL0 *
   *   d(dNdr[j])/d(Ys[j])   = -(1-Ys[p])*HX[p]/m[p]*rho*RHO0*Rp/NCOL0 *
   *   d(dNdr[j])/d(Ys[p])   = +Ys[j]*HX[p]/m[p]*rho*RHO0*Rp/NCOL0     *
   *   d(dNdr[j])/d(Ys[i])   = 0  for any other i                        *
   *----------------------------------------------------------------------*/
  get_dNcoldr(dNdr, &avgvars);   /* fills dNdr[j] = the unperturbed value  */

  for (j=0; j<NSPECIES; j++){
      parent = parameters.chain_parent[j];
      if (parent < 0) {
          /* Independent species */
          double nj_over_rho = parameters.HX[j]/parameters.atomic_mass[j];
          dNdr_drho[j] = -avgvars.Ys[j] * nj_over_rho * RHO0 * parameters.Rp/NCOL0;
          for (i=0; i<NSPECIES; i++)
              dNdr_dYs[j][i] = 0.0;
          dNdr_dYs[j][j] = -nj_over_rho * avgvars.rho * RHO0 * parameters.Rp/NCOL0;
      } else {
          /* Child species: effective density pool scales with (1-Ys[parent]) */
          double np_over_rho = parameters.HX[parent]/parameters.atomic_mass[parent];
          double fac = (1.0 - avgvars.Ys[parent]) * np_over_rho;
          dNdr_drho[j] = -avgvars.Ys[j] * fac * RHO0 * parameters.Rp/NCOL0;
          for (i=0; i<NSPECIES; i++)
              dNdr_dYs[j][i] = 0.0;
          /* d/d(Ys[j]): coefficient of Ys[j] in n_abs is (1-Ys[p])*np_over_rho*rho */
          dNdr_dYs[j][j]      = -fac * avgvars.rho * RHO0 * parameters.Rp/NCOL0;
          /* d/d(Ys[parent]): coefficient of (1-Ys[p]) differentiated gives -Ys[j]*np_over_rho*rho */
          dNdr_dYs[j][parent] =  avgvars.Ys[j] * np_over_rho
                                  * avgvars.rho * RHO0 * parameters.Rp/NCOL0;
      }
  }

  /*----------------------------------------------------------------------*
   * Macro: given col_pl[] and col_mi[] from a column perturbation of     *
   * step h, assign the numerical derivatives to all NON-NCOL rows of s   *
   * at column position col_s.                                             *
   *----------------------------------------------------------------------*/
#define ASSIGN_NUM_COL(col_s) \
  do { \
    s[1][col_s]    = (col_pl[0] - col_mi[0]) / h; \
    s[2][col_s]    = (col_pl[1] - col_mi[1]) / h; \
    for (j=0; j<NSPECIES; j++) { \
      s[4+j+NSPECIES][col_s] = (col_pl[2+NSPECIES+j] - col_mi[2+NSPECIES+j]) / h; \
    } \
    s[last-1][col_s] = (col_pl[2+2*NSPECIES] - col_mi[2+2*NSPECIES]) / h; \
    s[last][col_s]   = (col_pl[3+2*NSPECIES] - col_mi[3+2*NSPECIES]) / h; \
  } while(0)

  /*----------------------------------------------------------------------*
   * Residuals (jsf column) — computed from the already-cached res0[]    *
   *----------------------------------------------------------------------*/
  s[1][jsf]    = res0[0];
  s[2][jsf]    = res0[1];
  for (j=0; j<NSPECIES; j++) {
      s[4+j][jsf]          = res0[2+j];
      s[4+j+NSPECIES][jsf] = res0[2+NSPECIES+j];
  }
  s[last-1][jsf] = res0[2+2*NSPECIES];  /* T_eqn residual */
  s[last][jsf]   = res0[3+2*NSPECIES];  /* F_eqn residual */

  /*----------------------------------------------------------------------*
   * z-equation: fully analytic, no perturbation needed                   *
   *   z[k] - z[k-1] = 0  =>  s[3][indexv[2]] = -1, s[3][ne+indexv[2]] = 1
   *----------------------------------------------------------------------*/
  for (i=1; i<=2*ne+1; i++) s[3][i] = 0.0;
  s[3][indexv[2]]    = -1.0;
  s[3][ne+indexv[2]] =  1.0;
  s[3][jsf]          =  y[2][k]-y[2][k-1];

  /*----------------------------------------------------------------------*
   * Column by column: numerical derivatives for v, rho, ion, T rows     *
   * Column order: v[k-1], z[k-1], rho[k-1], T[k-1],                    *
   *               Ys[0..N-1][k-1], Ys[0..N-1][k],                       *
   *               v[k], z[k], rho[k], T[k]                              *
   * Ncol columns are handled analytically below (no eval calls).         *
   *----------------------------------------------------------------------*/

  /* --- v[k-1] column --- */
  h = avgvars.v / DERIVDIV;
  zero_ymod(ymod); ymod[1][0] =  h/2.0; eval_interior_eqns(k,x,y,ymod,col_pl);
  zero_ymod(ymod); ymod[1][0] = -h/2.0; eval_interior_eqns(k,x,y,ymod,col_mi);
  ASSIGN_NUM_COL(indexv[1]);
  /* analytic Ncol: dNdr does not depend on v */
  for (j=0; j<NSPECIES; j++) s[4+j][indexv[1]] = 0.0;

  /* --- z[k-1] column ---
   * z affects delta_r = z_avg*(q[k]-q[k-1]), but the Ncol equation's
   * dependence on z comes only through delta_r.
   * d(Ncol_eqn[j])/d(z[k-1]) = -dNdr[j] * delta_z * 0.5
   * (chain rule: d(delta_r)/d(z[k-1]) = 0.5*(q[k]-q[k-1]) = 0.5*delta_z) */
  h = avgvars.z / DERIVDIV;
  zero_ymod(ymod); ymod[2][0] =  h/2.0; eval_interior_eqns(k,x,y,ymod,col_pl);
  zero_ymod(ymod); ymod[2][0] = -h/2.0; eval_interior_eqns(k,x,y,ymod,col_mi);
  ASSIGN_NUM_COL(indexv[2]);
  for (j=0; j<NSPECIES; j++)
      s[4+j][indexv[2]] = -dNdr[j] * 0.5 * delta_z;

  /* --- rho[k-1] column ---
   * d(Ncol_eqn[j])/d(rho[k-1]) = -d(dNdr[j])/d(rho_avg) * delta_r * 0.5 */
  h = avgvars.rho / DERIVDIV;
  zero_ymod(ymod); ymod[3][0] =  h/2.0; eval_interior_eqns(k,x,y,ymod,col_pl);
  zero_ymod(ymod); ymod[3][0] = -h/2.0; eval_interior_eqns(k,x,y,ymod,col_mi);
  ASSIGN_NUM_COL(indexv[3]);
  for (j=0; j<NSPECIES; j++)
      s[4+j][indexv[3]] = -dNdr_drho[j] * delta_r * 0.5;

  /* --- T[k-1] column --- (dNdr has no T dependence) */
  h = avgvars.T / DERIVDIV;
  zero_ymod(ymod); ymod[4][0] =  h/2.0; eval_interior_eqns(k,x,y,ymod,col_pl);
  zero_ymod(ymod); ymod[4][0] = -h/2.0; eval_interior_eqns(k,x,y,ymod,col_mi);
  ASSIGN_NUM_COL(indexv[4]);
  for (j=0; j<NSPECIES; j++) s[4+j][indexv[4]] = 0.0;

  /* --- Ys[i][k-1] columns (i = 0..NSPECIES-1) --- */
  for (i=0; i<NSPECIES; i++) {
      h = avgvars.Ys[i] / DERIVDIV;
      zero_ymod(ymod); ymod[i+6][0] =  h/2.0; eval_interior_eqns(k,x,y,ymod,col_pl);
      zero_ymod(ymod); ymod[i+6][0] = -h/2.0; eval_interior_eqns(k,x,y,ymod,col_mi);
      ASSIGN_NUM_COL(indexv[i+6]);
      /* analytic Ncol: d(Ncol_eqn[j])/d(Ys[i][k-1]) = -dNdr_dYs[j][i]*delta_r*0.5 */
      for (j=0; j<NSPECIES; j++)
          s[4+j][indexv[i+6]] = -dNdr_dYs[j][i] * delta_r * 0.5;
  }

  /* --- v[k] column --- */
  h = avgvars.v / DERIVDIV;
  zero_ymod(ymod); ymod[1][1] =  h/2.0; eval_interior_eqns(k,x,y,ymod,col_pl);
  zero_ymod(ymod); ymod[1][1] = -h/2.0; eval_interior_eqns(k,x,y,ymod,col_mi);
  ASSIGN_NUM_COL(ne+indexv[1]);
  for (j=0; j<NSPECIES; j++) s[4+j][ne+indexv[1]] = 0.0;

  /* --- z[k] column --- */
  h = avgvars.z / DERIVDIV;
  zero_ymod(ymod); ymod[2][1] =  h/2.0; eval_interior_eqns(k,x,y,ymod,col_pl);
  zero_ymod(ymod); ymod[2][1] = -h/2.0; eval_interior_eqns(k,x,y,ymod,col_mi);
  ASSIGN_NUM_COL(ne+indexv[2]);
  for (j=0; j<NSPECIES; j++)
      s[4+j][ne+indexv[2]] = -dNdr[j] * 0.5 * delta_z;

  /* --- rho[k] column --- */
  h = avgvars.rho / DERIVDIV;
  zero_ymod(ymod); ymod[3][1] =  h/2.0; eval_interior_eqns(k,x,y,ymod,col_pl);
  zero_ymod(ymod); ymod[3][1] = -h/2.0; eval_interior_eqns(k,x,y,ymod,col_mi);
  ASSIGN_NUM_COL(ne+indexv[3]);
  for (j=0; j<NSPECIES; j++)
      s[4+j][ne+indexv[3]] = -dNdr_drho[j] * delta_r * 0.5;

  /* --- T[k] column --- */
  h = avgvars.T / DERIVDIV;
  zero_ymod(ymod); ymod[4][1] =  h/2.0; eval_interior_eqns(k,x,y,ymod,col_pl);
  zero_ymod(ymod); ymod[4][1] = -h/2.0; eval_interior_eqns(k,x,y,ymod,col_mi);
  ASSIGN_NUM_COL(ne+indexv[4]);
  for (j=0; j<NSPECIES; j++) s[4+j][ne+indexv[4]] = 0.0;

  /* --- Ys[i][k] columns --- */
  for (i=0; i<NSPECIES; i++) {
      h = avgvars.Ys[i] / DERIVDIV;
      zero_ymod(ymod); ymod[i+6][1] =  h/2.0; eval_interior_eqns(k,x,y,ymod,col_pl);
      zero_ymod(ymod); ymod[i+6][1] = -h/2.0; eval_interior_eqns(k,x,y,ymod,col_mi);
      ASSIGN_NUM_COL(ne+indexv[i+6]);
      for (j=0; j<NSPECIES; j++)
          s[4+j][ne+indexv[i+6]] = -dNdr_dYs[j][i] * delta_r * 0.5;
  }

  /*----------------------------------------------------------------------*
   * Ncol[i][k-1] and Ncol[i][k] columns                                 *
   *                                                                       *
   * dNdr[j] does not depend on Ncol at all, so the Ncol row is always 0  *
   * for the Ncol Jacobian columns (off-diagonal blocks).                  *
   * For all non-Ncol rows (v, rho, ion, T): these also depend on Ncol    *
   * through glq_rates, so we still need numerical derivatives there.      *
   * We call eval_interior_eqns here; the GLQ cache fires because a Ncol  *
   * perturbation IS a cache miss — but with column-by-column we pay this  *
   * cost once per perturbation, shared across all equations.              *
   *----------------------------------------------------------------------*/
  for (i=0; i<NSPECIES; i++) {
      /* Ncol[i][k-1] */
      h = avgvars.Ncol[i] / DERIVDIV;
      zero_ymod(ymod); ymod[i+6+NSPECIES][0] =  h/2.0; eval_interior_eqns(k,x,y,ymod,col_pl);
      zero_ymod(ymod); ymod[i+6+NSPECIES][0] = -h/2.0; eval_interior_eqns(k,x,y,ymod,col_mi);
      ASSIGN_NUM_COL(indexv[i+6+NSPECIES]);
      /* Ncol rows: Ncol_eqn[j] has no dependence on Ncol at all */
      for (j=0; j<NSPECIES; j++) s[4+j][indexv[i+6+NSPECIES]] = 0.0;

      /* Ncol[i][k] */
      zero_ymod(ymod); ymod[i+6+NSPECIES][1] =  h/2.0; eval_interior_eqns(k,x,y,ymod,col_pl);
      zero_ymod(ymod); ymod[i+6+NSPECIES][1] = -h/2.0; eval_interior_eqns(k,x,y,ymod,col_mi);
      ASSIGN_NUM_COL(ne+indexv[i+6+NSPECIES]);
      for (j=0; j<NSPECIES; j++) s[4+j][ne+indexv[i+6+NSPECIES]] = 0.0;
  }

  /* Diagonal of Ncol rows: d(Ncol_eqn[j])/d(Ncol[j][k-1]) = -1,
   *                        d(Ncol_eqn[j])/d(Ncol[j][k])   = +1  */
  for (j=0; j<NSPECIES; j++) {
      s[4+j][indexv[j+6+NSPECIES]]    = -1.0;
      s[4+j][ne+indexv[j+6+NSPECIES]] =  1.0;
  }

  /* --- F[k-1] column --- */
  h = (fabs(y[5][k-1]) > 1e-20) ? fabs(y[5][k-1]) / DERIVDIV : FSCALE / DERIVDIV;
  zero_ymod(ymod); ymod[5][0] =  h/2.0; eval_interior_eqns(k,x,y,ymod,col_pl);
  zero_ymod(ymod); ymod[5][0] = -h/2.0; eval_interior_eqns(k,x,y,ymod,col_mi);
  ASSIGN_NUM_COL(indexv[5]);
  for (j=0; j<NSPECIES; j++) s[4+j][indexv[5]] = 0.0;

  /* --- F[k] column --- */
  zero_ymod(ymod); ymod[5][1] =  h/2.0; eval_interior_eqns(k,x,y,ymod,col_pl);
  zero_ymod(ymod); ymod[5][1] = -h/2.0; eval_interior_eqns(k,x,y,ymod,col_mi);
  ASSIGN_NUM_COL(ne+indexv[5]);
  for (j=0; j<NSPECIES; j++) s[4+j][ne+indexv[5]] = 0.0;

#undef ASSIGN_NUM_COL

  return;
}

/*============================================================================*
 *! \fn static void get_bc_evals(VARLIST *varmod_p, double *eqn_eval_p,       *
 *                               VARLIST *eqn_pl_p, VARLIST *eqn_mi_p, int k, *
 *                               double *x, double **y, double (*ymod)[2], *
 *                               int eqnnum)                                  *
 *  \brief Calculates the numerical derivative of the given equation for each *
 *         variable at the boundaries                                         *
 *----------------------------------------------------------------------------*/

static void get_bc_evals(VARLIST *varmod_p, double *eqn_eval_p,
                         VARLIST *eqn_pl_p, VARLIST *eqn_mi_p, int k, double *x,
                         double **y, double (*ymod)[2], int eqnnum) {
  int j;
    
  zero_ymod(ymod);
  *eqn_eval_p = eval_eqn(k, x, y, ymod, eqnnum,1);

  /* ddrho */
  varmod_p->rho = y[3][k]/DERIVDIV;
  zero_ymod(ymod);
  ymod[3][1] = varmod_p->rho/2.0;
  eqn_pl_p->rho = eval_eqn(k, x, y, ymod, eqnnum,1); //Rho is not species dependent so only care about one value.
  zero_ymod(ymod);
  ymod[3][1] = -varmod_p->rho/2.0;
  eqn_mi_p->rho = eval_eqn(k, x, y, ymod, eqnnum,1);

  /* ddv */
  varmod_p->v = y[1][k]/DERIVDIV;
  zero_ymod(ymod);
  ymod[1][1]  = varmod_p->v/2.0;
  eqn_pl_p->v = eval_eqn(k, x, y, ymod, eqnnum,1);
  zero_ymod(ymod);
  ymod[1][1]  = -varmod_p->v/2.0;
  eqn_mi_p->v = eval_eqn( k, x, y, ymod, eqnnum,1);

  /* ddT */
  varmod_p->T = y[4][k]/DERIVDIV;
  zero_ymod(ymod);
  ymod[4][1]  = varmod_p->T/2.0;
  eqn_pl_p->T = eval_eqn(k, x, y, ymod, eqnnum,1);
  zero_ymod(ymod);
  ymod[4][1]  = -varmod_p->T/2.0;
  eqn_mi_p->T = eval_eqn(k, x, y, ymod, eqnnum,1);

  /* ddYs */
  for (j=0; j<NSPECIES; j++){
      varmod_p->Ys[j] = y[j+6][k]/DERIVDIV;
      zero_ymod(ymod);
      ymod[j+6][1]   = varmod_p->Ys[j]/2.0;
      eqn_pl_p->Ys[j] = eval_eqn( k, x, y, ymod, eqnnum,1); //FIX, this should maybe not be j+1, because all that species # controls is 
      //whether or not the EQNNUM plugged into is, e.g., IONEQN1 or IONEQN2 (a.k.a., as a function of Ys1 or Ys2). Since this function 
      //has no species dependence, we don't need to change the species number.
      //Technically, it shouldn't make a difference since IONEQN and NCOLEQN are never used with this function and it is meaningless to
      //have species set for the other EQNNUMs, since they have no species dependence.
      zero_ymod(ymod);
      ymod[j+6][1]   = -varmod_p->Ys[j]/2.0; //what is really making the difference for the species 1 and 2 here is calling Ys[j]
      eqn_mi_p->Ys[j] = eval_eqn( k, x, y, ymod, eqnnum,1);
//         printf("eqn_mi, j=%d: %.2e\n", j, eval_eqn( k, x, y, ymod, eqnnum,j+1));

  }

  /* ddNcol */
  for (j=0; j<NSPECIES; j++){
      varmod_p->Ncol[j] = y[j+6+NSPECIES][k]/DERIVDIV;
      zero_ymod(ymod);
      ymod[j+6+NSPECIES][1]   = varmod_p->Ncol[j]/2.0;
      eqn_pl_p->Ncol[j] = eval_eqn( k, x, y, ymod, eqnnum, 1);
      zero_ymod(ymod);
      ymod[j+6+NSPECIES][1]   = -varmod_p->Ncol[j]/2.0;
      eqn_mi_p->Ncol[j] = eval_eqn( k, x, y, ymod, eqnnum, 1);
  } 
    
  /* ddz */
  varmod_p->z = y[2][k]/DERIVDIV;
  zero_ymod(ymod);
  ymod[2][1]  = varmod_p->z/2.0;
  eqn_pl_p->z = eval_eqn( k, x, y, ymod, eqnnum,1);
  zero_ymod(ymod);
  ymod[2][1]  = -varmod_p->z/2.0;
  eqn_mi_p->z = eval_eqn( k, x, y, ymod, eqnnum,1);
    
  return;
}

/*============================================================================*
 *! \fn static void get_col_evals(int var_idx, int side, double h,            *
 *                                int k, double *x, double **y,               *
 *                                double (*ymod)[2], double res[])         *
 *  \brief Perturbs variable ymod[var_idx][side] by h/2 (plus and minus),    *
 *         evaluates ALL interior equations in each perturbed state, and      *
 *         returns the signed difference (pl - mi) in res[].                 *
 *                                                                            *
 *  This is the column-by-column analogue of the old get_eqn_evals.  One     *
 *  call here replaces (3+2*NSPECIES) separate get_eqn_evals invocations for *
 *  the same variable perturbation, and crucially shares the GLQ cache result *
 *  across all equations.                                                     *
 *                                                                            *
 *  Caller is responsible for dividing res[] by h to obtain the derivative.  *
 *----------------------------------------------------------------------------*/
static void get_col_evals(int var_idx, int side, double h,
                          int k, double *x, double **y,
                          double (*ymod)[2], double res[]) {
  /* This function is now unused since set_interior_eqns inlines the calls
   * to eval_interior_eqns directly.  It is retained here as a reference
   * implementation and in case future refactoring wants a helper. */
  double pl[3+2*NSPECIES_MAX+1], mi[3+2*NSPECIES_MAX+1];
  int r;

  zero_ymod(ymod);
  ymod[var_idx][side] =  h/2.0;
  eval_interior_eqns(k, x, y, ymod, pl);

  zero_ymod(ymod);
  ymod[var_idx][side] = -h/2.0;
  eval_interior_eqns(k, x, y, ymod, mi);

  for (r = 0; r < 2+2*NSPECIES+1; r++)
    res[r] = pl[r] - mi[r];
}


/*============================================================================*
 *! \fn static void zero_ymod(double (*ymod)[2])                           *
 *  \brief Zeroes the ymod matrix                                             *
 *----------------------------------------------------------------------------*/

static void zero_ymod(double (*ymod)[2]) {
  int i, j;

  for (i = 0; i <= NE; i++) {
    for (j = 0; j < 2; j++) {
      ymod[i][j] = 0;
    }
  }

  return;
}
