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
                         double **y, double ymod[NE+1][2], int eqnnum);
static void get_eqn_evals(VARLIST *varmod_p, double *eqn_eval_p,
                          VARLIST *eqn_km1_pl_p, VARLIST *eqn_km1_mi_p,
                          VARLIST *eqn_k_pl_p, VARLIST *eqn_k_mi_p,
                          I_EQNVARS avgvars, int k, double *x, double **y,
                          double ymod[NE+1][2], int eqnnum, int species);
static void zero_ymod(double ymod[NE+1][2]);

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
    char rownames[4+2*NSPECIES][20] = {             //each row of the matrix corresponds to one
                         "     v_eqn()",            //of the finite difference equations
                         "     rho_eqn()",          //ORDER MATTERS (see Numerical Recipes in C)
                         "     z_eqn()",
                         "     Ncol_eqn()",
                         "     Ncol2_eqn()",
                         "     Ys_eqn()",
                         "     Ys2_eqn()",
                         "     T_eqn()"
                     };
          char base_rownames[4+2*NSPECIES][15] = {
                     "",
                     "",
                     "",
                     "",
                     "     T_BC",
                     "     Ys_BC",
                     "     Ys2_BC",
                     "     rho_BC"
                 };
          char sp_rownames[4+2*NSPECIES][20] = {
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
 * set_interior_eqns - sets s at k-th interior point                          *
 * get_bc_evals      - calculates the boundary numerical jacobian             *
 * get_eqn_evals     - calculates the interior numerical jacobian             *
 * eval_eqn          - calls calculation of the residual                      *
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
  int j,i,last=2*NSPECIES+4;
    
  /* the density at the base is fixed */
  s[last][ne+indexv[1]] = 0.0;
  s[last][ne+indexv[2]] = 0.0;
  s[last][ne+indexv[3]] = 1.0;
  s[last][ne+indexv[4]] = 0.0;
  for (j=0; j<NSPECIES; j++){
      s[last][ne+indexv[j+5]] = 0.0; /*Iterating through Ys for all species*/
      s[last][ne+indexv[j+5+NSPECIES]] = 0.0; /*Iterating through Ncol for all species*/
  }
  s[last][jsf] = y[3][1]-parameters.rho_rmin; /*first y[#] index should be the same as indexv[#]*/

  /* the ionization fraction is fixed at the base for all species*/
  for (j=0;j<NSPECIES; j++){
      s[last-j-1][ne+indexv[1]] = 0.0;
      s[last-j-1][ne+indexv[2]] = 0.0;
      s[last-j-1][ne+indexv[3]] = 0.0;
      s[last-j-1][ne+indexv[4]] = 0.0;
      for (i=0; i<NSPECIES; i++){
          s[last-j-1][ne+indexv[i+5+NSPECIES]] = 0.0;
          if (i==j){
              s[last-j-1][ne+indexv[i+5]] = 1.0; /*setting s(Ys)=1 for all species*/
          }
          else{
              s[last-j-1][ne+indexv[i+5]] = 0.0;
          }
      }
      s[last-j-1][jsf] = y[j+5][1]-parameters.Ys_rmin[j];
  }


  /* the temperature is fixed at the base */
  s[last-NSPECIES-1][ne+indexv[1]] = 0.0;
  s[last-NSPECIES-1][ne+indexv[2]] = 0.0;
  s[last-NSPECIES-1][ne+indexv[3]] = 0.0;
  s[last-NSPECIES-1][ne+indexv[4]] = 1.0;
  for (j=0; j<NSPECIES; j++){
      s[last-NSPECIES-1][ne+indexv[j+5]] = 0.0;
      s[last-NSPECIES-1][ne+indexv[j+5+NSPECIES]] = 0.0;
  }
  s[last-NSPECIES-1][jsf] = y[4][1]-parameters.T_rmin;

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
  double ymod[NE+1][2];
  int j,i;

  /* Numerator (Bernoulli) condition at the sonic point */
  get_bc_evals(&varmod, &eqn, &eqn_pl, &eqn_mi, M, x, y, ymod, SPCRITEQN);

  s[1][ne+indexv[3]] = (eqn_pl.rho-eqn_mi.rho)/varmod.rho;
  s[1][ne+indexv[1]] = (eqn_pl.v-eqn_mi.v)/varmod.v;
  s[1][ne+indexv[2]] = (eqn_pl.z-eqn_mi.z)/varmod.z;
  s[1][ne+indexv[4]] = (eqn_pl.T-eqn_mi.T)/varmod.T;
  for (j=0; j<NSPECIES; j++){
      s[1][ne+indexv[j+5]] = (eqn_pl.Ys[j]-eqn_mi.Ys[j])/varmod.Ys[j];
      s[1][ne+indexv[j+5+NSPECIES]] = (eqn_pl.Ncol[j]-eqn_mi.Ncol[j])/varmod.Ncol[j];
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
  for (j=0; j<NSPECIES; j++){
      s[2][ne+indexv[j+5]] = (eqn_pl.Ys[j]-eqn_mi.Ys[j])/varmod.Ys[j];
      s[2][ne+indexv[j+5+NSPECIES]] = (eqn_pl.Ncol[j]-eqn_mi.Ncol[j])/varmod.Ncol[j];
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
      for (i=0; i<NSPECIES; i++){
          s[3+j][ne+indexv[i+5]] = 0.0;
          if (i==j){
              s[3+j][ne+indexv[i+5+NSPECIES]] = 1.0;
          }
          else{
              s[3+j][ne+indexv[i+5+NSPECIES]] = 0.0;
          }
      }
     s[3+j][jsf] = y[j+5+NSPECIES][M]-parameters.Ncol_sp[j];

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
  I_EQNVARS avgvars;
  double ravg;
  double eqn;
  VARLIST eqn_km1_pl, eqn_km1_mi, eqn_k_pl, eqn_k_mi;
  VARLIST varmod;
  double ymod[NE+1][2];
  int i,j,last=4+2*NSPECIES;

  avgvars.q     = 0.5*(x[k]+x[k-1]);
  avgvars.v     = 0.5*(y[1][k]+y[1][k-1]);
  avgvars.z     = 0.5*(y[2][k]+y[2][k-1]);
  avgvars.rho   = 0.5*(y[3][k]+y[3][k-1]);
  avgvars.T     = 0.5*(y[4][k]+y[4][k-1]);
  for (j=0; j<NSPECIES; j++){
      avgvars.Ys[j]    = 0.5*(y[j+5][k]+y[j+5][k-1]);
      avgvars.Ncol[j]  = 0.5*(y[j+5+NSPECIES][k]+y[j+5+NSPECIES][k-1]);
  }
  get_rad(&ravg, avgvars);
  avgvars.r    = ravg;

  /* Velocity Equation */
  get_eqn_evals(&varmod, &eqn, &eqn_km1_pl, &eqn_km1_mi, &eqn_k_pl, &eqn_k_mi,
                avgvars, k, x, y, ymod, VEQN, 1); //FIX trying 1 to avoid potential breakage when indexing in, e.g., ion_eqn
    
//   printf("k=%d: eqn %.2e, v %.2e \n",k,eqn,eqn_km1_pl.v); 

  s[1][indexv[1]] = (eqn_km1_pl.v-eqn_km1_mi.v)/varmod.v;
  s[1][indexv[2]] = (eqn_km1_pl.z-eqn_km1_mi.z)/varmod.z;
  s[1][indexv[3]] = (eqn_km1_pl.rho-eqn_km1_mi.rho)/varmod.rho;
  s[1][indexv[4]] = (eqn_km1_pl.T-eqn_km1_mi.T)/varmod.T;
  for (j=0; j<NSPECIES; j++){
//         get_eqn_evals(&varmod, &eqn, &eqn_km1_pl, &eqn_km1_mi, &eqn_k_pl, &eqn_k_mi,
//                 avgvars, k, x, y, ymod, VEQN, j+1);
      s[1][indexv[j+5]] = (eqn_km1_pl.Ys[j]-eqn_km1_mi.Ys[j])/varmod.Ys[j];
      s[1][indexv[j+5+NSPECIES]] = (eqn_km1_pl.Ncol[j]-eqn_km1_mi.Ncol[j])/varmod.Ncol[j];
  }
  s[1][ne+indexv[1]] = (eqn_k_pl.v-eqn_k_mi.v)/varmod.v;
  s[1][ne+indexv[2]] = (eqn_k_pl.z-eqn_k_mi.z)/varmod.z;
  s[1][ne+indexv[3]] = (eqn_k_pl.rho-eqn_k_mi.rho)/varmod.rho;
  s[1][ne+indexv[4]] = (eqn_k_pl.T-eqn_k_mi.T)/varmod.T;
  for (j=0; j<NSPECIES; j++){
      s[1][ne+indexv[j+5]] = (eqn_k_pl.Ys[j]-eqn_k_mi.Ys[j])/varmod.Ys[j];
      s[1][ne+indexv[j+5+NSPECIES]] = (eqn_k_pl.Ncol[j]-eqn_k_mi.Ncol[j])/varmod.Ncol[j];
  }
  s[1][jsf] = eqn;

  /* Density Equation */
  get_eqn_evals(&varmod, &eqn, &eqn_km1_pl, &eqn_km1_mi, &eqn_k_pl, &eqn_k_mi,
                avgvars, k, x, y, ymod, RHOEQN, 1);
    
//    printf("k=%d: eqn %.2e, v %.2e \n",k,eqn,eqn_km1_pl.rho); 

  s[2][indexv[1]] = (eqn_km1_pl.v-eqn_km1_mi.v)/varmod.v;
  s[2][indexv[2]] = (eqn_km1_pl.z-eqn_km1_mi.z)/varmod.z;
  s[2][indexv[3]] = (eqn_km1_pl.rho-eqn_km1_mi.rho)/varmod.rho;
  s[2][indexv[4]] = (eqn_km1_pl.T-eqn_km1_mi.T)/varmod.T;
  for (j=0; j<NSPECIES; j++){
      s[2][indexv[j+5]] = (eqn_km1_pl.Ys[j]-eqn_km1_mi.Ys[j])/varmod.Ys[j];
      s[2][indexv[j+5+NSPECIES]] = (eqn_km1_pl.Ncol[j]-eqn_km1_mi.Ncol[j])/varmod.Ncol[j];
  }
  s[2][ne+indexv[1]] = (eqn_k_pl.v-eqn_k_mi.v)/varmod.v;
//     if(k<5){printf("rho_eqn(v[%d]) = %.3e\n",k,(eqn_k_pl.v-eqn_k_mi.v)/varmod.v);}
//     if(k==1){printf("indexv[2]=%d\n",indexv[2]);}
  s[2][ne+indexv[2]] = (eqn_k_pl.z-eqn_k_mi.z)/varmod.z;
  s[2][ne+indexv[3]] = (eqn_k_pl.rho-eqn_k_mi.rho)/varmod.rho;
  s[2][ne+indexv[4]] = (eqn_k_pl.T-eqn_k_mi.T)/varmod.T;
  for (j=0; j<NSPECIES; j++){
      s[2][ne+indexv[j+5]] = (eqn_k_pl.Ys[j]-eqn_k_mi.Ys[j])/varmod.Ys[j];
      s[2][ne+indexv[j+5+NSPECIES]] = (eqn_k_pl.Ncol[j]-eqn_k_mi.Ncol[j])/varmod.Ncol[j];
  }
  s[2][jsf] = eqn;

  /* z parameter Equation */
  s[3][indexv[1]] = 0.0;
  s[3][indexv[2]] = -1.0;
  s[3][indexv[3]] = 0.0;
  s[3][indexv[4]] = 0.0;
  for (j=0; j<NSPECIES; j++){
      s[3][indexv[j+5]] = 0.0;
      s[3][indexv[j+5+NSPECIES]] = 0.0;
  }
  s[3][ne+indexv[1]] = 0.0;
  s[3][ne+indexv[2]] = 1.0;
  s[3][ne+indexv[3]] = 0.0;
  s[3][ne+indexv[4]] = 0.0;
  for (j=0; j<NSPECIES; j++){
      s[3][ne+indexv[j+5]] = 0.0;
      s[3][ne+indexv[j+5+NSPECIES]] = 0.0;
  }
  s[3][jsf] = y[2][k]-y[2][k-1];

  /* Column Density Equation */
   for (j=0; j<NSPECIES; j++){
     get_eqn_evals(&varmod, &eqn, &eqn_km1_pl, &eqn_km1_mi, &eqn_k_pl, &eqn_k_mi,
                avgvars, k, x, y, ymod, NCOLEQN, j+1); //j+1 is the species #. Species #s begin with 1. 
       //This acts to call the Ncol1 version of the Column Density Equation found in soe.c
     
      s[4+j][indexv[1]] = (eqn_km1_pl.v-eqn_km1_mi.v)/varmod.v;
      s[4+j][indexv[2]] = (eqn_km1_pl.z-eqn_km1_mi.z)/varmod.z;
      s[4+j][indexv[3]] = (eqn_km1_pl.rho-eqn_km1_mi.rho)/varmod.rho;
      s[4+j][indexv[4]] = (eqn_km1_pl.T-eqn_km1_mi.T)/varmod.T;
      for (i=0; i<NSPECIES; i++){
          s[4+j][indexv[i+5]] = (eqn_km1_pl.Ys[i]-eqn_km1_mi.Ys[i])/varmod.Ys[i];
          s[4+j][indexv[i+5+NSPECIES]] = (eqn_km1_pl.Ncol[i]-eqn_km1_mi.Ncol[i])/varmod.Ncol[i];
      }
      s[4+j][ne+indexv[1]] = (eqn_k_pl.v-eqn_k_mi.v)/varmod.v;
      s[4+j][ne+indexv[2]] = (eqn_k_pl.z-eqn_k_mi.z)/varmod.z;
      s[4+j][ne+indexv[3]] = (eqn_k_pl.rho-eqn_k_mi.rho)/varmod.rho;
      s[4+j][ne+indexv[4]] = (eqn_k_pl.T-eqn_k_mi.T)/varmod.T;
      for (i=0; i<NSPECIES; i++){
          s[4+j][ne+indexv[i+5]] = (eqn_k_pl.Ys[i]-eqn_k_mi.Ys[i])/varmod.Ys[i];
          s[4+j][ne+indexv[i+5+NSPECIES]] = (eqn_k_pl.Ncol[i]-eqn_k_mi.Ncol[i])/varmod.Ncol[i];
      }
      s[4+j][jsf] = eqn;
//         if(k>1500){ printf("k=%d, Ncol_eqn%d residual = %.17e\n\n", k, j+1, eqn);} //also check in 
//            if(k>1500){ printf("k=%d, km1 Ncol_eqn%d v=[%.17e - %.17e] / %.17e\n\n", k, j+1, eqn_km1_pl.v,eqn_km1_mi.v,varmod.v);}
//           if(k>1500){ printf("k=%d, k Ncol_eqn%d v=[%.17e - %.17e] / %.17e\n\n", k, j+1, eqn_k_pl.v,eqn_k_mi.v,varmod.v);}
//            if(k>1500){ printf("k=%d, Ncol_eqn%d residual = %.17e\n\n", k, j+1, eqn);}
//           if(k>1500){ printf("k=%d, Ncol_eqn%d residual = %.17e\n\n", k, j+1, eqn);}
//           if(k>1500){ printf("k=%d, Ncol_eqn%d residual = %.17e\n\n", k, j+1, eqn);}
   }

  /* Ionization Equation */
   for (j=0; j<NSPECIES; j++){
      get_eqn_evals(&varmod, &eqn, &eqn_km1_pl, &eqn_km1_mi, &eqn_k_pl, &eqn_k_mi,
                avgvars, k, x, y, ymod, IONEQN, j+1);
      s[4+j+NSPECIES][indexv[1]] = (eqn_km1_pl.v-eqn_km1_mi.v)/varmod.v;
      s[4+j+NSPECIES][indexv[2]] = (eqn_km1_pl.z-eqn_km1_mi.z)/varmod.z;
      s[4+j+NSPECIES][indexv[3]] = (eqn_km1_pl.rho-eqn_km1_mi.rho)/varmod.rho;
      s[4+j+NSPECIES][indexv[4]] = (eqn_km1_pl.T-eqn_km1_mi.T)/varmod.T;
      for (i=0; i<NSPECIES; i++){
          s[4+j+NSPECIES][indexv[i+5]] = (eqn_km1_pl.Ys[i]-eqn_km1_mi.Ys[i])/varmod.Ys[i];
//           if (k==2){
//           printf("j=%d\n",j);
//           printf("pl.Ys[%d] = %e\n", i, eqn_k_pl.Ys[i]); //pasted in the wrong place b/c interested in the second half
//           printf("mi.Ys[%d] = %e\n", i, eqn_k_mi.Ys[i]);
//           printf("varmod.Ys[%d] = %e\n", i, varmod.Ys[i]);
//           printf("s[%d][%d] = %f\n\n",4+j+NSPECIES, ne+indexv[i+5],(eqn_k_pl.Ys[i]-eqn_k_mi.Ys[i])/varmod.Ys[i]); 
//           }
          s[4+j+NSPECIES][indexv[i+5+NSPECIES]] = (eqn_km1_pl.Ncol[i]-eqn_km1_mi.Ncol[i])/varmod.Ncol[i];
      }
      s[4+j+NSPECIES][ne+indexv[1]] = (eqn_k_pl.v-eqn_k_mi.v)/varmod.v;
      s[4+j+NSPECIES][ne+indexv[2]] = (eqn_k_pl.z-eqn_k_mi.z)/varmod.z;
      s[4+j+NSPECIES][ne+indexv[3]] = (eqn_k_pl.rho-eqn_k_mi.rho)/varmod.rho;
      s[4+j+NSPECIES][ne+indexv[4]] = (eqn_k_pl.T-eqn_k_mi.T)/varmod.T;
      for (i=0; i<NSPECIES; i++){
          s[4+j+NSPECIES][ne+indexv[i+5]] = (eqn_k_pl.Ys[i]-eqn_k_mi.Ys[i])/varmod.Ys[i];
//           if(k==2){printf("Ys1: %.2e",(eqn_k_pl.Ys[1]-eqn_k_mi.Ys[1])/varmod.Ys[i]);}
          s[4+j+NSPECIES][ne+indexv[i+5+NSPECIES]] = (eqn_k_pl.Ncol[i]-eqn_k_mi.Ncol[i])/varmod.Ncol[i];
      }
      s[4+j+NSPECIES][jsf] = eqn;
//        if(k>1500){ printf("k=%d, Ys_eqn residual = %.17e\n\n", k, eqn);}
//               if(k>1500){ printf("k=%d, Ys_eqn%d Ncol2=[%.17e - %.17e] / %.17e\n\n", k, j+1, eqn_km1_pl.Ncol[1],eqn_km1_mi.Ncol[1],varmod.Ncol[1]);}
   }
    
      /* Ionization Equation 2 */
//   get_eqn_evals(&varmod, &eqn, &eqn_km1_pl, &eqn_km1_mi, &eqn_k_pl, &eqn_k_mi,
//                 avgvars, k, x, y, ymod, IONEQN2);
// 
//   s[6][indexv[3]] = (eqn_km1_pl.rho-eqn_km1_mi.rho)/varmod.rho;
//   s[6][indexv[1]] = (eqn_km1_pl.v-eqn_km1_mi.v)/varmod.v;
//   s[6][indexv[2]] = (eqn_km1_pl.z-eqn_km1_mi.z)/varmod.z;
//   s[6][indexv[5]] = (eqn_km1_pl.Ys-eqn_km1_mi.Ys)/varmod.Ys;
//   s[6][indexv[7]] = (eqn_km1_pl.Ncol-eqn_km1_mi.Ncol)/varmod.Ncol;
//   s[6][indexv[4]] = (eqn_km1_pl.T-eqn_km1_mi.T)/varmod.T;
//   s[6][indexv[6]] = (eqn_km1_pl.Ys2-eqn_km1_mi.Ys2)/varmod.Ys2;
//   s[6][indexv[8]] = (eqn_km1_pl.Ncol2-eqn_km1_mi.Ncol2)/varmod.Ncol2;
//   s[6][ne+indexv[3]] = (eqn_k_pl.rho-eqn_k_mi.rho)/varmod.rho;
//   s[6][ne+indexv[1]] = (eqn_k_pl.v-eqn_k_mi.v)/varmod.v;
//   s[6][ne+indexv[2]] = (eqn_k_pl.z-eqn_k_mi.z)/varmod.z;
//   s[6][ne+indexv[5]] = (eqn_k_pl.Ys-eqn_k_mi.Ys)/varmod.Ys;
//   s[6][ne+indexv[7]] = (eqn_k_pl.Ncol-eqn_k_mi.Ncol)/varmod.Ncol;
//   s[6][ne+indexv[4]] = (eqn_k_pl.T-eqn_k_mi.T)/varmod.T;
//   s[6][ne+indexv[6]] = (eqn_k_pl.Ys2-eqn_k_mi.Ys2)/varmod.Ys2;
// //     printf("Ys2 %.2e \n",(eqn_k_pl.Ys2-eqn_k_mi.Ys2)/varmod.Ys2);
//   s[6][ne+indexv[8]] = (eqn_k_pl.Ncol2-eqn_k_mi.Ncol2)/varmod.Ncol2;
// //     printf("IONEQN2 Ncol2 %.2e     k=%d\n",(eqn_k_pl.Ncol2-eqn_k_mi.Ncol2)/varmod.Ncol2,k);
//   s[6][jsf] = eqn;

  /* Temperature Equation */
  get_eqn_evals(&varmod, &eqn, &eqn_km1_pl, &eqn_km1_mi, &eqn_k_pl, &eqn_k_mi,
                avgvars, k, x, y, ymod, TEQN, 1);

  s[last][indexv[1]] = (eqn_km1_pl.v-eqn_km1_mi.v)/varmod.v;
  s[last][indexv[2]] = (eqn_km1_pl.z-eqn_km1_mi.z)/varmod.z;
  s[last][indexv[3]] = (eqn_km1_pl.rho-eqn_km1_mi.rho)/varmod.rho;
  s[last][indexv[4]] = (eqn_km1_pl.T-eqn_km1_mi.T)/varmod.T;
  for (j=0; j<NSPECIES; j++){
      s[last][indexv[j+5]] = (eqn_km1_pl.Ys[j]-eqn_km1_mi.Ys[j])/varmod.Ys[j];
      s[last][indexv[j+5+NSPECIES]] = (eqn_km1_pl.Ncol[j]-eqn_km1_mi.Ncol[j])/varmod.Ncol[j];
  }
  s[last][ne+indexv[1]] = (eqn_k_pl.v-eqn_k_mi.v)/varmod.v;
  s[last][ne+indexv[2]] = (eqn_k_pl.z-eqn_k_mi.z)/varmod.z;
  s[last][ne+indexv[3]] = (eqn_k_pl.rho-eqn_k_mi.rho)/varmod.rho;
  s[last][ne+indexv[4]] = (eqn_k_pl.T-eqn_k_mi.T)/varmod.T;
  for (j=0; j<NSPECIES; j++){
      s[last][ne+indexv[j+5]] = (eqn_k_pl.Ys[j]-eqn_k_mi.Ys[j])/varmod.Ys[j];
      s[last][ne+indexv[j+5+NSPECIES]] = (eqn_k_pl.Ncol[j]-eqn_k_mi.Ncol[j])/varmod.Ncol[j];
  }
  s[last][jsf] = eqn;
   

  return;
}

/*============================================================================*
 *! \fn static void get_bc_evals(VARLIST *varmod_p, double *eqn_eval_p,       *
 *                               VARLIST *eqn_pl_p, VARLIST *eqn_mi_p, int k, *
 *                               double *x, double **y, double ymod[NE+1][2], *
 *                               int eqnnum)                                  *
 *  \brief Calculates the numerical derivative of the given equation for each *
 *         variable at the boundaries                                         *
 *----------------------------------------------------------------------------*/

static void get_bc_evals(VARLIST *varmod_p, double *eqn_eval_p,
                         VARLIST *eqn_pl_p, VARLIST *eqn_mi_p, int k, double *x,
                         double **y, double ymod[NE+1][2], int eqnnum) {
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
      varmod_p->Ys[j] = y[j+5][k]/DERIVDIV;
      zero_ymod(ymod);
      ymod[j+5][1]   = varmod_p->Ys[j]/2.0;
      eqn_pl_p->Ys[j] = eval_eqn( k, x, y, ymod, eqnnum,1); //FIX, this should maybe not be j+1, because all that species # controls is 
      //whether or not the EQNNUM plugged into is, e.g., IONEQN1 or IONEQN2 (a.k.a., as a function of Ys1 or Ys2). Since this function 
      //has no species dependence, we don't need to change the species number.
      //Technically, it shouldn't make a difference since IONEQN and NCOLEQN are never used with this function and it is meaningless to
      //have species set for the other EQNNUMs, since they have no species dependence.
      zero_ymod(ymod);
      ymod[j+5][1]   = -varmod_p->Ys[j]/2.0; //what is really making the difference for the species 1 and 2 here is calling Ys[j]
      eqn_mi_p->Ys[j] = eval_eqn( k, x, y, ymod, eqnnum,1);
//         printf("eqn_mi, j=%d: %.2e\n", j, eval_eqn( k, x, y, ymod, eqnnum,j+1));

  }

  /* ddNcol */
  for (j=0; j<NSPECIES; j++){
      varmod_p->Ncol[j] = y[j+5+NSPECIES][k]/DERIVDIV;
      zero_ymod(ymod);
      ymod[j+5+NSPECIES][1]   = varmod_p->Ncol[j]/2.0;
      eqn_pl_p->Ncol[j] = eval_eqn( k, x, y, ymod, eqnnum, 1);
      zero_ymod(ymod);
      ymod[j+5+NSPECIES][1]   = -varmod_p->Ncol[j]/2.0;
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
 *! \fn static void get_eqn_evals(VARLIST *varmod_p, double *eqn_eval_p,      *
 *                                VARLIST *eqn_km1_pl_p,                      *
 *                                VARLIST *eqn_km1_mi_p, VARLIST *eqn_k_pl_p, *
 *                                VARLIST *eqn_k_mi_p, I_EQNVARS avgvars,     *
 *                                int k, double *x, double **y,               *
 *                                double ymod[NE+1][2], int eqnnum)           *
 *  \brief Calculates the numerical derivative of the given equation for each *
 *         variable at (k-1)-th and k-th interior point                       *
 *----------------------------------------------------------------------------*/

static void get_eqn_evals(VARLIST *varmod_p, double *eqn_eval_p,
                          VARLIST *eqn_km1_pl_p, VARLIST *eqn_km1_mi_p,
                          VARLIST *eqn_k_pl_p, VARLIST *eqn_k_mi_p,
                          I_EQNVARS avgvars, int k, double *x, double **y,
                          double ymod[NE+1][2], int eqnnum, int species) {
  int j;
    
  zero_ymod(ymod);
  *eqn_eval_p = eval_eqn( k, x, y, ymod, eqnnum, species);

  /* ddrho */
  varmod_p->rho = avgvars.rho/DERIVDIV;
  /* k-1 */
  zero_ymod(ymod);
  ymod[3][0] = varmod_p->rho/2.0;
  eqn_km1_pl_p->rho = eval_eqn( k, x, y, ymod, eqnnum, species);
  zero_ymod(ymod);
  ymod[3][0] = -varmod_p->rho/2.0;
  eqn_km1_mi_p->rho = eval_eqn( k, x, y, ymod, eqnnum, species);
  /* k */
  zero_ymod(ymod);
  ymod[3][1] = varmod_p->rho/2.0;
  eqn_k_pl_p->rho = eval_eqn( k, x, y, ymod, eqnnum, species);
  zero_ymod(ymod);
  ymod[3][1] = -varmod_p->rho/2.0;
  eqn_k_mi_p->rho = eval_eqn( k, x, y, ymod, eqnnum, species);

  /* ddv */
  varmod_p->v = avgvars.v/DERIVDIV;
  /* k-1 */
  zero_ymod(ymod);
  ymod[1][0] = varmod_p->v/2.0;
  eqn_km1_pl_p->v = eval_eqn( k, x, y, ymod, eqnnum, species);
  zero_ymod(ymod);
  ymod[1][0] = -varmod_p->v/2.0;
  eqn_km1_mi_p->v = eval_eqn( k, x, y, ymod, eqnnum, species);
  /* k */
  zero_ymod(ymod);
  ymod[1][1] = varmod_p->v/2.0;
  eqn_k_pl_p->v = eval_eqn( k, x, y, ymod, eqnnum, species);
  zero_ymod(ymod);
  ymod[1][1] = -varmod_p->v/2.0;
  eqn_k_mi_p->v = eval_eqn( k, x, y, ymod, eqnnum, species);

  /* ddT */
  varmod_p->T = avgvars.T/DERIVDIV;
  /* k-1 */
  zero_ymod(ymod);
  ymod[4][0] = varmod_p->T/2.0;
  eqn_km1_pl_p->T = eval_eqn( k, x, y, ymod, eqnnum, species);
  zero_ymod(ymod);
  ymod[4][0] = -varmod_p->T/2.0;
  eqn_km1_mi_p->T = eval_eqn( k, x, y, ymod, eqnnum, species);
  /* k */
  zero_ymod(ymod);
  ymod[4][1] = varmod_p->T/2.0;
  eqn_k_pl_p->T = eval_eqn( k, x, y, ymod, eqnnum, species);
  zero_ymod(ymod);
  ymod[4][1] = -varmod_p->T/2.0;
  eqn_k_mi_p->T = eval_eqn( k, x, y, ymod, eqnnum, species);

  /* ddYs */
  for (j=0; j<NSPECIES; j++){ //FIX this might need to not loop, but just be one value of Ys corresponding to the species
//     j=species-1;
      varmod_p->Ys[j] = avgvars.Ys[j]/DERIVDIV;
      /* k-1 */
      zero_ymod(ymod);
      ymod[j+5][0] = varmod_p->Ys[j]/2.0;
      eqn_km1_pl_p->Ys[j] = eval_eqn( k, x, y, ymod, eqnnum, species);
//         printf("eqn_pl, j=%d: %.2e\n", j, eval_eqn( k, x, y, ymod, eqnnum, species));
      zero_ymod(ymod);
      ymod[j+5][0] = -varmod_p->Ys[j]/2.0;
      eqn_km1_mi_p->Ys[j] = eval_eqn( k, x, y, ymod, eqnnum, species);
//         printf("eqn_mi, j=%d: %.2e\n", j, eval_eqn( k, x, y, ymod, eqnnum, species));
      /* k */
      zero_ymod(ymod);
      ymod[j+5][1] = varmod_p->Ys[j]/2.0;
      eqn_k_pl_p->Ys[j] = eval_eqn( k, x, y, ymod, eqnnum, species);
      zero_ymod(ymod);
      ymod[j+5][1] = -varmod_p->Ys[j]/2.0;
      eqn_k_mi_p->Ys[j] = eval_eqn( k, x, y, ymod, eqnnum, species);
  }
    

  /* ddNcol */
  for (j=0; j<NSPECIES; j++){
//     j=species-1;
      varmod_p->Ncol[j] = avgvars.Ncol[j]/DERIVDIV;
      /* k-1 */
      zero_ymod(ymod);
      ymod[j+5+NSPECIES][0] = varmod_p->Ncol[j]/2.0;
      eqn_km1_pl_p->Ncol[j] = eval_eqn( k, x, y, ymod, eqnnum, species);
      zero_ymod(ymod);
      ymod[j+5+NSPECIES][0] = -varmod_p->Ncol[j]/2.0;
      eqn_km1_mi_p->Ncol[j] = eval_eqn( k, x, y, ymod, eqnnum, species);
      /* k */
      zero_ymod(ymod);
      ymod[j+5+NSPECIES][1] = varmod_p->Ncol[j]/2.0;
      eqn_k_pl_p->Ncol[j] = eval_eqn( k, x, y, ymod, eqnnum, species);
      zero_ymod(ymod);
      ymod[j+5+NSPECIES][1] = -varmod_p->Ncol[j]/2.0;
      eqn_k_mi_p->Ncol[j] = eval_eqn( k, x, y, ymod, eqnnum, species);
  }

  /* ddz */
  varmod_p->z = avgvars.z/DERIVDIV;
  /* k-1 */
  zero_ymod(ymod);
  ymod[2][0] = varmod_p->z/2.0;
  eqn_km1_pl_p->z = eval_eqn( k, x, y, ymod, eqnnum, species);
  zero_ymod(ymod);
  ymod[2][0] = -varmod_p->z/2.0;
  eqn_km1_mi_p->z = eval_eqn( k, x, y, ymod, eqnnum, species);
  /* k */
  zero_ymod(ymod);
  ymod[2][1] = varmod_p->z/2.0;
  eqn_k_pl_p->z = eval_eqn( k, x, y, ymod, eqnnum, species);
  zero_ymod(ymod);
  ymod[2][1] = -varmod_p->z/2.0;
  eqn_k_mi_p->z = eval_eqn( k, x, y, ymod, eqnnum, species);

  return;
}

/*============================================================================*
 *! \fn static void zero_ymod(double ymod[NE+1][2])                           *
 *  \brief Zeroes the ymod matrix                                             *
 *----------------------------------------------------------------------------*/

static void zero_ymod(double ymod[NE+1][2]) {
  int i, j;

  for (i = 0; i <= NE; i++) {
    for (j = 0; j < 2; j++) {
      ymod[i][j] = 0;
    }
  }

  return;
}
