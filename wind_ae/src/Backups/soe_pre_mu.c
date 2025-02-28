/*============================================================================*
 *! \file soe.c (system of equations)                                         *
 *  \brief The explict declarations of the system of equations of the problem *
 *         both for the residuals of the relaxation method and the derivatives*
 *         of all the equations, as well as the system's numerical Jacobian   *
 *============================================================================*/

/* Standard C */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
/* My Headers */
#include "defs.h"
#include "wind.h"
#include "globals.h"
#include "prototypes.h"
#include "recombo_vars.h"
#include "line_cooling_coeffs.h"

/* static global variables */
static double dvdr_slope, dvdr_last, q_last;

/*----------------------------------------------------------------------------*
 *======================== PRIVATE FUNCTION PROTOTYPES =======================*
 *----------------------------------------------------------------------------*/
static void get_alpharec(double *alpharec_p, I_EQNVARS kvars);
static void get_mu(double *mu_p, I_EQNVARS kvars);
static void get_gamma(double *gamma_p);
static void set_vars(int k, double *x, double **y, I_EQNVARS *kvars,
                     I_EQNVARS *km1vars, I_EQNVARS *avgvars,
                     double ymod[NE+1][2]);
/*========================== EQUATIONS FUNCTIONS =============================*/
static double rho_eqn(int k, double *x, double **y, double ymod[NE+1][2], int species);
static double v_eqn(int k, double *x, double **y, double ymod[NE+1][2], int species);
static double T_eqn(int k, double *x, double **y, double ymod[NE+1][2], int species);
static double ion_eqn(int k, double *x, double **y, double ymod[NE+1][2], int species);
static double Ncol_eqn(int k, double *x, double **y, double ymod[NE+1][2], int species);
static double spv_eqn(int k, double *x, double **y, double ymod[NE+1][2], int species);
static double spcrit_eqn(int k, double *x, double **y, double ymod[NE+1][2], int species);
/*========================== DERIVATIVES FUNCTIONS ===========================*/
static double get_component(I_EQNVARS vars, int compnum);
static void get_component_derivs(VARLIST *derivs, I_EQNVARS vars, int compnum);
    
/*============================================================================*
 *=========================== PUBLIC FUNCTIONS ===============================*
 *============================================================================*
 * eval_eqn      - calls calculation of the residual                          *
 * get_dvdr      - calculates the derivative of velocity                      *
 * get_drhodr    - calculates the derivative of density                       *
 * get_dTdr      - calculates the derivative of temperature                   *
 * get_dYsdr     - calculates the derivative of ionization fraction           *
 * get_dNcoldr   - calculates the derivative of column density                *
 * get_spQ       - calculates the total specific heating                      *
 * get_rad       - get radius                                                 *
 *----------------------------------------------------------------------------*/

/*============================================================================*
 *! \fn double eval_eqn(int k, double *x, double **y, double ymod[NE+1][2]    *
 *                      int eqnnum)                                           *
 *  \brief Calculates the k-th residual of the relaxation method for given    *
 *         eqnnum                                                             *
 *----------------------------------------------------------------------------*/

double eval_eqn(int k, double *x, double **y, double ymod[NE+1][2],
                int eqnnum, int species) {
  double residual;

  if (eqnnum == VEQN) {
    residual = v_eqn(k, x, y, ymod, 1); /*Species Independent*/
//       if(k==200){printf("residual: %.2e    k=%d    species=%d \n", residual,k,1);}
  }
  else if (eqnnum == RHOEQN) {
    residual = rho_eqn(k, x, y, ymod, 1);
  }
  else if (eqnnum == TEQN) {
    residual = T_eqn(k, x, y, ymod, 1);
  }
  else if (eqnnum == SPVEQN) {
    residual = spv_eqn(k, x, y, ymod, 1);
  }
  else if (eqnnum == SPCRITEQN) {
    residual = spcrit_eqn(k, x, y, ymod, 1);
  }
  else if (eqnnum == IONEQN) {
    residual = ion_eqn(k, x, y, ymod, species); /*Species dependent*/
//       if(k==M){if(species=1){printf("residual: %.2e    k=%d    species=%d \n", residual,k,species);}}
  }
  else if (eqnnum == NCOLEQN) {
    residual = Ncol_eqn(k, x, y, ymod, species);
//         if(k==1501){if(species=1){printf("eval_eqn: residual: %.2e    k=%d    species=%d \n", residual,k,species);}}
//       printf("%.2e \n",residual);
  }
  else {
    fprintf(stderr, "ERROR: [relax_difeq: eval_eqn] Unrecognized eqnnum: %d\n",
            eqnnum);
    exit(701);
  }

  return residual;
}

/*============================================================================*
 *! \fn static void get_dvdr(double *dvdr, I_EQNVARS vars)                    *
 *  \brief Calculates the derivative of velocity with respects to r, both     *
 *         analytically and using L'Hopital's rule to approximate it near the *
 *         sonic point where it analytically tends towards 0/0.               *
 *----------------------------------------------------------------------------*/

void get_dvdr(double *dvdr, I_EQNVARS vars) {
  double mu, rad, gamma, Mach;
  double spQ, d1_phi, norm_semimajor, qinv;
  double weight, r, n, s;

  /* Initalize to zero */
  *dvdr = 0.;

  get_gamma(&gamma);
  get_mu(&mu, vars);
  Mach = vars.v/sqrt(gamma*vars.T/mu);
  if (Mach != Mach) {
    Mach = 0.0;
  }
  /* Calculate analytic dvdr weighting */
  n = parameters.erfn;
  if (Mach <= 1.) {
    r = parameters.rapidity;
    s = r/parameters.mach_limit-r+n;
  }
  else {
    /* The dvdr calculation seems more forgiving on the backend of Mach 1, so
       we use a larger rapiditiy and mach limit, but still not too large */
    r = 80.;
    s = r/(0.999)-r+n;
  }

  if (Mach < r/(r+s+n) || Mach > (r+2.*(s+n))/(r+s+n)) {
    /* In unweighted analytic dvdr region */
    weight = 1.;
  }
  else if (Mach <= 1.) {
    if (Mach >= parameters.mach_limit) { /* Mach >= r/(r+s-n) */
      /* In unweighted linearized dvdr region */
      weight = 0.;
    }
    else {
      /* In weighted dvdr region */
      weight = 0.5*erf(r*(1./Mach-1.)-s)+0.5*erf(n);
    }
  }
  else { /* Mach > 1. */
    if (Mach <= (r+2.*(s-n))/(r+s-n)) {
      /* In unweighted linearized dvdr region */
      weight = 0.0;
    }
    else {
      /* In weighted dvdr region */
      weight = 0.5*erf(r*(1./(2.-Mach)-1.)-s)+0.5*erf(n);
    }
  }

  if (weight != weight || weight < 0. || weight > 1.) {
    fprintf(stderr, "ERROR: erroneous weight:%e at Mach:%e\n", weight, Mach);
    exit(702);
  }

  /* If not within unweighted linearized dvdr region, calculate analytic dvdr */
  if (weight != 0.) {
    get_rad(&rad, vars);
    get_spQ(&spQ, vars,1e4);
    norm_semimajor = parameters.semimajor/parameters.Rp;

    /* Potential derivative terms */
    d1_phi = 1./SQR(rad);
    if (parameters.tidalforce == ON) {
      /* Inverse of tidal q parameter */
      qinv    = parameters.Mstar/parameters.Mp;
      d1_phi += -qinv/SQR(norm_semimajor-rad)
                +(1.+qinv)*(norm_semimajor-rad)/pow(norm_semimajor, 3.);
    }
    /* Scale potential terms */
    d1_phi *= (parameters.Rp/parameters.H0);

    /* Calculate the analytic velocity derivative */
    *dvdr = 2.*gamma*vars.T/(mu*rad)
            -(gamma-1.)*spQ/vars.v
            -d1_phi;
    /* Finish by common factors */
    *dvdr *= vars.v/(SQR(vars.v)-gamma*vars.T/mu);
    *dvdr *= weight;
  }
  /* If not in unweighted analytic dvdr region, add linearized dvdr */
  if (weight != 1.) {
    *dvdr += (1.-weight)*(dvdr_slope*(vars.q-q_last)+dvdr_last);
  }

  return;
}

void linearize_dvdr_crit(double *x, double **y) {
  int k,j;
  double Mach, rad, gamma, mu, last_Mach;
  double dYsdr[NSPECIES], dNdr[NSPECIES], dmudr, dgdr=0;
  double d1_phi, d2_phi, norm_semimajor, qinv;
  double spQ, spQ1, spQ2;
  VARLIST dQ;
  double a, b, c, sqrtarg, dvdr_crit;
  I_EQNVARS vars;
//           printf("linearize called\n");

  get_gamma(&gamma);
//       printf("gamma\n");
  /* The last purely analytic dvdr (used for linearization) */
  last_Mach = parameters.mach_limit
              /(1.+2.*(parameters.erfn*parameters.mach_limit
                       /parameters.rapidity));

  /* Find last analytic dvdr before interpolating */
  k = M-1;
  do {
      vars.q    = x[k];
      vars.v    = y[1][k];
      vars.z    = y[2][k];
      vars.rho  = y[3][k];
      vars.T    = y[4][k];
      for (j=0; j<NSPECIES; j++){
//           printf("linearize dvdr crit: Ys[%d]=%.2e \n",j,y[j+5][k]);
          vars.Ys[j]   = y[j+5][k];
          vars.Ncol[j] = y[j+5+NSPECIES][k];
      }

    get_mu(&mu, vars);
    Mach = vars.v/sqrt(gamma*vars.T/mu);
    if (Mach < last_Mach) {
      vars.q    = x[k+1];
      vars.v    = y[1][k+1];
      vars.z    = y[2][k+1];
      vars.rho  = y[3][k+1];
      vars.T    = y[4][k+1];
      for (j=0; j<NSPECIES; j++){
          vars.Ys[j]   = y[j+5][k+1];
          vars.Ncol[j] = y[j+5+NSPECIES][k+1];
      }

      break;
    }
    k--;
  } while(k > 0);

  /* Set dvdr_last and q_last */
  get_rad(&rad, vars);
//       printf("get_rad done\n");
  get_spQ(&spQ, vars,1e4);
//       printf("get_spQ done\n");
  norm_semimajor = parameters.semimajor/parameters.Rp;
//       printf("all get_ called\n");

  /* Potential derivative terms */
  d1_phi = 1./SQR(rad);
  if (parameters.tidalforce == ON) {
    /* Inverse of tidal q parameter */
    qinv    = parameters.Mstar/parameters.Mp;
    d1_phi += -qinv/SQR(norm_semimajor-rad)
              +(1.+qinv)*(norm_semimajor-rad)/pow(norm_semimajor, 3.);
  }
  /* Scale potential terms */
  d1_phi *= (parameters.Rp/parameters.H0);

  /* Calculate the analytic velocity derivative */
  dvdr_last = 2.*gamma*vars.T/(mu*rad)
              -(gamma-1.)*spQ/vars.v
              -d1_phi;
  /* Finish by common factors */
  dvdr_last *= vars.v/(SQR(vars.v)-gamma*vars.T/mu);
  q_last = vars.q;

  /* Calculate the velocity derivative at critical point via L'Hoptial rule */
  vars.q    = x[M];
  vars.v    = y[1][M];
  vars.z    = y[2][M];
  vars.rho  = y[3][M];
  vars.T    = y[4][M];
  for (j=0; j<NSPECIES; j++){
      vars.Ys[j]   = y[j+5][M];
      vars.Ncol[j] = y[j+5+NSPECIES][M];
  }
    
  get_rad(&rad, vars);
  get_mu(&mu, vars);
  get_spQ(&spQ, vars,1e4);

  get_component_derivs(&dQ, vars, QCOMP);
  get_dYsdr(dYsdr, vars,1e4);
  get_dNcoldr(dNdr, vars);
  for (j=0; j<NSPECIES; j++){
        dgdr+=-parameters.HX[j]*dYsdr[j];
  }
  dmudr = -pow(mu, 2)*dgdr;

  /* Potential derivative terms */
  norm_semimajor = parameters.semimajor/parameters.Rp;
  d2_phi = -2./pow(rad, 3.);
  if (parameters.tidalforce == ON) {
    /* Inverse of tidal q parameter */
    qinv    = parameters.Mstar/parameters.Mp;
    d2_phi += -2.*qinv/pow(norm_semimajor-rad, 3.)
              -(1.+qinv)/pow(norm_semimajor, 3.);
  }
  /* Scale potential terms */
  d2_phi *= (parameters.Rp/parameters.H0);

  /* Define dQ/dr = Q1 * dv/dr + Q2 */
  spQ1 = dQ.v
         +dQ.rho*(-vars.rho/vars.v)
         +dQ.T*(-(gamma-1.)*vars.T/vars.v);
  spQ1 /= vars.rho;
  
  spQ2 = dQ.rho*(-2.*vars.rho/rad)
     +dQ.T*(vars.T/mu*dmudr
            +(gamma-1.)*spQ*mu/vars.v
            -2.*(gamma-1.)*vars.T/rad);
 
  for (j=0; j<NSPECIES; j++){
      spQ2+= dQ.Ys[j]*fmin(dYsdr[j], 0.) /* min(dYsdr,0) hack powers through when Ys  */
                                        /* needs to increase at the critical point,  */
                                        /* i.e., rec_rate < ion_rate. Senstivity due */
                                        /* to |dQ.Ys| >> 1.                          */
             + dQ.Ncol[j]*dNdr[j]; 
      
        if (dYsdr[j] > 0.) {
            printf("WARNING: dYsdr_%s=%e > 0\n", parameters.species[j],dYsdr[j]);
          }
  }
  spQ2 /= vars.rho;

  /* Using L'Hospital's rule, calculate velocity derivative at the
     critical point. */
  /* End up with quadratic: a*(dv/dr)^2 + b*(dv/dr) + c = 0, which
     we can now solve for dv/dr */
  a = 2.*vars.v+gamma*(gamma-1.)*vars.T/(mu*vars.v);
  b = (gamma-1.)*(4.*gamma*vars.T/(mu*rad)-gamma*spQ/vars.v+spQ1);
  c = SQR(gamma-1.)*(vars.v*d2_phi/SQR(gamma-1.)
                     +spQ2/(gamma-1.)
                     +2.*gamma*(2.*gamma-1.)*vars.T*vars.v/
                     (mu*SQR((gamma-1.)*rad))
                     -2.*spQ/rad);
  sqrtarg = SQR(b)-4.*a*c;
  /* if the argument of the square root is negative, we have problems */
  if (sqrtarg < 0) {
    fprintf(stderr, "ERROR: negative square root argument\n");
    fprintf(stderr, "       arg: %e = %e + %e\n", sqrtarg, SQR(b), -4.*a*c);
    exit(703);
  }

  /* The positive quadratic root corresponds to a wind
     and the negative quadratic root corresponds to accrection */
  dvdr_crit  = (-b+sqrt(sqrtarg))/(2.*a);
  dvdr_slope = (dvdr_crit-dvdr_last)/(vars.q-q_last);

  return;
}

/*============================================================================*
 *! \fn static void get_drhodr(double *drhodr, I_EQNVARS vars, double dvdr)   *
 *  \brief Calculates the derivative of density with respects to r            *
 *----------------------------------------------------------------------------*/

void get_drhodr(double *drhodr, I_EQNVARS vars, double dvdr) {
  double rad;

  get_rad(&rad, vars);

  *drhodr = -vars.rho*(2./rad+dvdr/vars.v);
    
//    printf("drhodr = %.2e \n",drhodr);

  return;
}

/*============================================================================*
 *! \fn static void get_dTdr(double *dTdr, I_EQNVARS vars, double drhodr      *
 *                           double dYsdr)                                    *
 *  \brief Calculates the derivative of temperature with respects to r, both  *
 *         analytically and using L'Hopital's rule for the density derivative *
 *         term near the critial point where one term analytically tends      *
 *         towards 0/0                                                        *
 *----------------------------------------------------------------------------*/

void get_dTdr(double *dTdr, I_EQNVARS vars, double drhodr, double *dYsdr, double k) {
  double mu, gamma, dmudr, spQ, dgdr=0;
  int j;
    
// #if defined(ISOTHERMAL) //FIX
//   /* Not very well fleshed out, but hack to make atmosphere isothermal below
//      some kind of ionizing energy deposition criteria, e.g., tau > cutoff */
// /*     CONCERN */
//   if (1.e-20*NCOL0*vars.Ncol > 1.e20) {
//     *dTdr = 0.;
//     return;
//   }
// #endif /* ISOTHERMAL */
  /* Isothermal below Rmin, only should be used for inward integration */
//  if (vars.q <= 0.) {
//    *dTdr = 0.;
//    return;
//  }

  /* For mu(X) of the form 1/g(X), then set dgdr = d(g(X))/dr */
  get_mu(&mu, vars);
  
  for (j=0; j<NSPECIES; j++){
        dgdr += -parameters.HX[j]*dYsdr[j];
  }
  dmudr = -pow(mu, 2)*dgdr;

  get_gamma(&gamma);
  get_spQ(&spQ, vars,k);

  *dTdr = (gamma-1.0)*(spQ/vars.v*mu
                       +vars.T/vars.rho*drhodr)
          +vars.T/mu*dmudr; //also, none of these in right units maybe?

  return;
}

/*============================================================================*
 *! \fn static void get_dYsdr(double *dYsdr, I_EQNVARS vars)                  *
 *  \brief Calculates the derivative of ionization fraction with respects to r*
 *----------------------------------------------------------------------------*/

void get_dYsdr(double *dYsdr, I_EQNVARS vars, double k) {
  int j;
  double ion_rate, rec_rate, mu, ne=0;
  double alpharec[NSPECIES];
//     FILE *fptr;
    
  get_mu(&mu, vars);
  
  for (j=0; j<NSPECIES; j++){
        ne+=(RHO0/parameters.atomic_mass[j])*vars.rho*((1-vars.Ys[j])*parameters.HX[j]);
  }

  /* Calculate derived quantities */
  get_alpharec(alpharec, vars);
   
  
  for (j=0; j<NSPECIES; j++){
      ion_rate = vars.Ys[j]*parameters.Ftot*glq_ionization(vars.Ncol,vars.Ys,vars.rho,k)[j]; 
//             fptr = fopen("ion_rate.txt","a+");
//           fprintf(fptr,"%.0f %d %.2e\n",k,j,glq_ionization(vars.Ncol,vars.Ys,vars.rho,k)[j]);
//           fclose(fptr); 

//       if(k!=1e4){
//        fprintf(f,"%.0f,%d,%.2e\n",k,j,ion_rate);
//       }
//         printf("glq_ionization run\n");
//            printf("glq_ionization, species %d = %.2e \n",j,glq_ionization(Ncols)[0]);

      rec_rate = (alpharec[j])*ne*(1.-vars.Ys[j]);
      /* Sources minus sinks */
      dYsdr[j]  = (rec_rate-ion_rate)/vars.v * parameters.Rp/CS0;      /* Convert to dimensionless derivative */
  }
//                fclose(f);

  return;
}

/*============================================================================*
 *! \fn static void get_dNcoldr(double *dNdr, I_EQNVARS vars)                 *
 *  \brief Calculates the derivative of column density (optical depth) with   *
 *         respects to r                                                      *
 *----------------------------------------------------------------------------*/
void get_dNcoldr(double *dNdr, I_EQNVARS vars) {
  double nH;
  int j;
    
    for (j=0; j<NSPECIES; j++){
      /* Calculate derived quantities */
       nH = (RHO0/parameters.atomic_mass[j])*parameters.HX[j]*vars.rho;
//       if(j==0){printf("nH = %.5e\n",nH);}
//       printf("MH = %.5e\n",parameters.atomic_mass[j]);
      /* Calculate averaged derivative */
      dNdr[j] = -vars.Ys[j]*nH;// * parameters.Rp/NCOL0; /* Convert to dimensionless derivative */
      dNdr[j] *= parameters.Rp/NCOL0; //separated this out
//        printf("get_dNcoldr: dNdr: %.17e       species=%d \n", dNdr[j],j);
    }

  return;
}

/*============================================================================*
 *! \fn static void get_spQ(double *spQ, I_EQNVARS vars, double k)            *
 *  \brief Calculates the total specific heating,                             *
 *          k included for diagnostic print-outs                              *
 *----------------------------------------------------------------------------*/

void get_spQ(double *spQ, I_EQNVARS vars,double k) {
  double n0overrho, nIONoverrho=0,nIONjoverrho=0, ne=0, n_tot=0;
  double spQ_photoheat=0.0, spQ_lyacool=0.0, spQ_reccool=0.0, spQ_linecool=0.0;
  double spQ_boloheat=0.0, spQ_bolocool=0.0, smoothing_erf;
  double kappa_opt = 4e-3, kappa_IR = 1e-2, P;
  int j,sp,line;
//   double scale=0.05;
    
        FILE *fptr;

  /* The heating and cooling terms */
  *spQ = 0.0;
  /* photo-ionization heating */
  for (j=0; j<NSPECIES; j++){
      n0overrho = vars.Ys[j]*parameters.HX[j]/(parameters.atomic_mass[j]);   // number density of neutrals divided by density
      spQ_photoheat  += n0overrho*parameters.Ftot*glq_heating(vars.Ncol,vars.Ys,vars.rho,k)[j];
//           fptr = fopen("heating_rate.txt","a+");
//           fprintf(fptr,"%.0f %d %.2e\n",k,j,glq_heating(vars.Ncol,vars.Ys,vars.rho,k)[j]);
//           fclose(fptr); 
    } 
  spQ_photoheat *= (parameters.Rp/pow(CS0, 3));
  *spQ += spQ_photoheat;
  
    /* Bolometric Heating */
    /* Kappa optical and infrared are not valid at low pressures, so we send them 0 
    (over an order of mag in pressure) via an error function */
    for (j=0; j<NSPECIES; j++){
        n_tot+=(RHO0/parameters.atomic_mass[j])*vars.rho*parameters.HX[j];
    }
    P = n_tot*K*vars.T*T0; //Pressure in barye (g/cm/s2)
    smoothing_erf = erf(P/1); //transition pressure = 100 nanobars = 0.1 barye (g/cm2/s)
    spQ_boloheat = parameters.Lstar/(4*PI*pow(parameters.semimajor,2)) *(kappa_opt*smoothing_erf+0.25*kappa_IR*smoothing_erf); // divided by rho
  *spQ += RAMPFACTOR2*spQ_boloheat*(parameters.Rp/pow(CS0, 3)); 
    
    /* Bolometric cooling */
    spQ_bolocool = -2*SIG_SB*pow(vars.T*T0,4)*kappa_IR*smoothing_erf; // divided by rho
//     if (k==462){
//         printf("-------\n");
//         printf("R %.2e\n",vars.z+parameters.Rmin/parameters.Rp);
//         printf("boloheat  %.2e\n",spQ_bolocool);
//         printf("T.        %.2e\n",vars.T);
// //         printf("semimajor %.2e\n",parameters.semimajor);
// //         printf("kappa_opt %.2e\n",kappa_opt);
//         printf("kappa_IR  %.2e\n",kappa_IR);
//         printf("erf       %.2e\n",smoothing_erf);
        
//     }
  *spQ += RAMPFACTOR*spQ_bolocool*(parameters.Rp/pow(CS0, 3)); 

  
  /* Lyman-alpha cooling */
  if (parameters.lyacool == ON) {
    for (j=0; j<NSPECIES; j++){
            ne+=(RHO0/parameters.atomic_mass[j])*vars.rho*((1-vars.Ys[j])*parameters.HX[j]);
    }
//       fptr = fopen("ne.txt","a+");
//       fprintf(fptr,"%.0f %d %.2e\n",k,j,ne);
//       fclose(fptr); 
    n0overrho = vars.Ys[0]*parameters.HX[0]/(parameters.atomic_mass[0]); //This assumes H is the first species in list 
    spQ_lyacool  = LYACOOL_COEFF*exp(-LYACOOL_TEMP/vars.T)*(n0overrho)*ne; //LYACOOL_COEFF is negative
    spQ_lyacool *= (parameters.Rp/pow(CS0, 3));
//        if(k==700){printf("ly %.2e\n",spQ_lyacool);}
    /* Note that LYACOOL_COEFF < 0 */
    *spQ += spQ_lyacool;
      
//       if (k==50){printf("spQ before %.2e\n",spQ);}
//   }
    
  /* Line Cooling (OII, OIII, CII, CIII are relevant; James Owen, private correspondence) */
  /* OII, OIII, CII, CIII have 4,4,3,2 cooling lines respectively*/
//   if (parameters.lyacool == ON) {
    for (j=0; j<NSPECIES; j++){
        // CII Line Cooling
        if (strcmp(parameters.species[j],"CI") == 0){
            sp=j;
//             if (k==50){printf("-------------CII line cooling triggered\n");}
            nIONjoverrho = (1-vars.Ys[sp])*parameters.HX[sp]/(parameters.atomic_mass[sp]);
            for (line=0; line<3; line++){
                //CII are coefficients for the line cooling eqn. First index is which emission line, 2nd is which index
                spQ_linecool = nIONjoverrho*ne* CII[line][1]*exp(-CII[line][2]/(vars.T*T0)) / (ne*(1+CII[line][3]/ne)); 
//                 if (line==1){if(k==50){printf("%.2e\n",CII[line][2]/(vars.T*T0));}}
                spQ_linecool *= (parameters.Rp/pow(CS0, 3))/1.014;
//                 if (line==1){if(k==700){printf("%.2e\n",spQ_linecool);}}
                *spQ += spQ_linecool;
            }
        }
        spQ_linecool=0.0;
        // CIII Line Cooling
        if (strcmp(parameters.species[j],"CII") == 0){
            sp=j;
//             if (k==1){printf("CIII line cooling triggered\n");}
            nIONjoverrho = (1-vars.Ys[sp])*parameters.HX[sp]/(parameters.atomic_mass[sp]);
            for (line=0; line<2; line++){
                //CII are coefficients for the line cooling eqn. First index is which emission line, 2nd is which index
                spQ_linecool = nIONjoverrho*ne* CIII[line][1]*exp(-CIII[line][2]/(vars.T*T0)) / (ne*(1+CIII[line][3]/ne)); 
//                 if (line==1){if(k==50){printf("%.2e\n",spQ_linecool);}}
                spQ_linecool *= (parameters.Rp/pow(CS0, 3));
                *spQ += spQ_linecool;
            }
        }
        spQ_linecool=0.0;
        // OII Line Cooling
        if (strcmp(parameters.species[j],"OI") == 0){
            sp=j;
//             if (k==1){printf("OII line cooling triggered\n");}
            nIONjoverrho = (1-vars.Ys[sp])*parameters.HX[sp]/(parameters.atomic_mass[sp]);
            for (line=0; line<4; line++){
                //CII are coefficients for the line cooling eqn. First index is which emission line, 2nd is which index
                spQ_linecool = nIONjoverrho*ne* OII[line][1]*exp(-OII[line][2]/(vars.T*T0)) / (ne*(1+OII[line][3]/ne)); 
//                 if (line==1){if(k==50){printf("%.2e\n",spQ_linecool);}}
                spQ_linecool *= (parameters.Rp/pow(CS0, 3));
                *spQ += spQ_linecool;
            }
        }
        spQ_linecool=0.0;
        // OIII Line Cooling
        if (strcmp(parameters.species[j],"OII") == 0){
            sp=j;
//             if (k==1){printf("OIII line cooling triggered\n");}
            nIONjoverrho = (1-vars.Ys[sp])*parameters.HX[sp]/(parameters.atomic_mass[sp]);
            for (line=0; line<4; line++){
                //CII are coefficients for the line cooling eqn. First index is which emission line, 2nd is which index
                spQ_linecool = nIONjoverrho*ne* OIII[line][1]*exp(-OIII[line][2]/(vars.T*T0)) / (ne*(1+OIII[line][3]/ne)); 
                spQ_linecool *= (parameters.Rp/pow(CS0, 3));
                *spQ += spQ_linecool;
            }
        }
    }
  }
  /* Collisions with neutral H are relevant when e- density low => CI and OI line cooling matter */
    
  /* Recombination cooling */
    for (j=0; j<NSPECIES; j++){
        nIONoverrho  += (1.-vars.Ys[j])*parameters.HX[j]/(parameters.atomic_mass[j]); //TOTAL number density of ionized species
      }
    spQ_reccool  = -RECCOOL_COEFF*pow(vars.T, 0.11)*(nIONoverrho)*ne; //changed 5/4/2023 to be negative
    spQ_reccool *= (parameters.Rp/pow(CS0, 3));
    *spQ += spQ_reccool;

  return;
}

/*============================================================================*
 *! \fn static void get_rad(double *rad, I_EQNVARS vars)                      *
 *  \brief Calculates the radius given q and z                                *
 *----------------------------------------------------------------------------*/

void get_rad(double *rad, I_EQNVARS vars) {
  /* Note to self (not to or from John): is q vs. x[k] causing a problem? */
  *rad = parameters.Rmin+vars.q*vars.z;

  return;
}

/*============================================================================*
 *=========================== PRIVATE FUNCTIONS ==============================*
 *============================================================================*
 * get_alpharec         - gets the recombination coefficient                  *
 * get_mu               - get the mean molecular weight                       *
 * get_gamma            - get adiabatic index                                 *
 * set_vars             - sets variable values in around current neighborhood *
 * rho_eqn              - density residual equation                           *
 * v_eqn                - velocity residual equation                          *
 * T_eqn                - temperature residual equation                       *
 * ion_eqn              - ionization fraction residual equation               *
 * Ncol_eqn             - column density residual equation                    *
 * spv_eqn              - numerator residual equation                         *
 * spcrit_eqn           - denominator residual equation                       *
 * get_drhodr           - calculates the derivative of density                *
 * get_dTdr             - calculates the derivative of temperature            *
 * get_dYsdr            - calculates the derivative of ionization fraction    *
 * get_dNcoldr          - calculates the derivative of column density         *
 * get_component        - gets the value of requested component               *
 * get_component_derivs - gets derivative of component wrt relaxed variables  *
 *----------------------------------------------------------------------------*/

/*============================================================================*
 *! \fn static void get_alpharec(double *alpharec_p, I_EQNVARS kvars)         *
 *  \brief Calculates the temperature dependence of the recombination         *
 *----------------------------------------------------------------------------*/

static void get_alpharec(double *alpharec_p, I_EQNVARS kvars) {
    int j,i;
    int iz, in;
    double tt,T,r;
    
    //CLOUDY algorithm uses unscaled T, so mulitply by TEMPSCALE=1e4
    T = kvars.T*T0;                     

    for (j=0; j<NSPECIES; j++){
      iz = parameters.Z[j];
      in = parameters.N_e[j];
      //Adapted from the CLOUDY rrfit fortran agorithm
        if(in<3 || in==11 || (iz>5 && iz<9) || iz==10){ 
            for(i=0;i<126;i++){
                if(rnew[i][0]==iz && rnew[i][1]==in){
                     tt = sqrt(T/rnew[i][4]); //make kvar.T*1e4
                     r = rnew[i][2] / ( tt*pow(tt+1.0,1.0-rnew[i][3]) * pow(1.0+sqrt(T/rnew[i][5]),1.0+rnew[i][3]) );
                }
            }
        }
        else{
                tt=T*1.0e-04;
                if(iz==26 && in<13){
                    for(i=0;i<10;i++){
                          if(fe[i][0]==in){
                               r=fe[i][1]/pow(tt,(fe[i][2]+fe[i][3]*log10(tt)));
                          }
                    }
                }
                else{
                    for(i=0;i<349;i++){
                      if(rrec[i][0]==iz && rrec[i][1]==in){
                            r=rrec[i][2]/pow(tt,rrec[i][3]);
                        }
                    }
                }
            }
    
    alpharec_p[j] = r;
//     printf("alp %.2e \n",r);
    }
    
//     // Original H-only version
//     for (j=0; j<NSPECIES; j++){
//         alpharec_p[j] = pow(kvars.T, -0.7);
//     }

  return;
}



/*============================================================================*
 *! \fn static void get_mu(double *mu_p, I_EQNVARS kvars)                     *
 *  \brief Calculates the mean molecular weight, includes electrons           *
 *----------------------------------------------------------------------------*/

static void get_mu(double *mu_p, I_EQNVARS kvars) {
  /* Assumes only a hydrogen/helium mixture, where only hydrogen can ionize */
    int j;
    double denominator;
    
    if (strcmp(parameters.species[0],"HI") != 0){
        printf("WARNING: Mean molecular weight calculation assumes that the first species is HI.\n");
    }
        
    denominator = 0;
    for(j=0; j<NSPECIES; j++){
        denominator += (parameters.atomic_mass[0]/parameters.atomic_mass[j])*parameters.HX[j]*(2-kvars.Ys[j]);
    }
    *mu_p = 1./denominator; //FIX
    
//       *mu_p = 1./(MHOVERMHE+1*(2.-MHOVERMHE-kvars.Ys[0])); /*edited to work with double H*/ //KILL
//       *mu_p = 1./(MHOVERMHE+parameters.HX*(2.-MHOVERMHE-kvars.Ys)); /*original*/


  return;
}



/*============================================================================*
 *! \fn static void get_gamma(double *gamma_p)                                *
 *  \brief Sets the adiabatic index of the gas; in principle, this could vary *
 *----------------------------------------------------------------------------*/

static void get_gamma(double *gamma_p) {
  *gamma_p = GAMMA_ATOMIC;

  return;
}

/*============================================================================*
 *! \fn static void set_vars(int k, double *x, double **y, I_EQNVARS *kvars,  *
 *                           I_EQNVARS *km1vars, I_EQNVARS *avgvars,          *
 *                           double ymod[NE+1][2])                            *
 *  \brief Calculates the neighboring and average variables given ymod update *
 *  \note ymod[][0] modifies cell k-1 and ymod[][1] modifies cell k           *
 *----------------------------------------------------------------------------*/

static void set_vars(int k, double *x, double **y, I_EQNVARS *kvars,
                     I_EQNVARS *km1vars, I_EQNVARS *avgvars,
                     double ymod[NE+1][2]) {
  int j;
    
  kvars->q      = x[k];
  kvars->v      = y[1][k]+ymod[1][1];
  kvars->z      = y[2][k]+ymod[2][1];
  kvars->rho    = y[3][k]+ymod[3][1];
  kvars->T      = y[4][k]+ymod[4][1];
  for (j=0; j<NSPECIES; j++){
     kvars->Ys[j]     = y[j+5][k]+ymod[j+5][1];
     kvars->Ncol[j]   = y[j+5+NSPECIES][k]+ymod[j+5+NSPECIES][1];
  }

  km1vars->q    = x[k-1];
  km1vars->v    = y[1][k-1]+ymod[1][0];
  km1vars->z    = y[2][k-1]+ymod[2][0];
  km1vars->rho  = y[3][k-1]+ymod[3][0];
  km1vars->T    = y[4][k-1]+ymod[4][0];
  for (j=0; j<NSPECIES; j++){
      km1vars->Ys[j]   = y[j+5][k-1]+ymod[j+5][0];
      km1vars->Ncol[j] = y[j+5+NSPECIES][k-1]+ymod[j+5+NSPECIES][0];
  }

  avgvars->q    = 0.5*(kvars->q+km1vars->q);
  avgvars->rho  = 0.5*(kvars->rho+km1vars->rho);
  avgvars->v    = 0.5*(kvars->v+km1vars->v);
  avgvars->z    = 0.5*(kvars->z+km1vars->z);
  avgvars->T    = 0.5*(kvars->T+km1vars->T);
  for (j=0; j<NSPECIES; j++){
      avgvars->Ys[j]   = 0.5*(kvars->Ys[j]+km1vars->Ys[j]);
      avgvars->Ncol[j] = 0.5*(kvars->Ncol[j]+km1vars->Ncol[j]);
  }

  return;
}

/*============================================================================*
 *! \fn static double spv_eqn(int k, double *x, double **y,                   *
 *                            double ymod[NE+1][2])                           *
 *  \brief Calculates the residual of the denominator for the critical point  *
 *         criteria, i.e., that the outflow is Mach 1 at the critical point   *
 *----------------------------------------------------------------------------*/

static double spv_eqn(int k, double *x, double **y, double ymod[NE+1][2], int species) {
  double mu, gamma, residual; 
  I_EQNVARS kvars;
  int j;

  kvars.q    = x[k];
  kvars.v    = y[1][k]+ymod[1][1];
  kvars.z    = y[2][k]+ymod[2][1];
  kvars.rho  = y[3][k]+ymod[3][1];
  kvars.T    = y[4][k]+ymod[4][1];
  for (j=0; j<NSPECIES; j++){
      kvars.Ys[j]   = y[j+5][k]+ymod[j+5][1];
      kvars.Ncol[j] = y[j+5+NSPECIES][k]+ymod[j+5+NSPECIES][1];
  }
    
  get_mu(&mu, kvars);
//     printf("mu = %.2e \n",mu);
  get_gamma(&gamma);

  residual = kvars.v-sqrt(kvars.T*gamma/mu)*parameters.breezeparam;

  return residual;
}

/*============================================================================*
 *! \fn static double spcrit_eqn(int k, double *x, double **y,                *
 *                               double ymod[NE+1][2])                        *
 *  \brief Calculates the residual of the numerator for the critical point    *
 *         criteria, i.e., the Bernoulli criteria                             *
 *----------------------------------------------------------------------------*/

static double spcrit_eqn(int k, double *x, double **y, double ymod[NE+1][2], int species) {
  double T_term, grav_term, Q_term, spQ;
  double r, mu, gamma, residual;
  I_EQNVARS kvars;
  int j;

  kvars.q    = x[k];
  kvars.v    = y[1][k]+ymod[1][1];
  kvars.z    = y[2][k]+ymod[2][1];
  kvars.rho  = y[3][k]+ymod[3][1];
  kvars.T    = y[4][k]+ymod[4][1];
  for (j=0; j<NSPECIES; j++){
      kvars.Ys[j]   = y[j+5][k]+ymod[j+5][1];
      kvars.Ncol[j] = y[j+5+NSPECIES][k]+ymod[j+5+NSPECIES][1];
  }

  /* Conservative work terms */
  r = parameters.Rmin+kvars.q*kvars.z;
  grav_term = -1.0/r;
  if (parameters.tidalforce == ON) {
    double norm_semimajor = parameters.semimajor/parameters.Rp;
    grav_term += ((parameters.Mstar/parameters.Mp)*r*
                  (-(norm_semimajor-r)/pow(norm_semimajor, 3.)
                   +1./SQR(norm_semimajor-r)));
  }
  /* Scale gravitational term */
  grav_term *= parameters.Rp/parameters.H0;

  /* Enthalpy term */
  get_mu(&mu, kvars);
  get_gamma(&gamma);
  T_term = 2.0*gamma*kvars.T/mu;

  /* Heating terms */
  get_spQ(&spQ, kvars,k);
  Q_term = -(gamma-1.0)*spQ*r/kvars.v;

  /* Terms should sum to zero */
  residual = T_term+grav_term+Q_term;

  return residual;
}

/*============================================================================*
 *! \fn static double ion_eqn(int k, double *x, double **y,                   *
 *                          double ymod[NE+1][2])                             *
 *  \brief Calculates the residual of the ionization fraction equation        *
 *----------------------------------------------------------------------------*/

static double ion_eqn(int k, double *x, double **y, double ymod[NE+1][2], int species) {
  double delta_Ys_k, delta_r_k, dYsdr_avg[NSPECIES], residual;
  I_EQNVARS kvars, km1vars, avgvars;
  int spec_index = species-1;

  set_vars(k, x, y, &kvars, &km1vars, &avgvars, ymod);

  /* Calculate the finite difference */
  delta_r_k  = avgvars.z*(kvars.q-km1vars.q);
    
 /* Calculate averaged derivative */
//     printf("k = %d\n",k);
  get_dYsdr(dYsdr_avg, avgvars,k); 
//   print_eta(kvars.Ncol,kvars.Ys,kvars.rho);
//   printf('k=%d   eta=%.2f\n',k,eta);
    
  delta_Ys_k = kvars.Ys[spec_index]-km1vars.Ys[spec_index];
    
  delta_r_k  = avgvars.z*(kvars.q-km1vars.q);

  /* Calculate finite difference residual */
  residual = delta_Ys_k-dYsdr_avg[spec_index]*delta_r_k;

  return residual;
}

/*============================================================================*
 *! \fn static double Ncol_eqn(int k, double *x, double **y,                  *
 *                             double ymod[NE+1][2])                          *
 *  \brief Calculates the residual of the column density (optical depth)      *
 *         equation                                                           *
 *----------------------------------------------------------------------------*/

static double Ncol_eqn(int k, double *x, double **y, double ymod[NE+1][2], int species) {
  double delta_N_k, delta_r_k, dNdr_avg[NSPECIES], residual;
  I_EQNVARS kvars, km1vars, avgvars;
  int spec_index = species-1;

  set_vars(k, x, y, &kvars, &km1vars, &avgvars, ymod);

  delta_r_k = avgvars.z*(kvars.q-km1vars.q);
  get_dNcoldr(dNdr_avg, avgvars);
    
  /* Calculate the finite difference */
  delta_N_k = kvars.Ncol[spec_index]-km1vars.Ncol[spec_index];
  residual = delta_N_k-(dNdr_avg[spec_index])*delta_r_k;
//             if(k==1501){printf("Ncol_eqn: delta_N_k: %.17e    k=%d    species=%d \n", delta_N_k,k,species);}
//             if(k==1501){printf("Ncol_eqn: dNdr_avg: %.17e    k=%d    species=%d \n", dNdr_avg[spec_index],k,species);}
//             if(k==1501){printf("Ncol_eqn: delta_r_k: %.17e    k=%d    species=%d \n", delta_r_k,k,species);}
//     if(k==1501){printf("Ncol_eqn: residual: %.17e    k=%d    species=%d \n", residual,k,species);}
//        printf("k=%d",k);

  return residual;
}

/*============================================================================*
 *! \fn static double v_eqn(int k, double *x, double **y,                     *
 *                          double ymod[NE+1][2])                             *
 *  \brief Calculates the residual of the velocity equation                   *
 *----------------------------------------------------------------------------*/

static double v_eqn(int k, double *x, double **y, double ymod[NE+1][2], int species) {
  double delta_v_k, delta_r_k, dvdr_avg, residual;
  I_EQNVARS kvars, km1vars, avgvars, critvars;
  int j;

  set_vars(k, x, y, &kvars, &km1vars, &avgvars, ymod);
  critvars.q    = x[M];
  critvars.v    = y[1][M];
  critvars.z    = y[2][M];
  critvars.rho  = y[3][M];
  critvars.T    = y[4][M];
  for (j=0; j<NSPECIES; j++){
      critvars.Ys[j]   = y[j+5][M];
      critvars.Ncol[j] = y[j+5+NSPECIES][M];  
  }
    
  /* Calculate the finite difference */
  delta_v_k = kvars.v-km1vars.v;
  delta_r_k = avgvars.z*(kvars.q-km1vars.q);

  /* Calculate averaged derivative */
  get_dvdr(&dvdr_avg, avgvars);

  /* Calculate finite difference residue */
  residual = delta_v_k-dvdr_avg*delta_r_k;

  return residual;
}

/*============================================================================*
 *! \fn static double rho_eqn(int k, double *x, double **y,                   *
 *                            double ymod[NE+1][2])                           *
 *  \brief Calculates the residual of the density equation                    *
 *----------------------------------------------------------------------------*/

static double rho_eqn(int k, double *x, double **y, double ymod[NE+1][2], int species) {
  double delta_rho_k, delta_r_k, drhodr_avg, residual;
  double dvdr_avg;
  I_EQNVARS kvars, km1vars, avgvars, critvars;
  int j;

  set_vars(k, x, y, &kvars, &km1vars, &avgvars, ymod);
  critvars.q    = x[M];
  critvars.v    = y[1][M];
  critvars.z    = y[2][M];
  critvars.rho  = y[3][M];
  critvars.T    = y[4][M];
  for (j=0; j<NSPECIES; j++){
      critvars.Ys[j]   = y[j+5][M];
      critvars.Ncol[j] = y[j+5+NSPECIES][M];
  }
    
  /* Calculate the finite difference */
  delta_rho_k = kvars.rho-km1vars.rho;
  delta_r_k   = avgvars.z*(kvars.q-km1vars.q);

  /* Calculate averaged derivatives */
  get_dvdr(&dvdr_avg, avgvars);
  get_drhodr(&drhodr_avg, avgvars, dvdr_avg);

  /* Calculate finite difference residue */
  residual = delta_rho_k-drhodr_avg*delta_r_k;
    
//     if(k<5){printf("delta_rho_k[%d] = %.17e\n",k,delta_rho_k);}
    
//     if(k<5){printf("rho_eqn residual[%d] = %e\n",k,residual);}
   
//     if(k<5){printf("delta_r_k[%d] = %.17e \n",k, delta_r_k);}
//     if(k<5){printf("drhodr_avg[%d] = %.17e \n", k,drhodr_avg);}
//     if(k<5){printf("dvdr_avg[%d] = %e \n",k, dvdr_avg);}


  return residual;
}

/*============================================================================*
 *! \fn static double T_eqn(int k, double *x, double **y,                     *
 *                          double ymod[NE+1][2])                             *
 *  \brief Calculates the residual of the temperature equation                *
 *----------------------------------------------------------------------------*/

static double T_eqn(int k, double *x, double **y, double ymod[NE+1][2], int species) {
  double delta_T_k, delta_r_k, dTdr_avg, residual;
  double dYsdr_avg[NSPECIES], drhodr_avg, dvdr_avg;
  I_EQNVARS kvars, km1vars, avgvars, critvars;
  int j;

  set_vars(k, x, y, &kvars, &km1vars, &avgvars, ymod);
  critvars.q    = x[M];
  critvars.v    = y[1][M];
  critvars.z    = y[2][M];
  critvars.rho  = y[3][M];
  critvars.T    = y[4][M];
  for (j=0; j<NSPECIES; j++){
      critvars.Ys[j]   = y[j+5][M];
      critvars.Ncol[j] = y[j+5+NSPECIES][M];
  }

  /* Calculate the finite difference */
  delta_T_k = kvars.T-km1vars.T;
  delta_r_k = avgvars.z*(kvars.q-km1vars.q);

  /* Calculate averaged derivatives */
  get_dvdr(&dvdr_avg, avgvars);
  get_drhodr(&drhodr_avg, avgvars, dvdr_avg);
  get_dYsdr(dYsdr_avg, avgvars,k); 
  get_dTdr(&dTdr_avg, avgvars, drhodr_avg, dYsdr_avg,k);

  /* Calculate finite difference residue */
  residual = delta_T_k-dTdr_avg*delta_r_k;

  return residual;
}

/*============================================================================*
 *! \fn static double Energy_Conservation(int k, double *x, double **y,       *
 *                          double ymod[NE+1][2])                             *
 *  \brief Tracks the energy as a function of radius.                         *
 *----------------------------------------------------------------------------*/

// static double Energy_Conservation(int k, double *x, double **y, double ymod[NE+1][2], int species) {
//   double delta_rho_k, delta_r_k, dTdr_avg, residual,gamma, mu;
//   double dYsdr_avg,dYsdr_avg2, drhodr_avg, dvdr_avg, dgdr, dmudr;
//   I_EQNVARS kvars, km1vars, avgvars, critvars;

//   set_vars(k, x, y, &kvars, &km1vars, &avgvars, ymod);
//   critvars.q    = x[M];
//   critvars.rho  = y[1][M];
//   critvars.v    = y[2][M];
//   critvars.z    = y[3][M];
//   critvars.Ys   = y[4][M];
//   critvars.Ncol = y[5][M];
//   critvars.T    = y[6][M];
//   critvars.Ys2  = y[7][M];
//   critvars.Ncol2= y[8][M];


//    get_gamma(&gamma);
//    get_mu(&mu,kvars);
//   /* Calculate the finite difference */
// //   delta_rho_k = kvars.rho-km1vars.rho;
// //   delta_r_k = avgvars.z*(kvars.q-km1vars.q);

//   /* Calculate averaged derivatives */
//   get_dvdr(&dvdr_avg, avgvars);
//   get_drhodr(&drhodr_avg, avgvars, dvdr_avg);
//   get_dYsdr(&dYsdr_avg, avgvars, 1); 
//   get_dYsdr(&dYsdr_avg2, avgvars, 2); 
//   dgdr  = -(parameters.HX*dYsdr_avg+parameters.HX2*dYsdr_avg2); 
//   dmudr = -pow(mu, 2)*dgdr;
//   get_dTdr(&dTdr_avg, avgvars, drhodr_avg, dYsdr_avg, dYsdr_avg2);

//   /* Calculate finite difference residue */
// //   residual = delta_T_k-dTdr_avg*delta_r_k;

//     printf("E %.2e \n",kvars.rho*kvars.v);
//   return;
// }                              
                 
/*============================================================================*
 *! \fn static double get_component(I_EQNVARS vars, int compnum)              *
 *  \brief Returns the value of the specified component                       *
 *----------------------------------------------------------------------------*/

static double get_component(I_EQNVARS vars, int compnum) {
  double component;

  if (compnum == QCOMP) {
    get_spQ(&component, vars,1e4);
    component *= vars.rho;
  }
  else {
    fprintf(stderr, "ERROR: unknown component\n");
    exit(704);
  }

  return component;
}

/*============================================================================*
 *! \fn static void get_component_derivs(VARLIST *derivs, I_EQNVARS vars,     *
 *                                       int compnum)                         *
 *  \brief Calculates the numerical derivative of the given component with    *
 *         respects to all of the relax variables                             *
 *----------------------------------------------------------------------------*/

static void get_component_derivs(VARLIST *derivs, I_EQNVARS vars, int compnum) {
  double delta_var, comp_pl, comp_mi;
  I_EQNVARS varmod;
  int j;

  /* finite difference wrt rho */
  delta_var   = vars.rho/DERIVDIV;
  varmod      = vars;
  varmod.rho  = vars.rho+delta_var/2.;
  comp_pl     = get_component(varmod, compnum);
  varmod.rho  = vars.rho-delta_var/2.;
  comp_mi     = get_component(varmod, compnum);
  derivs->rho = (comp_pl-comp_mi)/delta_var;

  /* finite difference wrt v */
  delta_var = vars.v/DERIVDIV;
  varmod    = vars;
  varmod.v  = vars.v+delta_var/2.;
  comp_pl   = get_component(varmod, compnum);
  varmod.v  = vars.v-delta_var/2.;
  comp_mi   = get_component(varmod, compnum);
  derivs->v = (comp_pl-comp_mi)/delta_var;

  /* finite difference wrt T */
  delta_var = vars.T/DERIVDIV;
  varmod    = vars;
  varmod.T  = vars.T+delta_var/2.;
  comp_pl   = get_component(varmod, compnum);
  varmod.T  = vars.T-delta_var/2.;
  comp_mi   = get_component(varmod, compnum);
  derivs->T = (comp_pl-comp_mi)/delta_var;

  /* finite difference wrt Ys */
  for (j=0; j<NSPECIES; j++){
      delta_var  = vars.Ys[j]/DERIVDIV;
      varmod     = vars;
      varmod.Ys[j]  = vars.Ys[j]+delta_var/2.;
      comp_pl    = get_component(varmod, compnum);
      varmod.Ys[j]  = vars.Ys[j]-delta_var/2.;
      comp_mi    = get_component(varmod, compnum);
      derivs->Ys[j] = (comp_pl-comp_mi)/delta_var;
  }

  /* finite difference wrt Ncol */
  for (j=0; j<NSPECIES; j++){
      delta_var  = vars.Ncol[j]/DERIVDIV;
      varmod     = vars;
      varmod.Ncol[j]  = vars.Ncol[j]+delta_var/2.;
      comp_pl    = get_component(varmod, compnum);
      varmod.Ncol[j]  = vars.Ncol[j]-delta_var/2.;
      comp_mi    = get_component(varmod, compnum);
      derivs->Ncol[j] = (comp_pl-comp_mi)/delta_var;
  }
    
  /* finite difference wrt z */
  delta_var = vars.z/DERIVDIV;
  varmod    = vars;
  varmod.z  = vars.z+delta_var/2.;
  comp_pl   = get_component(varmod, compnum);
  varmod.z  = vars.z-delta_var/2.;
  comp_mi   = get_component(varmod, compnum);
  derivs->z = (comp_pl-comp_mi)/delta_var;

  return;
}
