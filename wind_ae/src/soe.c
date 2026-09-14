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
#include "ion_pots.h"
#include "recombo_vars.h"
#include "line_cooling_coeffs.h"

/* static global variables */
static double dvdr_slope, dvdr_last, q_last;
double erf_norm = 1; //need to normalize erf(velocity) s.t. it equals 1.
double rate = 1e2;

/*============================================================================*
 *! \fn static double get_erf_smooth(const I_EQNVARS *vars)                         *
 *  \brief Returns the normalised complementary-error-function smoothing       *
 *         weight that drives the transition between molecular (lower          *
 *         atmosphere) and atomic (wind) regimes.  Result is in [0,1]:        *
 *           ~1  in the molecular region  (low velocity / deep atmosphere)    *
 *           ~0  in the wind             (high velocity / upper atmosphere)   *
 *         The result is NOT multiplied by any flag; callers apply their own  *
 *         scaling (bolo_heat_cool, molec_layer, etc.) as appropriate.        *
 *----------------------------------------------------------------------------*/
static double get_erf_smooth(const I_EQNVARS *vars) {
  double s = 1.0 - erf((vars->v - parameters.erf_drop[0]) / parameters.erf_drop[1]);
  if (erf_norm > 0.0)
    s /= erf_norm;
  return s;
}

/*============================================================================*
 * Precomputed alpharec lookup tables — filled ONCE by init_soe() after       *
 * set_parameters().  Avoids O(126–349) linear scans per species per call.    *
 *                                                                             *
 * kind: AR_RNEW=0  rnew table (126 rows)                                     *
 *       AR_FE  =1  fe   table (10  rows, iron only)                          *
 *       AR_RREC=2  rrec table (349 rows)                                     *
 *       AR_ZERO=3  sentinel — lower-state lookup when N_e already == Z       *
 * idx : row index in the chosen table.                                        *
 *                                                                             *
 * T-keyed cache: get_alpharec depends only on T (species params are fixed).  *
 * In the column-by-column Jacobian only 2 of (8+4N) perturbation columns     *
 * change T; all others share the same T and return the cached result.        *
 * For N=8 this skips ~95% of all get_alpharec calls entirely.                *
 *============================================================================*/
typedef enum { AR_RNEW=0, AR_FE=1, AR_RREC=2, AR_ZERO=3 } ArKind;
typedef struct { ArKind kind; int idx; } ArLookup;

static ArLookup ar_hi[NSPECIES_MAX];      /* calc_lower_ion_state = 1 (higher state) */
static ArLookup ar_lo[NSPECIES_MAX];      /* calc_lower_ion_state = 0 (lower  state) */

static double ar_cache_T_hi = -1.0;  /* T at which ar_cache_hi was last filled  */
static double ar_cache_T_lo = -1.0;
static double ar_cache_hi[NSPECIES_MAX];
static double ar_cache_lo[NSPECIES_MAX];

/* Line-cooling species index lookup: set to -1 if species absent in run.     *
 * Replaces per-call strcmp loops in get_spQ.                                 */
static int lc_CI   = -1, lc_CII  = -1, lc_OI   = -1, lc_OII  = -1;
static int lc_FeI  = -1, lc_MgI  = -1, lc_CaI  = -1, lc_NeII = -1;

/*============================================================================*
 * Chain ionization helpers                                                   *
 *                                                                            *
 * When multiple ionization states of one element are included (e.g. CI, CII *
 * and CIII), species j is flagged as a "child" of parent p if               *
 *   Z[j]==Z[p]  and  N_e[j]==N_e[p]-1  (one more electron stripped).       *
 *                                                                            *
 * The three carbon populations at any point are then:                        *
 *   n_CI   = Ys_CI  * n_C                                                   *
 *   n_CII  = Ys_CII * (1-Ys_CI) * n_C                                      *
 *   n_CIII = (1-Ys_CII)*(1-Ys_CI) * n_C                                    *
 * so they sum to n_C exactly.                                                *
 *                                                                            *
 * get_child_species(j)       - returns child index of j, or -1              *
 * get_eff_ntot_over_rho(j,Y) - n_tot[j]/rho, scaled by (1-Ys[parent])      *
 *                              for child species                             *
 *============================================================================*/
static int get_child_species(int j) {
  int c;
  for (c = 0; c < NSPECIES; c++)
    if (parameters.chain_parent[c] == j) return c;
  return -1;
}

static double get_eff_ntot_over_rho(int j, const double *Ys) {
  /* Walk the full chain to the root, accumulating (1-Ys[ancestor]) at every
   * level.  For a depth-1 child (e.g. CII with parent CI) this gives
   *   (1-Ys[CI]) * HX[CI]/m[CI]
   * and for depth-2 (CIII) it correctly gives
   *   (1-Ys[CII])*(1-Ys[CI]) * HX[CI]/m[CI].
   * Independent species (chain_parent < 0) return HX[j]/m[j] unchanged. */
  if (parameters.chain_parent[j] < 0)
    return parameters.HX[j] / parameters.atomic_mass[j];
  double prod = 1.0;
  int cur = j, root = j;
  while (parameters.chain_parent[cur] >= 0) {
    int anc = parameters.chain_parent[cur];
    prod   *= (1.0 - Ys[anc]);
    root    = anc;
    cur     = anc;
  }
  return prod * parameters.HX[root] / parameters.atomic_mass[root];
}
//actual value is computed below for k=2, so this is guess for k=1

/*----------------------------------------------------------------------------*
 *======================== PRIVATE FUNCTION PROTOTYPES =======================*
 *----------------------------------------------------------------------------*/
static void get_alpharec(double *alpharec_p, const I_EQNVARS *kvars, int calc_lower_ion_state);
static void get_mu(double *mu_p, const I_EQNVARS *kvars);
static void get_gamma(double *gamma_p);
static void set_vars(int k, double *x, double **y, I_EQNVARS *kvars,
                     I_EQNVARS *km1vars, I_EQNVARS *avgvars,
                     double (*ymod)[2]);
/*========================== EQUATIONS FUNCTIONS =============================*/
static double rho_eqn(int k, double *x, double **y, double (*ymod)[2], int species);
static double v_eqn(int k, double *x, double **y, double (*ymod)[2], int species);
static double T_eqn(int k, double *x, double **y, double (*ymod)[2], int species);
static double F_eqn(int k, double *x, double **y, double (*ymod)[2], int species);
static double ion_eqn(int k, double *x, double **y, double (*ymod)[2], int species);
static double Ncol_eqn(int k, double *x, double **y, double (*ymod)[2], int species);
static double spv_eqn(int k, double *x, double **y, double (*ymod)[2], int species);
static double spcrit_eqn(int k, double *x, double **y, double (*ymod)[2], int species);
/*========================== DERIVATIVES FUNCTIONS ===========================*/
static double get_component(const I_EQNVARS *vars, int compnum);
static void get_component_derivs(VARLIST *derivs, const I_EQNVARS *vars, int compnum);
    
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
 *! \fn double eval_eqn(int k, double *x, double **y, double (*ymod)[2]    *
 *                      int eqnnum)                                           *
 *  \brief Calculates the k-th residual of the relaxation method for given    *
 *         eqnnum                                                             *
 *----------------------------------------------------------------------------*/

double eval_eqn(int k, double *x, double **y, double (*ymod)[2],
                int eqnnum, int species) {
  double residual;

  if (eqnnum == VEQN) {
    residual = v_eqn(k, x, y, ymod, 1); /*Species Independent*/
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
  }
  else if (eqnnum == NCOLEQN) {
    residual = Ncol_eqn(k, x, y, ymod, species);
  }
  else {
    fprintf(stderr, "ERROR: [relax_difeq: eval_eqn] Unrecognized eqnnum: %d\n",
            eqnnum);
    exit(701);
  }

  return residual;
}

/*============================================================================*
 *! \fn void eval_interior_eqns(int k, double *x, double **y,                 *
 *                              double (*ymod)[2], double res[])           *
 *  \brief Evaluates ALL (NE-1) interior equation residuals in a single pass. *
 *                                                                            *
 *  Shared-computation restructuring (Patch 2):                              *
 *    Previously v_eqn/rho_eqn/T_eqn/ion_eqn/Ncol_eqn each independently   *
 *    called set_vars, get_dvdr, and get_drhodr on the same avgvars.         *
 *    Now set_vars is called once, then each derivative is computed once      *
 *    and the residuals are assembled inline.  This eliminates:               *
 *      - 4 redundant set_vars calls (5 -> 1)                                *
 *      - 2 redundant get_dvdr calls (3 -> 1) — each was get_spQ+get_mu     *
 *      - 1 redundant get_drhodr call (2 -> 1)                               *
 *      - 1 redundant get_dYsdr call (2 -> 1) — each was glq_ionization     *
 *      - 4 redundant erf_norm updates and delta_r_k computations            *
 *    The GLQ cache-sharing from the column-by-column strategy is preserved: *
 *    get_dvdr triggers the GLQ cache via get_spQ; get_dYsdr/get_dTdr reuse  *
 *    the same cached result at no extra cost.                               *
 *                                                                            *
 *  res[] layout (z-equation is analytic; excluded here):                    *
 *    res[0]            = v_eqn   residual        (-> s row 1)               *
 *    res[1]            = rho_eqn residual        (-> s row 2)               *
 *    res[2+j]          = Ncol_eqn[j] residual    (-> s row 4+j)             *
 *    res[2+NSPECIES+j] = ion_eqn[j]  residual    (-> s row 4+NSPECIES+j)  *
 *    res[2+2*NSPECIES] = T_eqn   residual        (-> s row 4+2*NSPECIES)   *
 *    res[3+2*NSPECIES] = F_eqn   residual        (-> s row 5+2*NSPECIES)   *
 *----------------------------------------------------------------------------*/
void eval_interior_eqns(int k, double *x, double **y, double (*ymod)[2],
                        double res[]) {
  int j;
  I_EQNVARS kvars, km1vars, avgvars;
  double delta_r_k;
  double dvdr_avg, drhodr_avg;
  double dNdr_avg[NSPECIES_MAX], dYsdr_avg[NSPECIES];
  double dTdr_avg;

  /*------------------------------------------------------------------------*
   * Shared setup — computed ONCE for all five equation types               *
   *------------------------------------------------------------------------*/
  set_vars(k, x, y, &kvars, &km1vars, &avgvars, ymod);

  /* erf_norm normalisation (previously updated redundantly in each eqn fn) */
  if (k == 2) {
      erf_norm = (1.0 - erf((kvars.v - parameters.erf_drop[0])
                            / parameters.erf_drop[1]));
  }

  /* Spatial step shared by all finite-difference residuals */
  delta_r_k = avgvars.z * (kvars.q - km1vars.q);

  /*------------------------------------------------------------------------*
   * Derivative chain — each quantity computed exactly once                 *
   * get_dvdr -> get_spQ -> glq_heating sets the GLQ cache; everything      *
   * below reuses that cached result via the last_N guard in glq_rates.c.   *
   *------------------------------------------------------------------------*/
  get_dvdr   (&dvdr_avg,   &avgvars);
  get_drhodr (&drhodr_avg, &avgvars, dvdr_avg);
  get_dNcoldr(dNdr_avg,    &avgvars);
  get_dYsdr  (dYsdr_avg,   &avgvars, k, 1);
  get_dTdr   (&dTdr_avg,   &avgvars, drhodr_avg, dYsdr_avg, k);

  /*------------------------------------------------------------------------*
   * Assemble residuals inline (identical formulae to the original eqn fns) *
   *------------------------------------------------------------------------*/

  /* v_eqn: delta_v - dvdr * delta_r = 0 */
  res[0] = (kvars.v   - km1vars.v)   - dvdr_avg   * delta_r_k;

  /* rho_eqn: delta_rho - drhodr * delta_r = 0 */
  res[1] = (kvars.rho - km1vars.rho) - drhodr_avg * delta_r_k;

  /* Ncol_eqn[j]: delta_Ncol[j] - dNdr[j] * delta_r = 0 */
  for (j = 0; j < NSPECIES; j++) {
      res[2+j] = (kvars.Ncol[j] - km1vars.Ncol[j]) - dNdr_avg[j] * delta_r_k;
  }

  /* ion_eqn[j]: delta_Ys[j] - dYsdr[j] * delta_r = 0 */
  for (j = 0; j < NSPECIES; j++) {
      res[2+NSPECIES+j] = (kvars.Ys[j] - km1vars.Ys[j]) - dYsdr_avg[j] * delta_r_k;
  }

  /* T_eqn: delta_T - dTdr * delta_r = 0  (or F/kappa_norm if conduction ON) */
  res[2+2*NSPECIES] = T_eqn(k, x, y, ymod, 1);

  /* F_eqn: heat-flux evolution equation (trivially F=0 when conduction OFF) */
  res[3+2*NSPECIES] = F_eqn(k, x, y, ymod, 1);
}

/*============================================================================*
 *! \fn double get_kappa_norm(const I_EQNVARS *vars)                                 *
 *  \brief Returns dimensionless conductivity kappa_norm = kappa_phys         *
 *         * T0 / (RHO0 * CS0^3 * Rp).                                       *
 *                                                                            *
 *  This is the factor relating the code-unit heat flux F_code and the        *
 *  temperature gradient: dT_code/dr_code = F_code / kappa_norm,             *
 *  where F_code = F_phys / (RHO0 * CS0^3) is the dimensionless heat flux.   *
 *  Conductivities follow Banks & Kockarts (1973) eqs. 8-11.                 *
 *----------------------------------------------------------------------------*/
double get_kappa_norm(const I_EQNVARS *vars) {
  double T_phys  = vars->T * T0;
  double rho_phys = vars->rho * RHO0;
  int jc;
  double ne_cond = 0.0, nn_cond = 0.0;

  for (jc = 0; jc < NSPECIES; jc++) {
    int   childq  = get_child_species(jc);
    int   eps_j   = parameters.Z[jc] - parameters.N_e[jc];
    double n_tot_j = rho_phys * get_eff_ntot_over_rho(jc, vars->Ys);
    double n0_j   = vars->Ys[jc] * n_tot_j;
    double n_ion_j = (1.0 - vars->Ys[jc]) * n_tot_j;
    if (eps_j == 0){
        nn_cond += n0_j;}
    if (childq >= 0)  ne_cond += n0_j * eps_j;
    else              ne_cond += n0_j * eps_j + n_ion_j * (eps_j + 1);
  }
  double n_tot_cond = ne_cond + nn_cond;
  if (n_tot_cond <= 0.0 || T_phys <= 0.0) return 1.0; /* guard */

  /* Electron conductivity (Banks & Kockarts 1973, eq. 9) */
  double kappa_ei = 1.2e-6 * pow(T_phys, 2.5);
  double sum_inv_ken = 0.0;
  for (jc = 0; jc < NSPECIES; jc++) {
    double n0_j = vars->Ys[jc] * rho_phys * get_eff_ntot_over_rho(jc, vars->Ys);
    if (n0_j > 0.0)
      sum_inv_ken += 1.0 / (6.0e4 * sqrt(T_phys) * (ne_cond / n0_j));
  }
  double kappa_e = kappa_ei / (1.0 + kappa_ei * sum_inv_ken);

  /* Neutral conductivity (Banks & Kockarts 1973, eq. 11) */
  double kappa_n_num = 0.0;
  double Aj_HI = 1.22e-6, Aj_HeI = 3.84e-6, Aj_OI = 3.9e-6;
  for (jc = 0; jc < NSPECIES; jc++) {
    double mj   = parameters.atomic_mass[jc];
    double n0_j = vars->Ys[jc] * rho_phys * get_eff_ntot_over_rho(jc, vars->Ys);
    double Aj;
    if      (strcmp(parameters.species[jc], "HI")  == 0) Aj = Aj_HI;
    else if (strcmp(parameters.species[jc], "HeI") == 0) Aj = Aj_HeI;
    else if (strcmp(parameters.species[jc], "OI")  == 0) Aj = Aj_OI;
    else                                                   Aj = Aj_HI;
    kappa_n_num += (15.0/4.0) * (K / mj) * Aj * pow(T_phys, 0.69) * n0_j;
  }
  double kappa_n   = (nn_cond > 0.0) ? kappa_n_num / nn_cond : 0.0;
  double kappa_tot = (ne_cond / n_tot_cond) * kappa_e
                   + (nn_cond / n_tot_cond) * kappa_n;

  /* kappa_norm = kappa_phys * T0 / (RHO0 * CS0^3 * Rp) */
  double kappa_norm = kappa_tot * T0 / (RHO0 * pow(CS0, 3) * parameters.Rp);
  return (kappa_norm > 0.0) ? kappa_norm : 1.0;
}

/*============================================================================*
 *! \fn static void get_dvdr(double *dvdr, const I_EQNVARS *vars)                    *
 *  \brief Calculates the derivative of velocity with respects to r, both     *
 *         analytically and using L'Hopital's rule to approximate it near the *
 *         sonic point where it analytically tends towards 0/0.               *
 *----------------------------------------------------------------------------*/

void get_dvdr(double *dvdr, const I_EQNVARS *vars) {
  double mu, rad, gamma, Mach;
  double spQ, d1_phi, norm_semimajor, qinv;
  double weight, r, n, s;

  /* Initalize to zero */
  *dvdr = 0.;

  get_gamma(&gamma);
  get_mu(&mu, vars);
  Mach = vars->v/sqrt(gamma*vars->T/mu);
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
    get_spQ(&spQ, vars,1e4,1);
    norm_semimajor = parameters.semimajor/parameters.Rp;

    /* Potential derivative terms */
    d1_phi = 1./SQR(rad);
//     if (parameters.tidalforce == ON) {
      /* Inverse of tidal q parameter */
      qinv    = parameters.Mstar/parameters.Mp;
      d1_phi += ( -qinv/SQR(norm_semimajor-rad)
                 +(qinv*(norm_semimajor-rad) - rad)/pow(norm_semimajor, 3.) )*parameters.tidalforce;
//     }
    /* Scale potential terms */
    d1_phi *= (parameters.Rp/parameters.H0);

    /* Calculate the analytic velocity derivative */
    *dvdr = 2.*gamma*vars->T/(mu*rad)
            -(gamma-1.)*spQ/vars->v
            -d1_phi;
    /* Finish by common factors */
    *dvdr *= vars->v/(SQR(vars->v)-gamma*vars->T/mu);     
    *dvdr *= weight;
  }
  /* If not in unweighted analytic dvdr region, add linearized dvdr */
  if (weight != 1.) {
    *dvdr += (1.-weight)*(dvdr_slope*(vars->q-q_last)+dvdr_last);
  }

  return;
}

void linearize_dvdr_crit(double *x, double **y) {
  int k,j;
  double Mach, rad, gamma, mu, last_Mach;
  double dYsdr[NSPECIES_MAX], dNdr[NSPECIES], dmudr, dgdr=0;
  double d1_phi, d2_phi, norm_semimajor, qinv;
  double spQ, spQ1, spQ2;
  VARLIST dQ;
  double a, b, c, sqrtarg, dvdr_crit;
  I_EQNVARS vars;

  get_gamma(&gamma);
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
          vars.Ys[j]   = y[6+j][k];
          vars.Ncol[j] = y[6+NSPECIES+j][k];
      }

    get_mu(&mu, &vars);
    Mach = vars.v/sqrt(gamma*vars.T/mu);
    if (Mach < last_Mach) {
      vars.q    = x[k+1];
      vars.v    = y[1][k+1];
      vars.z    = y[2][k+1];
      vars.rho  = y[3][k+1];
      vars.T    = y[4][k+1];
      for (j=0; j<NSPECIES; j++){
          vars.Ys[j]   = y[6+j][k+1];
          vars.Ncol[j] = y[6+NSPECIES+j][k+1];
      }

      break;
    }
    k--;
  } while(k > 0);

  /* Set dvdr_last and q_last */
  get_rad(&rad, &vars);
  get_spQ(&spQ, &vars,1e4,1);
  norm_semimajor = parameters.semimajor/parameters.Rp;

  /* Potential derivative terms */
  d1_phi = 1./SQR(rad);
//   if (parameters.tidalforce == ON) {
    /* Inverse of tidal q parameter */
    qinv    = parameters.Mstar/parameters.Mp;
    d1_phi += ( -qinv/SQR(norm_semimajor-rad)
                 +(qinv*(norm_semimajor-rad) - rad)/pow(norm_semimajor, 3.) )*parameters.tidalforce;
//   }
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
      vars.Ys[j]   = y[6+j][M];
      vars.Ncol[j] = y[6+NSPECIES+j][M];
  }
    
  get_rad(&rad, &vars);
  get_mu(&mu, &vars);
  get_spQ(&spQ, &vars,1e4,1);

  get_component_derivs(&dQ, &vars, QCOMP);
  get_dYsdr(dYsdr, &vars,1e4,1);
  get_dNcoldr(dNdr, &vars);
  /* Chain-aware dgdr — same product-rule as get_dTdr; see that function for derivation */
  {
    int ki, depth, chain[NSPECIES+1], cur, root_j;
    double mr_root, pfx[NSPECIES+2], sfx;
    dgdr = 0.0;
    //Chain rule of dgdr = d(1/mu)dr. (see get_dTdr for derivation)
    for (j=0; j<NSPECIES; j++){
      depth = 0; cur = j;
      //making an array of the indices of all ancestors of j, including itself, up to the root
      do { chain[depth++] = cur; cur = parameters.chain_parent[cur]; } while (cur >= 0); 
      root_j  = chain[depth-1];
      mr_root = parameters.atomic_mass[0] / parameters.atomic_mass[root_j]; //this is the mass ratio used in mu
      pfx[0] = 1.0;
      for (ki = 0; ki < depth; ki++)
        pfx[ki+1] = pfx[ki] * (1.0 - vars.Ys[chain[ki]]);
      sfx = 1.0;
      for (ki = depth-1; ki >= 0; ki--) {
        dgdr -= mr_root * parameters.HX[root_j] * dYsdr[chain[ki]] * pfx[ki] * sfx;
        sfx  *= (1.0 - vars.Ys[chain[ki]]);
      }
    }
  }
  dmudr = -pow(mu, 2)*dgdr;

  /* Potential derivative terms */
  norm_semimajor = parameters.semimajor/parameters.Rp;
  d2_phi = -2./pow(rad, 3.);
//   if (parameters.tidalforce == ON) {
    /* Inverse of tidal q parameter */
    qinv    = parameters.Mstar/parameters.Mp;
    d2_phi += (-2.*qinv/pow(norm_semimajor-rad, 3.) - (1.+qinv)/pow(norm_semimajor, 3.)) *parameters.tidalforce;
//   }
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
 *! \fn static void get_drhodr(double *drhodr, const I_EQNVARS *vars, double dvdr)   *
 *  \brief Calculates the derivative of density with respects to r            *
 *----------------------------------------------------------------------------*/

void get_drhodr(double *drhodr, const I_EQNVARS *vars, double dvdr) {
  double rad;

  get_rad(&rad, vars);

  *drhodr = -vars->rho*(2./rad+dvdr/vars->v);
    
  return;
}

/*============================================================================*
 *! \fn static void get_dTdr(double *dTdr, const I_EQNVARS *vars, double drhodr      *
 *                           double dYsdr)                                    *
 *  \brief Calculates the derivative of temperature with respects to r, both  *
 *         analytically and using L'Hopital's rule for the density derivative *
 *         term near the critial point where one term analytically tends      *
 *         towards 0/0                                                        *
 *----------------------------------------------------------------------------*/

void get_dTdr(double *dTdr, const I_EQNVARS *vars, double drhodr, double *dYsdr, double k) {
  double mu, gamma, dmudr, spQ, dgdr=0;
  int j;

  /* For mu(X) of the form 1/g(X), then set dgdr = d(g(X))/dr */
  get_mu(&mu, vars);
  
  /* dgdr = d/dr of the denominator in get_mu.
   * For each species j, the denominator contribution is:
   *   root: mr_j * HX[j] * (2+eps-Ys[j])  =>  d/dr = -mr_j*HX[j]*dYsdr[j]
   *   child at any depth: mr_root*HX[root] * Π_anc(1-Ys[anc]) * (1-Ys[j])
   *   d/dr = -mr_root*HX[root] * Σ_k [ dYsdr[chain[k]] * pfx[k] * sfx[k] ]
   * where chain[] lists j first then ancestors up to root,
   * pfx[k] = Π_{i<k}(1-Ys[chain[i]]) and sfx[k] = Π_{i>k}(1-Ys[chain[i]]).
   * Independent species (depth 1) reduces to the -mr*HX*dYsdr term. */

   /* dYsdr_CI ==> 1 + (1−Ys_CII) + (1−Ys_CII)(1−Ys_CIII) 
    * dYsdr_CII ==> (1−Ys_CI) + (1−Ys_CI)(1−Ys_CIII) 
    * dYsdr_CIII ==> (1−Ys_CI)(1−Ys_CII)                    */
  {
    int ki, depth, chain[NSPECIES+1], cur, root_j;
    double mr_root, pfx[NSPECIES+2], sfx;
    dgdr = 0.0;
    for (j=0; j<NSPECIES; j++){
      /* Build chain: chain[0]=j, chain[1]=parent(j), ..., chain[depth-1]=root */
      depth = 0; cur = j;
      do { chain[depth++] = cur; cur = parameters.chain_parent[cur]; } while (cur >= 0);
      root_j  = chain[depth-1];
      mr_root = parameters.atomic_mass[0] / parameters.atomic_mass[root_j];
      /* Prefix products */
      pfx[0] = 1.0;
      for (ki = 0; ki < depth; ki++)
        pfx[ki+1] = pfx[ki] * (1.0 - vars->Ys[chain[ki]]);
      /* Accumulate dgdr with running suffix (backwards pass) */
      sfx = 1.0;
      for (ki = depth-1; ki >= 0; ki--) {
        dgdr -= mr_root * parameters.HX[root_j] * dYsdr[chain[ki]] * pfx[ki] * sfx;
        sfx  *= (1.0 - vars->Ys[chain[ki]]);
      }
    }
  }
  dmudr = -pow(mu, 2)*dgdr;

  get_gamma(&gamma);
  get_spQ(&spQ, vars,k,0);

  *dTdr = (gamma-1.0)*(spQ/vars->v*mu
                       +vars->T/vars->rho*drhodr)
          +vars->T/mu*dmudr; //also, none of these in right units maybe?

  /* --- physical (volumetric, erg/cm^3/s) PdV and advection terms --- */
  // double rho_phys = RHO0*vars->rho;
  // double v_phys   = CS0*vars->v;
  // double T_phys   = T0*vars->T;
  // double P_phys   = rho_phys*K*T_phys/(mu*MH);

  // double drhodr_phys = (RHO0/parameters.Rp)*drhodr;      /* physical d(rho)/dr */
  // double dTdr_phys   = (T0/parameters.Rp)*(*dTdr);        /* physical d(T)/dr   */
  // double dmudr_phys  = dmudr/parameters.Rp;               /* physical d(mu)/dr  */

  // double cool_PdV_phys = P_phys*v_phys/rho_phys*drhodr_phys;

  // double de_therm_dr = K/((gamma-1.0)*MH)
  //                     * (dTdr_phys/mu - T_phys*dmudr_phys/(mu*mu));
  // double heat_advect_phys = -v_phys*rho_phys*de_therm_dr;

  // FILE *fp = fopen("outputs/ion_rates.txt", "a"); // Open for appending
  // if (fp == NULL) {
  //     perror("Failed to open ion_rates.txt");
  //     exit(1);
  // }
  // fprintf(fp, "%.0f,%.5e,%.5e\n",k,cool_PdV_phys,heat_advect_phys);
  // fclose(fp);

  return;
}

/*============================================================================*
 *! \fn static void get_dYsdr(double *dYsdr, const I_EQNVARS *vars)                  *
 *  \brief Calculates the derivative of ionization fraction with respects to r*
 *----------------------------------------------------------------------------*/

void get_dYsdr(double *dYsdr, const I_EQNVARS *vars, double k, int print_rates) {
  int q,qq;
  double ion_rate, rec_rate, mu, ne=0;
  double alpharec[NSPECIES_MAX];
  double n_tot_q, n0_q, n0_qq, n_ion_q;
    
  get_mu(&mu, vars);
  
  for (q=0; q<NSPECIES; q++){
    int childq = get_child_species(q);
    int electrons_per_species = parameters.Z[q] - parameters.N_e[q];
    double n_tot_q = RHO0*vars->rho * get_eff_ntot_over_rho(q, vars->Ys);
    double n0_q    = vars->Ys[q]*n_tot_q;
    double n_ion_q = (1.0-vars->Ys[q])*n_tot_q;
    if (childq >= 0) {
      /* This species has a child (e.g. CI->CII): the "ionized" pool of this
       * species is the neutral pool of the child. Electrons from that pool are
       * already counted when the child is processed, so only count electrons
       * from the neutral state of this species here. */
      ne += n0_q * electrons_per_species;
    } else {
      ne += n0_q*electrons_per_species + n_ion_q*(electrons_per_species+1);
    }
  }

  /* Calculate derived quantities */
  get_alpharec(alpharec, vars, 1);
   
  /*Method for direct printout (minimizes duplications in printout)*/
  // if (print_rates==0){
  //   FILE *fp = fopen("outputs/ion_rates.txt", "a"); // Open for appending
  //   if (fp == NULL) {
  //       perror("Failed to open ion_rates.txt");
  //       exit(1);
  //   }
  //   fprintf(fp, "%.0f,%.5e,%.5e\n",
  //     k,
  //     vars->Ys[0]*parameters.HX[0]/parameters.atomic_mass[0]*vars->rho*RHO0*glq_ionization(vars->Ncol,vars->Ys,vars->rho,k)[0]+vars->Ys[0]*parameters.HX[0]/parameters.atomic_mass[0]*glq_secondary_ionization(vars->Ncol,vars->Ys,vars->rho,k,0)[0]+vars->Ys[1]*parameters.HX[1]/parameters.atomic_mass[1]*glq_secondary_ionization(vars->Ncol,vars->Ys,vars->rho,k,0)[1],
  //     vars->Ys[1]*parameters.HX[1]/parameters.atomic_mass[1]*vars->rho*RHO0*glq_ionization(vars->Ncol,vars->Ys,vars->rho,k)[1]+vars->Ys[0]*parameters.HX[0]/parameters.atomic_mass[0]*glq_secondary_ionization(vars->Ncol,vars->Ys,vars->rho,k,1)[0]+vars->Ys[1]*parameters.HX[1]/parameters.atomic_mass[1]*glq_secondary_ionization(vars->Ncol,vars->Ys,vars->rho,k,1)[1]);
  //     // glq_ionization(vars->Ncol,vars->Ys,vars->rho,k)[0],
  //     // glq_ionization(vars->Ncol,vars->Ys,vars->rho,k)[1]);
  //   fclose(fp);
  // }
    double *ion_rate_array = glq_ionization((double*)vars->Ncol, (double*)vars->Ys, vars->rho, k);

    ion_rate = 0.0;
    for (q=0; q<NSPECIES; q++){
      int childq = get_child_species(q);
      /* Effective total density for this species sub-pool (divided by rho).
       * For child species j with parent p: n_tot[j]/rho = (1-Ys[p])*HX[p]/m[p]
       * For independent species:          n_tot[j]/rho = HX[j]/m[j]          */
      n_tot_q  = get_eff_ntot_over_rho(q, vars->Ys);
      n0_q     = vars->Ys[q]*n_tot_q;
      n_ion_q  = (1.0-vars->Ys[q])*n_tot_q;

      /*Primary ionization rate of species q*/
      ion_rate = n0_q*ion_rate_array[q*(NSPECIES+1)];

      /*Adding the secondary ionization rate of species q*/
      /*ion_rate_array[q+1+(NSPECIES+1)*qq] returns primary_ion_rate[qq] multiplied
        by the number of secondary ionizations of species q per primary ionization
        of species qq: eta[q] = (E-IP[qq])*f_ion[q] / IP[q] */
      for (qq=0; qq<NSPECIES; qq++){
        n0_qq = vars->Ys[qq] * get_eff_ntot_over_rho(qq, vars->Ys);
        ion_rate += n0_qq*ion_rate_array[q+1+(NSPECIES+1)*qq];
      }

      /* Recombination: for a parent species that has a child (e.g. CI->CII),
       * the population recombining back into q is the lower (neutral) state of
       * the child: n_CII_neutral = Ys[child]*(1-Ys[q])*(HX[q]/m[q]).
       * For an independent species (no child), it is the usual n_ion_q.       */
      double n_recomb;
      if (childq >= 0) {
        n_recomb = vars->Ys[childq] * (1.0-vars->Ys[q])
                   * (parameters.HX[q]/parameters.atomic_mass[q]);
      } else {
        n_recomb = n_ion_q;
      }

      rec_rate = (alpharec[q])*ne*n_recomb;
      /* Sources minus sinks */
      dYsdr[q]  = (rec_rate-ion_rate)/(vars->v*n_tot_q) * parameters.Rp/CS0;
      /* Convert to dimensionless derivative */
    }

  /*==========================================================================*
   * Cascade dilution correction                                               *
   *                                                                           *
   * Ys_CII is defined as n_CII / ((1-Ys_CI)*n_C).  When CI ionises, new CII  *
   * atoms enter the denominator (1-Ys_CI)*n_C, diluting (or concentrating)   *
   * Ys_CII.  The main loop above evolves Ys_CII as if the denominator is      *
   * fixed; this post-processing step adds the missing term:                   *
   *                                                                           *
   *   d(Ys_child)/dr += -(1-Ys_child)/(1-Ys_parent) * d(Ys_parent)/dr       *
   *                                                                           *
   * Derivation: Ys_child = n_child / D, D = (1-Ys_parent)*n_C               *
   *   dYs_child/dr = [d(n_child)/dr * D - n_child * d(D)/dr] / D^2           *
   *   d(D)/dr = -n_C*dYs_parent/dr  =>  extra term = -(1-Ys_child)*dYs_parent/dr/(1-Ys_parent)
   *                                                                           *
   * Without this, fraction (1-Ys_child) of each ionised-parent atom is       *
   * incorrectly credited to the grandchild state instead of the child.        *
   *==========================================================================*/
  for (q=0; q<NSPECIES; q++) {
    int p = parameters.chain_parent[q];
    if (p >= 0) {
      double denom = 1.0 - vars->Ys[p];
      if (denom > 1.0e-30)   /* guard against fully-ionised parent */
        dYsdr[q] -= (1.0 - vars->Ys[q]) / denom * dYsdr[p];
    }
  }

  return;
}

/*============================================================================*
 *! \fn static void get_dNcoldr(double *dNdr, const I_EQNVARS *vars)                 *
 *  \brief Calculates the derivative of column density (optical depth) with   *
 *         respects to r                                                      *
 *----------------------------------------------------------------------------*/
void get_dNcoldr(double *dNdr, const I_EQNVARS *vars) {
  double n_abs;
  int j;
  int parent;

  for (j=0; j<NSPECIES; j++){
    parent = parameters.chain_parent[j];
    /* The column density tracks the photon-absorbing (lower) state of species j.
     * For an independent species this is simply Ys[j]*n_j.
     * For a child species (e.g. CII), the full sub-pool is (1-Ys[parent])*n_C,
     * and the absorbing fraction within that pool is Ys[j], giving
     *   n_abs = Ys[j]*(1-Ys[parent])*n_C  */
    if (parent >= 0) {
      n_abs = (RHO0/parameters.atomic_mass[parent])
              * parameters.HX[parent] * vars->rho
              * vars->Ys[j] * (1.0 - vars->Ys[parent]);
    } else {
      n_abs = (RHO0/parameters.atomic_mass[j])
              * parameters.HX[j] * vars->rho
              * vars->Ys[j];
    }
    dNdr[j] = -n_abs * parameters.Rp/NCOL0;
  }

  return;
}

/*============================================================================*
 *! \fn static void get_spQ(double *spQ, const I_EQNVARS *vars, double k)            *
 *  \brief Calculates the total specific heating,                             *
 *          k included for diagnostic print-outs                              *
 *----------------------------------------------------------------------------*/

void get_spQ(double *spQ, const I_EQNVARS *vars,double k,int printout) {
  double n0overrho, nIONoverrho=0,nIONjoverrho=0;
  double ne=0,n_ion_j=0,n0_j=0,n_tot_j=0;
  double spQ_photoheat=0.0, spQ_lyacool=0.0, spQ_reccool=0.0, spQ_linecool=0.0;
  double spQ_boloheat=0.0, spQ_bolocool=0.0, smoothing_erf;
  double kappa_opt = parameters.kappa_opt, kappa_IR = parameters.kappa_IR;
  int j,sp,line,electrons_per_species;

  /* The heating and cooling terms */
  *spQ = 0.0;
  /* photo-ionization heating */
  for (j=0; j<NSPECIES; j++){
      /* n0 = lower (absorbing) state density / rho — chain-aware */
      n0overrho = vars->Ys[j] * get_eff_ntot_over_rho(j, vars->Ys);
      spQ_photoheat  += n0overrho*glq_heating((double*)vars->Ncol,(double*)vars->Ys,vars->rho,k)[j];
    } 
  // if (printout==0){
  //     FILE *fp = fopen("outputs/ion_rates.txt", "a"); // Open for appending
  //     if (fp == NULL) {
  //         perror("Failed to open ion_rates.txt");
  //         exit(1);
  //     }
  //     fprintf(fp, "%.0f,%.5e\n",k,spQ_photoheat*RHO0*vars->rho);
  //     fclose(fp);
  //   }
  spQ_photoheat *= (parameters.Rp/pow(CS0, 3));
  *spQ += spQ_photoheat;
  
   if (parameters.bolo_heat_cool > 0) {
        /* Bolometric Heating */
        /* Kappa optical and infrared are not valid at low pressures, so we send them 0 
        (over an order of mag in pressure) via an error function */
       smoothing_erf = get_erf_smooth(vars);
       
        spQ_boloheat = parameters.Lstar/(4*PI*pow(parameters.semimajor,2)) *(kappa_opt*smoothing_erf+0.25*kappa_IR*smoothing_erf); // divided by rho
        *spQ += parameters.bolo_heat_cool*spQ_boloheat*(parameters.Rp/pow(CS0, 3)); 

        /* Bolometric cooling */
        spQ_bolocool = -2*SIG_SB*pow(vars->T*T0,4)*kappa_IR*smoothing_erf; // divided by rho
        *spQ += parameters.bolo_heat_cool*spQ_bolocool*(parameters.Rp/pow(CS0, 3)); 
   }
  
  /* === Electron number density — needed by lyacool, reccool, AND freefree ===
   * Accumulated unconditionally so all three terms have access to ne.         */
  ne = 0.0;
  for (j=0; j<NSPECIES; j++){
    int childj = get_child_species(j);
    electrons_per_species = parameters.Z[j] - parameters.N_e[j];
    n_tot_j = RHO0*vars->rho * get_eff_ntot_over_rho(j, vars->Ys);
    n0_j    = vars->Ys[j]  * n_tot_j;
    n_ion_j = (1.0-vars->Ys[j]) * n_tot_j;
    if (childj >= 0) {
      ne += n0_j * electrons_per_species;
    } else {
      ne += n0_j*electrons_per_species + n_ion_j*(electrons_per_species+1);
    }
  }

  /* Lyman-alpha cooling */
  if (parameters.linecool == ON) {
    double nHIoverrho = vars->Ys[0]*parameters.HX[0]/(parameters.atomic_mass[0]);
    spQ_lyacool  = LYACOOL_COEFF*exp(-LYACOOL_TEMP/vars->T)*(nHIoverrho)*ne; //LYACOOL_COEFF is negative
    // if (printout==0){
    //   FILE *fp = fopen("outputs/ion_rates.txt", "a"); // Open for appending
    //   if (fp == NULL) {
    //       perror("Failed to open ion_rates.txt");
    //       exit(1);
    //   }
    //   fprintf(fp, "%.0f,%.5e\n",k,spQ_lyacool*RHO0*vars->rho);
    //   fclose(fp);
    // }
    spQ_lyacool *= (parameters.Rp/pow(CS0, 3));
    /* Note that LYACOOL_COEFF < 0 */
    *spQ += spQ_lyacool;
 
    
  /* Line Cooling*/
  /* OII, OIII, CII, CIII are relevant; James Owen, private correspondence */
  /* FeII, MgII, CaII, NeIII are relevant at higher metallicities (Linssen et al. 2024) */
  /* nIONjoverrho = density of the ionized (emitting) state divided by rho.
   * For chain species (e.g. CII with parent CI):
   *   n_CII_ionized / rho = (1-Ys[CII]) * (1-Ys[CI]) * HX[CI]/m[CI]
   * For independent species (e.g. OI with no parent):
   *   n_OII_ionized / rho = (1-Ys[OI]) * HX[OI]/m[OI]
   *
   * lc_* indices are precomputed by init_soe(); each block runs only when
   * the species is present, with no strcmp in the hot path.                  */

    /* CII line cooling (emitter = ionized state of CI slot) */
    if (lc_CI >= 0) {
        sp = lc_CI;
        nIONjoverrho = (1.0-vars->Ys[sp]) * get_eff_ntot_over_rho(sp, vars->Ys);
        for (line=0; line<3; line++) {
            spQ_linecool = nIONjoverrho*ne * CII[line][1]*exp(-CII[line][2]/(vars->T*T0))
                           / (ne*(1+CII[line][3]/ne));
            spQ_linecool *= (parameters.Rp/pow(CS0, 3));
            *spQ += spQ_linecool;
        }
    }

    /* CIII line cooling (emitter = ionized state of CII slot) */
    if (lc_CII >= 0) {
        sp = lc_CII;
        nIONjoverrho = (1.0-vars->Ys[sp]) * get_eff_ntot_over_rho(sp, vars->Ys);
        for (line=0; line<2; line++) {
            spQ_linecool = nIONjoverrho*ne * CIII[line][1]*exp(-CIII[line][2]/(vars->T*T0))
                           / (ne*(1+CIII[line][3]/ne));
            spQ_linecool *= (parameters.Rp/pow(CS0, 3));
            *spQ += spQ_linecool;
        }
    }

    /* OII line cooling */
    if (lc_OI >= 0) {
        sp = lc_OI;
        nIONjoverrho = (1.0-vars->Ys[sp]) * get_eff_ntot_over_rho(sp, vars->Ys);
        for (line=0; line<4; line++) {
            spQ_linecool = nIONjoverrho*ne * OII[line][1]*exp(-OII[line][2]/(vars->T*T0))
                           / (ne*(1+OII[line][3]/ne));
            spQ_linecool *= (parameters.Rp/pow(CS0, 3));
            *spQ += spQ_linecool;
        }
    }

    /* OIII line cooling */
    if (lc_OII >= 0) {
        sp = lc_OII;
        nIONjoverrho = (1.0-vars->Ys[sp]) * get_eff_ntot_over_rho(sp, vars->Ys);
        for (line=0; line<4; line++) {
            spQ_linecool = nIONjoverrho*ne * OIII[line][1]*exp(-OIII[line][2]/(vars->T*T0))
                           / (ne*(1+OIII[line][3]/ne));
            spQ_linecool *= (parameters.Rp/pow(CS0, 3));
            *spQ += spQ_linecool;
        }
    }

    /* FeII line cooling */
    if (lc_FeI >= 0) {
        sp = lc_FeI;
        nIONjoverrho = (1.0-vars->Ys[sp]) * get_eff_ntot_over_rho(sp, vars->Ys);
        for (line=0; line<3; line++) {
            spQ_linecool = nIONjoverrho*ne * FeII[line][1]*exp(-FeII[line][2]/(vars->T*T0))
                           / (ne*(1+FeII[line][3]/ne));
            spQ_linecool *= (parameters.Rp/pow(CS0, 3));
            *spQ += spQ_linecool;
        }
    }

    /* MgII line cooling */
    if (lc_MgI >= 0) {
        sp = lc_MgI;
        nIONjoverrho = (1.0-vars->Ys[sp]) * get_eff_ntot_over_rho(sp, vars->Ys);
        for (line=0; line<1; line++) {
            spQ_linecool = nIONjoverrho*ne * MgII[line][1]*exp(-MgII[line][2]/(vars->T*T0))
                           / (ne*(1+MgII[line][3]/ne));
            spQ_linecool *= (parameters.Rp/pow(CS0, 3));
            *spQ += spQ_linecool;
        }
    }

    /* CaII line cooling */
    if (lc_CaI >= 0) {
        sp = lc_CaI;
        nIONjoverrho = (1.0-vars->Ys[sp]) * get_eff_ntot_over_rho(sp, vars->Ys);
        for (line=0; line<4; line++) {
            spQ_linecool = nIONjoverrho*ne * CaII[line][1]*exp(-CaII[line][2]/(vars->T*T0))
                           / (ne*(1+CaII[line][3]/ne));
            spQ_linecool *= (parameters.Rp/pow(CS0, 3));
            *spQ += spQ_linecool;
        }
    }

    /* NeIII line cooling */
    if (lc_NeII >= 0) {
        sp = lc_NeII;
        nIONjoverrho = (1.0-vars->Ys[sp]) * get_eff_ntot_over_rho(sp, vars->Ys);
        for (line=0; line<3; line++) {
            spQ_linecool = nIONjoverrho*ne * NeIII[line][1]*exp(-NeIII[line][2]/(vars->T*T0))
                           / (ne*(1+NeIII[line][3]/ne));
            spQ_linecool *= (parameters.Rp/pow(CS0, 3));
            *spQ += spQ_linecool;
        }
    }
  }
  /* Collisions with neutral H are relevant when e- density low => CI and OI line cooling matter */

    /* Recombination cooling: chain-aware.
     * For a parent species p with a child c (e.g. CI->CII):
     *   - Recombination of CII->CI uses n_CII_ion = (1-Ys[c])*(1-Ys[p])*n_C
     *     with the alpharec for the CII slot (child).
     *   - The original CI slot still has its own alpharec for CI->C0 if needed,
     *     but since CI is never truly neutral here, we skip its ion contribution
     *     and let the child's entry handle it.
     * For an independent species (Z==1 hydrogen, or species with no chain link):
     *   uses the standard n_ion/n_neutral logic as before.                    */
    if (parameters.recombo_cool) {
    double alpharec_higher_ion_state[NSPECIES_MAX];
    double alpharec_lower_ion_state[NSPECIES_MAX];

    get_alpharec(alpharec_higher_ion_state,vars,1);
    get_alpharec(alpharec_lower_ion_state,vars,0);

    spQ_reccool = 0.0;
    // /*OLD METHOD: Using H type B recombination for all species*/
    //   for (j=0; j<NSPECIES; j++){
    //   nIONoverrho  += (1.-vars->Ys[j])*parameters.HX[j]/(parameters.atomic_mass[j]); 
    //   //TOTAL number density of ionized species
    // }
    // nIONoverrho = 0.0; // Reset for next use
    // //RECCOOL_COEFF has a factor of 1/T0^0.11 rolled in
    // spQ_reccool  = -RECCOOL_COEFF*pow(vars->T, 0.11)*(nIONoverrho)*ne; //changed 5/4/2023 to be negative
    // spQ_reccool *= (parameters.Rp/pow(CS0, 3));
    // *spQ += spQ_reccool;
    /*New method*/
    for (j=0; j<NSPECIES; j++){
        int parent_j = parameters.chain_parent[j];
        double eff_ntot_over_rho = get_eff_ntot_over_rho(j, vars->Ys);
        nIONoverrho = 0.0; n0overrho = 0.0;
        nIONoverrho = (1.0-vars->Ys[j]) * eff_ntot_over_rho;
        n0overrho   =      vars->Ys[j]  * eff_ntot_over_rho;

        if (parameters.Z[j] == 1){
          /* Case B recombination for HII — no chain link for hydrogen */
          spQ_reccool += -2.85e-27 * ne * nIONoverrho * sqrt(T0 * vars->T) * (5.914 - 0.5 * log(T0 * vars->T) + 0.01184 * pow(T0 * vars->T, 1.0 / 3.0));
        }
        else if (parameters.Z[j] != parameters.N_e[j]){
          /* Species whose lowest tracked state is already ionized (e.g. CII in
           * a chain where CI is also present): cool from the lower (neutral) state
           * of this slot back to the parent's ionized pool. */
          spQ_reccool += - alpharec_lower_ion_state[j] * n0overrho * ne * (1.5 * K * T0 * vars->T); 
        }
        else{
          /* Species whose lowest tracked state is neutral (e.g. CI, HeI, OI):
           * cool from the ionized pool. */
          spQ_reccool += - alpharec_higher_ion_state[j] * nIONoverrho * ne * (1.5 * K * T0 * vars->T); 
        }
        (void)parent_j; /* available for future per-chain adjustments */
    }
    spQ_reccool *= (parameters.Rp/pow(CS0, 3));
    *spQ += spQ_reccool;
    }

  /* === Free-free (bremsstrahlung) cooling ===
   * Λ_ff = -1.426e-27 * g_ff * T^{1/2} * n_e * Σ_j Z_j^2 * n_j
   * where the sum runs over all charged states of all species (Rybicki & Lightman §5.3).
   * Each species slot j contributes:
   *   ionized state (charge eps+1): (eps+1)^2 * n_ion_j
   *   neutral state (charge eps):   eps^2     * n_0_j     (zero unless standalone pre-ionized)
   * where eps = Z[j] - N_e[j] is the charge of the lower (neutral) state.
   * Reduces to the rmc2009 eq. A3 form for a pure-H plasma.                  */
  if (parameters.free_free_cool) {  /* free-free (bremsstrahlung) cooling */
    double Z2n_over_rho = 0.0;
    double spQ_freecool = 0.0;
    for (j=0; j<NSPECIES; j++){
      int eps = parameters.Z[j] - parameters.N_e[j];
      double eff = get_eff_ntot_over_rho(j, vars->Ys);
      double nIONj_over_rho = (1.0 - vars->Ys[j]) * eff;
      double n0j_over_rho   =        vars->Ys[j]  * eff;
      Z2n_over_rho += (double)(eps+1)*(eps+1) * nIONj_over_rho;
      if (eps > 0)
        Z2n_over_rho += (double)eps*eps * n0j_over_rho;
    }
    /* ne [cm^{-3}], Z2n_over_rho [cm^{-3} per (RHO0*vars->rho g/cm^3)]
     * T in code units (T0=1e4 K), so physical T = vars->T * T0               */
    spQ_freecool = -1.426e-27 * 1.3 * sqrt(vars->T * T0) * ne * Z2n_over_rho;
    spQ_freecool *= (parameters.Rp/pow(CS0, 3));
    *spQ += spQ_freecool;
  }

  return;
}

// /*===========================================================================*
//  * Returns the ionization potential from ion_pots array for given Z and N_e.
//  *---------------------------------------------------------------------------*/
// double get_ion_pot(int Z, int N_e, int num_rows) {
//     for (int i = 0; i < num_rows; i++) {
//         if ((int)ion_pots[i][0] == Z && (int)ion_pots[i][1] == N_e) {
//             return ion_pots[i][2];
//         }
//     }
//     // Not found: return 0 or error value
//     return 0.0;
// }

/*============================================================================*
 *! \fn static void get_rad(double *rad, const I_EQNVARS *vars)                      *
 *  \brief Calculates the radius given q and z                                *
 *----------------------------------------------------------------------------*/

void get_rad(double *rad, const I_EQNVARS *vars) {
  *rad = parameters.Rmin+vars->q*vars->z;

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
 *! \fn static ArLookup find_ar_lookup(int iz, int in)                        *
 *  \brief Scans the recombination tables ONCE at init time to find the row   *
 *         for (iz, in).  Never called in the hot path.                       *
 *----------------------------------------------------------------------------*/
static ArLookup find_ar_lookup(int iz, int in) {
    ArLookup lut;
    int i;

    if (in < 3 || in == 11 || (iz > 5 && iz < 9) || iz == 10) {
        lut.kind = AR_RNEW;
        lut.idx  = 0;
        for (i = 0; i < 126; i++) {
            if ((int)rnew[i][0] == iz && (int)rnew[i][1] == in) {
                lut.idx = i;
                break;
            }
        }
    } else if (iz == 26 && in < 13) {
        lut.kind = AR_FE;
        lut.idx  = 0;
        for (i = 0; i < 10; i++) {
            if ((int)fe[i][0] == in) {
                lut.idx = i;
                break;
            }
        }
    } else {
        lut.kind = AR_RREC;
        lut.idx  = 0;
        for (i = 0; i < 349; i++) {
            if ((int)rrec[i][0] == iz && (int)rrec[i][1] == in) {
                lut.idx = i;
                break;
            }
        }
    }
    return lut;
}

/*============================================================================*
 *! \fn void init_soe(void)                                                   *
 *  \brief One-time initialisation called from main() after set_parameters(). *
 *                                                                             *
 *  (1) Precomputes the alpharec table row indices for each species so that   *
 *      get_alpharec uses direct array access instead of linear scans.        *
 *  (2) Records integer indices for line-cooling species so that get_spQ      *
 *      avoids strcmp loops on every call.                                    *
 *----------------------------------------------------------------------------*/
void init_soe(void) {
    int j, iz, in;

    printf("Precomputing alpharec and line-cooling indices...\n");

    /* --- alpharec lookup --- */
    for (j = 0; j < NSPECIES; j++) {
        iz = parameters.Z[j];
        in = parameters.N_e[j];

        ar_hi[j] = find_ar_lookup(iz, in);       /* higher-ion state (calc_lower=1) */

        if (in < iz) {
            ar_lo[j] = find_ar_lookup(iz, in+1); /* lower-ion  state (calc_lower=0) */
        } else {
            ar_lo[j].kind = AR_ZERO;             /* already neutral — zero output    */
            ar_lo[j].idx  = 0;
        }
    }

    /* Invalidate T-caches so first call always computes fresh */
    ar_cache_T_hi = -1.0;
    ar_cache_T_lo = -1.0;

    /* --- line-cooling species indices --- */
    lc_CI = lc_CII = lc_OI = lc_OII = lc_FeI = lc_MgI = lc_CaI = lc_NeII = -1;
    for (j = 0; j < NSPECIES; j++) {
        if      (strcmp(parameters.species[j], "CI"  ) == 0) lc_CI   = j;
        else if (strcmp(parameters.species[j], "CII" ) == 0) lc_CII  = j;
        else if (strcmp(parameters.species[j], "OI"  ) == 0) lc_OI   = j;
        else if (strcmp(parameters.species[j], "OII" ) == 0) lc_OII  = j;
        else if (strcmp(parameters.species[j], "FeI" ) == 0) lc_FeI  = j;
        else if (strcmp(parameters.species[j], "MgI" ) == 0) lc_MgI  = j;
        else if (strcmp(parameters.species[j], "CaI" ) == 0) lc_CaI  = j;
        else if (strcmp(parameters.species[j], "NeII") == 0) lc_NeII = j;
    }

    printf("  Done.\n");
    return;
}

/*============================================================================*
 *! \fn static void get_alpharec(double *alpharec_p, const I_EQNVARS *kvars)         *
 *  \brief Calculates the temperature dependence of the recombination         *
 *----------------------------------------------------------------------------*/

static void get_alpharec(double *alpharec_p, const I_EQNVARS *kvars, int calc_lower_ion_state) {
    int j, i;
    double T, tt, r;
    ArLookup *lut;

    /* CLOUDY algorithm uses unscaled T */
    T = kvars->T * T0;

    /* ---------------------------------------------------------------------- *
     * T-keyed cache: get_alpharec depends only on T (species params fixed).   *
     * In the column-by-column Jacobian ~95% of calls share the same T and    *
     * return here without any arithmetic.                                     *
     * ---------------------------------------------------------------------- */
    if (calc_lower_ion_state == 1 && T == ar_cache_T_hi) {
        memcpy(alpharec_p, ar_cache_hi, NSPECIES * sizeof(double));
        return;
    }
    if (calc_lower_ion_state == 0 && T == ar_cache_T_lo) {
        memcpy(alpharec_p, ar_cache_lo, NSPECIES * sizeof(double));
        return;
    }

    /* Use precomputed row indices — direct array access, no linear scan */
    lut = (calc_lower_ion_state == 1) ? ar_hi : ar_lo;

    for (j = 0; j < NSPECIES; j++) {
        i = lut[j].idx;
        switch (lut[j].kind) {

            case AR_RNEW:
                /* Adapted from the CLOUDY rrfit fortran algorithm */
                tt = sqrt(T / rnew[i][4]);
                r  = rnew[i][2] / ( tt * pow(tt + 1.0, 1.0 - rnew[i][3])
                                    * pow(1.0 + sqrt(T / rnew[i][5]),
                                          1.0 + rnew[i][3]) );
                break;

            case AR_FE: {
                double tt2 = T * 1.0e-4;
                r = fe[i][1] / pow(tt2, fe[i][2] + fe[i][3] * log10(tt2));
                break;
            }

            case AR_RREC: {
                double tt2 = T * 1.0e-4;
                r = rrec[i][2] / pow(tt2, rrec[i][3]);
                break;
            }

            case AR_ZERO:
            default:
                /* Lower-state lookup when species is already neutral */
                r = 0.0;
                break;
        }
        alpharec_p[j] = r;
    }

    /* Store result in T-keyed cache for reuse by subsequent calls with same T */
    if (calc_lower_ion_state == 1) {
        ar_cache_T_hi = T;
        memcpy(ar_cache_hi, alpharec_p, NSPECIES * sizeof(double));
    } else {
        ar_cache_T_lo = T;
        memcpy(ar_cache_lo, alpharec_p, NSPECIES * sizeof(double));
    }

    return;
}



/*============================================================================*
 *! \fn static void get_mu(double *mu_p, const I_EQNVARS *kvars)                     *
 *  \brief Calculates the mean molecular weight, includes electrons           *
 *----------------------------------------------------------------------------*/

static void get_mu(double *mu_p, const I_EQNVARS *kvars) {
  /* 1/mu = sum_roots[mr_root * HX_root]   (one atom per atom, roots only)
   *       + n_e * m_H / rho               (one term per free electron)
   *
   * The electron term uses the same accumulation as get_spQ: for a species
   * with a child, only its neutral state contributes eps electrons (the
   * ionized pool is the child's sub-pool and counted there); for a leaf
   * species, both neutral (eps) and ionized (eps+1) states contribute.
   * This handles all cases -- HI, HeI, standalone CII, full chains -- without
   * any special-casing on eps or parent/child status.                       */
  int j, eps, childj;
  double ne_over_rho, n_tot_j, n0_j, n_ion_j;
  double denominator, smoothing_erf, mu_atom, mu_mol;

  if (strcmp(parameters.species[0],"HI") != 0)
    printf("WARNING: Mean molecular weight assumes first species is HI.\n");

  /* get_erf_smooth returns the normalised erfc weight; scale by molec_layer
   * (independent flag for the mu transition, separate from bolo_heat_cool). */
  smoothing_erf  = get_erf_smooth(kvars);
  smoothing_erf *= parameters.molec_layer;

  /* Atoms: one per atom, sum roots only (chains share the root's HX/m) */
  denominator = 0.0;
  for (j = 0; j < NSPECIES; j++) {
    if (parameters.chain_parent[j] < 0)
      denominator += (parameters.atomic_mass[0]/parameters.atomic_mass[j])
                   * parameters.HX[j];
  }

  /* Electrons: n_e * m_H / rho (dimensionless, same accumulation as get_spQ) */
  ne_over_rho = 0.0;
  for (j = 0; j < NSPECIES; j++) {
    childj  = get_child_species(j);
    eps     = parameters.Z[j] - parameters.N_e[j];
    n_tot_j = get_eff_ntot_over_rho(j, kvars->Ys);
    n0_j    =        kvars->Ys[j]  * n_tot_j;
    n_ion_j = (1.0 - kvars->Ys[j]) * n_tot_j;
    if (childj >= 0)
      ne_over_rho += n0_j * eps; //if you have a child, only count yourself
    else
      ne_over_rho += n0_j * eps + n_ion_j * (eps + 1); //if no child, then go ahead and count ion also
  }
  denominator += ne_over_rho * parameters.atomic_mass[0];

  mu_atom = (1.0 - smoothing_erf) / denominator;
  mu_mol  = smoothing_erf * parameters.molec_adjust;
  *mu_p   = mu_atom + mu_mol;
  return;
}

static void get_gamma(double *gamma_p) {
  *gamma_p = parameters.gamma;

  return;
}

/*============================================================================*
 *! \fn static void set_vars(int k, double *x, double **y, I_EQNVARS *kvars,  *
 *                           I_EQNVARS *km1vars, I_EQNVARS *avgvars,          *
 *                           double (*ymod)[2])                            *
 *  \brief Calculates the neighboring and average variables given ymod update *
 *  \note ymod[][0] modifies cell k-1 and ymod[][1] modifies cell k           *
 *----------------------------------------------------------------------------*/

static void set_vars(int k, double *x, double **y, I_EQNVARS *kvars,
                     I_EQNVARS *km1vars, I_EQNVARS *avgvars,
                     double (*ymod)[2]) {
  int j;
    
  kvars->q      = x[k];
  kvars->v      = y[1][k]+ymod[1][1];
  kvars->z      = y[2][k]+ymod[2][1];
  kvars->rho    = y[3][k]+ymod[3][1];
  kvars->T      = y[4][k]+ymod[4][1];
  for (j=0; j<NSPECIES; j++){
     kvars->Ys[j]     = y[6+j][k]+ymod[6+j][1];
     kvars->Ncol[j]   = y[6+NSPECIES+j][k]+ymod[6+NSPECIES+j][1];
  }

  km1vars->q    = x[k-1];
  km1vars->v    = y[1][k-1]+ymod[1][0];
  km1vars->z    = y[2][k-1]+ymod[2][0];
  km1vars->rho  = y[3][k-1]+ymod[3][0];
  km1vars->T    = y[4][k-1]+ymod[4][0];
  for (j=0; j<NSPECIES; j++){
      km1vars->Ys[j]   = y[6+j][k-1]+ymod[6+j][0];
      km1vars->Ncol[j] = y[6+NSPECIES+j][k-1]+ymod[6+NSPECIES+j][0];
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
 *                            double (*ymod)[2])                           *
 *  \brief Calculates the residual of the denominator for the critical point  *
 *         criteria, i.e., that the outflow is Mach 1 at the critical point   *
 *----------------------------------------------------------------------------*/

static double spv_eqn(int k, double *x, double **y, double (*ymod)[2], int species) {
  double mu, gamma, residual; 
  I_EQNVARS kvars;
  int j;

  kvars.q    = x[k];
  kvars.v    = y[1][k]+ymod[1][1];
  kvars.z    = y[2][k]+ymod[2][1];
  kvars.rho  = y[3][k]+ymod[3][1];
  kvars.T    = y[4][k]+ymod[4][1];
  for (j=0; j<NSPECIES; j++){
      kvars.Ys[j]   = y[6+j][k]+ymod[6+j][1];
      kvars.Ncol[j] = y[6+NSPECIES+j][k]+ymod[6+NSPECIES+j][1];
  }
  

  get_mu(&mu, &kvars);
  get_gamma(&gamma);

  residual = kvars.v-sqrt(kvars.T*gamma/mu)*parameters.breezeparam;

  return residual;
}

/*============================================================================*
 *! \fn static double spcrit_eqn(int k, double *x, double **y,                *
 *                               double (*ymod)[2])                        *
 *  \brief Calculates the residual of the numerator for the critical point    *
 *         criteria, i.e., the Bernoulli criteria                             *
 *----------------------------------------------------------------------------*/

static double spcrit_eqn(int k, double *x, double **y, double (*ymod)[2], int species) {
  double T_term, grav_term, Q_term, spQ;
  double r, mu, gamma, residual, qinv;
  I_EQNVARS kvars;
  int j;

  kvars.q    = x[k];
  kvars.v    = y[1][k]+ymod[1][1];
  kvars.z    = y[2][k]+ymod[2][1];
  kvars.rho  = y[3][k]+ymod[3][1];
  kvars.T    = y[4][k]+ymod[4][1];
  for (j=0; j<NSPECIES; j++){
      kvars.Ys[j]   = y[6+j][k]+ymod[6+j][1];
      kvars.Ncol[j] = y[6+NSPECIES+j][k]+ymod[6+NSPECIES+j][1];
  }

  /* Conservative work terms */
  r = parameters.Rmin+kvars.q*kvars.z;
  grav_term = -1.0/r;
//   if (parameters.tidalforce == ON) {
    double norm_semimajor = parameters.semimajor/parameters.Rp;
    qinv = parameters.Mstar/parameters.Mp;
//     grav_term += (qinv/(norm_semimajor-r) + (qinv+1)/(2*pow(norm_semimajor,3.))*term) *parameters.tidalforce;
    grav_term += (qinv*r*
                  (-(norm_semimajor-r+r/qinv)/pow(norm_semimajor, 3.)
                   +1./SQR(norm_semimajor-r)))*parameters.tidalforce;
//   }
  /* Scale gravitational term */
  grav_term *= parameters.Rp/parameters.H0;

  /* Enthalpy term */
  get_mu(&mu, &kvars);
  get_gamma(&gamma);
  T_term = 2.0*gamma*kvars.T/mu;

  /* Heating terms */
  get_spQ(&spQ, &kvars,k,1);
  Q_term = -(gamma-1.0)*spQ*r/kvars.v;

  /* Terms should sum to zero */
  residual = T_term+grav_term+Q_term;

  return residual;
}

/*============================================================================*
 *! \fn static double ion_eqn(int k, double *x, double **y,                   *
 *                          double (*ymod)[2])                             *
 *  \brief Calculates the residual of the ionization fraction equation        *
 *----------------------------------------------------------------------------*/

static double ion_eqn(int k, double *x, double **y, double (*ymod)[2], int species) {
  double delta_Ys_k, delta_r_k, dYsdr_avg[NSPECIES], residual;//,Hsc0;
  I_EQNVARS kvars, km1vars, avgvars;
  int spec_index = species-1;

  set_vars(k, x, y, &kvars, &km1vars, &avgvars, ymod);
  if (k==2){   
      erf_norm = (1-erf((kvars.v-parameters.erf_drop[0])/parameters.erf_drop[1]));
  }
  /* Calculate the finite difference */
  delta_r_k  = avgvars.z*(kvars.q-km1vars.q);
    
 /* Calculate averaged derivative */
  get_dYsdr(dYsdr_avg, &avgvars,k,0); 
    
  delta_Ys_k = kvars.Ys[spec_index]-km1vars.Ys[spec_index];
    
  delta_r_k  = avgvars.z*(kvars.q-km1vars.q);

  /* Calculate finite difference residual */
  residual = delta_Ys_k-dYsdr_avg[spec_index]*delta_r_k;

  return residual;
}

/*============================================================================*
 *! \fn static double Ncol_eqn(int k, double *x, double **y,                  *
 *                             double (*ymod)[2])                          *
 *  \brief Calculates the residual of the column density (optical depth)      *
 *         equation                                                           *
 *----------------------------------------------------------------------------*/

static double Ncol_eqn(int k, double *x, double **y, double (*ymod)[2], int species) {
  double delta_N_k, delta_r_k, dNdr_avg[NSPECIES], residual;//, Hsc0;
  I_EQNVARS kvars, km1vars, avgvars;
  int spec_index = species-1;

  set_vars(k, x, y, &kvars, &km1vars, &avgvars, ymod);
  if (k==2){
      erf_norm = (1-erf((kvars.v-parameters.erf_drop[0])/parameters.erf_drop[1]));
  }
  delta_r_k = avgvars.z*(kvars.q-km1vars.q);
  get_dNcoldr(dNdr_avg, &avgvars);
    
  /* Calculate the finite difference */
  delta_N_k = kvars.Ncol[spec_index]-km1vars.Ncol[spec_index];
  residual = delta_N_k-(dNdr_avg[spec_index])*delta_r_k;

  return residual;
}

/*============================================================================*
 *! \fn static double v_eqn(int k, double *x, double **y,                     *
 *                          double (*ymod)[2])                             *
 *  \brief Calculates the residual of the velocity equation                   *
 *----------------------------------------------------------------------------*/

static double v_eqn(int k, double *x, double **y, double (*ymod)[2], int species) {
  double delta_v_k, delta_r_k, dvdr_avg, residual;//,Hsc0;
  I_EQNVARS kvars, km1vars, avgvars, critvars;
  int j;

  set_vars(k, x, y, &kvars, &km1vars, &avgvars, ymod);
  if (k==2){
      erf_norm = (1-erf((kvars.v-parameters.erf_drop[0])/parameters.erf_drop[1]));
  }
  critvars.q    = x[M];
  critvars.v    = y[1][M];
  critvars.z    = y[2][M];
  critvars.rho  = y[3][M];
  critvars.T    = y[4][M];
  for (j=0; j<NSPECIES; j++){
      critvars.Ys[j]   = y[6+j][M];
      critvars.Ncol[j] = y[6+NSPECIES+j][M];  
  }
    
  /* Calculate the finite difference */
  delta_v_k = kvars.v-km1vars.v;
  delta_r_k = avgvars.z*(kvars.q-km1vars.q);

  /* Calculate averaged derivative */
  get_dvdr(&dvdr_avg, &avgvars);

  /* Calculate finite difference residue */
  residual = delta_v_k-dvdr_avg*delta_r_k;

  return residual;
}

/*============================================================================*
 *! \fn static double rho_eqn(int k, double *x, double **y,                   *
 *                            double (*ymod)[2])                           *
 *  \brief Calculates the residual of the density equation                    *
 *----------------------------------------------------------------------------*/

static double rho_eqn(int k, double *x, double **y, double (*ymod)[2], int species) {
  double delta_rho_k, delta_r_k, drhodr_avg, residual;//,Hsc0;
  double dvdr_avg;
  I_EQNVARS kvars, km1vars, avgvars, critvars;
  int j;

  set_vars(k, x, y, &kvars, &km1vars, &avgvars, ymod);
  if (k==2){
      erf_norm = (1-erf((kvars.v-parameters.erf_drop[0])/parameters.erf_drop[1]));
  }
  critvars.q    = x[M];
  critvars.v    = y[1][M];
  critvars.z    = y[2][M];
  critvars.rho  = y[3][M];
  critvars.T    = y[4][M];
  for (j=0; j<NSPECIES; j++){
      critvars.Ys[j]   = y[6+j][M];
      critvars.Ncol[j] = y[6+NSPECIES+j][M];
  }
    
  /* Calculate the finite difference */
  delta_rho_k = kvars.rho-km1vars.rho;
  delta_r_k   = avgvars.z*(kvars.q-km1vars.q);

  /* Calculate averaged derivatives */
  get_dvdr(&dvdr_avg, &avgvars);
  get_drhodr(&drhodr_avg, &avgvars, dvdr_avg);

  /* Calculate finite difference residue */
  residual = delta_rho_k-drhodr_avg*delta_r_k;

  return residual;
}

/*============================================================================*
 *! \fn static double T_eqn(int k, double *x, double **y,
 *                          double (*ymod)[2])
 *  \brief Calculates the residual of the temperature equation.
 *         When conduction is ON:  dT/dr = F/kappa_norm.
 *         When conduction is OFF: standard energy equation dT/dr.
 *----------------------------------------------------------------------------*/

static double T_eqn(int k, double *x, double **y, double (*ymod)[2], int species) {
  I_EQNVARS kvars, km1vars, avgvars;
  double delta_T_k, delta_r_k;

  set_vars(k, x, y, &kvars, &km1vars, &avgvars, ymod);
  if (k == 2)
    erf_norm = (1.0 - erf((kvars.v - parameters.erf_drop[0]) / parameters.erf_drop[1]));

  delta_T_k = kvars.T - km1vars.T;
  delta_r_k = avgvars.z * (kvars.q - km1vars.q);

  if (parameters.conduction) {
    /* Conduction ON: dT/dr = F/kappa_norm */
    double F_avg = 0.5 * (y[5][k] + ymod[5][1] + y[5][k-1] + ymod[5][0]);
    double kappa_norm_avg = get_kappa_norm(&avgvars);
    return delta_T_k - F_avg / kappa_norm_avg * delta_r_k;
  } else {
    /* Conduction OFF: standard energy equation */
    double dYsdr_avg[NSPECIES_MAX], drhodr_avg, dvdr_avg, dTdr_avg;
    get_dvdr(&dvdr_avg, &avgvars);
    get_drhodr(&drhodr_avg, &avgvars, dvdr_avg);
    get_dYsdr(dYsdr_avg, &avgvars, k, 1);
    get_dTdr(&dTdr_avg, &avgvars, drhodr_avg, dYsdr_avg, k);
    return delta_T_k - dTdr_avg * delta_r_k;
  }
}

/*============================================================================*
 *! \fn static double F_eqn(int k, ...)
 *  \brief Heat-flux evolution equation; closes the conduction system.
 *
 *  When conduction is OFF: trivially enforce F[k] = 0.
 *  When conduction is ON the full evolution equation is solved. After
 *  substituting dT/dr = F/kappa_norm into the energy equation, dF/dr
 *  can be written without any d2T/dr2 term:
 *    dF/dr = rho*spQ + v*T/mu*drhodr + v*rho*T*dmudr/(mu^2*(g-1))
 *            - v*rho*F/(mu*(g-1)*kappa_norm) - 2F/r
 *  FDE: F[k] - F[k-1] - dF/dr_avg * dr = 0
 *----------------------------------------------------------------------------*/
static double F_eqn(int k, double *x, double **y, double (*ymod)[2], int species) {
  if (!parameters.conduction) {
    return y[5][k] + ymod[5][1];   /* Pin F = 0 when conduction is off */
  }

  I_EQNVARS kvars, km1vars, avgvars;
  double spQ, dvdr, drhodr, dYsdr[NSPECIES];
  double mu, gamma, dgdr, dmudr;
  double F_k, F_km1, F_avg, kappa_norm_avg;
  double r_avg, delta_r_k, dFdr;
  int j;

  set_vars(k, x, y, &kvars, &km1vars, &avgvars, ymod);
  delta_r_k = avgvars.z * (kvars.q - km1vars.q);

  F_k   = y[5][k]   + ymod[5][1];
  F_km1 = y[5][k-1] + ymod[5][0];
  F_avg = 0.5 * (F_k + F_km1);

  kappa_norm_avg = get_kappa_norm(&avgvars);
  get_spQ(&spQ, &avgvars, k, 1);
  get_dvdr(&dvdr, &avgvars);
  get_drhodr(&drhodr, &avgvars, dvdr);
  get_dYsdr(dYsdr, &avgvars, k, 1);
  get_mu(&mu, &avgvars);
  get_gamma(&gamma);

  /* dmu/dr: from chain rule on the 1/mu denominator */
  dgdr = 0.0;
  for (j = 0; j < NSPECIES; j++) {
    int parent_j = parameters.chain_parent[j];
    double mr = parameters.atomic_mass[0] / parameters.atomic_mass[j];
    if (parent_j >= 0)
      dgdr += mr * parameters.HX[j]
              * ((1.0 - avgvars.Ys[j])        * dYsdr[parent_j]
               + (1.0 - avgvars.Ys[parent_j]) * dYsdr[j]);
    else
      dgdr -= mr * parameters.HX[j] * dYsdr[j];
  }
  dmudr = -mu * mu * dgdr;

  get_rad(&r_avg, &avgvars);

  dFdr =   avgvars.rho * spQ
         + avgvars.v * avgvars.T / mu * drhodr
         + avgvars.v * avgvars.rho * avgvars.T / (mu * mu * (gamma - 1.0)) * dmudr
         - avgvars.v * avgvars.rho * F_avg / (mu * (gamma - 1.0) * kappa_norm_avg)
         - 2.0 * F_avg / r_avg;

  return (F_k - F_km1) - dFdr * delta_r_k;
}


/*============================================================================*
 *! \fn static double get_component(const I_EQNVARS *vars, int compnum)              *
 *  \brief Returns the value of the specified component                       *
 *----------------------------------------------------------------------------*/

static double get_component(const I_EQNVARS *vars, int compnum) {
  double component;

  if (compnum == QCOMP) {
    get_spQ(&component, vars,1e4,1);
    component *= vars->rho;
  }
  else {
    fprintf(stderr, "ERROR: unknown component\n");
    exit(704);
  }

  return component;
}

/*============================================================================*
 *! \fn static void get_component_derivs(VARLIST *derivs, const I_EQNVARS *vars,     *
 *                                       int compnum)                         *
 *  \brief Calculates the numerical derivative of the given component with    *
 *         respects to all of the relax variables                             *
 *----------------------------------------------------------------------------*/

static void get_component_derivs(VARLIST *derivs, const I_EQNVARS *vars, int compnum) {
  double delta_var, comp_pl, comp_mi;
  I_EQNVARS varmod;
  int j;

  /* finite difference wrt rho */
  delta_var   = vars->rho/DERIVDIV;
  varmod      = *vars;
  varmod.rho  = vars->rho+delta_var/2.;
  comp_pl     = get_component(&varmod, compnum);
  varmod.rho  = vars->rho-delta_var/2.;
  comp_mi     = get_component(&varmod, compnum);
  derivs->rho = (comp_pl-comp_mi)/delta_var;

  /* finite difference wrt v */
  delta_var = vars->v/DERIVDIV;
  varmod    = *vars;
  varmod.v  = vars->v+delta_var/2.;
  comp_pl   = get_component(&varmod, compnum);
  varmod.v  = vars->v-delta_var/2.;
  comp_mi   = get_component(&varmod, compnum);
  derivs->v = (comp_pl-comp_mi)/delta_var;

  /* finite difference wrt T */
  delta_var = vars->T/DERIVDIV;
  varmod    = *vars;
  varmod.T  = vars->T+delta_var/2.;
  comp_pl   = get_component(&varmod, compnum);
  varmod.T  = vars->T-delta_var/2.;
  comp_mi   = get_component(&varmod, compnum);
  derivs->T = (comp_pl-comp_mi)/delta_var;

  /* finite difference wrt Ys */
  for (j=0; j<NSPECIES; j++){
      delta_var  = vars->Ys[j]/DERIVDIV;
      varmod     = *vars;
      varmod.Ys[j]  = vars->Ys[j]+delta_var/2.;
      comp_pl    = get_component(&varmod, compnum);
      varmod.Ys[j]  = vars->Ys[j]-delta_var/2.;
      comp_mi    = get_component(&varmod, compnum);
      derivs->Ys[j] = (comp_pl-comp_mi)/delta_var;
  }

  /* finite difference wrt Ncol */
  for (j=0; j<NSPECIES; j++){
      delta_var  = vars->Ncol[j]/DERIVDIV;
      varmod     = *vars;
      varmod.Ncol[j]  = vars->Ncol[j]+delta_var/2.;
      comp_pl    = get_component(&varmod, compnum);
      varmod.Ncol[j]  = vars->Ncol[j]-delta_var/2.;
      comp_mi    = get_component(&varmod, compnum);
      derivs->Ncol[j] = (comp_pl-comp_mi)/delta_var;
  }
    
  /* finite difference wrt z */
  delta_var = vars->z/DERIVDIV;
  varmod    = *vars;
  varmod.z  = vars->z+delta_var/2.;
  comp_pl   = get_component(&varmod, compnum);
  varmod.z  = vars->z-delta_var/2.;
  comp_mi   = get_component(&varmod, compnum);
  derivs->z = (comp_pl-comp_mi)/delta_var;

  return;
}
