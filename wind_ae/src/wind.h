#ifndef _WIND_H_
#define _WIND_H_
/*============================================================================*
 *! \file wind.h                                                              *
 *  \brief Contains typedefs unquie to the relax_ae code                      *
 *============================================================================*/

/* Structure for simulation input parameters */
#ifndef _PARAMLIST_DECLARE_T_
typedef struct _paramlist {
  /* planetary system properties */
  double Mp;
  double Rp;
  double Mstar;
  double semimajor;
  double Ftot;
  double Lstar;
  double H0;
  /* boundary conditions */
  double Rmin;
  double Rmax;
  double T_rmin;
  double Ys_rmin[NSPECIES_MAX];
  double rho_rmin;
  double Ncol_sp[NSPECIES_MAX];
  double erf_drop[2];
  /* term indicators */
  int    integrate_outward;
  double    tidalforce;
  int    linecool;
  double    bolo_heat_cool;
  int    conduction;       /* 1 = solve F-variable self-consistently; 0 = F=0 (dTdr from energy eqn) */
  int    recombo_cool;     /* 1 = recombination cooling on; 0 = off */
  int    free_free_cool;   /* 1 = free-free (bremsstrahlung) cooling on; 0 = off */
  /* technical parameters */
  double breezeparam;
  double rapidity;
  double erfn;
  double mach_limit;
  /* physics parameters */
  double HX[NSPECIES_MAX]; // Mass fraction in each species
  char species[NSPECIES_MAX][10]; // Species name in Roman numeral format
  double atomic_mass[NSPECIES_MAX]; // Atomic mass in g
  int Z[NSPECIES_MAX]; // Atomic number
  int N_e[NSPECIES_MAX]; // Number of electrons in atom
  double molec_adjust;
  double kappa_opt;      /* optical opacity for bolometric heating/cooling (default 4e-3) */
  double kappa_IR;       /* IR opacity for bolometric heating/cooling (default 1e-2) */
  double gamma;          /* adiabatic index (default 5.0/3.0) */
  double molec_layer;    /* erfc multiplier for mu adjustment, independent of bolo_heat_cool (default = bolo_heat_cool) */
  double add_params[N_ADD_PARAMS]; // This allows users to add additional vars to the soln files for future changes
  /* Chain ionization linkage: chain_parent[j] = index of parent species, or -1 if independent.
   * If chain_parent[j] = p, species j is a higher ionization state of the same element as p,
   * with N_e[j] = N_e[p]-1 (one more electron removed). The effective total number density
   * for species j is scaled by (1 - Ys[p]) relative to the base HX[p]/m[p] value.
   * This is auto-detected from Z[] and N_e[] in io.c after reading phys_params.inp. */
  int chain_parent[NSPECIES_MAX];
} PARAMLIST;
#define _PARAMLIST_DECLARE_T_
#endif /* _PARAMLIST_DECLARE_T_ */

/* Structure for the equation variables */
#ifndef _EQNVARS_DECLARE_T_
typedef struct _eqnvars {
  double r[TOTALPTS_MAX];    /* radius (units: Rp); r = z*q + Rmin */
  double rho[TOTALPTS_MAX];  /* density (units: RHO0), multiplied times fraction to get density for each species */
  double v[TOTALPTS_MAX];    /* velocity (units: CS0) */
  double z[TOTALPTS_MAX];    /* rs-rmin (units: Rp) */
  double **Ys;   /* ionization fraction:     Ys[0..TOTALPTS-1][0..NSPECIES-1]
                 * heap-allocated in main() after g_nspecies is known;
                 * freed in main() before exit. */
  double **Ncol; /* optical depth to photoionization: same layout as Ys */
  double T[TOTALPTS_MAX];    /* temperature (units: T0) */
  double q[TOTALPTS_MAX];    /* independent variable; 0 <= q <= 1 */
  double dTdr[TOTALPTS_MAX]; /* dT/dr = Fp/kappa_norm (code units: T0/Rp); populated after solve */
  double Fp[TOTALPTS_MAX];   /* conductive heat flux F = kappa*dT/dr, code units: F_phys/(RHO0*CS0^3) */
} EQNVARS;
#define _EQNVARS_DECLARE_T_
#endif /* _EQNVARS_DECLARE_T_ */

/* Structure for an individual set of equation variables */
#ifndef _I_EQNVARS_DECLARE_T_
typedef struct _i_eqnvars {
  double r;    /* radius (units: Rp); r = z*q + Rmin */
  double rho;  /* density (units: RHO0) */
  double v;    /* velocity (units: CS0) */
  double z;    /* rs-rmin (units: Rp) */
  double Ys[NSPECIES_MAX];   /* ionization fraction */
  double Ncol[NSPECIES_MAX]; /* optical depth to photoionization */
  double T;    /* temperature (units: T0) */
  double q;    /* independent variable; 0 <= q <= 1 */
} I_EQNVARS;
#define _I_EQNVARS_DECLARE_T_
#endif /* _I_EQNVARS_DECLARE_T_ */

/* An expression with a value for each variable */
#ifndef _VARLIST_DECLARE_T_
typedef struct _varlist {
  double rho;
  double v;
  double T;
  double Ys[NSPECIES_MAX];
  double Ncol[NSPECIES_MAX];
  double z;
} VARLIST;
#define _VARLIST_DECLARE_T_
#endif /* _VARLIST_DECLARE_T_ */
#endif /* _WIND_H_ */