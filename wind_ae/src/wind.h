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
  double Ys_rmin[NSPECIES];
  double rho_rmin;
  double Ncol_sp[NSPECIES];
  double erf_drop[2];
  /* term indicators */
  int    lyacool;
  double    tidalforce;
  double    bolo_heat_cool;
  int    integrate_outward;
  /* technical parameters */
  double breezeparam;
  double rapidity;
  double erfn;
  double mach_limit;
  /* physics parameters */
  double HX[NSPECIES]; // Mass fraction in each species
  char species[NSPECIES][10]; // Species name in Roman numeral format
  double atomic_mass[NSPECIES]; // Atomic mass in g
  int Z[NSPECIES]; // Atomic number
  int N_e[NSPECIES]; // Number of electrons in atom
  double molec_adjust;
  double add_params[N_ADD_PARAMS]; // This allows users to add additional vars to the soln files for future changes
} PARAMLIST;
#define _PARAMLIST_DECLARE_T_
#endif /* _PARAMLIST_DECLARE_T_ */

/* Structure for the equation variables */
#ifndef _EQNVARS_DECLARE_T_
typedef struct _eqnvars {
  double r[TOTALPTS];    /* radius (units: Rp); r = z*q + Rmin */
  double rho[TOTALPTS];  /* density (units: RHO0), multiplied times fraction to get density for each species */
  double v[TOTALPTS];    /* velocity (units: CS0) */
  double z[TOTALPTS];    /* rs-rmin (units: Rp) */
  double Ys[TOTALPTS][NSPECIES];   /* ionization fraction */
  double Ncol[TOTALPTS][NSPECIES]; /* optical depth to photoionization */
  double T[TOTALPTS];    /* temperature (units: T0) */
  double q[TOTALPTS];    /* independent variable; 0 <= q <= 1 */
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
  double Ys[NSPECIES];   /* ionization fraction */
  double Ncol[NSPECIES]; /* optical depth to photoionization */
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
  double Ys[NSPECIES];
  double Ncol[NSPECIES];
  double z;
} VARLIST;
#define _VARLIST_DECLARE_T_
#endif /* _VARLIST_DECLARE_T_ */
#endif /* _WIND_H_ */