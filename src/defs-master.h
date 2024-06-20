#ifndef DEFINITIONS_H
#define DEFINITIONS_H
/*============================================================================*
 *! \file defs.h                                                              *
 *  \brief Contains preprocessor directives defintions used throughout        *
 *============================================================================*/

#include <math.h>

#define ON 1
#define OFF 0

#define SQR(x) ((x)*(x))

/* Relaxation code parameters */
#define NSPECIES 1
#define NE (4+2*NSPECIES)         /* number of equations */
#define M 1501            /* number of points */
#define NB (2+NSPECIES)          /* number of boundary conditions at first point */
#define ITMAX 100     /* max number of iterations */
#define CONV 1.0E-10  /* convergence criterion */
#define SLOWC 1.0E-2  /* slow down when far from convergence */
#define N_ADD_PARAMS 0  /* Number of additional parameters. Allows future users to easily add params */
// #define NADD 0.0  /*Number of additional parameters Allows future users to easily add params */

/* odeint continuation parameters */
#define RIN 0.9     /* migrate this parameter in bc inputs */
#define INPTS 0  /* number of extra points inwards of Rmin */
#define ADDPTS 1000 /* number of extra points past the relaxation range */
#define NUM_EQNS (3+2*NSPECIES) /*formerly 5 in single species version*/
#define ACCURACY 1.0E-13
#define STEP_GUESS 1.0E-1
#define STEP_MIN 0.0
#define MAX_STORED_STEPS 10

/* total number of grid points */
#define TOTALPTS (INPTS+M+ADDPTS)

/* Constants */
#define PI 3.141592653589793

/* Physical constants */
#define MH 1.6733E-24             /* g */
#define K 1.380658E-16            /* erg/K */
#define SIG_SB 5.6705E-5          /* erg/cm^2 */
#define BAR 1.0E6                 /* dyne/cm^2 */
#define G 6.67259E-8              /* cm^3/g/s^2 */
#define H 6.6260755E-27           /* erg s */
#define EV 1.6021772E-12          /* erg */
#define MHEOVERMH 4.0             /* mass of He over H */
#define MHOVERMHE (1.0/MHEOVERMH)
#define AU 1.496E13               /* cm */
#define MSUN 1.99E33              /* g */

/* Scale values for non-dimensionalizing */
#define T0 1.0E4           /* K approx wind temp */
#define RHO0 1.0E-15       /* ~nBAR*MH/K/T0 */
#define CS0 sqrt(K*T0/MH)
#define NCOL0 1.0E+17      /* tau~few at ionization edge */

/* Coefficients for equations */
#define ALPHAREC_B_COEFF 2.59E-13       /* cm^3/s */
#define ALPHAREC_COEFF ALPHAREC_B_COEFF
#define BETA_B_COEFF 1.73E-13           /* cm^3/s */
#define BETA_A_COEFF 3.23E-13           /* cm^3/s */
#define GAMMA_ATOMIC (5.0/3.0)
#define LYACOOL_COEFF -7.5E-19          /* erg/cm^3/s */
#define LYACOOL_TEMP (118348.0/T0)      /* units of T0 */
#define RECCOOL_COEFF 8.44E-26          /* erg/cm^3/s/K^0.11 */

#define RHOSCALE 100.0
#define VSCALE 2.0                      /* units of CS0 */
#define TEMPSCALE 1.0                   /* units of T0 */
#define NCOLSCALE 100.0
#define FPSCALE 1.0
#define ZSCALE 3.0                      /* z = rs-rmin; units of RP */

#define PLANET_PARAM_FILE "inputs/planet_params.inp"
#define BC_PARAM_FILE "inputs/bcs.inp"
#define TERM_PARAM_FILE "inputs/term_ind.inp"
#define TECH_PARAM_FILE "inputs/tech_params.inp"
#define PHYS_PARAM_FILE "inputs/phys_params.inp"
#define SPECTRUM_FILE "inputs/spectrum.inp"
#define GUESS_FILE "inputs/guess.inp"
#define ADD_PARAM_FILE "inputs/additional_params.inp"


#define HEADER_DIFF_VERBOSE ON
#define GET_PARAMS_VERBOSE ON

#define NOCLOBBER OFF
#define SOLNFILE "saves/windsoln.csv"

#define DERIVDIV 1.0E7 /* divisor for derivative steps */

/* Give each equation a number */
#define VEQN 1
#define RHOEQN 2
#define IONEQN 3
#define NCOLEQN 4
#define TEQN 5
#define SPVEQN 6
#define SPCRITEQN 7

/* And each component as well */
#define QCOMP 100

#endif /* DEFINITIONS_H */