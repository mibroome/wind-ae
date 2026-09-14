#if !defined(GLOBALS_H)
#define GLOBALS_H
/*============================================================================*
 *! \file globals.h                                                           *
 *  \brief Contains global variables. If main.c define globals, else declare. *
 *============================================================================*/
/* Runtime species count and derived equation counts.
 * Set in set_phys_params() from inputs/phys_params.inp before any solver call. */
#if defined(MAIN_C)
int g_nspecies = 0;  /* NSPECIES — read from phys_params.inp at startup */
int g_ne       = 0;  /* 5 + 2*g_nspecies */
int g_nb       = 0;  /* 3 + g_nspecies   */
int g_num_eqns = 0;  /* 3 + 2*g_nspecies */
int    g_m       = 0;     /* M        — read from tech_params.inp at startup */
int    g_itmax   = 100;   /* ITMAX    — read from tech_params.inp at startup */
double g_rhoscale= 10.0;  /* RHOSCALE — read from tech_params.inp at startup */
PARAMLIST parameters;
double x[M_MAX+1]; /* used for communication between relax and difeq */
double erf_norm;
#else
extern int g_nspecies, g_ne, g_nb, g_num_eqns;
extern int    g_m, g_itmax;
extern double g_rhoscale;
extern PARAMLIST parameters;
extern double x[M_MAX+1];
extern double erf_norm;
#endif /* MAIN_C */

/* Convenience aliases: use these throughout the code instead of hardcoded values.
 * They resolve to runtime ints, so they work in loops and expressions.
 * Do NOT use them as struct member array dimensions — use NSPECIES_MAX there. */
#ifndef NSPECIES
#define NSPECIES g_nspecies
#endif
#ifndef M
#define M        g_m
#endif
/* TOTALPTS as a runtime expression (safe in loops/conditions; use TOTALPTS_MAX for arrays) */
#ifndef TOTALPTS
#define TOTALPTS (INPTS+g_m+ADDPTS)
#endif
#ifndef NE
#define NE       g_ne
#endif
#ifndef NB
#define NB       g_nb
#endif
#ifndef NUM_EQNS
#define NUM_EQNS g_num_eqns
#endif
#endif /* GLOBALS_H */
