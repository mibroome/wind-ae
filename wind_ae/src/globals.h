#if !defined(GLOBALS_H)
#define GLOBALS_H
/*============================================================================*
 *! \file globals.h                                                           *
 *  \brief Contains global variables. If main.c define globals, else declare. *
 *============================================================================*/
#if defined(MAIN_C)
PARAMLIST parameters;
double x[M+1]; /* used for communication between relax and difeq */
#else
extern PARAMLIST parameters;
extern double x[M+1];
extern double erf_norm;
#endif /* MAIN_C */
#endif /* GLOBALS_H */
