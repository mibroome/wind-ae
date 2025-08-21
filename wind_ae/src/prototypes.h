#ifndef PROTOTYPES_H
#define PROTOTYPES_H

#include "nr_prototypes.h"

/* glq_rates */
void init_glq(void);
void free_glq(void);
double *glq_ionization(double *N, double *Ys, double rho, double k);
// double *glq_secondary_ionization(double *N, double *Ys, double rho, double k, int m);
double *glq_heating(double *N, double *Ys, double rho, double k);
double *ionization_potentials(void);
void save_spec(FILE *file);
void check_spec(char *hline);

/* io.c */
void set_parameters(void);
void initial_guess(EQNVARS *equationvars);
void save_solution(EQNVARS equationvars);

/* intode.c */
void integrate_ode(EQNVARS *equationvars_p);

/* relax.c */
void relax(EQNVARS *equationvars_p);

/* soe.c */
double eval_eqn(int k, double *x, double **y, double ymod[NE+1][2], int eqnnum, int species);
void get_dvdr(double *dvdr, I_EQNVARS vars);
void linearize_dvdr_crit(double *x, double **y);
void get_drhodr(double *drhodr, I_EQNVARS vars, double dvdr);
void get_dTdr(double *dTdr, I_EQNVARS vars, double drhodr, double dXdr[NSPECIES],double k);
void get_dYsdr(double *dXdr, I_EQNVARS vars,double k, int print_rates);
void get_dNcoldr(double *dNdr, I_EQNVARS vars);
void get_spQ(double *spQ, I_EQNVARS vars,double k,int printout);
// void Energy_Conservation(int k, double *x, double **y, double ymod[NE+1][2], int species);
void get_rad(double *rad, I_EQNVARS vars);

/* utils.c */
void *calloc_1d_array(size_t nc, size_t size);
void free_1d_array(void *array);
double *calloc_1d_array_gross(long is, long ie);
double **calloc_2d_array_gross(long is, long ie, long js, long je);
double ***calloc_3d_array_gross(long is, long ie, long js, long je, long ks,
                                long ke);
void free_1d_array_gross(double *array, long is);
void free_2d_array_gross(double **array, long is, long js);
void free_3d_array_gross(double ***array, long is, long js, long ks);

#endif /* PROTOTYPES_H */
