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
void save_solution(EQNVARS *equationvars_p);

/* intode.c */
void integrate_ode(EQNVARS *equationvars_p);

/* relax.c */
void relax(EQNVARS *equationvars_p);

/* soe.c */
void init_soe(void);
double get_kappa_norm(const I_EQNVARS *vars);
double eval_eqn(int k, double *x, double **y, double (*ymod)[2], int eqnnum, int species);
void eval_interior_eqns(int k, double *x, double **y, double (*ymod)[2], double res[]);
void get_dvdr(double *dvdr, const I_EQNVARS *vars);
void linearize_dvdr_crit(double *x, double **y);
void get_drhodr(double *drhodr, const I_EQNVARS *vars, double dvdr);
void get_dTdr(double *dTdr, const I_EQNVARS *vars, double drhodr, double *dXdr,double k);
void get_dYsdr(double *dXdr, const I_EQNVARS *vars,double k, int print_rates);
void get_dNcoldr(double *dNdr, const I_EQNVARS *vars);
void get_spQ(double *spQ, const I_EQNVARS *vars,double k,int printout);
// void Energy_Conservation(int k, double *x, double **y, double (*ymod)[2], int species);
void get_rad(double *rad, const I_EQNVARS *vars);
// void print_cond_cool(EQNVARS *ev);

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
