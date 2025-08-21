#define MAIN_C
/*============================================================================*
 *! \file main.c                                                              *
 *  \brief Code for solving the ordinary differential equations that describe *
 *         atmospheres undergoing hydrodynamic escape driven by ionization    *
 *         heating.                                                           *
 *  \refs 1. Murray-Clay et al. 2009                                          *
 *        2. McCann & Murray-Clay 2020                                        *
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

int main(void){
  /* the equation variables */
  EQNVARS equationvars;

  /* load Gauss-Legendre quadrature rates */
  init_glq();
  
    /* read in the changeable parameters */
  set_parameters();

  /* read in the initial guess */
  initial_guess(&equationvars);
        
  /* solve from the base to the sonic point using a relaxation code */
  relax(&equationvars);

  /* integrate system of odes outside relaxation domain */
  integrate_ode(&equationvars);

  /* print the solution to file */
  save_solution(equationvars);

  /* free leftover allocated arrays */
  free_glq();

  return 0;
}
