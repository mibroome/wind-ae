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

  /* read in the changeable parameters — MUST come before alloc so g_nspecies/g_m are set */
  set_parameters();

  /* Allocate Ys and Ncol in EQNVARS now that g_nspecies and g_m are known.
   * Sized exactly to the runtime species count to avoid NSPECIES_MAX overhead. */
  {
    int totalpts = INPTS + g_m + ADDPTS;
    equationvars.Ys   = calloc_2d_array_gross(0, totalpts-1, 0, g_nspecies-1);
    equationvars.Ncol = calloc_2d_array_gross(0, totalpts-1, 0, g_nspecies-1);
    if (!equationvars.Ys || !equationvars.Ncol) {
      fprintf(stderr, "ERROR: Failed to allocate EQNVARS.Ys/Ncol\n");
      exit(1);
    }
  }

  /* precompute alpharec table row indices and line-cooling species indices */
  init_soe();

  /* read in the initial guess */
  initial_guess(&equationvars);

  /* solve from the base to the sonic point using a relaxation code */
  relax(&equationvars);

  /* integrate system of odes outside relaxation domain */
  integrate_ode(&equationvars);

  /* print the solution to file */
  save_solution(&equationvars);

  /* free Ys and Ncol heap arrays */
  free_2d_array_gross(equationvars.Ys,   0, 0);
  free_2d_array_gross(equationvars.Ncol, 0, 0);

  /* free leftover allocated arrays */
  free_glq();

  return 0;
}
