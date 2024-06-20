/*============================================================================*
 *! \file glq_rates.c (Gauss-Legendre quadrature ionizing rates)              *
 *  \brief Handles all that has to do with the multifrequency ionization.     *
 *         Does the full integration over the spectrum via GLQ                *
 *============================================================================*/

/* Standard C */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
/* My Headers */
#include "defs.h"
#include "wind.h"
#include "prototypes.h"
#include "globals.h"

/* static global variables */
static char *glq_table_header, *date, *kind, *src_file;
static int npts, nspecies;
static double *ion_pot;
static double wndw_lob, wndw_upb, rslv_lob, rslv_upb, norm_lob, norm_upb;
static double *Phi, *ionization_rate, *heating_rate;
static double *last_N, *hc_over_wl, *tau_array, *wPhi_wl, **sigma_wl;

/*----------------------------------------------------------------------------*
 *======================== PRIVATE FUNCTION PROTOTYPES =======================*
 *----------------------------------------------------------------------------*/
static void calc_gql_rates(double *N, double *Ys, double rho, double k);

/*============================================================================*
 *=========================== PUBLIC FUNCTIONS ===============================*
 *============================================================================*
 * init_glq       - initalize the gauss-legendre quadrature spectrum          *
 * free_glq       - free heap memory allocated in the glq scheme              *
 * glq_ionization - calculate multifrequency ionization rate                  *
 * glq_heating    - calculate multifrequency heating rate                     *
 * save_spec      - saves the glq table to the output header                  *
 * check_spec     - checks if spectrum parameters are different than guess's  *
 *----------------------------------------------------------------------------*/

/*============================================================================*
 *! \fn void init_glq(void)                                                   *
 *  \brief Initalize all the variables needed for the glq multifrequency algo *
 *----------------------------------------------------------------------------*/

#define BUFSIZE 1024
void init_glq(void) {
  int i, j;
  char *line, *tok;
  FILE *filep;
  char *filename, dline[BUFSIZE];
  long start;
  size_t len, nheader;
  ssize_t read;

  printf("Loading Gauss-Legendre quadrature data...\n");
  filename = SPECTRUM_FILE;

  /* open the guess file for reading */
  if ((filep = fopen(filename, "r")) == NULL) {
    fprintf(stderr, "Error: Cannot open %s\n", filename);
    exit(201);
  }

  /* Read header information */
  line = NULL;
  len  = 0;
  nheader = 0;
  do {
    start = ftell(filep);
    read  = getline(&line, &len, filep);
    nheader++;
    tok   = strtok(line, ":");
    if (strstr(line, "NPTS") != NULL) {
      tok  = strtok(NULL, ",");
      npts = atoi(tok);
    }
    else if (strstr(line, "NSPECIES") != NULL) {
      tok = strtok(NULL, ",");
      nspecies = atoi(tok);
//       printf("nspecies = %d", nspecies);
      ion_pot  = (double*)calloc_1d_array(nspecies, sizeof(double));
      last_N   = (double*)calloc_1d_array(nspecies, sizeof(double));
      ionization_rate = (double*)calloc_1d_array(nspecies, sizeof(double));
      heating_rate = (double*)calloc_1d_array(nspecies, sizeof(double));
    }
    else if (strstr(line, "DATE") != NULL) {
      tok  = strtok(NULL, ",");
      date = malloc((strlen(tok)+1)*sizeof(char));
      strcpy(date, tok);
      i = 0;
      j = 0;
      while (date[i]) {
        if (date[i] != ' ' && date[i] != '\n') {
          date[j++] = date[i];
        }
        i++;
      }
      date[j++] = '\0';
      date = (char*)realloc(date, j);
    }
    else if (strstr(line, "FILE") != NULL) {
      tok  = strtok(NULL, ",");
//         printf('--------%s',tok);
      src_file = malloc((strlen(tok)+1)*sizeof(char));
      strcpy(src_file, tok);
      i = 0;
      j = 0;
      while (src_file[i]) {
        if (src_file[i] != ' ' && src_file[i] != '\n') {
          src_file[j++] = src_file[i];
        }
        i++;
      }
      src_file[j++] = '\0';
      src_file = (char*)realloc(src_file, j);
//       printf('%s',src_file);
    }
    else if (strstr(line, "KIND") != NULL) {
      tok  = strtok(NULL, ",");
      kind = malloc((strlen(tok)+1)*sizeof(char));
      strcpy(kind, tok);
      i = 0;
      j = 0;
      while (kind[i]) {
        if (kind[i] != ' ' && kind[i] != '\n') {
          kind[j++] = kind[i];
        }
        i++;
      }
      kind[j++] = '\0';
      kind = (char*)realloc(kind, j);
    }
    else if (strstr(line, "WINDOW") != NULL) {
      tok = strtok(NULL, ",");
      wndw_lob = atof(tok);
      tok = strtok(NULL, ",");
      wndw_upb = atof(tok);
    }
    else if (strstr(line, "RESOLVED") != NULL) {
      tok = strtok(NULL, ",");
      rslv_lob = atof(tok);
      tok = strtok(NULL, ",");
      rslv_upb = atof(tok);
    }
    else if (strstr(line, "NORMALIZED") != NULL) {
      tok = strtok(NULL, ",");
      norm_lob = atof(tok);
      tok = strtok(NULL, ",");
      norm_upb = atof(tok);
    }
    else if (strstr(line, "IONPOTS") != NULL) {
      for (i = 0; i < nspecies; i++) {
        tok = strtok(NULL, ",");
        ion_pot[i] = atof(tok);
      }
    }
    else if ((line[0]) == '#') {
      /* Assuming only remaining (or last) comment is tabler header */
      glq_table_header = malloc((strlen(line)+1)*sizeof(char));
      strcpy(glq_table_header, line);
    }
  } while ((line[0]) == '#');
  sigma_wl = (double**)calloc_1d_array(npts, sizeof(double*));
  for (i = 0; i < npts; i++) {
    sigma_wl[i] = (double*)calloc_1d_array(nspecies, sizeof(double));
    }
  hc_over_wl = (double*)calloc_1d_array(npts, sizeof(double));
  tau_array  = (double*)calloc_1d_array(npts, sizeof(double));
  wPhi_wl = (double*)calloc_1d_array(npts, sizeof(double));
  Phi = (double*)calloc_1d_array(npts, sizeof(double));

  nheader--;
  printf("  Skipping %zu lines of header in %s\n", nheader, filename);
  fseek(filep, start, SEEK_SET);

  /* Read in the glq weights and variables evaluated at the abscissas */
  /* Each line of the input file has the form: */
  /* hc_over_wl, wPhi_wl, sigma_wl (for each species)... */
  i = 0;
  while (fgets(dline, BUFSIZE, filep) && i < npts) {
    tok = strtok(dline, ",");
    hc_over_wl[i] = atof(tok);
    tok = strtok(NULL, ",");
    wPhi_wl[i] = atof(tok);
    for (j = 0; j < nspecies; j++) {
      tok = strtok(NULL, ",");
      sigma_wl[i][j] = atof(tok);
    }
    i++;
  }

  /* close the file */
  fclose(filep);
  printf("  npts:%d, nspecies:%d\n", npts, nspecies);
  printf("Successfully loaded Gauss-Legendre quadrature data.\n");

  return;
}
#undef BUFSIZE

/*============================================================================*
 *! \fn void free_glq(void)                                                   *
 *  \brief Free all heap memory allocated by the glq scheme                   *
 *----------------------------------------------------------------------------*/

void free_glq(void) {
  int i;

  free_1d_array(ion_pot);
  free_1d_array(last_N);
  free_1d_array(hc_over_wl);
  free_1d_array(wPhi_wl);
  free_1d_array(Phi);
  free_1d_array(ionization_rate);
  free_1d_array(heating_rate);
  for (i = 0; i < npts; i++) {
    free_1d_array(sigma_wl[i]);
  }
  free_1d_array(sigma_wl);
  free(glq_table_header);
//   free(date);
  free(kind);
  free(src_file);

  return;
}

/*============================================================================*
 *! \fn double glq_ionization(double *N, double *Ys, double rho)                                      *
 *  \brief Given the column densities of each species, calculates the         *
 *         ionization rate of each species due to the multifrequency spectrum *
 *----------------------------------------------------------------------------*/

double *glq_ionization(double *N, double *Ys, double rho, double k) {
  calc_gql_rates(N, Ys, rho, k);

  return ionization_rate;
}

/*============================================================================*
 *! \fn double glq_heating(double *N, double *Ys, double rho)                                         *
 *  \brief Given the column densities of each species, calculates the total   *
 *         heating rate due to the multifrequency spectrum                    *
 *----------------------------------------------------------------------------*/

double *glq_heating(double *N, double *Ys, double rho, double k) {
//     FILE *fptr;
  calc_gql_rates(N, Ys, rho, k);
//     fptr = fopen("heating_rate.txt","a+");
//   fprintf(fptr,"%.0f %.2e\n",k,*heating_rate);
//   fclose(fptr);  
  return heating_rate;
}

/*============================================================================*
 *! \fn void save_spec(FILE *file)                                            *
 *  \brief Saves the spectrum data table to the header when writing the save  *
 *         file.                                                              *
 *----------------------------------------------------------------------------*/

void save_spec(FILE *file) {
  int i, j;
  char *specline = "%d,%d,%s,%s,%s,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e";

  /* Print spectrum parameters */
  fprintf(file, "#spec_prms: ");
  fprintf(file, specline, npts, nspecies, date, src_file, kind,
          wndw_lob, wndw_upb, rslv_lob, rslv_upb, norm_lob, norm_upb);
  for (j = 0; j < nspecies; j++) {
    fprintf(file, ",%.17e", ion_pot[j]);
  }
  fprintf(file, "\n");

  /* Print GLQ table */
  fprintf(file, "%s", glq_table_header);
  for (i = 0; i < npts; i++) {
    fprintf(file, "## %.17e,%.17e", hc_over_wl[i], wPhi_wl[i]);
    for (j = 0; j < nspecies; j++) {
      fprintf(file, ",%.17e", sigma_wl[i][j]);
    }
    fprintf(file, "\n");
  }

  return;
}

/*============================================================================*
 *! \fn void check_spec(char *hline)                                          *
 *  \brief Checks the header of the guess.inp to see if the spectrum has been *
 *         altered. Prints a warning to let the user know of any changes.     *
 *----------------------------------------------------------------------------*/

void check_spec(char *hline) {
  int i, g_nspecies;
  char *tok;

  tok = strtok(hline, ":");
  if (atoi(tok = strtok(NULL, ",")) != npts) {
    fprintf(stdout, "WARNING: npts differs from guess's:\n"
            "         guess:%d, run:%d\n", atoi(tok), npts);
  }
  if (atoi(tok = strtok(NULL, ",")) != nspecies) {
    fprintf(stdout, "WARNING: nspecies differs from guess's:\n"
            "         guess:%d, run:%d\n", atoi(tok), nspecies);
  }

//     printf("NSPECIES spec_check done");
  g_nspecies = atoi(tok);
  tok = strtok(NULL, ",");
  if (strcmp(tok, date) != 0) {
    fprintf(stdout, "WARNING: date differs from guess's:\n"
            "         guess:%s, run:%s\n", tok, date);
  }
    
  tok = strtok(NULL, ",");
  if (strcmp(tok, src_file) != 0) {
    fprintf(stdout, "WARNING: Spectrum source file differs from guess's:\n"
            "         guess:%s, run:%s\n", tok, src_file);
  }
  tok = strtok(NULL, ",");
  if (strcmp(tok, kind) != 0) {
    fprintf(stdout, "WARNING: kind differs from guess's:\n"
            "         guess:%s, run:%s\n", tok, kind);
  }
  if (atof(tok = strtok(NULL, ",")) != wndw_lob) {
    fprintf(stdout, "WARNING: wndw_lob differs from guess's:\n"
            "         guess:%e, run:%e\n", atof(tok), wndw_lob);
  }
  if (atof(tok = strtok(NULL, ",")) != wndw_upb) {
    fprintf(stdout, "WARNING: window_upb differs from guess's:\n"
            "         guess:%e, run:%e\n", atof(tok), wndw_upb);
  }
  if (atof(tok = strtok(NULL, ",")) != rslv_lob) {
    fprintf(stdout, "WARNING: rslv_lob differs from guess's:\n"
            "         guess:%e, run:%e\n", atof(tok), rslv_lob);
  }
  if (atof(tok = strtok(NULL, ",")) != rslv_upb) {
    fprintf(stdout, "WARNING: rslv_upb differs from guess's:\n"
            "         guess:%e, run:%e\n", atof(tok), rslv_upb);
  }
  if (atof(tok = strtok(NULL, ",")) != norm_lob) {
    fprintf(stdout, "WARNING: norm_lob differs from guess's:\n"
            "         guess:%e, run:%e\n", atof(tok), norm_lob);
  }
  if (atof(tok = strtok(NULL, ",")) != norm_upb) {
    fprintf(stdout, "WARNING: norm_upb differs from guess's:\n"
            "         guess:%e, run:%e\n", atof(tok), norm_upb);
  }
  if (g_nspecies == nspecies) {
    for (i = 0; i < nspecies; i++) {
      if (atof(tok = strtok(NULL, ",")) != ion_pot[i]) {
        fprintf(stdout, "WARNING: IONPOTS[%d] differs from guess's:\n"
                "         guess:%e, run:%e\n", i, atof(tok), ion_pot[i]);
      }
    }
  }

  return;
}

/*============================================================================*
 *=========================== PRIVATE FUNCTIONS ==============================*
 *============================================================================*
 * calc_gql_rates - Calculates the heating and ionization rates               *
 *----------------------------------------------------------------------------*/

/*============================================================================*
 *! \fn static void calc_gql_rates(double *N)                                 *
 *  \brief Given the column densities of each species, does the calculation   *
 *         of the ionizationg and heating rates by first calculating the      *
 *         attenuated spectrum at the given column densities                  *
 *----------------------------------------------------------------------------*/

static void calc_gql_rates(double *N, double *Ys, double rho, double k) {
  int i, j, m, recalc = 0;
  double sigma, sigma_m, tau, eta=0, secondary_H_ion_rate;
  double background_ioniz_frac, E_0, frac_in_ion, frac_in_heat=1, n_tot=0, n_H=0, n_HII=0;
  double C, b, a, f;
  double Ncol;

  /* Check if prior rates can be reused */
  for (j = 0; j < nspecies; j++) {
    if (N[j] != last_N[j]) {
      recalc++;
      break;
    }
  }
   
  /* Need total column density across species to calculate background ionization fraction */
  for (j = 0; j < nspecies; j++) {
         n_tot += (RHO0/parameters.atomic_mass[j])*parameters.HX[j]*rho;
  }  
   n_HII = (1-Ys[0])*(RHO0/parameters.atomic_mass[0])*parameters.HX[0]*rho;
   n_H = (RHO0/parameters.atomic_mass[0])*parameters.HX[0]*rho;
  background_ioniz_frac = n_HII/n_H; //background ionization fraction in keeping with (Shull & Steenberg 1985)
    //They define background ionization fraction as n_HII/n_H, even for an H-He atmosphere. 
    //Coulomb collisions with thermal electrons heat and collision excitation/ionization of neutrals ionize,
    //so, technically, this should be n_HII/n_tot, but, for H dominated atmos, the two yield very similar results
    //To keep with the fitting of Shull & Steenberg (1985), we elect to use n_HII/n_H...
    //***HOWEVER, this is likely change for high Z atmospheres.***
    
  /* If we do not need to update rates return. else, calculate updated rates */
  if (recalc == 0) {
    return;
  }
  else {
    /* Calculate Phi */
    for (i = 0; i < npts; i++) {
      tau = 0.;
      for (j = 0; j < nspecies; j++) {
        tau += sigma_wl[i][j]*N[j]*NCOL0;
      }
      tau_array[i] = tau;
      Phi[i] = wPhi_wl[i]*exp(-tau);
    }
      
    /* Calculate ionization and heating rate */
    for (j = 0; j < nspecies; j++) {
      ionization_rate[j] = 0.;
      heating_rate[j] = 0.;
      for (i = 0; i < npts; i++) {
        sigma = sigma_wl[i][j];
        if (N[j] < 0){Ncol = 0;} //temporarily enforce this because going negative 
        else {Ncol = N[j];}          
        f = sigma*Ncol*NCOL0 / tau_array[i];
        //f is fraction of energy into each species per Osterbrock Ch 2.4: f = Ncol[j]*sigma[j] / sum,j(Ncol[j]*sigma[j])
        E_0 = hc_over_wl[i] - ion_pot[j]; //Energy carried by primary electron
//          printf("%.2e\n",E_0);
//         if (E_0 < 0){ E_0 = 0; } //this should be take care of by the sigma condition below, but just being safe
        if (sigma == 0) {
          break;
        }
      /*Adapted from Mocassin (Shull & Steenberg 1985)*/
       //Their solution for fraction into heat and ionization heating as a function of background ionization frac
       //tends to that solution as E_0 approaches 100 eV. However, there should also be leftover energy for secondary e-
       //with E_0>E_ionization (e.g., E_0>13.6 for pure H), so, while the 100eV function understimates the fraction
       //deposited in heat for lower E_0 (see Fig. 3 (Shull & Steenberg 1985)), it is still closer than assuming all
       //energy is deposited in heat, which would overestimate the heating. So, we select an artificial boundary of 40 eV
        if (E_0 > 6.408707e-11){ //Lower limit of x-ray energy range in ergs
//             printf("%.4f\n",background_ioniz_frac);
            if (background_ioniz_frac < 0.000){ background_ioniz_frac = 1e-10;} //sus
            /*fraction of energy of primary electron deposited in heat*/
            frac_in_heat = 0.9971 * (1 - pow(1-pow(background_ioniz_frac,0.2663),1.3163));            
            for (m=0; m < nspecies; m++){
                sigma_m = sigma_wl[i][m]; 
                if (strcmp(parameters.species[m],"HI") == 0){
                    C=0.3908; a=0.4092; b=1.7592;
                    frac_in_ion = C * pow(1-pow(background_ioniz_frac,a),b);
                    eta = (frac_in_ion*E_0/ion_pot[m]); //# of secondary ionizations (as a func of background ioniz. frac)
                    secondary_H_ion_rate = eta*sigma_m*Phi[i]*f;
                    ionization_rate[m] += secondary_H_ion_rate;  //secondary ionization rate         
                }
                else if (strcmp(parameters.species[m],"HeI") == 0){
                    C=0.0554; a=0.4614; b=1.6660;
                    frac_in_ion = C * pow(1-pow(background_ioniz_frac,a),b);
                    eta = (frac_in_ion*E_0/ion_pot[m]); 
                    ionization_rate[m] += eta*sigma_m*f*Phi[i];
                }
                else { 
                    /*To lowest order, can estimate metal ioniz. rate as H secondary ioniz. rate x sigma[metal]/sigma[H]*/
                    ionization_rate[m] += secondary_H_ion_rate*sigma_m/sigma_wl[i][0];
                }
            }
        }
        else{ frac_in_heat = 1.0; }
//           if (k==200){if (j==1){printf("k=%.1f   %.1f   %.4f\n",k,E_0/1.602177e-12,frac_in_heat);}}
//           printf("k=%.1f   %.1f   %.4f\n",k,E_0/1.602177e-12,frac_in_heat);
        ionization_rate[j] += sigma*f*Phi[i]; //initial ionziation (secondary calculated above)
        heating_rate[j] += frac_in_heat*E_0*sigma*f*Phi[i]; 
      }
    }
//    fclose(f);
  }
    

  /* Save state used for rates for next check */
  for (j = 0; j < nspecies; j++) {
    last_N[j] = N[j];
  }

  return;
    
}