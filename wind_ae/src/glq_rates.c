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
#include "rate_coeffs.h"


/* static global variables */
static char *glq_table_header, *date, *kind, *src_file;
static int npts, nspecies;
static double *ion_pot;
static double wndw_lob, wndw_upb, rslv_lob, rslv_upb, norm_lob, norm_upb;
static double *Phi, *ionization_rate,  *heating_rate;
static double *last_N, *hc_over_wl, *tau_array, *f_denom, *wPhi_wl, **sigma_wl;

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
      ion_pot  = (double*)calloc_1d_array(nspecies, sizeof(double));
      last_N   = (double*)calloc_1d_array(nspecies, sizeof(double));
      ionization_rate = (double*)calloc_1d_array(nspecies+nspecies*nspecies, sizeof(double));
      // secondary_ionization_rate = (double*)calloc_1d_array(nspecies, sizeof(double));
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
  //     fgets(line, sizeof(line), filep);
  // sscanf(line, "%s %lf%n %[^,]", temp, &(params->atomic_mass[0]), &n, temp);
  // for (i = 1; i<NSPECIES; i=i+1){
  //     sscanf(line+n+1+i, "%lf%n %[^,]", &(params->atomic_mass[i]), &m, temp);
  //     n=n+m;
  // }
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
  f_denom  = (double*)calloc_1d_array(npts, sizeof(double));
  wPhi_wl = (double*)calloc_1d_array(npts, sizeof(double));
  Phi = (double*)calloc_1d_array(npts, sizeof(double));

  nheader--;
  (void)read;
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
  free_1d_array(tau_array);
  free_1d_array(f_denom);
  free_1d_array(Phi);
  free_1d_array(ionization_rate);
  free_1d_array(heating_rate);
  for (i = 0; i < npts; i++) {
    free_1d_array(sigma_wl[i]);
  }
  free_1d_array(sigma_wl);
  free(glq_table_header);
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
  calc_gql_rates(N, Ys, rho, k);

  return heating_rate;
}

double *ionization_potentials(void) { 
  return ion_pot;
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
  int i, j, jj, m, recalc = 0;
  double sigma, tau, f, eta_m, denom;
  double background_ioniz_frac, E_0;
  double frac_in_ion_tot=0,frac_in_ion_m=0, frac_in_heat, frac_in_excite=0;
  double Ncol, primary_ion_rate_j,n_ion_tot=0, n_tot=0, n_H=0, n_HII=0, n_j=0, n0_m=0;
  // int e, e_index, r, r_index;
  double R_tot;//,frac_in_ion_H,frac_in_ion_He;


  /* Check if prior rates can be reused */
  for (j = 0; j < nspecies; j++) {
    if (N[j] != last_N[j]) {
      recalc++;
      break;
    }
  }

  // (void)n_tot;

   /* Background ionization fraction in keeping with (Shull & van Steenberg 1985)
    * They define background ionization fraction as n_HII/n_H, even for an H-He atmosphere. 
    * Coulomb collisions with thermal electrons heat and collision excitation/ionization of neutrals ionize,
    * so, technically, this should be n_HII/n_tot, but, for H dominated atmos, the two yield very similar results
    * To keep with the fitting of Shull & Steenberg (1985), we elect to use n_HII/n_H...
    * ***HOWEVER, this is likely change for high Z atmospheres.***
    */
    
  /* If we do not need to update rates return. else, calculate updated rates */
  if (recalc == 0) {    
    return;
  }
  else {
    /* Need total column density across species to calculate background ionization fraction */
    // for (e=0;e<100;e++){ R_tot[e] = 0;}
    for (j = 0; j < nspecies; j++) {
      n_j = parameters.HX[j]/parameters.atomic_mass[j]; //over rho
      n_tot += n_j;
      n_ion_tot += n_j*(1-Ys[j]);

      /*Calculating the number density weighted sum of the rate coefficients for all species j*/
      /*Rate coefficients from Dere (2007) Table 29*/
      // for (r=0;r<465;r++){
      //   if (R[r][0]==parameters.Z[j] && R[r][1]==parameters.N_e[j]){
      //     for (e=0;e<100;e++){
      //       R_tot[e] += R[r][e+2]*n_j*Ys[j]; //times neutral number density
      //     }
      //   }
      // }
    }
    n_H = parameters.HX[0]/parameters.atomic_mass[0]; //over rho
    n_HII = (1-Ys[0])*n_HII; 
    (void)Ncol;
    (void)n_H;
    //  background_ioniz_frac = n_HII/n_H; 
    /* Calculate Phi */
    for (i = 0; i < npts; i++) {
      tau = 0.;
      denom = 0.;
      for (j = 0; j < nspecies; j++) {
        tau += sigma_wl[i][j]*N[j]*NCOL0;
        denom += sigma_wl[i][j]*parameters.HX[j]*Ys[j]/parameters.atomic_mass[j]; //over rho
      }
      tau_array[i] = tau;
      f_denom[i] = denom; //populating denominator of on-the-spot opacity fraction
      Phi[i] = wPhi_wl[i]*exp(-tau);
    }
      
    /* Calculate ionization and heating rate */
    for (j = 0; j < nspecies; j++) {
      ionization_rate[j*(nspecies+1)] = 0.; //primary ionization rate
      for (m=0;m<nspecies;m++){
        /*secondary ionization rate induced in species m by phototelectrons released by ionizing species j*/
        ionization_rate[j*(nspecies+1)+1+m] = 0.; 
      }
      heating_rate[j] = 0.;

      for (i = 0; i < npts; i++) {
        background_ioniz_frac = n_ion_tot/n_tot;
        frac_in_excite = (0.4766 * pow(1-pow(background_ioniz_frac,0.2735),1.5221));
        frac_in_heat = 0.9971 * (1 - pow(1-pow(background_ioniz_frac,0.2663),1.3163)); 
        frac_in_ion_tot = 1-frac_in_heat-frac_in_excite;
        sigma = sigma_wl[i][j];       
        f = sigma*(parameters.HX[j]*Ys[j]/parameters.atomic_mass[j]) / f_denom[i];
        //f is fraction of energy into each species per Osterbrock Ch 2.4: f = n0[j]*sigma[j] / sum,j(n0[j]*sigma[j])
        E_0 = hc_over_wl[i] - ion_pot[j]; //Energy carried by primary photoelectron
        if (sigma == 0) {
          break;
        }
      /*Adapted from Mocassin (Shull & Steenberg 1985)*/
       //Their solution for fraction into heat and ionization heating as a function of background ionization frac
       //tends to that solution as E_0 approaches 100 eV. However, there should also be leftover energy for secondary e-
       //with E_0>E_ionization (e.g., E_0>13.6 for pure H), so, while the 100eV function understimates the fraction
       //deposited in heat for lower E_0 (see Fig. 3 (Shull & Steenberg 1985)), it is still closer than assuming all
       //energy is deposited in heat, which would overestimate the heating. So, we select an artificial boundary of 40 eV
        // if(k!=1e4){printf("k=%.0f %s %.2feV %.6f\n", k ,parameters.species[j],E_0*6.242e+11,frac_in_heat);}
        primary_ion_rate_j = parameters.Ftot*sigma*f*Phi[i];
        if (E_0 > 6.408707e-11){ //Lower limit of x-ray energy range in ergs
          if (background_ioniz_frac < 0.000){ background_ioniz_frac = 1e-10;} //b/c sometimes Ys>1
          /*fraction of energy of primary electron deposited in heat*/

          /*m is the species from whence the electron is emitted to ionize species j*/
          
          for (m=0;m<nspecies;m++){          
            R_tot = 0;
            for (jj=0; jj<nspecies; jj++){
              R_tot += R[nspecies*jj+m][i]*Ys[jj]*parameters.HX[jj]/parameters.atomic_mass[jj];
            }
            n0_m = Ys[m]*parameters.HX[m]/parameters.atomic_mass[m]; //over rho
            frac_in_ion_m = R[nspecies*j+m][i]*n0_m/R_tot;

            // frac_in_ion_m = R[r_index][e_index+2]*n0_m / R_tot[e_index];
            eta_m = E_0*(frac_in_ion_tot*frac_in_ion_m)/ion_pot[m];
            // if (k==200){printf("Ion of species %s by a %.1f photon yields %.1f %s ionizations\n",parameters.species[j],hc_over_wl[i]*6.24e+11,eta_m,parameters.species[m]);}

            //for unknown numerical reasons, have to multiply by n0_j in soe.c
            ionization_rate[j*(nspecies+1)+1+m] += eta_m*primary_ion_rate_j;
            
          }
          /*Analytically we are doing: sum_E(n[j]*primary_ion[j]) + sum_m(sum_E(n[m]*secondary_ion[m]*eta[j])) *
          * Which is equivalent to: n[j]*sum_E(primary_ion[j]) + sum_mm(n[mm]*sum_E(secondary_ion[mm]*eta[j])) *
          * so, confusingly, when we call for mm in range(nspecies): glq_secondary_ionization(j)[mm]           */
        }  
        else{ frac_in_heat = 1.0; }

        /*Yes, it would be lovely and sensible to multiply by n0_j here*/
        /*but for unknown and unknowable numerical reasons that does not work*/
        ionization_rate[j*(nspecies+1)] += primary_ion_rate_j; //initial ionziation (secondary calculated above)
        // if(k!=1e4){printf("k=%.0f %s %.2feV %.6f\n", k ,parameters.species[j],E_0*6.242e+11,frac_in_heat);}
        // printf("%s ",)

        heating_rate[j] += frac_in_heat*E_0*sigma*f*Phi[i]*parameters.Ftot; 
        // if (k==10){
        //   FILE *fp = fopen("outputs/ion_rates.txt", "a"); // Open for appending
        //   if (fp == NULL) {
        //       perror("Failed to open ion_rates.txt");
        //       exit(1);
        //   }
        //   fprintf(fp, "%.0f,%d,%.2f,%.6f,%.2e,%.6f,%f\n",k,j,E_0*6.242e+11,frac_in_heat,sigma,f,Phi[i]);
        //   fclose(fp);
        // }
      }
    }
  }
    

  /* Save state used for rates for next check */
  for (j = 0; j < nspecies; j++) {
    last_N[j] = N[j];
  }

  return;
    
}