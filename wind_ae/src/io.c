/*============================================================================*
 *! \file io.c (inputs/outputs)                                               *
 *  \brief Handles the input and output of the model                          *
 *============================================================================*/

/* Standard C */
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
/* For mkdir */
#include <sys/stat.h>
#include <sys/types.h>
#include <errno.h>
/* My Headers */
#include "defs.h"
#include "wind.h"
#include "globals.h"
#include "prototypes.h"

/*----------------------------------------------------------------------------*
 *======================== PRIVATE FUNCTION PROTOTYPES =======================*
 *----------------------------------------------------------------------------*/
/*============================== INPUT FUNCTIONS =============================*/
static int mkpath(char* path, mode_t mode);
static void set_bcs(PARAMLIST *params);
static void set_term_ind(PARAMLIST *params);
static void set_tech_params(PARAMLIST *params);
static void set_phys_params(PARAMLIST *params);
static void set_plnt_params(PARAMLIST *params);
/*============================= OUTPUT FUNCTIONS =============================*/
static void print_params(void);

/*============================================================================*
 *=========================== PUBLIC FUNCTIONS ===============================*
 *============================================================================*
 * set_parameters - reads in the parameters                                   *
 * initial_guess  - reads in the inital guess                                 *
 *----------------------------------------------------------------------------*/

/*============================================================================*
 *! \fn void set_parameters(void)                                                 *
 *  \brief Reads in and sets the parameters for the problem                   *
 *----------------------------------------------------------------------------*/

void set_parameters(void) {
  PARAMLIST params;
  char path[] = SOLNFILE; /* Must be allocated to heap so mkpath can alter */

  /* Before trying relaxation code make sure path for solution exist */
  if (mkpath(path, 0700)) {
    fprintf(stderr, "ERROR: Failed to make path for solution.\n");
    exit(401);
  }

  /* read and set the inputs for the model */
  set_bcs(&params);
  set_term_ind(&params);
  set_tech_params(&params);
  set_phys_params(&params);
  set_plnt_params(&params);

  /* Initalize global */
  parameters = params;

  if (GET_PARAMS_VERBOSE == ON) {
    print_params();
  }

  return;
}

/*============================================================================*
 *! \fn int initial_guess(EQNVARS *equationvars_p)                            *
 *  \brief Reads in and sets the initial guess for the relaxation domain      *
 *----------------------------------------------------------------------------*/

#define BUFSIZE 1024
void initial_guess(EQNVARS *equationvars_p) {
  int i,j,k,m;
  double radius, density, velocity, temperature, fraction[NSPECIES], column[NSPECIES], qew, zed;
  char *filename = GUESS_FILE;
  FILE *filep;
  char dline[BUFSIZE], *tok, *hline = NULL;
  long start;
  size_t len = 0, nheader = 0;
  ssize_t read;
    
  /* open the guess file for reading */
  if ((filep = fopen(filename, "r")) == NULL) {
    fprintf(stderr, "ERROR: Cannot open %s\n", filename);
    exit(402);
  }

  /* parse header comments */
  do {
    start = ftell(filep);
    read  = getline(&hline, &len, filep);
    nheader++;
#if HEADER_DIFF_VERBOSE
    strcpy(dline, hline);
    tok = strtok(hline, ":");
    /* planet parameters */
    if (strstr(hline, "nspecies") != NULL) {
      if (atof(tok = strtok(NULL, "\n")) != NSPECIES) {
        fprintf(stdout, "WARNING: NSPECIES defined in def.h differs from guess's:\n"
                "         guess:%d, run:%d\n",
                atoi(tok), NSPECIES);
      }
    }
    if (strstr(hline, "plnt_prms") != NULL) {
      if (atof(tok = strtok(NULL, ",")) != parameters.Mp) {
        fprintf(stdout, "WARNING: Mp differs from guess's:\n"
                "         guess:%e, run:%e\n",
                atof(tok), parameters.Mp);
      }
      if (atof(tok = strtok(NULL, ",")) != parameters.Rp) {
        fprintf(stdout, "WARNING: Rp differs from guess's:\n"
                "         guess:%e, run:%e\n",
                atof(tok), parameters.Rp);
      }
      if (atof(tok = strtok(NULL, ",")) != parameters.Mstar) {
        fprintf(stdout, "WARNING: Mstar differs from guess's:\n"
                "         guess:%e, run:%e\n",
                atof(tok), parameters.Mstar);
      }
      if (atof(tok = strtok(NULL, ",")) != parameters.semimajor) {
        fprintf(stdout, "WARNING: semimajor differs from guess's:\n"
                "         guess:%e, run:%e\n",
                atof(tok), parameters.semimajor);
      }
      if (atof(tok = strtok(NULL, ",")) != parameters.Ftot) {
        fprintf(stdout, "WARNING: Ftot differs from guess's:\n"
                "         guess:%e, run:%e\n",
                atof(tok), parameters.Ftot);
      }
      if (atof(tok = strtok(NULL, ",")) != parameters.Lstar) {
        fprintf(stdout, "WARNING: Lstar differs from guess's:\n"
                "         guess:%e, run:%e\n",
                atof(tok), parameters.Lstar);
      }
    }
    /* physics parameters */
    if (strstr(hline, "phys_prms") != NULL) {
      for (j=0;j<NSPECIES; j++){
          if (atof(tok = strtok(NULL, ",")) != parameters.HX[j]) {
            fprintf(stdout, "WARNING: Mass Fraction differs from guess's:\n"
                    "         guess:%e, run:%e\n",
                    atof(tok), parameters.HX[j]);
          }
      }
//         printf("past most warnings \n");
      for (j=0;j<NSPECIES; j++){
          tok = strtok(NULL, ",");
          if (strcmp(tok, parameters.species[j]) != 0) {
            fprintf(stdout, "WARNING: Species name differs from guess's:\n"
                    "         guess:%s, run:%s\n", tok, parameters.species[j]);
          }
      }
      for (j=0;j<NSPECIES; j++){
          if (atof(tok = strtok(NULL, ",")) != parameters.atomic_mass[j]) {
            fprintf(stdout, "WARNING: Atomic mass differs from guess's:\n"
                    "         guess:%e, run:%e\n",
                    atof(tok), parameters.atomic_mass[j]);
          }
      }
      if (atof(tok = strtok(NULL, ",")) != parameters.molec_adjust) {
        fprintf(stdout, "WARNING: Molecular adjustment factor differs from guess's:\n"
                "         guess:%e, run:%e\n",
                atof(tok), parameters.molec_adjust);
      }
    }
    /* boundary conditions */
    if (strstr(hline, "bcs") != NULL) {
      if (atof(tok = strtok(NULL, ",")) != parameters.Rmin) {
        fprintf(stdout, "WARNING: Rmin differs from guess's:\n"
                "         guess:%e, run:%e\n",
                atof(tok), parameters.Rmin);
      }
      if (atof(tok = strtok(NULL, ",")) != parameters.Rmax) {
        fprintf(stdout, "WARNING: Rmax differs from guess's:\n"
                "         guess:%e, run:%e\n",
                atof(tok), parameters.Rmax);
      }
      if (atof(tok = strtok(NULL, ",")) != parameters.rho_rmin) {
        fprintf(stdout, "WARNING: rho_rmin differs from guess's:\n"
                "         guess:%e, run:%e\n",
                atof(tok), parameters.rho_rmin);
      }
      if (atof(tok = strtok(NULL, ",")) != parameters.T_rmin) {
        fprintf(stdout, "WARNING: T_rmin differs from guess's:\n"
                "         guess:%e, run:%e\n",
                atof(tok), parameters.T_rmin);
      }
      for (j=0;j<NSPECIES; j++){
          if (atof(tok = strtok(NULL, ",")) != parameters.Ys_rmin[j]) {
            fprintf(stdout, "WARNING: Ys_rmin_%s differs from guess's:\n"
                    "         guess:%e, run:%e\n",
                    parameters.species[j], atof(tok), parameters.Ys_rmin[j]);
          }
      }
      for (j=0;j<NSPECIES; j++){
          if (atof(tok = strtok(NULL, ",")) != parameters.Ncol_sp[j]) {
            fprintf(stdout, "WARNING: Ncol_sp_%s differs from guess's:\n"
                    "         guess:%e, run:%e\n",
                    parameters.species[j], atof(tok), parameters.Ncol_sp[j]);
          }
      }
      for (m=0;m<2; m++){
        if (atof(tok = strtok(NULL, ",")) != parameters.erf_drop[m]) {
         fprintf(stdout, "WARNING: erf_drop[%d] differs from guess's:\n"
                "         guess:%e run:%e\n",
                m, atof(tok), parameters.erf_drop[m]);
        }
      }
    }
    /* techincal parameters */
    if (strstr(hline, "tech") != NULL) {
      if (atof(tok = strtok(NULL, ",")) != parameters.breezeparam) {
        fprintf(stdout, "WARNING: breezeparam differs from guess's:\n"
                "         guess:%e, run:%e\n",
                atof(tok), parameters.breezeparam);
      }
      if (atof(tok = strtok(NULL, ",")) != parameters.rapidity) {
        fprintf(stdout, "WARNING: rapidity differs from guess's:\n"
                "         guess:%e, run:%e\n",
                atof(tok), parameters.rapidity);
      }
      if (atof(tok = strtok(NULL, ",")) != parameters.erfn) {
        fprintf(stdout, "WARNING: erfn differs from guess's:\n"
                "         guess:%e, run:%e\n",
                atof(tok), parameters.erfn);
      }
      if (atof(tok = strtok(NULL, ",")) != parameters.mach_limit) {
        fprintf(stdout, "WARNING: mach_limit differs from guess's:\n"
                "         guess:%e, run:%e\n",
                atof(tok), parameters.mach_limit);
      }
    }
    /* indicator flags */
    if (strstr(hline, "flags") != NULL) {
      if (atoi(tok = strtok(NULL, ",")) != parameters.lyacool) {
        fprintf(stdout, "WARNING: lyacool differs from guess's:\n"
                "         guess:%d, run:%d\n",
                atoi(tok), parameters.lyacool);
      }
      if (atoi(tok = strtok(NULL, ",")) != parameters.tidalforce) {
        fprintf(stdout, "WARNING: tidalforce differs from guess's:\n"
                "         guess:%lf, run:%lf\n",
                atof(tok), parameters.tidalforce);
      }
      if (atof(tok = strtok(NULL, ",")) != parameters.bolo_heat_cool) {
        fprintf(stdout, "WARNING: bolo_heat_cool differs from guess's:\n"
                "         guess:%lf, run:%lf\n",
                atof(tok), parameters.bolo_heat_cool);
      }
      if (atoi(tok = strtok(NULL, ",")) != parameters.integrate_outward) {
        fprintf(stdout, "WARNING: integrate_outward differs from guess's:\n"
                "         guess:%d, run:%d\n",
                atoi(tok), parameters.integrate_outward);
      }
    }
//       /* Additional parameters */
//     if (strstr(hline, "add_prms") != NULL) {
//       if (atoi(tok = strtok(NULL, ",")) != parameters.lyacool) {
//         fprintf(stdout, "WARNING: lyacool differs from guess's:\n"
//                 "         guess:%d, run:%d\n",
//                 atoi(tok), parameters.lyacool);
//       }
//     }
//     #add_prms:
//       printf("past all warnings \n");
    /* spectrum */
    if (strstr(hline, "spec_prms") != NULL) {
      check_spec(dline);
    }
#endif
  } while ((hline[0]) == '#');
  nheader--;
  (void)read;
  printf("Skipping %zu lines of header\n", nheader);
  fseek(filep, start, SEEK_SET);

  /* Read in and assign the guesses */
  /* Each line of the input file has the form: */
  /* r rho v T Ys Ys2 Ncol Ncol2 q z */
  i = INPTS;
  while (fgets(dline, BUFSIZE, filep) && i < INPTS+M) { 
    tok = strtok(dline, ",");
    radius = atof(tok);
    tok = strtok(NULL, ",");
    density = atof(tok);
    tok = strtok(NULL, ",");
    velocity = atof(tok);
    tok = strtok(NULL, ",");
    temperature = atof(tok);
    for (k=0;k<NSPECIES;k++){
        tok = strtok(NULL, ",");
        fraction[k] = atof(tok);
    }
    for (k=0;k<NSPECIES;k++){
        tok = strtok(NULL, ",");
        column[k] = atof(tok);
    }
    tok = strtok(NULL, ",");
    qew = atof(tok);
    tok = strtok(NULL, ",");
    zed = atof(tok);
    if (qew >= 0.) {
      equationvars_p->r[i]    = radius;
      equationvars_p->rho[i]  = density;
      equationvars_p->v[i]    = velocity;
      equationvars_p->T[i]    = temperature;
      for (k=0;k<NSPECIES;k++){
          equationvars_p->Ys[i][k]   = fraction[k];
          equationvars_p->Ncol[i][k] = column[k];
      }
      equationvars_p->q[i]    = qew;
      equationvars_p->z[i]    = zed;
      if (equationvars_p->q[i++] == 1.0) {
        break;
      }
    }
  }
  /* close the file */
  fclose(filep);

  /* sanity check results read in */
  if (i-INPTS != M) {
    fprintf(stderr, "ERROR: Number of relaxtion points does not match "
            "the guess's.\n");
    fprintf(stderr, "       #pts: %d, M: %d (increase M in defs.h)\n",
            i-INPTS, M);
    exit(403);
  }
  if (equationvars_p->q[i-1] != 1.) {
    fprintf(stderr, "ERROR: Last read in q:%e != 1.0\n",
            equationvars_p->q[i-1]);
    exit(404);
  }

  return;
}
#undef BUFSIZE

/*============================================================================*
 *! \fn int save_solution(EQNVARS equationvars)                               *
 *  \brief Saves the solution of relaxation and odeint if requested to a csv  *
 *----------------------------------------------------------------------------*/

void save_solution(EQNVARS equationvars) {
  int i,j,m;
  FILE *filep;
  char *filename = SOLNFILE;

  /* check whether the file is already there */
  if ((filep = fopen(filename, "r")) != NULL && NOCLOBBER == ON) {
    fprintf(stderr, "ERROR: %s already exists & NOCLOBBER flag is ON.\n",
            filename);
    exit(405);
  }

  /* open the output file for writing */
  if ((filep = fopen(filename, "w")) == NULL) {
    fprintf(stderr, "ERROR: Cannot open %s\n", filename);
    exit(406);
  }

      /* Print physics parameters */
    //   char *physline = "";
  fprintf(filep, "#nspecies: ");
  fprintf(filep,"%d\n",NSPECIES);
    
  /* Using %.17e instead of more efficient %.17g for human readability */
      /* Print variables */
    char varline[1024] = "r,rho,v,T";
    fprintf(filep, "#vars: ");
    for (i=0; i<NSPECIES; i++){
        char ysname[1024] = ",Ys_";
        strncat(ysname,parameters.species[i],strlen(parameters.species[i]));
        strncat(varline,ysname,strlen(ysname));
    }
    for (i=0; i<NSPECIES; i++){
        char ncolname[1024] = ",Ncol_";
        strncat(ncolname,parameters.species[i],strlen(parameters.species[i]));
        strncat(varline,ncolname,strlen(ncolname));
    }
    strncat(varline,",q,z\n",6);
    fprintf(filep,"%s",varline);
    
      /* Print scalings */
  char *scaleline = "%.17e,%.17e,%.17e,%.17e,";
  fprintf(filep, "#scales: ");
  fprintf(filep, scaleline, parameters.Rp, RHO0, CS0, T0);
  char *scaleline_end = "%.17e,%.17e\n";
  for (i=0; i<NSPECIES; i++){
    fprintf(filep,"%.17e,",1.);
  }
  for (i=0; i<NSPECIES; i++){
    fprintf(filep,"%.17e,",NCOL0);
  }
  fprintf(filep, scaleline_end, 1., parameters.Rp);
  
  /* Print planet parameters */
  char *plntline = "%.17e,%.17e,%.17e,%.17e,%.17e,%.17e\n";
  fprintf(filep, "#plnt_prms: ");
  fprintf(filep, plntline, parameters.Mp, parameters.Rp, parameters.Mstar,
          parameters.semimajor, parameters.Ftot, parameters.Lstar);
          
  /* Print physics parameters */
//   char *physline = "";
  fprintf(filep, "#phys_prms: %.17e",parameters.HX[0]);
  for (i=1; i<NSPECIES; i++){
    fprintf(filep,",%.17e",parameters.HX[i]);
  }
  for (i=0; i<NSPECIES; i++){
    fprintf(filep,",%s",parameters.species[i]);
  }
  for (i=0; i<NSPECIES; i++){
    fprintf(filep,",%.17e",parameters.atomic_mass[i]);
  }
  fprintf(filep,",%.17e",parameters.molec_adjust);
  fprintf(filep,"%s","\n");
  
    /* Print boundary conditions */
  char *bcsline = "%.17e,%.17e,%.17e,%.17e";
  fprintf(filep, "#bcs: ");
  fprintf(filep, bcsline, parameters.Rmin, parameters.Rmax, parameters.rho_rmin,
          parameters.T_rmin);
  for (i=0; i<NSPECIES; i++){
    fprintf(filep,",%.17e",parameters.Ys_rmin[i]);
  }
  for (i=0; i<NSPECIES; i++){
    fprintf(filep,",%.17e",parameters.Ncol_sp[i]);
  }
  for (m=0; m<2; m++){
      fprintf(filep, ",%.17e", parameters.erf_drop[m]);
  }
  fprintf(filep,"%s","\n");
  
  /* Print technical terms */
  char *techline = "%.17e,%.17e,%.17e,%.17e \n";
  fprintf(filep, "#tech: ");
  fprintf(filep, techline, parameters.breezeparam, parameters.rapidity,
          parameters.erfn, parameters.mach_limit);

  /* Print flags */
  char *flagline = "%d,%.5f,%.5f,%d \n";
  fprintf(filep, "#flags: ");
  fprintf(filep, flagline, parameters.lyacool, parameters.tidalforce,
          parameters.bolo_heat_cool, parameters.integrate_outward);
    
  /* Print additional parameters */
    char *addline = "%d";
    fprintf(filep, "#add_prms: ");
    fprintf(filep, addline, N_ADD_PARAMS);
    if (N_ADD_PARAMS > 0){
      for (i=0; i<N_ADD_PARAMS; i++){
        fprintf(filep,",%.17e",parameters.add_params[i]);
      }
    }
      fprintf(filep,"%s","\n");
//     else{
//       fprintf(filep,",%.0f",0);
//       fprintf(filep,"%s","\n");
//     }
    
  /* Print spectrum parameters */
  save_spec(filep);
    
  /* Print out the columns */
  int start, end;
//   start = parameters.integrate_inward ? 0 : INPTS;
  start = 0;
  end   = parameters.integrate_outward ? TOTALPTS : INPTS+M;
  printf("Printing %d lines to %s\n", end-start, filename);
  /* Each line of the output file has the form: */
  /* Each line of the output file has the form: */
  /* r rho v T Ys Ys2 Ncol Ncol2 q z */
  char *dataformat = "%.17e,%.17e,%.17e,%.17e";
  char *dataformat_end = ",%.17e,%.17e \n";
  for (i = start; i < end; i++) {
    fprintf(filep, dataformat, equationvars.r[i], equationvars.rho[i], 
            equationvars.v[i], equationvars.T[i]);
            for (j=0; j<NSPECIES; j++){
                fprintf(filep, ",%.17e",equationvars.Ys[i][j]);
            }
            for (j=0; j<NSPECIES; j++){
                fprintf(filep, ",%.17e",equationvars.Ncol[i][j]);
            }
    fprintf(filep, dataformat_end, equationvars.q[i], equationvars.z[i]);
  }
    
  /* close the file */
  fclose(filep);

  return;
}

/*============================================================================*
 *=========================== PRIVATE FUNCTIONS ==============================*
 *============================================================================*
 * set_bcs                 - reads and sets boundary conditions               *
 * set_term_ind            - reads and sets term indicators                   *
 * set_tech_params         - reads and sets techincal parameters              *
 * set_phys_params         - reads and sets physics parameters                *
 * set_planet_params       - reads and sets planet parameters                 *
 *----------------------------------------------------------------------------*/

static int mkpath(char* path, mode_t mode) {
  char *slash = strchr(path+1, '/');

  if (path == NULL || !(*path)) {
    fprintf(stderr, "ERROR: path not properly specified\n");
    exit(407);
  }

  /* Traverse down path directories, making directories as needed */
  for (; slash; slash = strchr(slash+1, '/')) {
    /* truncate path and make parent directory */
    *slash = '\0';
    if (mkdir(path, mode) == -1) {
      if (errno != EEXIST) {
        *slash = '/';
        return 1;
      }
    }
    /* restore path char to forward slash */
    *slash = '/';
  }

  return 0;
}

static void set_bcs(PARAMLIST *params) {
  char *filename = BC_PARAM_FILE;
  FILE *filep;
  char line[1024], temp[1024], *hline = NULL;
  long start;
  size_t len = 0, nheader = 0;
  ssize_t read;
  int i,n,m;

  /* open the planet parameters file for reading */
  if ((filep = fopen(filename, "r")) == NULL) {
    fprintf(stderr, "ERROR: Cannot open %s\n", filename);
    exit(408);
  }

  /* Skip the header if it exist */
  do {
    start = ftell(filep);
    read  = getline(&hline, &len, filep);
    nheader++;
  } while ((hline[0]) == '#');
  nheader--;
  (void)nheader;
  (void)read;
  fseek(filep, start, SEEK_SET);

  /* Read in and assign the parameters */
  /* the lowest height of the calculation in Rp */
  fgets(line, sizeof(line), filep);
  sscanf(line, "%s %lf %[^\n]", temp, &(params->Rmin), temp);
  /* the highest height in the calculation in Rp */
  fgets(line, sizeof(line), filep);
  sscanf(line, "%s %lf %[^\n]", temp, &(params->Rmax), temp);
  /* the density at Rmin in RHO0 */
  fgets(line, sizeof(line), filep);
  sscanf(line, "%s %lf %[^\n]", temp, &(params->rho_rmin), temp);
  /* the temperature at Rmin in T0 */
  fgets(line, sizeof(line), filep);
  sscanf(line, "%s %lf %[^\n]", temp, &(params->T_rmin), temp);
  /* the ionization fraction of H at Rmin */
  fgets(line, sizeof(line), filep);
  sscanf(line, "%s %lf%n %[^,]", temp, &(params->Ys_rmin[0]), &n, temp);
  for (i = 1; i<NSPECIES; i=i+1){
      sscanf(line+n+1+i, "%lf%n %[^,]", &(params->Ys_rmin[i]), &m, temp);
      n=n+m;
  }  
  /* the column density of H at the sonic point */
  fgets(line, sizeof(line), filep);
  sscanf(line, "%s %lf%n %[^,]", temp, &(params->Ncol_sp[0]), &n, temp);
  for (i = 1; i<NSPECIES; i=i+1){
      sscanf(line+n+1+i, "%lf%n %[^,]", &(params->Ncol_sp[i]), &m, temp);
      n=n+m;
  }
  /*  */
  fgets(line, sizeof(line), filep);
  sscanf(line, "%s %lf%n %[^,]", temp, &(params->erf_drop[0]), &n, temp);
  sscanf(line+n+1, "%lf%n %[^\n]", &(params->erf_drop[1]), &m, temp);
  /* close the file */
  fclose(filep);

  return;
}

static void set_term_ind(PARAMLIST *params) {
  char *filename = TERM_PARAM_FILE;
  FILE *filep;
  char line[1024], temp[1024], *hline = NULL;
  long start;
  size_t len = 0, nheader = 0;
  ssize_t read;

  /* open the planet parameters file for reading */
  if ((filep = fopen(filename, "r")) == NULL) {
    fprintf(stderr, "ERROR: Cannot open %s\n", filename);
    exit(409);
  }

  /* Skip the header if it exist */
  do {
    start = ftell(filep);
    read  = getline(&hline, &len, filep);
    nheader++;
  } while ((hline[0]) == '#');
  nheader--;
  (void)nheader;
  (void)read;
  fseek(filep, start, SEEK_SET);

  /* Read in and assign the parameters */
  /* Lyman alpha cooling on or off */
  fgets(line, sizeof(line), filep);
  sscanf(line, "%s %d %[^\n]", temp, &(params->lyacool), temp);
  /* Tidal force included or not */
  fgets(line, sizeof(line), filep);
  sscanf(line, "%s %lf %[^\n]", temp, &(params->tidalforce), temp);
  /* Inward integration included or not */
  fgets(line, sizeof(line), filep);
  sscanf(line, "%s %lf %[^\n]", temp, &(params->bolo_heat_cool), temp);
  /* Outward integration included or not */
  fgets(line, sizeof(line), filep);
  sscanf(line, "%s %d %[^\n]", temp, &(params->integrate_outward), temp);

  /* close the file */
  fclose(filep);

  return;
}

static void set_tech_params(PARAMLIST *params) {
  char *filename = TECH_PARAM_FILE;
  FILE *filep;
  char line[1024], temp[1024], *hline = NULL;
  long start;
  size_t len = 0, nheader = 0;
  ssize_t read;

  /* open the planet parameters file for reading */
  if ((filep = fopen(filename, "r")) == NULL) {
    fprintf(stderr, "ERROR: Cannot open %s\n", filename);
    exit(410);
  }

  /* Skip the header if it exist */
  do {
    start = ftell(filep);
    read  = getline(&hline, &len, filep);
    nheader++;
  } while ((hline[0]) == '#');
  nheader--;
  (void)read;
  (void)nheader;
  fseek(filep, start, SEEK_SET);

  /* Read in and assign the parameters */
  /* this is 1 for a wind, < 1 for a breeze */
  fgets(line, sizeof(line), filep);
  sscanf(line, "%s %lf %[^\n]", temp, &(params->breezeparam), temp);
  /* rapidity parameter for how fast weighting transitions */
  fgets(line, sizeof(line), filep);
  sscanf(line, "%s %lf %[^\n]", temp, &(params->rapidity), temp);
  /* erfn parameter controls how close to get to 1. in weighting */
  fgets(line, sizeof(line), filep);
  sscanf(line, "%s %lf %[^\n]", temp, &(params->erfn), temp);
  /* mach_limit parameter for limit after which dvdr purely linearized */
  fgets(line, sizeof(line), filep);
  sscanf(line, "%s %lf %[^\n]", temp, &(params->mach_limit), temp);

  /* close the file */
  fclose(filep);

  return;
}

static void set_phys_params(PARAMLIST *params) {
  char *filename = PHYS_PARAM_FILE;
  FILE *filep;
  char line[1024], temp[1024], *hline = NULL;
  long start;
  size_t len = 0, nheader = 0;
  ssize_t read;
  int i,n,m;

  /* open the planet parameters file for reading */
  if ((filep = fopen(filename, "r")) == NULL) {
    fprintf(stderr, "ERROR: Cannot open %s\n", filename);
    exit(411);
  }

  /* Skip the header if it exist */
  do {
    start = ftell(filep);
    read  = getline(&hline, &len, filep);
    nheader++;
  } while ((hline[0]) == '#');
  nheader--;
  (void)nheader;
  (void)read;
  fseek(filep, start, SEEK_SET);

  /* Read in and assign the parameters */
  /* the mass fraction */
  fgets(line, sizeof(line), filep);
  sscanf(line, "%s %lf%n %[^,]", temp, &(params->HX[0]), &n, temp);
  for (i = 1; i<NSPECIES; i=i+1){
      sscanf(line+n+1+i, "%lf%n %[^,]", &(params->HX[i]), &m, temp);
      n=n+m;
  }
  
  /*The species name in Arabic+Roman format*/
  fgets(line, sizeof(line), filep);
  char *tok; 
  tok = strtok(line, ":");
  for (i=1;i<=NSPECIES;i++){
      if (i==NSPECIES){
          tok = strtok(NULL, "#");
      }
      else{
          tok = strtok(NULL, ",");
      }
      sprintf((params->species[i-1]), "%s", tok);
  }
  
    /* the atomic mass in g */
  fgets(line, sizeof(line), filep);
  sscanf(line, "%s %lf%n %[^,]", temp, &(params->atomic_mass[0]), &n, temp);
  for (i = 1; i<NSPECIES; i=i+1){
      sscanf(line+n+1+i, "%lf%n %[^,]", &(params->atomic_mass[i]), &m, temp);
      n=n+m;
  }
    
  
/* Z and N_e, number of electrons in atom, for computing recombination coefficient*/
  fgets(line, sizeof(line), filep);
  sscanf(line, "%s %d%n %[^,]", temp, &(params->Z[0]), &n, temp);
  for (i = 1; i<NSPECIES; i=i+1){
      sscanf(line+n+1+i, "%d%n %[^,]", &(params->Z[i]), &m, temp);
      n=n+m;
  }
  
  fgets(line, sizeof(line), filep);
  sscanf(line, "%s %d%n %[^,]", temp, &(params->N_e[0]), &n, temp);
  for (i = 1; i<NSPECIES; i=i+1){
      sscanf(line+n+1+i, "%d%n %[^,]", &(params->N_e[i]), &m, temp);
      n=n+m;
  }

  /* Molecular adjustment (dim'less) multiplication factor to account for higher mean mol weight below wind */  
  fgets(line, sizeof(line), filep);
  sscanf(line, "%s %lf%n %[^,]", temp, &(params->molec_adjust), &n, temp);

  /* close the file */
  fclose(filep);

  return;
}

static void set_plnt_params(PARAMLIST *params) {
  char *filename = PLANET_PARAM_FILE;
  FILE *filep;
  char line[1024], temp[1024], *hline = NULL;
  long start;
  size_t len = 0, nheader = 0;
  ssize_t read;

  /* open the planet parameters file for reading */
  if ((filep = fopen(filename, "r")) == NULL) {
    fprintf(stderr, "ERROR: Cannot open %s\n", filename);
    exit(412);
  }

  /* Skip the header if it exist */
  do {
    start = ftell(filep);
    read  = getline(&hline, &len, filep);
    nheader++;
  } while ((hline[0]) == '#');
  nheader--;
  (void)nheader;
  (void)read;
  fseek(filep, start, SEEK_SET);

  /* Read in and assign the planet parameters */
  /* the planet's mass in grams */
  fgets(line, sizeof(line), filep);
  sscanf(line, "%s %lf %[^\n]", temp, &(params->Mp), temp);
  /* the planet's radius in cm */
  fgets(line, sizeof(line), filep);
  sscanf(line, "%s %lf %[^\n]", temp, &(params->Rp), temp);
  /* the star's mass: read in in Msun and convert to g */
  fgets(line, sizeof(line), filep);
  sscanf(line, "%s %lf %[^\n]", temp, &(params->Mstar), temp);
  /* the planet's semi-major axis: read in in AU and convert to cm */
  fgets(line, sizeof(line), filep);
  sscanf(line, "%s %lf %[^\n]", temp, &(params->semimajor), temp);
  /* the UV flux at the planet's surface in erg/cm^2/s */
  fgets(line, sizeof(line), filep);
  sscanf(line, "%s %lf %[^\n]", temp, &(params->Ftot), temp);
  /* the star's luminosity in ergs/s */
  fgets(line, sizeof(line), filep);
  sscanf(line, "%s %lf %[^\n]", temp, &(params->Lstar), temp);
  /* the scale height of the atmosphere */
  params->H0 = pow(CS0, 2)*pow(params->Rp, 2)/G/params->Mp;

  /* close the file */
  fclose(filep);

  return;
}
    
// static void add_params(PARAMLIST *params) {
//   char *filename = ADD_PARAM_FILE;
//   FILE *filep;
//   char line[1024], temp[1024], *hline = NULL;
//   long start;
//   size_t len = 0, nheader = 0;
//   ssize_t read;
//   int i,n,m;

//   /* open the planet parameters file for reading */
//   if ((filep = fopen(filename, "r")) == NULL) {
//     fprintf(stderr, "ERROR: Cannot open %s\n", filename);
//     exit(408);
//   }

//   /* Skip the header if it exist */
//   do {
//     start = ftell(filep);
//     read  = getline(&hline, &len, filep);
//     nheader++;
//   } while ((hline[0]) == '#');
//   nheader--;
//   fseek(filep, start, SEEK_SET);

//   /* Read in and assign the parameters */
//   /* Number of additional params */
//   fgets(line, sizeof(line), filep);
// //  sscanf(line, "%s %lf %[^\n]", temp, &(params->num_additional_params), temp);

//   for (i = 1; i<N_ADD_PARAMS; i=i+1){
//       sscanf(line+i, "%lf%n %[^,]", &(params->add_params[i]), &m, temp);
//       n=n+m;
//   }

//   /* close the file */
//   fclose(filep);

//   return;
// }

static void print_params(void) {
  int i;
  /* print out all of the parameters */
  printf("Full parameter list:\n");
  /* print out the planet parameters */
  printf("The planet parameters are:\n");
  printf("  Mp = %g g\n", parameters.Mp);
  printf("  Rp = %g cm\n", parameters.Rp);
  printf("  Mstar = %g g\n", parameters.Mstar);
  printf("  semimajor = %g cm\n", parameters.semimajor);
  printf("  Ftot = %g erg/cm^2/s\n", parameters.Ftot);
  printf("  Lstar = %g erg/s\n", parameters.Lstar);
  /* print out the physics parameters */
  printf("The physics parameters are:\n");
  printf("  HX = %g\n", parameters.HX[0]);
  for (i=1;i<NSPECIES; i=i+1){
  printf("       %g\n", parameters.HX[i]);
  }
  /* print out the boundary conditions */
  printf("The boundary conditions are:\n");
  printf("  Rmin = %g Rp\n", parameters.Rmin);
  printf("  Rmax = %g Rp\n", parameters.Rmax);
  printf("  rho_rmin = %g RHO0\n", parameters.rho_rmin);
  printf("  T_rmin = %g T0\n", parameters.T_rmin);
  printf("  Ys_rmin = %g\n", parameters.Ys_rmin[0]);
  for (i=1;i<NSPECIES; i=i+1){
  printf("            %g\n", parameters.Ys_rmin[i]);
  }
  printf("  N_sp = %g NCOL0\n", parameters.Ncol_sp[0]);
  for (i=1;i<NSPECIES; i=i+1){
      printf("         %g NCOL0\n", parameters.Ncol_sp[i]);
  }
  printf("  erf_drop (center)   = %g V0\n", parameters.erf_drop[0]);
  printf("  erf_drop (gradient) = %g V0\n", parameters.erf_drop[1]);
    
  /* print out which terms are ON/OFF */
  printf("The term flags are:\n");
  printf("  lyacool = %d\n", parameters.lyacool);
  printf("  tidalforce = %.5f\n", parameters.tidalforce);
  printf("  bolo_heat_cool = %.5f\n", parameters.bolo_heat_cool);
  printf("  integrate_outward = %d\n", parameters.integrate_outward);
  /* print out the technical parameters */
  printf("The technical parameters are:\n");
  printf("  breezeparam = %g\n", parameters.breezeparam);
  printf("  rapidity = %g\n", parameters.rapidity);
  printf("  erfn = %g\n", parameters.erfn);
  printf("  Mach limit = %g\n", parameters.mach_limit);

  return;
}
