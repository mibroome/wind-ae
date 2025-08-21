/*============================================================================*
 *! \file utils.c                                                             *
 *  \brief File with all the utilities used throughout the code               *
 *============================================================================*/

#include <stdio.h>
#include <stdlib.h>

/*============================================================================*
 *=========================== PUBLIC FUNCTIONS ===============================*
 *============================================================================*
 * calloc_1d_array       - callocs a 1d array with 0-based indexing           *
 * free_1d_array         - frees 1d array that had 0-based indexing           *
 * calloc_1d_array_gross - callocs a 1d array with 1-based indexing           *
 * calloc_2d_array_gross - callocs a 2d array with 1-based indexing           *
 * calloc_3d_array_gross - callocs a 3d array with 1-based indexing           *
 * free_1d_array_gross   - frees 1d array that had 1-based indexing           *
 * free_2d_array_gross   - frees 2d array that had 1-based indexing           *
 * free_3d_array_gross   - frees 3d array that had 1-based indexing           *
 *----------------------------------------------------------------------------*/

/*============================================================================*
 *! \fn void *calloc_1d_array(size_t nc, size_t size)                         *
 *  \brief Callocs a 1d array                                                 *
 *----------------------------------------------------------------------------*/

void *calloc_1d_array(size_t nc, size_t size){
  void *array;

  if ((array = (void*)calloc(nc, size)) == NULL) {
    fprintf(stderr, "ERROR: [calloc_1d] failed to allocate memory"
            "(%d of size %d)\n", (int)nc, (int)size);
    exit(801);
  }

  return array;
}

void **calloc_2d_array(size_t nrows, size_t ncols, size_t size) {
    void **array = malloc(nrows * sizeof(double *));
    if (!array) return NULL;
    double *data = calloc(nrows * ncols, size);
    if (!data) {
        free(array);
        return NULL;
    }
    for (size_t i = 0; i < nrows; i++) {
        array[i] = data + i * ncols;
    }
    return array;
}

/*============================================================================*
 *! \fn void free_1d_array(void *array)                                       *
 *  \brief Free memory used by 1D array                                       *
 *----------------------------------------------------------------------------*/

void free_1d_array(void *array){
  free(array);

  return;
}

void free_2d_array(void **array){
  free(array[0]);
  free(array);

  return;
}


/*============================================================================*
 *! \fn void *calloc_1d_array_gross(long is, long ie)                         *
 *  \brief Callocs a 1d array with 1-based indexing that spans array[is..ie]  *
 *----------------------------------------------------------------------------*/

double *calloc_1d_array_gross(long is, long ie){
  return (double*)calloc_1d_array(ie-is+2, sizeof(double))-is+1;
}

/*============================================================================*
 *! \fn void *calloc_2d_array_gross(long is, long ie, long js, long je)       *
 *  \brief Callocs a 2d array with 1-based indexing that spans                *
 *          array[is..ie][js..je]                                             *
 *----------------------------------------------------------------------------*/

double **calloc_2d_array_gross(long is, long ie, long js, long je){
  long i, nx, ny;
  double **array;

  nx = ie-is+1;
  ny = je-js+1;

  /* calloc pointers to rows */
  array  = (double**)calloc_1d_array(nx+1, sizeof(double*));
  array += 1;
  array -= is;

  /* calloc row and set pointers */
  array[is]  = (double*)calloc_1d_array(nx*ny+1, sizeof(double));
  array[is] += 1;
  array[is] -= js;
  for (i = is+1; i <= ie; i++) {
    array[i] = array[i-1]+ny;
  }

  return array;
}

/*============================================================================*
 *! \fn void *calloc_3d_array_gross(long is, long ie, long js, long je,       *
 *                                  long ks, long ke)                         *
 *  \brief Callocs a 3d array with 1-based indexing that spans                *
 *         array[is..ie][js..je][ks..ke]                                      *
 *----------------------------------------------------------------------------*/

double ***calloc_3d_array_gross(long is, long ie, long js, long je, long ks,
                                long ke){
  long i, j, nx, ny, nz;
  double ***array;

  nx = ie-is+1;
  ny = je-js+1;
  nz = ke-ks+1;

  /* calloc pointers to pointers to rows */
  array  = (double***)calloc_1d_array(nx+1, sizeof(double**));
  array += 1;
  array -= is;

  /* calloc pointers to rows and set pointers */
  array[is]  = (double**)calloc_1d_array(nx*ny+1, sizeof(double*));
  array[is] += 1;
  array[is] -= js;

  /* calloc rows and set pointers */
  array[is][js]  = (double*)calloc_1d_array(nx*ny*nz+1, sizeof(double**));
  array[is][js] += 1;
  array[is][js] -= ks;
  for (j = js+1; j <= je; j++) {
    array[is][j] = array[is][j-1]+nz;
  }
  for (i = is+1; i <= ie; i++) {
    array[i] = array[i-1]+ny;
    array[i][js] = array[i-1][js]+ny*nz;
    for (j = js+1; j <= je; j++) {
      array[i][j] = array[i][j-1]+nz;
    }
  }

  return array;
}

/*============================================================================*
 *! \fn void free_1d_array_gross(double *array, long is)                      *
 *  \brief Free memory used by a gross 1D array                               *
 *----------------------------------------------------------------------------*/

void free_1d_array_gross(double *array, long is){
  free((char*)(array+is-1));

  return;
}

/*============================================================================*
 *! \fn void free_2d_array_gross(double **array, long is, long js)            *
 *  \brief Free memory used by a gross 2D array                               *
 *----------------------------------------------------------------------------*/

void free_2d_array_gross(double **array, long is, long js){
  free((char*)(array[is]+js-1));
  free((char*)(array+is-1));

  return;
}

/*============================================================================*
 *! \fn void free_3d_array_gross(double ***array, long is, long js, long ks)  *
 *  \brief Free memory used by a gross 3D array                               *
 *----------------------------------------------------------------------------*/

void free_3d_array_gross(double ***array, long is, long js, long ks){
  free((char*)(array[is][js]+ks-1));
  free((char*)(array[is]+js-1));
  free((char*)(array+is-1));

  return;
}
