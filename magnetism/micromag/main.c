#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <stdbool.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <omp.h>
#include <fftw3.h>

#include "constants.h"
#include "memory.c"
#include "plan_and_store.c"
#include "savedata.c"
#include "initial.c"
#include "initialise_tensors.c"
#include "generate_N_ij.c"
#include "eigen_strain.c"
#include "generate_Xj.c"
#include "generate_strain.c"
#include "E_anisotropy.c"
#include "E_external.c"
#include "E_demag.c"
#include "E_mstr.c"
#include "h_eff.c"
#include "gauss_seidel.c"
#include "time_evolution.c"

int main(void)
{
  srand(time(NULL));
  clock_t start , end;
  double rt;
  start = clock();

  ALLOCATE_MEMORY();

  omp_set_num_threads(4);

  /* fftw_init_threads(); */

  ////////////////
  TimeEvolution();
  ////////////////

  FREE_MEMORY();

  /* fftw_cleanup_threads(); */
  
  end = clock();
  rt = (double)(end - start) / CLOCKS_PER_SEC;
  printf("Execution complete in %lf sec.\n",rt);

  return 0;
}
