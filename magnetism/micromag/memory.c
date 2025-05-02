void ALLOCATE_MEMORY()
{
  if((input = (fftw_complex*)fftw_malloc(N * sizeof(fftw_complex))) == NULL)
    {
      printf("input could not be allocated!"); exit(1);
    }
  if((output = (fftw_complex*)fftw_malloc(N * sizeof(fftw_complex))) == NULL)
    {
      printf("output could not be allocated!"); exit(1);
    }
  if((m1 = (fftw_complex*)fftw_malloc(N * sizeof(fftw_complex))) == NULL)
    {
      printf("m1 could not be allocated!"); exit(1);
    }
  if((m2 = (fftw_complex*)fftw_malloc(N * sizeof(fftw_complex))) == NULL)
    {
      printf("m2 could not be allocated!"); exit(1);
    }
  if((m3 = (fftw_complex*)fftw_malloc(N * sizeof(fftw_complex))) == NULL)
    {
      printf("m3 could not be allocated!"); exit(1);
    }
  if((h1 = (fftw_complex*)fftw_malloc(N * sizeof(fftw_complex))) == NULL)
    {
      printf("h1 could not be allocated!"); exit(1);
    }
  if((h2 = (fftw_complex*)fftw_malloc(N * sizeof(fftw_complex))) == NULL)
    {
      printf("h2 could not be allocated!"); exit(1);
    }
  if((h3 = (fftw_complex*)fftw_malloc(N * sizeof(fftw_complex))) == NULL)
    {
      printf("h3 could not be allocated!"); exit(1);
    }
  if((gs_coeff_1 = (fftw_complex*)fftw_malloc(N * sizeof(fftw_complex))) == NULL)
    {
      printf("gs_coeff_1 could not be allocated!"); exit(1);
    }
  if((gs_coeff_2 = (fftw_complex*)fftw_malloc(N * sizeof(fftw_complex))) == NULL)
    {
      printf("gs_coeff_2 could not be allocated!"); exit(1);
    }
  if((gs_coeff_3 = (fftw_complex*)fftw_malloc(N * sizeof(fftw_complex))) == NULL)
    {
      printf("gs_coeff_3 could not be allocated!"); exit(1);
    }
  if((sclr_pot = (fftw_complex*)fftw_malloc(N * sizeof(fftw_complex))) == NULL)
    {
      printf("sclr_pot could not be allocated!"); exit(1);
    }
  if((temp_gs = (fftw_complex*)fftw_malloc(N * sizeof(fftw_complex))) == NULL)
    {
      printf("temp_gs could not be allocated!"); exit(1);
    }
  if((fourier_m1 = (fftw_complex*)fftw_malloc(N * sizeof(fftw_complex))) == NULL)
    {
      printf("fourier_m1 could not be allocated!"); exit(1);
    }
  if((fourier_m2 = (fftw_complex*)fftw_malloc(N * sizeof(fftw_complex))) == NULL)
    {
      printf("fourier_m2 could not be allocated!"); exit(1);
    }
  if((Hdemag1 = (fftw_complex*)fftw_malloc(N * sizeof(fftw_complex))) == NULL)
    {
      printf("Hdemag1 could not be allocated!"); exit(1);
    }
  if((Hdemag2 = (fftw_complex*)fftw_malloc(N * sizeof(fftw_complex))) == NULL)
    {
      printf("Hdemag2 could not be allocated!"); exit(1);
    }
  for(p = 0; p < 6; p++)
    {
      if((SFS[p].eps0 = (fftw_complex*)fftw_malloc(N * sizeof(fftw_complex))) == NULL)
	{
	  printf("SFS[%d].eps0 could not be allocated!",p); exit(1);
	}
      if((strain[p].strn_fld = (double*)malloc(N * sizeof(double))) == NULL)
	{
	  printf("strain[%d].strn_fld could not be allocated!",p); exit(1);
	}
    }
  if((N_ij = (Cof_K*)malloc(N * sizeof(Cof_K))) == NULL)
    {
      printf("N_ij could not be allocated!"); exit(1);
    }
  if((X = (homog_stress_f*)malloc(N * sizeof(homog_stress_f))) == NULL)
    {
      printf("X could not be allocated!"); exit(1);
    }
  for(point = 0; point < N; point++)
    {
      if((X[point].comp = (fftw_complex*)malloc(3 * sizeof(fftw_complex))) == NULL)
	{
	  printf("X[%d].comp could not be allocated!",p); exit(1);
	}
    }
  if((strain_hetero = (fftw_complex*)fftw_malloc(N * sizeof(fftw_complex))) == NULL)
    {
      printf("strain_hetero could not be allocated!"); exit(1);
    }
  if((dist_fx = (double*)malloc(N * sizeof(double))) == NULL)
    {
      printf("dist_fx could not be allocated!"); exit(1);
    }
  if((dist_fy = (double*)malloc(N * sizeof(double))) == NULL)
    {
      printf("dist_fy could not be allocated!"); exit(1);
    }
  if((zeta_mag = (double*)malloc(N * sizeof(double))) == NULL)
    {
      printf("zeta_mag could not be allocated!"); exit(1);
    }
  if((inv_Det_K = (double*)malloc(N * sizeof(double))) == NULL)
    {
      printf("inv_Det_K could not be allocated!"); exit(1);
    }
}

void FREE_MEMORY()
{
  fftw_destroy_plan(do_fourier);
  fftw_destroy_plan(do_fourier2);
  fftw_destroy_plan(do_inv_fourier);
    
  fftw_free(input);
  fftw_free(output);
  fftw_free(m1);
  fftw_free(m2);
  fftw_free(m3);
  fftw_free(h1);
  fftw_free(h2);
  fftw_free(h3);
  fftw_free(gs_coeff_1);
  fftw_free(gs_coeff_2);
  fftw_free(gs_coeff_3);
  fftw_free(sclr_pot);
  fftw_free(fourier_m1);
  fftw_free(fourier_m2);
  fftw_free(temp_gs);
  fftw_free(Hdemag1);
  fftw_free(Hdemag2);
  for(p = 0; p < 6; p++)
    {
      fftw_free(SFS[p].eps0);
      free(strain[p].strn_fld);
    }
  for(point = 0; point < N; point++) free(X[point].comp);
  free(X);
  fftw_free(strain_hetero);
  free(N_ij);
  free(dist_fx);
  free(dist_fy);
  free(zeta_mag);
  free(inv_Det_K);
}
