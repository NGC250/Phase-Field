void Demagnetization(void)
{
#pragma omp parallel for if(N > N0)
  for(point = 0; point < N; point++)
    {
      m1[point][1] = 0.0;
      m2[point][1] = 0.0;    
    }
    
  run_fft(m1 , fourier_m1 , 1);
  run_fft(m2 , fourier_m2 , 1);

  sclr_pot[0][0] = 0.0;
  sclr_pot[0][1] = 0.0;
    
  double frac;
    
  for(point = 1; point < N; point++)
    {
      zeta_x = dist_fx[point];
      zeta_y = dist_fy[point];
        
      sclr_pot[point][0] = fourier_m1[point][0] * zeta_x + fourier_m2[point][0] * zeta_y;
      sclr_pot[point][1] = fourier_m1[point][1] * zeta_x + fourier_m2[point][1] * zeta_y;

      frac = 1.0 / zeta_mag[point];
        
      sclr_pot[point][0] *= frac;
      sclr_pot[point][1] *= frac;    
    }
    
#pragma omp parallel for if(N > N0)
  for(point = 0; point < N; point++)
    {
      Hdemag1[point][0] = sclr_pot[point][0] * dist_fx[point];
      Hdemag1[point][1] = sclr_pot[point][1] * dist_fx[point];
        
      Hdemag2[point][0] = sclr_pot[point][0] * dist_fy[point];
      Hdemag2[point][1] = sclr_pot[point][1] * dist_fy[point];
    }

  run_fft(Hdemag1 , Hdemag1 , -1);
  run_fft(Hdemag2 , Hdemag2 , -1);

#pragma omp parallel for if(N > N0)
  for(point = 0; point < N; point++)
    {
      h1[point][0] -= 0.5 * norm * Hdemag1[point][0];
      h2[point][0] -= 0.5 * norm * Hdemag2[point][0];
    }
}
