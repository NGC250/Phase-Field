void GaussSeidel(fftw_complex* mag_comp , fftw_complex* h_comp , fftw_complex* return_array , bool coeff)
{
  double const1 , const2;
    
  const1 = (coeff == true) ? tick : aux2;
  const2 = (coeff == true) ? aux1 : aux3;
    
#pragma omp parallel for if(N > N0)
  for(point = 0; point < N; point++)
    {        
      temp_gs[point][0] = mag_comp[point][0] + const1 * h_comp[point][0];
      temp_gs[point][1] = 0.0;
    }
  
  run_fft(temp_gs , temp_gs , 1);

  for(point = 0; point < N; point++)
    {        
      double frac = 1.0/(1.0 + zeta_mag[point] * const2);
        
      temp_gs[point][0] *= frac;
      temp_gs[point][1] *= frac;
    }
    
  run_fft(temp_gs , temp_gs , -1);
  
#pragma omp parallel for if(N > N0)
  for(point = 0; point < N; point++) return_array[point][0] = temp_gs[point][0] * norm;
}
