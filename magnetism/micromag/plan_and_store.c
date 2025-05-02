void plan_fft()
{
  /* fftw_plan_with_nthreads(4); */
    
  do_fourier = fftw_plan_dft_1d(N , input , output , FFTW_FORWARD , FFTW_MEASURE);
  do_inv_fourier = fftw_plan_dft_1d(N , input , output , FFTW_BACKWARD , FFTW_MEASURE);
}

void run_fft(fftw_complex* in_array , fftw_complex* out_array , int8_t forward)
{
  
#pragma omp parallel for if(N > N0)
  for(point = 0; point < N; point++)
    {   
      input[point][0] = in_array[point][0];
      input[point][1] = in_array[point][1];
    }
    
  if(forward == 1) fftw_execute(do_fourier);
  else if(forward == -1) fftw_execute(do_inv_fourier);
    
#pragma omp parallel for if(N > N0)
  for(point = 0; point < N; point++)
    {   
      out_array[point][0] = output[point][0];
      out_array[point][1] = output[point][1];
    }
}
