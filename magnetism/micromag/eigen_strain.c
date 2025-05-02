void gen_eigen_strain(void)
{
  
#pragma omp parallel for if(N > N0)
  for(point = 0; point < N; point++)
    {
      SFS[0].eps0[point][0] = Mstr_coeff[0][0] * (m1[point][0] * m1[point][0] - 1.0/3.0);
      SFS[1].eps0[point][0] = Mstr_coeff[1][1] * (m2[point][0] * m2[point][0] - 1.0/3.0);
      SFS[2].eps0[point][0] = Mstr_coeff[2][2] * (m3[point][0] * m3[point][0] - 1.0/3.0);
      SFS[3].eps0[point][0] = Mstr_coeff[0][1] * m1[point][0] * m2[point][0];
      SFS[4].eps0[point][0] = Mstr_coeff[0][2] * m1[point][0] * m3[point][0];
      SFS[5].eps0[point][0] = Mstr_coeff[1][2] * m2[point][0] * m3[point][0];
        
      SFS[0].eps0[point][1] = 0.0;
      SFS[1].eps0[point][1] = 0.0;
      SFS[2].eps0[point][1] = 0.0;
      SFS[3].eps0[point][1] = 0.0;
      SFS[4].eps0[point][1] = 0.0;
      SFS[5].eps0[point][1] = 0.0;
    }
    
  run_fft(SFS[0].eps0 , SFS[0].eps0 , 1);
  run_fft(SFS[1].eps0 , SFS[1].eps0 , 1);
  run_fft(SFS[2].eps0 , SFS[2].eps0 , 1);
  run_fft(SFS[3].eps0 , SFS[3].eps0 , 1);
  run_fft(SFS[4].eps0 , SFS[4].eps0 , 1);
  run_fft(SFS[5].eps0 , SFS[5].eps0 , 1);
}
