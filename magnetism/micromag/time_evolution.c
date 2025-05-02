void TimeEvolution(void)
{
  InitialCondition();
  // m(t) obtained
  
  plan_fft();
    
  initialise_tensors();
    
  generate_N();
    
  uint32_t t = 0;
  
  while(t < iterations)
    {
      H_eff();
      // h(m) obtained

      GaussSeidel(m2 , h2 , gs_coeff_2 , true);
      // g2 obtained
        
      GaussSeidel(m3 , h3 , gs_coeff_3 , true);
      // g3 obtained
		
#pragma omp parallel for if(N > N0)
      for(point = 0; point < N; point++)
       	{
	  m1[point][0] += gs_coeff_2[point][0] * m3[point][0] - gs_coeff_3[point][0] * m2[point][0];
       	}
      // m1* obtained
      
      GaussSeidel(m1 , h1 , gs_coeff_1 , true);
      // g1* obtained
		
#pragma omp parallel for if(N > N0)
      for(point = 0; point < N; point++)
       	{
	  m2[point][0] += gs_coeff_3[point][0] * m1[point][0] - gs_coeff_1[point][0] * m3[point][0];
       	}
      // m2* obtained
       	
      GaussSeidel(m2 , h2 , gs_coeff_2 , true);
      // g2* obtained
		
#pragma omp parallel for if(N > N0)
      for(point = 0; point < N; point++)
       	{
	  m3[point][0] += gs_coeff_1[point][0] * m2[point][0] - gs_coeff_2[point][0] * m1[point][0];
       	}
      // m3* obtained
		
      H_eff();
      // h(m*) obtained
      
      GaussSeidel(m1 , h1 , m1 , false);
      GaussSeidel(m2 , h2 , m2 , false);
      GaussSeidel(m3 , h3 , m3 , false);
      // m** obtained
      
#pragma omp parallel for if(N > N0)
      for(point = 0; point < N; point++)
       	{
	  double m_mag = sqrt(m1[point][0] * m1[point][0] + m2[point][0] * m2[point][0] + m3[point][0] * m3[point][0]);

	  m1[point][0] /= m_mag;
	  m2[point][0] /= m_mag;
	  m3[point][0] /= m_mag;
       	}
      // m(t+tick) obtained
      
      if((++t) % save_after == 0) savedata(t);
    }
}
