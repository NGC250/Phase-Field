void H_eff()
{
  
#pragma omp parallel for if(N > N0)
  for(point = 0; point < N; point++)
    {        
      h1[point][0] = 0.0;
      h2[point][0] = 0.0;
      h3[point][0] = 0.0;
        
      h1[point][1] = 0.0;
      h2[point][1] = 0.0;
      h3[point][1] = 0.0;
    }

  Anisotropy();
    
  /* External(); */

  /* Demagnetization(); */

  /* Magnetostriction(); */
}
