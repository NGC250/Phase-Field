void Magnetostriction(void)
{
  gen_eigen_strain();
    
  generate_Xj();
    
  //gen_homog_strain();
    
  gen_hetero_strain(0,0);
#pragma omp parallel for if(N > N0)
  for(point = 0; point < N; point++) strain[0].strn_fld[point] = strain_hetero[point][0];
    
  gen_hetero_strain(1,1);
#pragma omp parallel for if(N > N0)
  for(point = 0; point < N; point++) strain[1].strn_fld[point] = strain_hetero[point][0];
    
  //gen_hetero_strain(2,2);
  //#pragma omp parallel for if(N > N0)
  //for(point = 0; point < N; point++) strain[2].strn_fld[point] = strain_hetero[point][0];
    
  //gen_hetero_strain(0,1);
  //#pragma omp parallel for if(N > N0)
  //for(point = 0; point < N; point++) strain[3].strn_fld[point] = strain_hetero[point][0];
    
  //gen_hetero_strain(0,2);
  //#pragma omp parallel for if(N > N0)
  //for(point = 0; point < N; point++) strain[4].strn_fld[point] = strain_hetero[point][0];
    
  //gen_hetero_strain(1,2);
  //#pragma omp parallel for if(N > N0)
  //for(point = 0; point < N; point++) strain[5].strn_fld[point] = strain_hetero[point][0];
	
#pragma omp parallel for if(N > N0)
  for(point = 0; point < N; point++)
    {
      /*
      // First contribution is implicit in magnetization(hence not included)
        
      // Second contribution
      h1[point][0] -= aux4 * m1[point][0] * (1.0 - m1[point][0] * m1[point][0]);
      h2[point][0] -= aux4 * m2[point][0] * (1.0 - m2[point][0] * m2[point][0]);
      h3[point][0] -= aux4 * m3[point][0] * (1.0 - m3[point][0] * m3[point][0]);
        
      // Third contribution
      h1[point][0] += aux5 * strain[0].strn_fld[point] * m1[point][0];
      //+ aux6 * (strain[3].strn_fld[point] * m2[point][0] + strain[4].strn_fld[point] * m3[point][0]);
                     
      h2[point][0] += aux5 * strain[1].strn_fld[point] * m2[point][0];
      //+ aux6 * (strain[3].strn_fld[point] * m1[point][0] + strain[5].strn_fld[point] * m3[point][0]);
                     
      h3[point][0] += aux5 * strain[2].strn_fld[point] * m3[point][0];
      //+ aux6 * (strain[5].strn_fld[point] * m2[point][0] + strain[4].strn_fld[point] * m1[point][0]);
      */
      h1[point][0] += m1[point][0] * (aux5 * strain[0].strn_fld[point] - aux4 * (1.0 - m1[point][0] * m1[point][0]));
      h2[point][0] += m2[point][0] * (aux5 * strain[1].strn_fld[point] - aux4 * (1.0 - m2x[point][0] * m2[point][0]));
      h3[point][0] += m3[point][0] * (aux5 * strain[2].strn_fld[point] - aux4 * (1.0 - m3[point][0] * m3[point][0]));
    }
}
