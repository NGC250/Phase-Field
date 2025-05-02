void generate_Xj(void)
{
  double temp , zeta_j , fourier_eps0[2];
    
  for(point = 0; point < N; point++)
    {
      for(p = 0; p < 3; p++)
        {            
	  X[point].comp[p][0] = 0.0;
	  X[point].comp[p][1] = 0.0;
            
	  for(q = 0; q < 3; q++)
            {
	      zeta_j = (q == 2) ? 0.0 : ((q == 0) ? dist_fx[point] : dist_fy[point]);
                
	      for(r = 0; r < 3; r++)
                {
		  for(s = 0; s < 3; s++)
                    {
		      if(r == s) 
                        {                            
			  fourier_eps0[0] = SFS[r].eps0[point][0];
			  fourier_eps0[1] = SFS[r].eps0[point][1];
                        }
		      else
                        {
			  if(r + s == 1) //epsilon 12
                            {   
			      fourier_eps0[0] = SFS[3].eps0[point][0];
			      fourier_eps0[1] = SFS[3].eps0[point][1];
                            }
                            
			  else if(r + s == 2) //epsilon 13
                            {   
			      fourier_eps0[0] = SFS[4].eps0[point][0];
			      fourier_eps0[1] = SFS[4].eps0[point][1];
                            } 
                            
			  else //epsilon 23
                            {   
			      fourier_eps0[0] = SFS[5].eps0[point][0];
			      fourier_eps0[1] = SFS[5].eps0[point][1];
                            }
                        }
                        
		      temp = Stiffness[p][q][r][s] * zeta_j;
                        
		      X[point].comp[p][0] += temp * fourier_eps0[0];
		      X[point].comp[p][1] += temp * fourier_eps0[1];
                    }
                }
            }
        }
    }
}
