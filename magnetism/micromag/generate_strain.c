void gen_homog_strain(void)
{
  for(p = 0; p < 3; p++)
    {
      for(q = 0; q < 3; q++)
        {
	  m_avg.mimj[p][q] = 0.0;
            
	  for(point = 0; point < N; point++)
            {   
	      mzn[0] = m1[point][0];
	      mzn[1] = m2[point][0];
	      mzn[2] = m3[point][0];
                
	      m_avg.mimj[p][q] += mzn[p] * mzn[q];
            }
            
	  m_avg.mimj[p][q] *= norm;
            
	  Strain_homog[p][q] = Mstr_coeff[p][q] * (m_avg.mimj[p][q] - Delta[p][q] / 3.0);
        }
    }
}

void gen_hetero_strain(uint8_t comp1 , uint8_t comp2)
{
  double zeta_comp1 , zeta_comp2 , u_ij[2] , u_ji[2];
    
  for(point = 0; point < N; point++)
    {
      u_ij[0] = 0.0;
      u_ij[1] = 0.0;
      u_ji[0] = 0.0;
      u_ji[1] = 0.0;

      for(p = 0; p < 3; p++)
        {
	  u_ij[0] += X[point].comp[p][0] * N_ij[point].cell[comp1][p];
	  u_ij[1] += X[point].comp[p][1] * N_ij[point].cell[comp1][p];
        }
        
      for(q = 0; q < 3; q++)
        {
	  u_ji[0] += X[point].comp[q][0] * N_ij[point].cell[comp2][q];
	  u_ji[1] += X[point].comp[q][1] * N_ij[point].cell[comp2][q];
        }
        
      zeta_comp1 = (comp1 == 2) ? 0.0 : ((comp1 == 0) ? dist_fx[point] : dist_fy[point]);
      zeta_comp2 = (comp2 == 2) ? 0.0 : ((comp2 == 0) ? dist_fx[point] : dist_fy[point]);
        
      strain_hetero[point][0] = 0.5 * inv_Det_K[point] * (u_ij[0] * zeta_comp2 + u_ji[0] * zeta_comp1);
      strain_hetero[point][1] = 0.5 * inv_Det_K[point] * (u_ij[1] * zeta_comp2 + u_ji[1] * zeta_comp1);
    }
    
  run_fft(strain_hetero , strain_hetero , -1);
    
#pragma omp parallel for if(N > N0)
  for(point = 0; point < N; point++) strain_hetero[point][0] *= norm;
}
