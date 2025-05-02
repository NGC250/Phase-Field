void generate_N(void)
{
  for(point = 0; point < N; point++)
    {
      double zeta1_sq = dist_fx[point] * dist_fx[point];
      double zeta2_sq = dist_fy[point] * dist_fy[point];
        
      double sh_fc = c44 * zeta_mag[point]; // same dimensions as shear force
        
      N_ij[point].cell[0][0] = sh_fc * (sh_fc + aux7 * zeta2_sq);
      N_ij[point].cell[1][1] = sh_fc * (sh_fc + aux7 * zeta1_sq);
      N_ij[point].cell[2][2] = sh_fc * zeta_mag[point] * c11 + aux8 * zeta1_sq * zeta2_sq;
        
      N_ij[point].cell[0][1] = aux9 * dist_fx[point] * dist_fy[point] * zeta_mag[point];
      N_ij[point].cell[1][0] = N_ij[point].cell[0][1];
        
      N_ij[point].cell[0][2] = 0.0;
      N_ij[point].cell[2][0] = N_ij[point].cell[0][2];
        
      N_ij[point].cell[1][2] = 0.0;
      N_ij[point].cell[2][1] = N_ij[point].cell[1][2];
    }
}
