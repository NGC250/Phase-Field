void InitialCondition(void)
{    
  double x_comp , y_comp , z_comp , normal;
    
  for(point = 0; point < N; point++)
    {
      do
	{
	  x_comp = 2.0 * ( (double)(rand()) / (double)(RAND_MAX) ) - 1.0;
	  y_comp = 2.0 * ( (double)(rand()) / (double)(RAND_MAX) ) - 1.0;
	  z_comp = 2.0 * ( (double)(rand()) / (double)(RAND_MAX) ) - 1.0;

	  normal = sqrt(x_comp * x_comp + y_comp * y_comp + z_comp * z_comp);
	}
      while(normal == 0.0);

      x_comp /= normal;
      y_comp /= normal;
      z_comp /= normal;       

      m1[point][0] = x_comp;
      m2[point][0] = y_comp;
      m3[point][0] = z_comp;
      
      m1[point][1] = 0.0;
      m2[point][1] = 0.0;
      m3[point][1] = 0.0;
    }
    
  savedata(0);

  for(k = 0; k < H; k++)
    {
      zeta_x = (double) ( (k <= H/2) ? k : k - H );
      zeta_x *= 2.0 * pi/((double)(H));

      for(l = 0; l < W; l++)
	{
	  zeta_y = (double) ( (l <= W/2) ? l : l - W );
	  zeta_y *= 2.0 * pi/((double)(W));

	  point = k * W + l;

	  dist_fx[point] = zeta_x;
	  dist_fy[point] = zeta_y;
            
	  zeta_mag[point] = zeta_x * zeta_x + zeta_y * zeta_y;
	}
    }

  double zeta1_sq , zeta2_sq , Det_K;
    
  inv_Det_K[0] = 0.0;

  for(point = 1; point < N; point++)
    {
	zeta1_sq = dist_fx[point] * dist_fx[point];
	zeta2_sq = dist_fy[point] * dist_fy[point];
        
	Det_K = zeta_mag[point] * (theta1 * zeta_mag[point] * zeta_mag[point] + theta2 * zeta1_sq * zeta2_sq);
        
	inv_Det_K[point] = 1.0/Det_K;
  }
}
