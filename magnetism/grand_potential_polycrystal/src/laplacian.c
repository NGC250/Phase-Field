void Laplacian(void)
{
	double up , down , left , right , center;
	
	for(i = 0; i < H; i++)
	{
		for(j = 0; j < W; j++)
		{
			for(grain = 0; grain < max_eta; grain++)
			{
				if(i == 0) up = phase_field[(H-1) * W + j].eta[grain];
				else up = phase_field[(i-1) * W + j].eta[grain];
				
				if(i == H-1) down = phase_field[j].eta[grain];
				else down = phase_field[(i+1) * W + j].eta[grain];
				
				if(j == 0) left = phase_field[i * W + (W-1)].eta[grain];
				else left = phase_field[i * W + (j-1)].eta[grain];
				
				if(j == W-1) right = phase_field[i * W].eta[grain];
				else right = phase_field[i * W + (j+1)].eta[grain];
			
				center = phase_field[i * W + j].eta[grain];
				
				function[i * W + j].laplacian[grain] = (up + down - 2.0 * center) + (left + right - 2.0 * center);
			}			
		}
	}
}