void Laplacian_array(void)
{	
    double up, down, left, right, center;

	for(grain = 0; grain < num_grains; grain++)
	{		
		for(i = 0; i < H; i++)
		{
			for(j = 0; j < W; j++)
			{
				if(i == H-1) down = buffer_layer[grain].down_buffer[j];					
				else down = phase_field[(i+1) * W + j].vol_frac[grain];
				
				if(i == 0) up = buffer_layer[grain].up_buffer[j];
				else up = phase_field[(i-1) * W + j].vol_frac[grain];
				
				if(j == 0) left = buffer_layer[grain].left_buffer[i];
				else left = phase_field[i * W + (j-1)].vol_frac[grain];
				
				if(j == W-1) right = buffer_layer[grain].right_buffer[i];
				else right = phase_field[i * W + (j+1)].vol_frac[grain];

				center = phase_field[i * W + j].vol_frac[grain];
				
				L_vol_frac[i * W + j].laplacian[grain] = (right + left - 2.0 * center) + (up + down - 2.0 * center);
			}
		}
	}
}
