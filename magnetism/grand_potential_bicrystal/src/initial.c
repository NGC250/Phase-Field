void InitialCondition(void)
{
	uint64_t* matrix = malloc(N * sizeof(uint64_t));
	
	matrix = Voronoi(num_grains);
	
	for(i = 0; i < H; i++)
	{
		for(j = 0; j < W; j++)
		{
			point = i * W + j;
			
			for(grain = 0; grain < num_grains; grain++)
			{
				if(grain == matrix[point]) phase_field[point].vol_frac[grain] = 1.0;
				else phase_field[point].vol_frac[grain] = 0.0;
				
				if(i == 0) buffer_layer[grain].up_buffer[j] = phase_field[j].vol_frac[grain];
				
				if(i == H-1) buffer_layer[grain].down_buffer[j] = phase_field[N - W + j].vol_frac[grain];
				
				if(j == 0) buffer_layer[grain].left_buffer[i] = phase_field[point].vol_frac[grain];
				
				if(j == W-1) buffer_layer[grain].right_buffer[i] = phase_field[point].vol_frac[grain];
			}
		}
	}

	FILE *imagefile0;
	char *filename;
	filename = "../Data/images/structure0.pgm";
	
	imagefile0 = fopen(filename , "w");
	fprintf(imagefile0 , "P2\n%d %d\n255\n", W, H);
	
	for(i = 0; i < H; i++)
	{
		for(j = 0; j < W; j++)
		{
			fprintf(imagefile0 , "%d " , (uint8_t)(phase_field[i * W + j].boundary_reveal * 255));
		}
		fprintf(imagefile0 ,"\n");
	}
	
	fclose(imagefile0);
	
	for(point = 0; point < N; point++)
	{
		Parameters* pf = &phase_field[point];
		if((pf->vol_frac[0]) * (pf->vol_frac[1]) == 0)
		{
			if(pf->vol_frac[0] > pf->vol_frac[1]) pf->boundary_reveal = 0.5;
			else pf->boundary_reveal = 1;
		}
		else pf->boundary_reveal = 0;
	}
	free(matrix);
}