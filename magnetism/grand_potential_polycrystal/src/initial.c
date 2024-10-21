void InitialCondition(void)
{
	uint64_t* grain_index = Voronoi(num_grains);
	
	uint64_t gi;
	
	for(point = 0; point < N; point++)
	{
		for(grain = 0; grain < max_eta; grain++) phase_field[point].eta[grain] = 0.0;
	}
	
	for(point = 0; point < N; point++)
	{
		gi = grain_index[point];
		
		if(gi < num_eta1) phase_field[point].eta[0] = 1.0;
		else if((gi >= num_eta1) && (gi < num_eta1 + num_eta2)) phase_field[point].eta[1] = 1.0;
		
		else if(gi >= num_eta1 + num_eta2) 
		{
			phase_field[point].eta[gi - random_eta_index_start] = 1.0;
		}
	}
	
	Extract extremes;
	double eta_values[max_eta];
	
	for(point = 0; point < N; point++)
	{
		for(grain = 0; grain < max_eta; grain++)
		{
			eta_values[grain] = phase_field[point].eta[grain];
		}
		
		extremes = In_Array(eta_values , max_eta);

		if(extremes.max_value >= 0.5)
		{
			if(extremes.max_index == 0) phase_field[point].boundary_reveal = 1.0;
			else if(extremes.max_index == 1) phase_field[point].boundary_reveal = 0.5;
			else phase_field[point].boundary_reveal = 0.25;
		}
		else phase_field[point].boundary_reveal = 0.0;
	}
	
	FILE *image0;
	char filename0[256];
	sprintf(filename0 , "../Data/images/structure0.pgm");
	
	image0 = fopen(filename0 , "w");
	
	fprintf(image0 , "P2\n%d %d\n255\n", W, H);
	
	for(i = 0; i < H; i++)
	{
		for(j = 0; j < W; j++)
		{
			point = i * W + j;	
			fprintf(image0 , "%d " , (uint8_t)(phase_field[point].boundary_reveal * 255));
		}
		fprintf(image0 , "\n");
	}
	fclose(image0);
	
	free(grain_index);
}