void InitialCondition(void)
{
	uint64_t* grain_index = Voronoi(num_grains);
	
	for(point = 0; point < N; point++)
	{
		phase_field[point].boundary_reveal = 0.0;
		for(grain = 0; grain < num_grains; grain++)
		{
			if(grain_index[point] == grain) phase_field[point].eta[grain] = 1.0;
			else phase_field[point].eta[grain] = 0.0;
			
			phase_field[point].boundary_reveal += phase_field[point].eta[grain] * phase_field[point].eta[grain];
		}
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