Extract In_Array(double* array, uint64_t size)
{
	Extract extremes;
	
	extremes.min_value = array[0];
	extremes.max_value = array[0];
	extremes.min_index = 0;
	extremes.max_index = 0;

	for(i = 1; i < size; i++)
	{
		if(array[i] < extremes.min_value)
		{
			extremes.min_value = array[i];
			extremes.min_index = i;
		}
		if(array[i] > extremes.max_value)
		{
			extremes.max_value = array[i];
			extremes.max_index = i;
		}
	}
	return extremes;
}

uint64_t* Voronoi(uint64_t num_cells)
{
	uint64_t cell_centers[2 * num_cells];
	
	srand(time(NULL));
	
	uint64_t elements = 0;
	while(elements < num_cells){
		
		uint64_t rand_x = rand() % H;
		uint64_t rand_y = rand() % W;
		cell_centers[elements] = rand_x;
		cell_centers[elements + num_cells] = rand_y;

		elements++;
	}
	
	uint64_t* tessellation;
	if((tessellation = (uint64_t *)malloc(N*sizeof(uint64_t))) == NULL)
	{
		printf("Tessellation could not be allocated!"); exit(0);
	}
	
	Extract extremes;
	
	for(point = 0; point < N; point++)
	{
		double x = point/W;
		double y = point%W;
		double distances[num_cells];
		
		for(grain = 0; grain < num_cells; grain++) 
			distances[grain] = (x - cell_centers[grain] ) * (x - cell_centers[grain] ) + (y - cell_centers[num_cells + grain] ) * (y - cell_centers[num_cells + grain]);

		extremes = In_Array(distances, num_cells);

		tessellation[point] = extremes.min_index;
	}
	
	return tessellation;
}