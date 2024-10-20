MinResult minimumArray(double* array, uint64_t size)
{
	MinResult result;
	result.value = array[0];
	result.index = 0;

	for(i = 1; i < size; i++)
	{
		if (array[i] < result.value)
		{
			result.value = array[i];
			result.index = i;
		}
	}
	return result;
}

uint64_t* Voronoi(uint64_t num_cells){
	
	uint64_t cell_centers[2 * num_cells];
	
	srand(time(NULL));
	
	uint64_t elements = 0;
	while(elements < num_cells){
		
		int rand_x = rand() % H;
		int rand_y = rand() % W;
		cell_centers[elements] = rand_x;
		cell_centers[elements + num_cells] = rand_y;

		elements++;
	}
	
	uint64_t* tessellation = malloc(N*sizeof(uint64_t));
if(tessellation == NULL){ fprintf(stderr, "Memory allocation failed for Voronoi tessellation"); exit(EXIT_FAILURE); }
	
	
	for(point = 0; point < N; point++){
		double x = point/W;
		double y = point%W;
		double distances[num_cells];
		for(grain = 0; grain < num_cells; grain++) distances[grain] = (x - cell_centers[grain] ) * (x - cell_centers[grain] ) + (y - cell_centers[num_cells + grain] ) * (y - cell_centers[num_cells + grain]);

		MinResult result = minimumArray(distances, num_cells);
		double min_value = result.value;
		uint64_t min_index = result.index;

		tessellation[point] = min_index;
	}

	return tessellation;
}
