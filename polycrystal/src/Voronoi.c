double* Voronoi(int num_cells){
	
	int cell_centers[2 * num_cells];
	
	srand(time(NULL));
	
	int elements = 0;
	while(elements < num_cells){
		
		int rand_x = rand() % H + 1;
		int rand_y = rand() % W + 1;
		cell_centers[elements] = rand_x;
		cell_centers[elements + num_cells] = rand_y;

		elements++;
	}
	
	double* tessellation = malloc(N*sizeof(double));
	if(tessellation == NULL){ fprintf(stderr, "Memory allocation failed for Voronoi tessellation"); exit(EXIT_FAILURE); }
	
	for(int i = 0; i < N; i++){
		int x = i/W;
		int y = i%W;
		double distances[num_cells];
		for(int j = 0; j < num_cells;j++){ distances[j] = sqrt((x - cell_centers[j] ) * (x - cell_centers[j] ) + (y - cell_centers[num_cells+j] ) * (y - cell_centers[num_cells+j])); }

		MinResult result = minimumArray(distances, num_cells);
		double min_value = result.value;
		int min_index = result.index;

		for(int k = 0; k < num_cells;k++){
			if(k == min_index){tessellation[i] = k;}
		}
	}
	
	return tessellation;
}
