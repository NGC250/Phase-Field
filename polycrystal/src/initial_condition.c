void InitialCondition(orientation* grid_point, double* Matrix){
	
	long i; int j;

	for(int i = 0; i < N; i++){
		grid_point[i].boundary_reveal = 0;
	}

	Matrix = Voronoi(num_grains);

	for(i = 0; i < N; i++){
		for(j = 0; j < num_grains; j++){
			if(j == Matrix[i]){ grid_point[i].eta[j] = 1.0; }
			else{ grid_point[i].eta[j] = 0.0; }
		}
	}
}