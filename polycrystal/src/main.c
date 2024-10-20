#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>

#include "constants.h"
#include "Laplacian.c"
#include "saveImage.c"
#include "extract_minimum.c"
#include "Voronoi.c"
#include "initial_condition.c"

void PolycrystallineMatrix(int iterations){
	
	orientation* grid_point = malloc(N * sizeof(orientation));
	double* Matrix = malloc(N * sizeof(double));
	double* update_phi = malloc(N * sizeof(double));
	double* eta_layer = malloc(N * sizeof(double));
	double* dot_product = malloc(N * sizeof(double));
	double* L_phi = malloc(N * sizeof(double));
	
	if (update_phi == NULL){fprintf(stderr, "Memory allocation failed while defining 'update_phi' in 'PolycrystallineMatrix'\n"); exit(EXIT_FAILURE);}
	if (Matrix == NULL){fprintf(stderr, "Memory allocation failed while defining 'Matrix' in 'PolycrystallineMatrix'\n"); exit(EXIT_FAILURE);}
	if (eta_layer == NULL){fprintf(stderr, "Memory allocation failed while defining 'eta_layer' in 'PolycrystallineMatrix'\n"); exit(EXIT_FAILURE);}
	if (dot_product == NULL){ fprintf(stderr, "Memory allocation failed while defining 'dot_product' in 'PolycrystallineMatrix'\n"); exit(EXIT_FAILURE);}
	if (L_phi == NULL){fprintf(stderr, "Memory allocation failed while defining 'L_phi' in 'PolycrystallineMatrix'\n"); exit(EXIT_FAILURE);}
	
	InitialCondition(grid_point, Matrix);
	
	for(i = 0; i < N; i++){ dot_product[i] = grid_point[i].boundary_reveal; }
	Array2Image(dot_product, H , W , 0 , iterations/save_after, save_location);
	
	char SizeInfo[] = "./AvgGrainSize.csv";
	FILE* arrayCSV = fopen(SizeInfo,"w");
	if(!arrayCSV){perror("Failed to open file"); exit(EXIT_FAILURE);}
	fprintf(arrayCSV,"time step, avg grain size, grains left\n");
	fclose(arrayCSV);
	
	int t = 0;
	while(t < iterations){
		
		for(j = 0; j < num_grains; j++){
			
			for(i = 0; i < N; i++){ eta_layer[i] = grid_point[i].eta[j]; }
			
			Laplacian(eta_layer, L_phi);
			
			for(i = 0; i < N; i++){
				update_phi[i] = eta_layer[i] + tick * interface_mobility * (eta_layer[i] + (eta_layer[i] * eta_layer[i] * eta_layer[i]) - 2.0 * eta_layer[i] * grid_point[i].boundary_reveal + 2.0 * L_phi[i]);
				grid_point[i].eta[j] = update_phi[i];
			}
		}
		
		for(i = 0; i < N; i++){
			double sum = 0;
			for(j = 0 ; j < num_grains; j++){
				sum += (grid_point[i].eta[j]) * (grid_point[i].eta[j]);
			}
			grid_point[i].boundary_reveal = sum;
			dot_product[i] = 1.0 - grid_point[i].boundary_reveal;
		}
		
		if((t+1)%save_after == 0){ Array2Image(dot_product, H , W , (t+1)/save_after , iterations/save_after , save_location); }
		
		double grain_size = 0.0;
		int remaining_grains = num_grains;
		
		for(j = 0; j < num_grains; j++){
			int sum = 0;
			for(i = 0; i < N; i++){
				if(grid_point[i].eta[j] > 0.5) sum += 1;
			}
			if(sum == 0){ remaining_grains -= 1; }
			grain_size += sqrt(sum);
		}
		
		double avg_size = 2 * K * (grain_size/remaining_grains);
		
		double r_sq = avg_size * avg_size;
		
		arrayCSV = fopen(SizeInfo, "a");
        if (!arrayCSV) { perror("Failed to open file"); exit(EXIT_FAILURE); }
        fprintf(arrayCSV, "%le,%le,%d\n", t * tick, r_sq, remaining_grains);
        fclose(arrayCSV);
		
		t++;
	}
	
	fclose(arrayCSV);
	
	free(update_phi);
	free(eta_layer);
    free(grid_point);
	free(L_phi);
	free(dot_product);
}

int main(){
	
	PolycrystallineMatrix(5000);
	return 0;
}

// command to create animation: ffmpeg -framerate 10 -i polycrystal%03d.pgm -c:v libx264 -r 30 -pix_fmt yuv420p polycrystal.mp4