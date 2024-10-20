#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>

#include "constants.h"
#include "image_from_array.c"
#include "Laplacian.c"

//NOTE: The double well potential is taken to be [barrier_height * C^2 * (1 - C)^2]
void CahnHilliard(double* phi, int iterations){

	double* update_phi = malloc(N*sizeof(double));
	double* dF_dPhi = malloc(N * sizeof(double));
	double* laplacian = malloc(N * sizeof(double));
	double* L_dF_dPhi = malloc(N * sizeof(double));

	if ((laplacian == NULL) || (dF_dPhi == NULL) || (update_phi == NULL) || (L_dF_dPhi == NULL)){
        fprintf(stderr, "Memory allocation failed in Laplacian\n");
        exit(EXIT_FAILURE);
    }
	
	int t = 0;
	while(t < iterations){

		Laplacian(phi, laplacian);

		for(int j = 0; j < N; j++){
			double m = phi[j];
			dF_dPhi[j] = 2 * barrier_height * m * (1 - m) * (1-2*m) - gradient_coefficient * laplacian[j];}

		Laplacian(dF_dPhi , L_dF_dPhi);

		for(int k = 0; k < N; k++){ update_phi[k] = phi[k] + tick * diffusion_coefficient * L_dF_dPhi[k]; }

		for(int l = 0; l < N; l++){ phi[l] = update_phi[l]; }
		
		if((t+1)%save_after == 0){ ImagefromArray(phi, H , W , t+1 , iterations , save_location); }
		
		t++;
	}
	
	free(dF_dPhi);
	free(phi);
	free(laplacian);
}

int main(){

	double* matrix = malloc(N*sizeof(double));
	if (matrix == NULL) {
        fprintf(stderr, "Memory allocation failed in main\n");
        exit(EXIT_FAILURE);
    }

	for(int i = 0; i < N; i++){ matrix[i] = 0.5; }
	
	srand(time(NULL));
	for(int i = 0; i < N; i++){
		double random_num = (double)rand()/RAND_MAX;
		matrix[i] += noise_intensity * (2 * random_num - 1);
	}
	
	int iterations = 10000;
	ImagefromArray(matrix , H , W , 0 , iterations , save_location);
	
	CahnHilliard(matrix , iterations);

	return 0;
}

/* To create animation from saved images, open the terminal or command prompt in the directory of saved images.
   For windows and for linux: ffmpeg -framerate 10 -i spinodal%03d.pgm -c:v libx264 -r 30 -pix_fmt yuv420p spinodal_animation.mp4
*/