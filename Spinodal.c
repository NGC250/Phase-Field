#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>

#define H 100	// These represent the height(H) and the
#define W 100	// width(W) of your system. Preferentially,
#define N H*W	// the H and W should be equal.


/*This function creates and saves images in from the arrays. The function saves files in a .pgm format for
faster saving. To open a .pgm image, use GIMP or ImageJ or ImageMagick OR you can directly create the animation
from the saved images for visualisation. Rules to create animation are at the end of the file */
void ImagefromArray(double array[N], int index, int save_period, const char* filename, const char* location) {
    char files[128];
    char buffer[128];

    snprintf(buffer, sizeof(buffer), "%s%s", location, filename);

    int num_digits = snprintf(NULL, 0, "%d", save_period);

    if (snprintf(files, sizeof(files), "%s%0*d.pgm", buffer, num_digits, index) >= sizeof(files)) {
        fprintf(stderr, "Error: File path buffer too small\n");
        exit(EXIT_FAILURE);
    }

    FILE* pgmimg = fopen(files, "wb");
    if (!pgmimg) {
        perror("Failed to open file");
        exit(EXIT_FAILURE);
    }
    fprintf(pgmimg, "P2\n%d %d\n255\n", W, H);

    for (int i = 0; i < H; i++) {
        for (int j = 0; j < W; j++) {
            int temp = (int)(array[i * W + j] * 255);
            fprintf(pgmimg, "%d ", temp);
        }
        fprintf(pgmimg, "\n");
    }
    fclose(pgmimg);
}



/* This function is used to calculate the Laplacian of an array */
double* Laplacian(double array1D[N], int step){

	double* laplacian = (double*)malloc(N * sizeof(double));
	if (laplacian == NULL) {
        fprintf(stderr, "Memory allocation failed in Laplacian\n");
        exit(EXIT_FAILURE);
    }
	double up,down,left,right,center;

	for(int i = 0; i < N; i++){

		if(i >= H){up = array1D[i-H];}
		else{up = array1D[N-H+i];}

		if(N-H > i){down = array1D[i+H];}
		else{down = array1D[i+H-N];}

		if(i%W == 0){left = array1D[i+W-1];}
		else{left = array1D[i-1];}

		if((i+1)%W == 0){right = array1D[i+1-W];}
		else{right = array1D[i+1];}

		center = array1D[i];

		laplacian[i] = (right+left-2*center)/(step*step) + (up+down-2*center)/(step*step);}

	return laplacian;
}


/* This function adds compositional noise to initiate the spinodal decomposition */
double* Noise(double array1D[N], float noise_intensity){

	srand(time(NULL));
	for(int i = 0; i < N; i++){
		double random_num = (double)rand()/RAND_MAX;
		array1D[i] = array1D[i] + noise_intensity * (2 * random_num - 1);}

	return array1D;
}


/* The main function where the physics is programmed.
NOTE: The double well potential is taken to be [barrier_height * C^2 * (1 - C)^2] where C is the composition at a point.*/
double* CahnHilliard(double phi1D[N], float diffusion_coefficient, float barrier_height, float gradient_coefficient, float space_step, float time_step, int iterations, const char* save_location){

	double* update_phi = (double*)malloc(N*sizeof(double));
	double* dF_dPhi = (double*)malloc(N * sizeof(double));
	if (update_phi == NULL || dF_dPhi == NULL) {
        fprintf(stderr, "Memory allocation failed in CahnHilliard\n");
        exit(EXIT_FAILURE);
    }

	int t = 0;
	int image_index = 0;
	while(t < iterations){

		double* L_phi = Laplacian(phi1D,space_step);

		for(int j = 0; j < N; j++){
			double m = phi1D[j];
			dF_dPhi[j] = 2 * barrier_height * m * (1 - m) * (1-2*m) - gradient_coefficient * L_phi[j];} //Change the double well potential if you wish to.

		double* L_dF_dPhi = Laplacian(dF_dPhi,space_step);

		for(int k = 0; k < N; k++){

			update_phi[k] = phi1D[k] + time_step * diffusion_coefficient * L_dF_dPhi[k];}

		for(int l = 0; l < N; l++){phi1D[l] = update_phi[l];}
		
		int save_period = 100;
		float num_digits = t % save_period;
		if(num_digits == 0 || t == iterations-1){
            ImagefromArray(phi1D, image_index, save_period, "spinodal", save_location); //Change filename if you wish to but also change the name in .bat or .sh file to create animation.
            image_index++;
        }
		
		t++;
	}
	free(dF_dPhi);

	return phi1D;
}


int main(){

	double* matrix1D = (double*)malloc(N*sizeof(double));
	if (matrix1D == NULL) {
        fprintf(stderr, "Memory allocation failed in main\n");
        exit(EXIT_FAILURE);
    }

	for(int i = 0; i <N; i++){matrix1D[i] = 0.5;}
	double noise_intensity = 0.05; // Change noise_intensity accordingly
	double* phi = Noise(matrix1D, noise_intensity);
	
	const char* locationCH = "img/spinodal/"; // Change the location according to you.
	
	float diffusion_coefficient = 1, barrier_height = 1, gradient_coefficient = 2, time_step = 0.005; //Change accordingly
	int iterations = 10000; //Change accordingly
	double* finalCH = CahnHilliard(phi, diffusion_coefficient, barrier_height, gradient_coefficient, 1, time_step, iterations, locationCH);

	return 0;
}

/* To create animation from saved images, open the terminal or command prompt in the directory of saved images.
   For windows and for linux: ffmpeg -framerate 10 -i spinodal%02d.pgm -c:v libx264 -r 30 -pix_fmt yuv420p animation.mp4
*/