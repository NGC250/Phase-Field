#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <time.h>
#include <string.h>

#include "constants.h"
#include "Voronoi.c"
#include "initial.c"
#include "Laplacian.c"

void PolycrystallineMatrix(uint64_t iterations){
	
	if((phase_field = (Parameters *)malloc(N * sizeof(Parameters))) == NULL)
	{
		printf("phase_field variable could not be allocated!"); exit(0);
	}
	if((function = (Functions *)malloc(N * sizeof(Functions))) == NULL)
	{
		printf("function variable could not be allocated!"); exit(0);
	}
	
	InitialCondition();
	
	const char* SizeInfo = "./AvgGrainSize.csv";
	FILE* arrayCSV = fopen(SizeInfo,"w");
	
	fprintf(arrayCSV,"time step, avg grain size, grains left\n");
	
	uint64_t t = 0;
	while(t < iterations)
	{
		Laplacian();

		for(grain = 0; grain < num_grains; grain++)
		{	
			for(point = 0; point < N; point++)
			{
				double eta_layer = phase_field[point].eta[grain];
				
				phase_field[point].update_eta[grain] = eta_layer + tick * interface_mobility * (eta_layer + (eta_layer * eta_layer * eta_layer) - 
				2.0 * eta_layer * phase_field[point].boundary_reveal + 2.0 * function[point].laplacian[grain]);				
			}
		}
		
		for(point = 0; point < N; point++)
		{
			for(grain = 0; grain < num_grains; grain++) phase_field[point].eta[grain] = phase_field[point].update_eta[grain];
		}
		
		for(point = 0; point < N; point++)
		{
			double sum = 0.0;
			for(grain = 0 ; grain < num_grains; grain++)
			{
				sum += (phase_field[point].eta[grain]) * (phase_field[point].eta[grain]);
			}
			phase_field[point].boundary_reveal = sum;
		}
		
		t++;

		if(t%save_after == 0) // writing loop
		{
			FILE *imagefile;
			char filename[256];
			sprintf(filename , "../Data/images/structure%lu.pgm",t);
			
			imagefile = fopen(filename , "w");
			
			fprintf(imagefile , "P2\n%d %d\n255\n", W, H);
			
			for(i = 0; i < H; i++)
			{
				for(j = 0; j < W; j++)
				{
					point = i * W + j;	
					fprintf(imagefile , "%d " , (uint8_t)(phase_field[point].boundary_reveal * 255));
				}
				fprintf(imagefile , "\n");
			}
			fclose(imagefile);
		}
		
		double grain_size = 0.0;
		int remaining_grains = num_grains;
		
		for(grain = 0; grain < num_grains; grain++)
		{
			uint64_t sum = 0;
			for(point = 0; point < N; point++)
			{
				if(phase_field[point].eta[grain] > 0.5) sum += 1;
			}
			if(sum == 0) remaining_grains -= 1;
			grain_size += sqrt(sum);
		}
		
		double avg_size = 2.0 * K * (grain_size/remaining_grains);
		
		double r_sq = avg_size * avg_size;
		
		fprintf(arrayCSV, "%le,%le,%d\n", t * tick, r_sq, remaining_grains);
	}
	
	fclose(arrayCSV);
	
    free(phase_field);
}

int main()
{
	PolycrystallineMatrix(5000);
	return 0;
}