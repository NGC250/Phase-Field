#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>

//Simulation
#define step (1e-4)
#define tick (1e-3)
#define iterations 250
#define save_after 1

//system
#define H 100
#define W 100
#define N (H*W)
#define num_eta1 2
#define num_eta2 2
#define num_eta_random 6
#define num_grains (num_eta1 + num_eta2 + num_eta_random) //2(0,35,phi2) + 2(180,35,phi2) + 6 random
#define max_eta (2 + num_eta_random)
#define random_eta_index_start (num_grains - max_eta)
#define interacting_members 2.0
#define diff_area (1/(step * step))


// Physical Parameters
#define pi 3.141592653589793238
#define field_strength (13.5e6)
#define delta_chi (11.8e-6)
#define mu0 (1.26e-6)

#define bulk_energy_coeff (0.2e6)
#define grad_energy_coeff (1.46e-6)
#define mobility (0.7511)
#define interface_width (6e-4)
#define omega (1.0/(mobility * interface_width * interface_width))
#define interfacial_energy 0.30


uint64_t point , grain , other , i , j;

double cos_gamma2[max_eta] = {(0.8473291852) , (0.002739052) , (0.7191855734) , (0.2808144266) , (0.2808144266) , (0.0) , (0.7191855734) , (0.0)};

typedef struct
{
	double eta[max_eta];
	double dEta_dt[max_eta];
	double boundary_reveal;
}
Parameters;

typedef struct
{
	double laplacian[max_eta];
}
Functions;

typedef struct
{
	uint64_t max_index;
	uint64_t min_index;
	double max_value;
	double min_value;
}
Extract;

Parameters* phase_field;
Functions* function;

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
		
		int rand_x = rand() % H;
		int rand_y = rand() % W;
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
	
	for(point = 0; point < N; point++){
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
		
		// double sum = 0;
		// for(uint64_t elem = 0; elem < max_eta - 1; elem++)
		// {
			// sum += eta_values[elem] * eta_values[elem + 1];
		// }
		// phase_field[point].boundary_reveal = 1.0 - sum;
		
		extremes = In_Array(eta_values , max_eta);

		if(extremes.max_value >= 0.5)
		{
			if(extremes.max_index == 0) phase_field[point].boundary_reveal = 1.0;
			else if(extremes.max_index == 1) phase_field[point].boundary_reveal = 0.5;
			else phase_field[point].boundary_reveal = 0.25;
		}
		else phase_field[point].boundary_reveal = 0.0;
	}
	
	FILE *image0 , *file0;
	char filename0[256] , filename1[256];
	sprintf(filename0 , "../Data/images/structure0.pgm");
	sprintf(filename1 , "../Data/matrix/reveal_0.dat");
	
	image0 = fopen(filename0 , "w");
	file0 = fopen(filename1 , "w");
	
	fprintf(image0 , "P2\n%d %d\n255\n", W, H);
	
	for(i = 0; i < H; i++)
	{
		for(j = 0; j < W; j++)
		{
			point = i * W + j;	
			fprintf(image0 , "%d " , (uint8_t)(phase_field[point].boundary_reveal * 255));
			fprintf(file0 , "%le," , phase_field[point].boundary_reveal);
		}
		fprintf(image0 , "\n");
		fprintf(file0 , "\n");
	}
	fclose(image0);
	fclose(file0);
	
	free(grain_index);
}

void Laplacian(void)
{
	double up , down , left , right , center;
	
	for(i = 0; i < H; i++)
	{
		for(j = 0; j < W; j++)
		{
			for(grain = 0; grain < max_eta; grain++)
			{
				if(i == 0) up = phase_field[(H-1) * W + j].eta[grain];
				else up = phase_field[(i-1) * W + j].eta[grain];
				
				if(i == H-1) down = phase_field[j].eta[grain];
				else down = phase_field[(i+1) * W + j].eta[grain];
				
				if(j == 0) left = phase_field[i * W + (W-1)].eta[grain];
				else left = phase_field[i * W + (j-1)].eta[grain];
				
				if(j == W-1) right = phase_field[i * W].eta[grain];
				else right = phase_field[i * W + (j+1)].eta[grain];
			
				center = phase_field[i * W + j].eta[grain];
				
				function[i * W + j].laplacian[grain] = (up + down - 2.0 * center) + (left + right - 2.0 * center);
			}			
		}
	}
	// printf("%le , ",function[2500].laplacian[0]);
}

double g(double phi) //interpolation function (g(phi))
{
	double value = phi * phi * (3.0 - 2.0 * phi);
	return value;
}

void TimeEvolution(void)
{
	double* L_save ;
	
	uint64_t* save_index;
	double* save_value;
	
	if((phase_field = (Parameters *)malloc(N * sizeof(Parameters))) == NULL)
	{
		printf("phase_field variable could not be allocated!"); exit(0);
	}
	if((function = (Functions *)malloc(N * sizeof(Functions))) == NULL)
	{
		printf("function variable could not be allocated!"); exit(0);
	}
	if((L_save = (double *)malloc(N * sizeof(double))) == NULL)
	{
		printf("L_save could not be allocated!"); exit(0);
	}

	InitialCondition();
	
	uint64_t t = 0;
	while(t < iterations)
	{
		Laplacian();
		
		for(i = 0; i < H; i++)
		{
			for(j = 0; j < W; j++)
			{
				point = i * W + j;
				Parameters* pf = &phase_field[point];
				
				L_save[point] = function[point].laplacian[0];
				
				double lambda = 0.0 , sum1 = 0.0;
				
				for(grain = 0; grain < max_eta; grain++) sum1 += g(pf->eta[grain]);
				sum1 = sum1 * sum1;
				
				for(grain = 0; grain < max_eta; grain++)
				{
					
					double bulk_term = 0.0 , grad_term = 0.0 , mag_term = 0.0;
					
					for(other = 0; other < max_eta; other++)
					{
						if(other != grain)
						{
							bulk_term += pf->eta[other] * function[point].laplacian[grain] - pf->eta[grain] * function[point].laplacian[other];
							grad_term += pf->eta[other] - pf->eta[grain];
							
							double sum2 = 0.0;
							for(uint64_t k = 0; k < max_eta; k++)
							{
								if(k != grain) sum2 += g(pf->eta[k]) * (cos_gamma2[k] - cos_gamma2[other]);
							}
							mag_term += g(pf->eta[other]) * sum2;
						}
					}
					
					double sum3 = 0.0;
					for(other = 0; other < max_eta; other++)
					{
						if(other != grain) sum3 += g(pf->eta[other]) * (cos_gamma2[other] - cos_gamma2[grain]);
					}
					
					mag_term = 0.5 * mu0 * field_strength * field_strength * delta_chi * (g(pf->eta[grain]) * sum3 - mag_term) / (sum1) ;
					
					pf->dEta_dt[grain] = ( (interfacial_energy/omega) * (diff_area * bulk_term - (8.0/(pi*pi * interface_width*interface_width)) * grad_term) ) + mag_term;
					
					lambda += pf->dEta_dt[grain];
				}
				
				lambda /= max_eta;
				
				for(grain = 0; grain < max_eta; grain++) pf->dEta_dt[grain] -= lambda;
			}
		}
		
		for(point = 0; point < N; point++)
		{
			Parameters* pf = &phase_field[point];
			
			for(grain = 0; grain < max_eta; grain++)
			{
				pf->eta[grain] = pf->eta[grain] + pf->dEta_dt[grain] * tick;
				pf->eta[grain] = (pf->eta[grain] > 1.0) ? 1.0 : ((pf->eta[grain] < 0.0) ? 0.0 : pf->eta[grain]);
			}
		}
		
		Extract extremes;
		double eta_values[max_eta];
		
		for(point = 0; point < N; point++)
		{
			for(grain = 0; grain < max_eta; grain++) eta_values[grain] = phase_field[point].eta[grain];
			
			// double sum = 0;
			// for(grain = 0; grain < max_eta; grain++)
			// {
				// for(other = grain + 1; other < max_eta - 1; other++)
				// {
					// sum += eta_values[other] * eta_values[grain];
				// }
			// }
			// phase_field[point].boundary_reveal = 1.0 - sum;
			
			extremes = In_Array(eta_values , max_eta);
			
			if(extremes.max_value == 1.0)
			{
				if(extremes.max_index == 0) phase_field[point].boundary_reveal = 1.0;
				else if(extremes.max_index == 1) phase_field[point].boundary_reveal = 0.5;
				else phase_field[point].boundary_reveal = 0.25;
			}
			else phase_field[point].boundary_reveal = 0.0;
		}
		
		t++;

		uint64_t rem = t%save_after;
		if(rem == 0)
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
	}
	
	free(phase_field); free(function); free(L_save);
}

int main()
{
	TimeEvolution();
	return 0;
}