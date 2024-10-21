#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>

#include "constants.h"
#include "voronoi.c"
#include "initial.c"
#include "laplacian.c"

double g(double phi) //interpolation function (g(eta))
{
	double value = phi * phi * (3.0 - 2.0 * phi);
	return value;
}

void TimeEvolution(void)
{
	double* L_save;
	
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

	InitialCondition();
	
	uint64_t t = 0;
	while(t < iterations) //main time loop
	{
		Laplacian();
		
		for(i = 0; i < H; i++)
		{
			for(j = 0; j < W; j++)
			{
				point = i * W + j;
				Parameters* pf = &phase_field[point];
				
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
		
		for(point = 0; point < N; point++) //update loop
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
		
		for(point = 0; point < N; point++) //revealing loop
		{
			for(grain = 0; grain < max_eta; grain++) eta_values[grain] = phase_field[point].eta[grain];
			
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
		if(rem == 0) // writing loop
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
	free(phase_field); free(function);
}

int main()
{
	TimeEvolution();
	return 0;
}