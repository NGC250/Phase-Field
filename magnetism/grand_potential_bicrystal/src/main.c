#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#include "constants.h"
#include "voronoi.c"
#include "initial.c"
#include "laplacian.c"

void TimeEvolution(void)
{
	if((phase_field = (Parameters *)malloc(N * sizeof(Parameters))) == NULL)
	{
		printf("phase field cannot be allocated\n");
        exit(0);
	}
	if((L_vol_frac = (Laplacian *)malloc(N * sizeof(Laplacian))) == NULL)
	{
		printf("L_vol_frac could not be allocated\n");
		exit(0);
	}
	if((buffer_layer = (Buffer *)malloc(num_grains * sizeof(Buffer))) == NULL)
	{
		printf("Buffer cannot be allocated\n");
        exit(0);
	}
	
	InitialCondition();
	
	uint64_t t = 0;
	while(t < iterations) //main time loop
	{
		
		Laplacian_array();
		
		for(i = 0; i < H; i++)
		{
			for(j = 0; j < W; j++){
				
				point = i * W + j;
				
				Parameters* pf = &phase_field[point];
				
				double grad_term = diff_area * ( pf->vol_frac[1] * L_vol_frac[point].laplacian[0] - pf->vol_frac[0] * L_vol_frac[point].laplacian[1] );
				double bulk_term = -1.0 * (8.0/(pi * pi * interface_width * interface_width)) * (pf->vol_frac[1] - pf->vol_frac[0]);

				double magnetic_term_coeff = mu0 * field_strength * field_strength * delta_chi;
				double phi0 = pf->vol_frac[0];
				double interpolation = (phi0 - phi0 * phi0)*(phi0 - phi0 * phi0) * (3.0 + 4.0*phi0 - 4.0*phi0*phi0)/((1 + 4.0*phi0 - 4.0*phi0*phi0)*(1 + 4.0*phi0 - 4.0*phi0*phi0));
				double magnetic_term = -magnetic_term_coeff * interpolation;
				
				pf->update_vol_frac[0] = pf->vol_frac[0] + tick * ((interfacial_energy/omega) * (grad_term + bulk_term) + magnetic_term);
				pf->update_vol_frac[1] = 1.0 - pf->update_vol_frac[0];
			}
		}
		
		for(point = 0; point < N; point++) // updating and revealing the interafce
		{
			Parameters* pf = &phase_field[point];
			
			pf->vol_frac[0] = (pf->update_vol_frac[0] > 1.0) ? 1.0 : ((pf->update_vol_frac[0] < 0.0) ? 0.0 : pf->update_vol_frac[0]);
			pf->vol_frac[1] = (pf->update_vol_frac[1] > 1.0) ? 1.0 : ((pf->update_vol_frac[1] < 0.0) ? 0.0 : pf->update_vol_frac[1]);
			
			if((pf->vol_frac[0]) * (pf->vol_frac[1]) == 0)
			{
				if(pf->vol_frac[0] > pf->vol_frac[1]) pf->boundary_reveal = 0.5;
				else pf->boundary_reveal = 1;
			}
			else pf->boundary_reveal = 0;
		}
		
		for(grain = 0; grain < num_grains; grain++) //buffer layer update
		{
			for(i = 0; i < H; i++)
			{
				buffer_layer[grain].left_buffer[i] = phase_field[i * W].vol_frac[grain];
				buffer_layer[grain].right_buffer[i] = phase_field[i * W + W-1].vol_frac[grain];
			}
			for(j = 0;j < W; j++)
			{
				buffer_layer[grain].up_buffer[j] = phase_field[j].vol_frac[grain];
				buffer_layer[grain].down_buffer[j] = phase_field[(H-1) * W + j].vol_frac[grain];
			}
		}
		
		t++;
		uint64_t rem = t%save_after;
		if(rem == 0) // writing loop
		{
			FILE *imagefile;
			char filename[256];
			
			sprintf(filename,"../Data/images/structure%lu.pgm",t);
			imagefile = fopen(filename , "w");
			
			fprintf(imagefile , "P2\n%d %d\n255\n", W, H);
			
			for(i = 0; i < H; i++)
			{
				for(j = 0; j < W; j++)
				{
					fprintf(imagefile , "%d " , (uint8_t)(phase_field[i * W + j].boundary_reveal * 255));
				}
				fprintf(imagefile ,"\n");
			}
			fclose(imagefile);
		}
	}
	
	free(phase_field);
}

int main()
{
	TimeEvolution();
	return 0;
}
