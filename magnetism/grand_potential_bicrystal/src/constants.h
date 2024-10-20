//Simulation
#define step (1e-4)
#define tick (1e-3)

#define iterations 250
#define save_after 1
#define H 100
#define W 100
#define N (H*W)
#define num_grains 2

#define pi 3.141592653589793238

// Physical Parameters
#define field_strength (13.5e6)
#define delta_chi (4.01e-6)
#define mu0 (1.26e-6)

#define bulk_energy_coeff (0.3066e6)
#define grad_energy_coeff (2.2371e-6)
#define mobility (1809.42)
#define interface_width (6e-4)
#define omega (1.0/(mobility * interface_width * interface_width))
#define interfacial_energy 0.002

#define interacting_members 2.0

double diff_area = 1.0/(step * step);

typedef struct
{
	double vol_frac[num_grains];
	double update_vol_frac[num_grains];
	double boundary_reveal;
}
Parameters;

typedef struct
{
	double laplacian[num_grains];
}
Laplacian;

typedef struct
{
	double left_buffer[H];
	double right_buffer[H];
	double up_buffer[W];
	double down_buffer[W];
}
Buffer;

typedef struct
{
	uint64_t index;
	double value;
}
MinResult;

Parameters* phase_field;
Buffer* buffer_layer;
Laplacian* L_vol_frac;

uint64_t grain , point , i , j;