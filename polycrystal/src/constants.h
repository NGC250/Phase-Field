//Lattice
#define H 100
#define W 100
#define N (H*W)
#define step 1.0
#define tick 0.1
#define save_after 50
#define num_grains 10
#define diff_area 1.0/(step * step)

//physical parameters
#define interface_mobility 1.0
#define barrier_height 1.0
#define gradient_coefficient 2.0

#define K (1/sqrt(3.141592653589793238))

typedef struct
{
	double eta[num_grains];
	double boundary_reveal;
	double update_eta[num_grains];
}
Parameters;

typedef struct
{
	double laplacian[num_grains];
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

uint64_t grain , point , i , j;