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
#define diff_area (1.0/(step * step))

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