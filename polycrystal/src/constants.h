//Lattice
#define H 110
#define W 110
#define N (H*W)
#define step 1.0
#define tick 0.1
#define save_after 100
#define num_grains 15

//physical parameters
#define interface_mobility 1.0
#define barrier_height 1.0
#define gradient_coefficient 2.0

#define K (1/sqrt(3.141592653589793238))

typedef struct{
	double eta[num_grains];
	double boundary_reveal;
}orientation;

const char* save_location = "../images/polycrystal";
const char* save_location2 = "../images/layerFIVE";

long i;
int j;