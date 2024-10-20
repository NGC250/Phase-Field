//simulation
#define H 100
#define W 100
#define N (H * W)
#define step 1.0
#define tick 0.005
#define save_after 100

//lattice
#define diffusion_coefficient 1.0
#define barrier_height 1.0
#define gradient_coefficient 2.0
#define noise_intensity 0.05

const char* save_location = "../images/CH_spinodal";