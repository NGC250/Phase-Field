//Simulation
#define step 1.0 // dx,dy,dz
#define tick 0.1 // dt
#define iterations 100000// number of cycles
#define save_after 500 // save at every time step
#define N0 100000 // threshold for parallelization

//system
#define H 256 // System size along x(machine)
#define W 256 // System size along y(machine)
#define D 1 // System size along z(machine)
#define N (H * W * D) // System volume(machine)

// Physical Parameters
#define pi 3.141592653589793238
#define mu0 (1.256e-6) // Permeability of free space

#define Ms (1.432e6) // Saturation magnetization

#define K1 (2.0e4) // First-order anisotropy constant
#define K2 (-4.5e4) // Second-order anisotropy constant

#define H1 (5.0e3) // External field along x
#define H2 0.0 // External field along y
#define H3 0.0 // External field along z

#define alpha 0.5 // damping coeff

#define A_star 0.0625 // reduced exchange constant

#define l111 0.0
#define l100 (2.64e-4)

#define c11 (1.96e11)
#define c12 (1.56e11)
#define c44 (1.23e11)
#define chi (c11 - c12 - 2.0 * c44)

// auxiliary variables : To reduce computations

#define aux0 (1.0/(mu0 * Ms * Ms))
#define aux1 (A_star * tick)
#define aux2 (alpha * tick)
#define aux3 (A_star * alpha * tick)

#define K1_S (2.0 * aux0 * K1)
#define K2_S (2.0 * aux0 * K2)

#define H_EXT_1 (H1/Ms)
#define H_EXT_2 (H2/Ms)
#define H_EXT_3 (H3/Ms)

#define c11_S (c11 * aux0)
#define c12_S (c12 * aux0)
#define c44_S (c44 * aux0)
#define aux4 (4.5 * (2.0 * c44_S * l111 * l111 - (c11_S - c12_S) * l100 * l100))
#define aux5 (3.0 * l100 * (c11_S - c12_S))
#define aux6 (6.0 * l111 * c44_S)
#define aux7 (c11 - c44)
#define aux8 (chi * (c11 + c12))
#define aux9 (-c44 * (c12 + c44))
#define theta1 (c44 * c44 * c11)
#define theta2 (c44 * chi * (c11 + c12))

////////////////////////////////////////////////////////////////////////////////////////////

fftw_plan do_fourier , do_fourier2 , do_inv_fourier;

fftw_complex *m1 , *m2 , *m3;
fftw_complex *h1 , *h2 , *h3;
fftw_complex *gs_coeff_1 , *gs_coeff_2 , *gs_coeff_3;
fftw_complex *sclr_pot , *Hdemag1 , *Hdemag2;
fftw_complex *fourier_m1 , *fourier_m2 , *temp_gs;
fftw_complex *strain_hetero;

fftw_complex *input , *output;

double *dist_fx , *dist_fy , *zeta_mag , *inv_Det_K;

uint64_t point , grain , i , j;

int64_t k , l;

uint8_t p , q , r , s;

double zeta_x , zeta_y , mzn[3];
const double norm = 1.0/((double)(N));

double Delta[3][3] , Strain_homog[3][3] , Mstr_coeff[3][3] , Stiffness[3][3][3][3]; // tensors

typedef struct
{
  double mimj[3][3]; // average value of pair products
}
  Avg_m_vect;

typedef struct
{
  fftw_complex* eps0; // Stress free eigen strain
}
  stress_free_strain;

typedef struct
{
  double* strn_fld; // strain field
}
  Strain;

typedef struct
{
  fftw_complex* comp; // Components of X vector
}
  homog_stress_f;

typedef struct
{
  double cell[3][3]; // Components of Nij tensor
}
  Cof_K;

Avg_m_vect m_avg;

stress_free_strain SFS[6];

Strain strain[6];

homog_stress_f* X;

Cof_K* N_ij;
