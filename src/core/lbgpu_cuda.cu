/* 
   Copyright (C) 2010,2011,2012,2013,2014 The ESPResSo project

   This file is part of ESPResSo.
  
   ESPResSo is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.
   
   ESPResSo is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.
   
   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

/** \file lbgpu_cuda.cu
 *
 * Cuda (.cu) file for the Lattice Boltzmann implementation on GPUs.
 * Header file for \ref lbgpu.hpp.
 */

#include "config.hpp"

#ifdef LB_GPU
#include <stdio.h>
#include <cuda.h>
#include <stdlib.h>

#include "electrokinetics.hpp"
#include "electrokinetics_pdb_parse.hpp"
#include "lbgpu.hpp"
#include "cuda_interface.hpp"
#include "cuda_utils.hpp"

#if (!defined(FLATNOISE) && !defined(GAUSSRANDOMCUT) && !defined(GAUSSRANDOM))
#define FLATNOISE
#endif

#define checkCudaErrors(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true)
{
   if (code != cudaSuccess) 
   {
      fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
      if (abort) exit(code);
   }
}

int extended_values_flag=0; /* TODO: this has to be set to one by
                               appropriate functions if there is 
                               the need to compute pi at every 
                               step (e.g. moving boundaries)*/

/**defining structures residing in global memory */

/** device_rho_v: struct for hydrodynamic fields: this is for internal use 
    (i.e. stores values in LB units) and should not used for 
    printing values  */
static LB_rho_v_gpu *device_rho_v= NULL;

/** device_rho_v_pi: extended struct for hydrodynamic fields: this is the interface
    to tcl, and stores values in MD units. It should not be used
    as an input for any LB calculations. TODO: This structure is not yet 
    used, and it is here to allow access to the stress tensor at any
    timestep, e.g. for future implementations of moving boundary codes */
static LB_rho_v_pi_gpu *device_rho_v_pi= NULL;

/** print_rho_v_pi: struct for hydrodynamic fields: this is the interface
    to tcl, and stores values in MD units. It should not used
    as an input for any LB calculations. TODO: in the future,
    one might want to have several structures for printing 
    separately rho, v, pi without having to compute/store 
    the complete set. */
static LB_rho_v_pi_gpu *print_rho_v_pi= NULL;

/** structs for velocity densities */
static LB_nodes_gpu nodes_a = {.vd=NULL,.seed=NULL,.boundary=NULL};
static LB_nodes_gpu nodes_b = {.vd=NULL,.seed=NULL,.boundary=NULL};;
/** struct for node force */

LB_node_force_gpu node_f = {.force=NULL,.scforce=NULL} ;

static LB_extern_nodeforce_gpu *extern_nodeforces = NULL;

#ifdef LB_BOUNDARIES_GPU
static float* lb_boundary_force = NULL;

static float* lb_boundary_velocity = NULL;

/** pointer for bound index array*/
static int *boundary_node_list;
static int *boundary_index_list;

static __device__ __constant__ int n_lb_boundaries_gpu = 0;
static __device__ __constant__ int n_lb_moving_boundaries_gpu = 0;
static int host_n_lb_moving_boundaries = 0;
static size_t size_of_boundindex;
static double tau_half;
static double tau_sq;
/** pointer for boundary struct*/
static LB_moving_boundary *lb_moving_boundary;
static Particle *host_lb_moving_boundary;
static float *d_lab_anchors = NULL;
static float *h_lab_anchors = NULL;
#endif

#if defined(ELECTROKINETICS)
static __device__ __constant__ int ek_initialized_gpu = 0;
#endif

EK_parameters* lb_ek_parameters_gpu;

/** pointers for additional cuda check flag*/
static int *gpu_check = NULL;
static int *h_gpu_check = NULL;

static unsigned int intflag = 1;
LB_nodes_gpu *current_nodes = NULL;
LB_nodes_gpu *next_nodes = NULL;
/**defining size values for allocating global memory */
static size_t size_of_rho_v;
static size_t size_of_rho_v_pi;
static size_t size_of_extern_nodeforces;

/**parameters residing in constant memory */
static __device__ __constant__ LB_parameters_gpu para;
static const float c_sound_sq = 1.0f/3.0f;

//DEBUG
unsigned long int n_integration_steps = 0;
double add_track = 0.0f;
double hyd_track = 0.0f;
double torque_track[3] = {0.0f, 0.0f, 0.0f};
/*-------------------------------------------------------*/
/*********************************************************/
/** \name device functions called by kernel functions */
/*********************************************************/
/*-------------------------------------------------------*/

/*-------------------------------------------------------*/

/** atomic add function for sveral cuda architectures 
*/
__device__ inline void atomicadd(float* address, float value){
#if !defined __CUDA_ARCH__ || __CUDA_ARCH__ >= 200 // for Fermi, atomicAdd supports floats
  atomicAdd(address, value);
#elif __CUDA_ARCH__ >= 110
#warning Using slower atomicAdd emulation
// float-atomic-add from 
// [url="http://forums.nvidia.com/index.php?showtopic=158039&view=findpost&p=991561"]
  float old = value;
  while ((old = atomicExch(address, atomicExch(address, 0.0f)+old))!=0.0f);
#else
#error I need at least compute capability 1.1
#endif
}

__device__ inline void atomicadd (double* address, double value) {
  unsigned long long oldval, newval, readback;
  oldval = __double_as_longlong(*address);
  newval = __double_as_longlong(__longlong_as_double(oldval) + value);
  while ((readback=atomicCAS((unsigned long long *)address, oldval, newval)) != oldval)
  {
    oldval = readback;
    newval = __double_as_longlong(__longlong_as_double(oldval) + value);
  }
}

/**randomgenerator which generates numbers [0,1]
 * @param *rn Pointer to randomnumber array of the local node or particle 
*/
__device__ void random_01(LB_randomnr_gpu *rn){

  const float mxi = 1.0f/(float)(1ul<<31);
  unsigned int curr = rn->seed;

  curr = 1103515245 * curr + 12345;
  rn->randomnr[0] = (float)(curr & ((1ul<<31)-1))*mxi;
  curr = 1103515245 * curr + 12345;
  rn->randomnr[1] = (float)(curr & ((1ul<<31)-1))*mxi;
  rn->seed = curr;

}

/**randomgenerator which generates numbers between -2 sigma and 2 sigma in the form of a Gaussian with standard deviation sigma=1.118591404 resulting in 
 * an actual standard deviation of 1.
 * @param *rn Pointer to randomnumber array of the local node or particle 
*/
__device__ void gaussian_random_cut(LB_randomnr_gpu *rn){

  float x1, x2;
  float r2, fac;
  /** On every second call two gaussian random numbers are calculated
   via the Box-Muller transformation.*/
  /** draw two uniform random numbers in the unit circle */
  do {
    random_01(rn);
    x1 = 2.0f*rn->randomnr[0] - 1.0f;
    x2 = 2.0f*rn->randomnr[1] - 1.0f;
    r2 = x1*x1 + x2*x2;
  } while (r2 >= 1.0f || r2 == 0.0f);

  /** perform Box-Muller transformation and cutoff the ends and replace with flat noise */
  /*
  fac = sqrtf(-2.0f*__logf(r2)/r2)*1.118591404f;
  rn->randomnr[0] = x2*fac;
  rn->randomnr[1] = x1*fac;
  random_01(rn);
  if ( fabs(rn->randomnr[0]) > 2.0f*1.118591404f) {
    rn->randomnr[0] = (2.0f*rn->randomnr[0]-1.0f)*2.0f*1.118591404f;
  }
  if ( fabs(rn->randomnr[1]) > 2.0f*1.118591404f ) {
    rn->randomnr[0] = (2.0f*rn->randomnr[1]-1.0f)*2.0f*1.118591404f;
  }
  */
  
  fac = sqrtf(-2.0f*__logf(r2)/r2)*1.042267973f;
  rn->randomnr[0] = x2*fac;
  rn->randomnr[1] = x1*fac;
  if ( fabs(rn->randomnr[0]) > 2.0f*1.042267973f) {
    if ( rn->randomnr[0] > 0 ) rn->randomnr[0] = 2.0f*1.042267973f;
    else rn->randomnr[0] = -2.0f*1.042267973f;
  }
  if ( fabs(rn->randomnr[1]) > 2.0f*1.042267973f ) {
    if ( rn->randomnr[1] > 0 ) rn->randomnr[1] = 2.0f*1.042267973f;
    else rn->randomnr[1] = -2.0f*1.042267973f;
  }
}

/** gaussian random nummber generator for thermalisation
 * @param *rn Pointer to randomnumber array of the local node node or particle 
*/
__device__ void gaussian_random(LB_randomnr_gpu *rn){

  float x1, x2;
  float r2, fac;
  /** On every second call two gaussian random numbers are calculated
   via the Box-Muller transformation.*/
  /** draw two uniform random numbers in the unit circle */
  do {
    random_01(rn);
    x1 = 2.0f*rn->randomnr[0]-1.0f;
    x2 = 2.0f*rn->randomnr[1]-1.0f;
    r2 = x1*x1 + x2*x2;
  } while (r2 >= 1.0f || r2 == 0.0f);

  /** perform Box-Muller transformation */
  fac = sqrtf(-2.0f*__logf(r2)/r2);
  rn->randomnr[0] = x2*fac;
  rn->randomnr[1] = x1*fac;
  
}
/* wrapper */
__device__ void random_wrapper(LB_randomnr_gpu *rn) { 

#if defined(FLATNOISE)
#define sqrt12 3.46410161514f
  random_01(rn);
  rn->randomnr[0]-=0.5f;
  rn->randomnr[0]*=sqrt12;
  rn->randomnr[1]-=0.5f;
  rn->randomnr[1]*=sqrt12;
#elif defined(GAUSSRANDOMCUT)
  gaussian_random_cut(rn);
#elif defined(GAUSSRANDOM)
  gaussian_random(rn);
#else
#error No noise type defined for the GPU LB
#endif  
  
}


/**tranformation from 1d array-index to xyz
 * @param index   node index / thread index (Input)
 * @param xyz     Pointer to calculated xyz array (Output)
 */
__device__ void index_to_xyz(unsigned int index, unsigned int *xyz){

  xyz[0] = index%para.dim_x;
  index /= para.dim_x;
  xyz[1] = index%para.dim_y;
  index /= para.dim_y;
  xyz[2] = index;
} 
/** compute minimum image distance vector from node to boundary center
 * @param xyz pointer to node xyz
 * @param center pointer to center xyz
 * @param delta_xyz pointer to distance vector
 */
__device__ void xyz_center_to_delta(unsigned int *xyz, float *center, float *delta_xyz){

  delta_xyz[0] = center[0] - xyz[0];
  delta_xyz[0] -= rintf(delta_xyz[0]/para.dim_x)*para.dim_x;
  delta_xyz[1] = center[1] - xyz[1];
  delta_xyz[1] -= rintf(delta_xyz[1]/para.dim_y)*para.dim_y;
  delta_xyz[2] = center[2] - xyz[2];
  delta_xyz[2] -= rintf(delta_xyz[2]/para.dim_z)*para.dim_z;
}

/**calculation of the modes from the velocity densities (space-transform.)
 * @param n_a     Pointer to local node residing in array a (Input)
 * @param index   Node index / thread index (Input)
 * @param mode    Pointer to the local register values mode (Output)
*/
__device__ void calc_m_from_n(LB_nodes_gpu n_a, unsigned int index, float *mode){

  #pragma unroll
  for(int ii=0;ii<LB_COMPONENTS;++ii)
  {
    // The following convention is used:
    // The $\hat{c}_i$ form B. Duenweg's paper are given by:

    /* c_0  = { 0, 0, 0}
       c_1  = { 1, 0, 0}
       c_2  = {-1, 0, 0}
       c_3  = { 0, 1, 0}
       c_4  = { 0,-1, 0}
       c_5  = { 0, 0, 1}
       c_6  = { 0, 0,-1}
       c_7  = { 1, 1, 0}
       c_8  = {-1,-1, 0}
       c_9  = { 1,-1, 0}
       c_10 = {-1, 1, 0}
       c_11 = { 1, 0, 1}
       c_12 = {-1, 0,-1}
       c_13 = { 1, 0,-1}
       c_14 = {-1, 0, 1}
       c_15 = { 0, 1, 1}
       c_16 = { 0,-1,-1}
       c_17 = { 0, 1,-1}
       c_18 = { 0,-1, 1} */

    // The basis vectors (modes) are constructed as follows
    // $m_k = \sum_{i} e_{ki} n_{i}$, where the $e_{ki}$ form a 
    // linear transformation (matrix) that is given by

    /* $e{ 0,i} = 1$
       $e{ 1,i} = c_{i,x}$
       $e{ 2,i} = c_{i,y}$
       $e{ 3,i} = c_{i,z}$
       $e{ 4,i} = c_{i}^2 - 1$
       $e{ 5,i} = c_{i,x}^2 - c_{i,y}^2$
       $e{ 6,i} = c_{i}^2 - 3*c_{i,z}^2$
       $e{ 7,i} = c_{i,x}*c_{i,y}$
       $e{ 8,i} = c_{i,x}*c_{i,z}$
       $e{ 9,i} = c_{i,y}*c_{i,z}$
       $e{10,i} = (3*c_{i}^2 - 5)*c_{i,x}$
       $e{11,i} = (3*c_{i}^2 - 5)*c_{i,y}$
       $e{12,i} = (3*c_{i}^2 - 5)*c_{i,z}$
       $e{13,i} = (c_{i,y}^2 - c_{i,z}^2)*c_{i,x}$
       $e{14,i} = (c_{i,x}^2 - c_{i,z}^2)*c_{i,y}$
       $e{15,i} = (c_{i,x}^2 - c_{i,y}^2)*c_{i,z}$
       $e{16,i} = 3*c_{i}^2^2 - 6*c_{i}^2 + 1$
       $e{17,i} = (2*c_{i}^2 - 3)*(c_{i,x}^2 - c_{i,y}^2)$
       $e{18,i} = (2*c_{i}^2 - 3)*(c_{i}^2 - 3*c_{i,z}^2)$ */

    // Such that the transformation matrix is given by

    /* {{ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
        { 0, 1,-1, 0, 0, 0, 0, 1,-1, 1,-1, 1,-1, 1,-1, 0, 0, 0, 0}, 
        { 0, 0, 0, 1,-1, 0, 0, 1,-1,-1, 1, 0, 0, 0, 0, 1,-1, 1,-1}, 
        { 0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 1,-1,-1, 1, 1,-1,-1, 1}, 
        {-1, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
        { 0, 1, 1,-1,-1, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1,-1,-1,-1,-1}, 
        { 0, 1, 1, 1, 1,-2,-2, 2, 2, 2, 2,-1,-1,-1,-1,-1,-1,-1,-1}, 
        { 0, 0, 0, 0, 0, 0, 0, 1, 1,-1,-1, 0, 0, 0, 0, 0, 0, 0, 0}, 
        { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1,-1,-1, 0, 0, 0, 0}, 
        { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1,-1,-1}, 
        { 0,-2, 2, 0, 0, 0, 0, 1,-1, 1,-1, 1,-1, 1,-1, 0, 0, 0, 0}, 
        { 0, 0, 0,-2, 2, 0, 0, 1,-1,-1, 1, 0, 0, 0, 0, 1,-1, 1,-1}, 
        { 0, 0, 0, 0, 0,-2, 2, 0, 0, 0, 0, 1,-1,-1, 1, 1,-1,-1, 1}, 
        { 0, 0, 0, 0, 0, 0, 0, 1,-1, 1,-1,-1, 1,-1, 1, 0, 0, 0, 0}, 
        { 0, 0, 0, 0, 0, 0, 0, 1,-1,-1, 1, 0, 0, 0, 0,-1, 1,-1, 1}, 
        { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1,-1, 1,-1, 1, 1,-1}, 
        { 1,-2,-2,-2,-2,-2,-2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
        { 0,-1,-1, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1,-1,-1,-1,-1}, 
        { 0,-1,-1,-1,-1, 2, 2, 2, 2, 2, 2,-1,-1,-1,-1,-1,-1,-1,-1}} */

    // With weights 

    /* q^{c_{i}} = { 1/3, 1/18, 1/18, 1/18,
                    1/18, 1/18, 1/18, 1/36,
                    1/36, 1/36, 1/36, 1/36, 
                    1/36, 1/36, 1/36, 1/36, 
                    1/36, 1/36, 1/36 } */

    // Which makes the transformation satisfy the following
    // orthogonality condition:
    // \sum_{i} q^{c_{i}} e_{ki} e_{li} = w_{k} \delta_{kl},
    // where the weights are:

    /* w_{i} = {  1, 1/3, 1/3, 1/3,
                2/3, 4/9, 4/3, 1/9,
                1/9, 1/9, 2/3, 2/3,
                2/3, 2/9, 2/9, 2/9, 
                  2, 4/9, 4/3 } */

    // mass mode

    mode[0 + ii * LBQ] =   n_a.vd[( 0 + ii*LBQ ) * para.number_of_nodes + index]
                         + n_a.vd[( 1 + ii*LBQ ) * para.number_of_nodes + index] + n_a.vd[( 2 + ii*LBQ ) * para.number_of_nodes + index]
                         + n_a.vd[( 3 + ii*LBQ ) * para.number_of_nodes + index] + n_a.vd[( 4 + ii*LBQ ) * para.number_of_nodes + index]
                         + n_a.vd[( 5 + ii*LBQ ) * para.number_of_nodes + index] + n_a.vd[( 6 + ii*LBQ ) * para.number_of_nodes + index]
                         + n_a.vd[( 7 + ii*LBQ ) * para.number_of_nodes + index] + n_a.vd[( 8 + ii*LBQ ) * para.number_of_nodes + index]
                         + n_a.vd[( 9 + ii*LBQ ) * para.number_of_nodes + index] + n_a.vd[(10 + ii*LBQ ) * para.number_of_nodes + index]
                         + n_a.vd[(11 + ii*LBQ ) * para.number_of_nodes + index] + n_a.vd[(12 + ii*LBQ ) * para.number_of_nodes + index]
                         + n_a.vd[(13 + ii*LBQ ) * para.number_of_nodes + index] + n_a.vd[(14 + ii*LBQ ) * para.number_of_nodes + index]
                         + n_a.vd[(15 + ii*LBQ ) * para.number_of_nodes + index] + n_a.vd[(16 + ii*LBQ ) * para.number_of_nodes + index]
                         + n_a.vd[(17 + ii*LBQ ) * para.number_of_nodes + index] + n_a.vd[(18 + ii*LBQ ) * para.number_of_nodes + index];

    // momentum modes

    mode[1 + ii * LBQ] =   (n_a.vd[( 1 + ii*LBQ ) * para.number_of_nodes + index] - n_a.vd[( 2 + ii*LBQ ) * para.number_of_nodes + index])
                         + (n_a.vd[( 7 + ii*LBQ ) * para.number_of_nodes + index] - n_a.vd[( 8 + ii*LBQ ) * para.number_of_nodes + index])
                         + (n_a.vd[( 9 + ii*LBQ ) * para.number_of_nodes + index] - n_a.vd[(10 + ii*LBQ ) * para.number_of_nodes + index])
                         + (n_a.vd[(11 + ii*LBQ ) * para.number_of_nodes + index] - n_a.vd[(12 + ii*LBQ ) * para.number_of_nodes + index])
                         + (n_a.vd[(13 + ii*LBQ ) * para.number_of_nodes + index] - n_a.vd[(14 + ii*LBQ ) * para.number_of_nodes + index]);

    mode[2 + ii * LBQ] =   (n_a.vd[( 3 + ii*LBQ ) * para.number_of_nodes + index] - n_a.vd[( 4 + ii*LBQ ) * para.number_of_nodes + index])
                         + (n_a.vd[( 7 + ii*LBQ ) * para.number_of_nodes + index] - n_a.vd[( 8 + ii*LBQ ) * para.number_of_nodes + index])
                         - (n_a.vd[( 9 + ii*LBQ ) * para.number_of_nodes + index] - n_a.vd[(10 + ii*LBQ ) * para.number_of_nodes + index])
                         + (n_a.vd[(15 + ii*LBQ ) * para.number_of_nodes + index] - n_a.vd[(16 + ii*LBQ ) * para.number_of_nodes + index])
                         + (n_a.vd[(17 + ii*LBQ ) * para.number_of_nodes + index] - n_a.vd[(18 + ii*LBQ ) * para.number_of_nodes + index]);

    mode[3 + ii * LBQ] =   (n_a.vd[( 5 + ii*LBQ ) * para.number_of_nodes + index] - n_a.vd[( 6 + ii*LBQ ) * para.number_of_nodes + index])
                         + (n_a.vd[(11 + ii*LBQ ) * para.number_of_nodes + index] - n_a.vd[(12 + ii*LBQ ) * para.number_of_nodes + index])
                         - (n_a.vd[(13 + ii*LBQ ) * para.number_of_nodes + index] - n_a.vd[(14 + ii*LBQ ) * para.number_of_nodes + index])
                         + (n_a.vd[(15 + ii*LBQ ) * para.number_of_nodes + index] - n_a.vd[(16 + ii*LBQ ) * para.number_of_nodes + index])
                         - (n_a.vd[(17 + ii*LBQ ) * para.number_of_nodes + index] - n_a.vd[(18 + ii*LBQ ) * para.number_of_nodes + index]);

    // stress modes
    mode[4 + ii * LBQ] = - n_a.vd[( 0 + ii*LBQ ) * para.number_of_nodes + index]
                         + n_a.vd[( 7 + ii*LBQ ) * para.number_of_nodes + index] + n_a.vd[( 8 + ii*LBQ ) * para.number_of_nodes + index]
                         + n_a.vd[( 9 + ii*LBQ ) * para.number_of_nodes + index] + n_a.vd[(10 + ii*LBQ ) * para.number_of_nodes + index]
                         + n_a.vd[(11 + ii*LBQ ) * para.number_of_nodes + index] + n_a.vd[(12 + ii*LBQ ) * para.number_of_nodes + index]
                         + n_a.vd[(13 + ii*LBQ ) * para.number_of_nodes + index] + n_a.vd[(14 + ii*LBQ ) * para.number_of_nodes + index]
                         + n_a.vd[(15 + ii*LBQ ) * para.number_of_nodes + index] + n_a.vd[(16 + ii*LBQ ) * para.number_of_nodes + index]
                         + n_a.vd[(17 + ii*LBQ ) * para.number_of_nodes + index] + n_a.vd[(18 + ii*LBQ ) * para.number_of_nodes + index];

    mode[5 + ii * LBQ] =   (n_a.vd[( 1 + ii*LBQ ) * para.number_of_nodes + index] + n_a.vd[( 2 + ii*LBQ ) * para.number_of_nodes + index])
                         - (n_a.vd[( 3 + ii*LBQ ) * para.number_of_nodes + index] + n_a.vd[( 4 + ii*LBQ ) * para.number_of_nodes + index])
                         + (n_a.vd[(11 + ii*LBQ ) * para.number_of_nodes + index] + n_a.vd[(12 + ii*LBQ ) * para.number_of_nodes + index])
                         + (n_a.vd[(13 + ii*LBQ ) * para.number_of_nodes + index] + n_a.vd[(14 + ii*LBQ ) * para.number_of_nodes + index])
                         - (n_a.vd[(15 + ii*LBQ ) * para.number_of_nodes + index] + n_a.vd[(16 + ii*LBQ ) * para.number_of_nodes + index])
                         - (n_a.vd[(17 + ii*LBQ ) * para.number_of_nodes + index] + n_a.vd[(18 + ii*LBQ ) * para.number_of_nodes + index]);

    mode[6 + ii * LBQ] =   (n_a.vd[( 1 + ii*LBQ ) * para.number_of_nodes + index] + n_a.vd[( 2 + ii*LBQ ) * para.number_of_nodes + index])
                         + (n_a.vd[( 3 + ii*LBQ ) * para.number_of_nodes + index] + n_a.vd[( 4 + ii*LBQ ) * para.number_of_nodes + index])
                         - (n_a.vd[(11 + ii*LBQ ) * para.number_of_nodes + index] + n_a.vd[(12 + ii*LBQ ) * para.number_of_nodes + index])
                         - (n_a.vd[(13 + ii*LBQ ) * para.number_of_nodes + index] + n_a.vd[(14 + ii*LBQ ) * para.number_of_nodes + index])
                         - (n_a.vd[(15 + ii*LBQ ) * para.number_of_nodes + index] + n_a.vd[(16 + ii*LBQ ) * para.number_of_nodes + index])
                         - (n_a.vd[(17 + ii*LBQ ) * para.number_of_nodes + index] + n_a.vd[(18 + ii*LBQ ) * para.number_of_nodes + index])
                         - 2.0f*( 
                                    (n_a.vd[( 5 + ii*LBQ ) * para.number_of_nodes + index] + n_a.vd[( 6 + ii*LBQ ) * para.number_of_nodes + index])
                                  - (n_a.vd[( 7 + ii*LBQ ) * para.number_of_nodes + index] + n_a.vd[( 8 + ii*LBQ ) * para.number_of_nodes + index])
                                  - (n_a.vd[( 9 + ii*LBQ ) * para.number_of_nodes + index] + n_a.vd[(10 + ii*LBQ ) * para.number_of_nodes + index])
                                );

    mode[7 + ii * LBQ] =   (n_a.vd[( 7 + ii*LBQ ) * para.number_of_nodes + index] + n_a.vd[( 8 + ii*LBQ ) * para.number_of_nodes + index])
                         - (n_a.vd[( 9 + ii*LBQ ) * para.number_of_nodes + index] + n_a.vd[(10 + ii*LBQ ) * para.number_of_nodes + index]);

    mode[8 + ii * LBQ] =   (n_a.vd[(11 + ii*LBQ ) * para.number_of_nodes + index] + n_a.vd[(12 + ii*LBQ ) * para.number_of_nodes + index])
                         - (n_a.vd[(13 + ii*LBQ ) * para.number_of_nodes + index] + n_a.vd[(14 + ii*LBQ ) * para.number_of_nodes + index]);

    mode[9 + ii * LBQ] =   (n_a.vd[(15 + ii*LBQ ) * para.number_of_nodes + index] + n_a.vd[(16 + ii*LBQ ) * para.number_of_nodes + index])
                         - (n_a.vd[(17 + ii*LBQ ) * para.number_of_nodes + index] + n_a.vd[(18 + ii*LBQ ) * para.number_of_nodes + index]);

    // kinetic modes

    mode[10 + ii * LBQ] = - 2.0f*(n_a.vd[( 1 + ii*LBQ ) * para.number_of_nodes + index] - n_a.vd[( 2 + ii*LBQ ) * para.number_of_nodes + index])
                               + (n_a.vd[( 7 + ii*LBQ ) * para.number_of_nodes + index] - n_a.vd[( 8 + ii*LBQ ) * para.number_of_nodes + index])
                               + (n_a.vd[( 9 + ii*LBQ ) * para.number_of_nodes + index] - n_a.vd[(10 + ii*LBQ ) * para.number_of_nodes + index])
                               + (n_a.vd[(11 + ii*LBQ ) * para.number_of_nodes + index] - n_a.vd[(12 + ii*LBQ ) * para.number_of_nodes + index])
                               + (n_a.vd[(13 + ii*LBQ ) * para.number_of_nodes + index] - n_a.vd[(14 + ii*LBQ ) * para.number_of_nodes + index]);

    mode[11 + ii * LBQ] = - 2.0f*(n_a.vd[( 3 + ii*LBQ ) * para.number_of_nodes + index] - n_a.vd[( 4 + ii*LBQ ) * para.number_of_nodes + index])
                               + (n_a.vd[( 7 + ii*LBQ ) * para.number_of_nodes + index] - n_a.vd[( 8 + ii*LBQ ) * para.number_of_nodes + index])
                               - (n_a.vd[( 9 + ii*LBQ ) * para.number_of_nodes + index] - n_a.vd[(10 + ii*LBQ ) * para.number_of_nodes + index])
                               + (n_a.vd[(15 + ii*LBQ ) * para.number_of_nodes + index] - n_a.vd[(16 + ii*LBQ ) * para.number_of_nodes + index])
                               + (n_a.vd[(17 + ii*LBQ ) * para.number_of_nodes + index] - n_a.vd[(18 + ii*LBQ ) * para.number_of_nodes + index]);

    mode[12 + ii * LBQ] = - 2.0f*(n_a.vd[( 5 + ii*LBQ ) * para.number_of_nodes + index] - n_a.vd[( 6 + ii*LBQ ) * para.number_of_nodes + index])
                               + (n_a.vd[(11 + ii*LBQ ) * para.number_of_nodes + index] - n_a.vd[(12 + ii*LBQ ) * para.number_of_nodes + index])
                               - (n_a.vd[(13 + ii*LBQ ) * para.number_of_nodes + index] - n_a.vd[(14 + ii*LBQ ) * para.number_of_nodes + index])
                               + (n_a.vd[(15 + ii*LBQ ) * para.number_of_nodes + index] - n_a.vd[(16 + ii*LBQ ) * para.number_of_nodes + index])
                               - (n_a.vd[(17 + ii*LBQ ) * para.number_of_nodes + index] - n_a.vd[(18 + ii*LBQ ) * para.number_of_nodes + index]);

    mode[13 + ii * LBQ] =   (n_a.vd[( 7 + ii*LBQ ) * para.number_of_nodes + index] - n_a.vd[( 8 + ii*LBQ ) * para.number_of_nodes + index])
                          + (n_a.vd[( 9 + ii*LBQ ) * para.number_of_nodes + index] - n_a.vd[(10 + ii*LBQ ) * para.number_of_nodes + index])
                          - (n_a.vd[(11 + ii*LBQ ) * para.number_of_nodes + index] - n_a.vd[(12 + ii*LBQ ) * para.number_of_nodes + index])
                          - (n_a.vd[(13 + ii*LBQ ) * para.number_of_nodes + index] - n_a.vd[(14 + ii*LBQ ) * para.number_of_nodes + index]);

    mode[14 + ii * LBQ] =   (n_a.vd[( 7 + ii*LBQ ) * para.number_of_nodes + index] - n_a.vd[( 8 + ii*LBQ ) * para.number_of_nodes + index])
                          - (n_a.vd[( 9 + ii*LBQ ) * para.number_of_nodes + index] - n_a.vd[(10 + ii*LBQ ) * para.number_of_nodes + index])
                          - (n_a.vd[(15 + ii*LBQ ) * para.number_of_nodes + index] - n_a.vd[(16 + ii*LBQ ) * para.number_of_nodes + index])
                          - (n_a.vd[(17 + ii*LBQ ) * para.number_of_nodes + index] - n_a.vd[(18 + ii*LBQ ) * para.number_of_nodes + index]);

    mode[15 + ii * LBQ] =   (n_a.vd[(11 + ii*LBQ ) * para.number_of_nodes + index] - n_a.vd[(12 + ii*LBQ ) * para.number_of_nodes + index])
                          - (n_a.vd[(13 + ii*LBQ ) * para.number_of_nodes + index] - n_a.vd[(14 + ii*LBQ ) * para.number_of_nodes + index])
                          - (n_a.vd[(15 + ii*LBQ ) * para.number_of_nodes + index] - n_a.vd[(16 + ii*LBQ ) * para.number_of_nodes + index])
                          + (n_a.vd[(17 + ii*LBQ ) * para.number_of_nodes + index] - n_a.vd[(18 + ii*LBQ ) * para.number_of_nodes + index]);

    mode[16 + ii * LBQ] =   n_a.vd[( 0 + ii*LBQ ) * para.number_of_nodes + index]
                          + n_a.vd[( 7 + ii*LBQ ) * para.number_of_nodes + index] + n_a.vd[( 8 + ii*LBQ ) * para.number_of_nodes + index]
                          + n_a.vd[( 9 + ii*LBQ ) * para.number_of_nodes + index] + n_a.vd[(10 + ii*LBQ ) * para.number_of_nodes + index]
                          + n_a.vd[(11 + ii*LBQ ) * para.number_of_nodes + index] + n_a.vd[(12 + ii*LBQ ) * para.number_of_nodes + index]
                          + n_a.vd[(13 + ii*LBQ ) * para.number_of_nodes + index] + n_a.vd[(14 + ii*LBQ ) * para.number_of_nodes + index]
                          + n_a.vd[(15 + ii*LBQ ) * para.number_of_nodes + index] + n_a.vd[(16 + ii*LBQ ) * para.number_of_nodes + index]
                          + n_a.vd[(17 + ii*LBQ ) * para.number_of_nodes + index] + n_a.vd[(18 + ii*LBQ ) * para.number_of_nodes + index]
                          - 2.0f*(
                                     (n_a.vd[( 1 + ii*LBQ ) * para.number_of_nodes + index] + n_a.vd[( 2 + ii*LBQ ) * para.number_of_nodes + index])
                                   + (n_a.vd[( 3 + ii*LBQ ) * para.number_of_nodes + index] + n_a.vd[( 4 + ii*LBQ ) * para.number_of_nodes + index])
                                   + (n_a.vd[( 5 + ii*LBQ ) * para.number_of_nodes + index] + n_a.vd[( 6 + ii*LBQ ) * para.number_of_nodes + index])
                                 );

    mode[17 + ii * LBQ] = - (n_a.vd[( 1 + ii*LBQ ) * para.number_of_nodes + index] + n_a.vd[( 2 + ii*LBQ ) * para.number_of_nodes + index])
                          + (n_a.vd[( 3 + ii*LBQ ) * para.number_of_nodes + index] + n_a.vd[( 4 + ii*LBQ ) * para.number_of_nodes + index])
                          + (n_a.vd[(11 + ii*LBQ ) * para.number_of_nodes + index] + n_a.vd[(12 + ii*LBQ ) * para.number_of_nodes + index])
                          + (n_a.vd[(13 + ii*LBQ ) * para.number_of_nodes + index] + n_a.vd[(14 + ii*LBQ ) * para.number_of_nodes + index])
                          - (n_a.vd[(15 + ii*LBQ ) * para.number_of_nodes + index] + n_a.vd[(16 + ii*LBQ ) * para.number_of_nodes + index])
                          - (n_a.vd[(17 + ii*LBQ ) * para.number_of_nodes + index] + n_a.vd[(18 + ii*LBQ ) * para.number_of_nodes + index]);

    mode[18 + ii * LBQ] = - (n_a.vd[( 1 + ii*LBQ ) * para.number_of_nodes + index] + n_a.vd[( 2 + ii*LBQ ) * para.number_of_nodes + index])
                          - (n_a.vd[( 3 + ii*LBQ ) * para.number_of_nodes + index] + n_a.vd[( 4 + ii*LBQ ) * para.number_of_nodes + index])
                          - (n_a.vd[(11 + ii*LBQ ) * para.number_of_nodes + index] + n_a.vd[(12 + ii*LBQ ) * para.number_of_nodes + index])
                          - (n_a.vd[(13 + ii*LBQ ) * para.number_of_nodes + index] + n_a.vd[(14 + ii*LBQ ) * para.number_of_nodes + index])
                          - (n_a.vd[(15 + ii*LBQ ) * para.number_of_nodes + index] + n_a.vd[(16 + ii*LBQ ) * para.number_of_nodes + index])
                          - (n_a.vd[(17 + ii*LBQ ) * para.number_of_nodes + index] + n_a.vd[(18 + ii*LBQ ) * para.number_of_nodes + index])
                          + 2.0f*(
                                     (n_a.vd[( 5 + ii*LBQ ) * para.number_of_nodes + index] + n_a.vd[( 6 + ii*LBQ ) * para.number_of_nodes + index])
                                   + (n_a.vd[( 7 + ii*LBQ ) * para.number_of_nodes + index] + n_a.vd[( 8 + ii*LBQ ) * para.number_of_nodes + index])
                                   + (n_a.vd[( 9 + ii*LBQ ) * para.number_of_nodes + index] + n_a.vd[(10 + ii*LBQ ) * para.number_of_nodes + index])
                                 );
  }
}

__device__ void reset_LB_forces(unsigned int index, LB_node_force_gpu node_f, bool buffer = true) {

  float force_factor=powf(para.agrid,2)*para.tau*para.tau;
  for(int ii=0;ii<LB_COMPONENTS;++ii)
  {

#if defined(IMMERSED_BOUNDARY) || defined(EK_DEBUG)
// Store backup of the node forces
    if (buffer)
    {
      node_f.force_buf[(0 + ii*3 ) * para.number_of_nodes + index] = node_f.force[(0 + ii*3 ) * para.number_of_nodes + index];
      node_f.force_buf[(1 + ii*3 ) * para.number_of_nodes + index] = node_f.force[(1 + ii*3 ) * para.number_of_nodes + index];
      node_f.force_buf[(2 + ii*3 ) * para.number_of_nodes + index] = node_f.force[(2 + ii*3 ) * para.number_of_nodes + index];
    }
#endif

#ifdef EXTERNAL_FORCES
      if(para.external_force)
      {
        node_f.force[(0 + ii*3 ) * para.number_of_nodes + index] = para.ext_force[0 + ii*3 ]*force_factor;
        node_f.force[(1 + ii*3 ) * para.number_of_nodes + index] = para.ext_force[1 + ii*3 ]*force_factor;
        node_f.force[(2 + ii*3 ) * para.number_of_nodes + index] = para.ext_force[2 + ii*3 ]*force_factor;
      }
      else
      {
        node_f.force[(0 + ii*3 ) * para.number_of_nodes + index] = 0.0f;
        node_f.force[(1 + ii*3 ) * para.number_of_nodes + index] = 0.0f;
        node_f.force[(2 + ii*3 ) * para.number_of_nodes + index] = 0.0f;
      }
#else
      /** reset force */
      node_f.force[(0 + ii*3 ) * para.number_of_nodes + index] = 0.0f;
      node_f.force[(1 + ii*3 ) * para.number_of_nodes + index] = 0.0f;
      node_f.force[(2 + ii*3 ) * para.number_of_nodes + index] = 0.0f;
#endif
  }
}

__global__ void reset_LB_forces_kernel(LB_node_force_gpu node_f, bool buffer = true) {
  unsigned int index = blockIdx.y * gridDim.x * blockDim.x + blockDim.x * blockIdx.x + threadIdx.x;

  if( index < para.number_of_nodes )
    reset_LB_forces(index, node_f, buffer);
}

void reset_LB_forces_GPU(bool buffer) {
  int threads_per_block = 64;
  int blocks_per_grid_y = 4;
  int blocks_per_grid_x = (lbpar_gpu.number_of_nodes + threads_per_block * blocks_per_grid_y - 1) /(threads_per_block * blocks_per_grid_y);
  dim3 dim_grid = make_uint3(blocks_per_grid_x, blocks_per_grid_y, 1);

  KERNELCALL(reset_LB_forces_kernel, dim_grid, threads_per_block, (node_f, buffer));
}


__device__ void update_rho_v(float *mode, unsigned int index, LB_node_force_gpu node_f, LB_rho_v_gpu *d_v){

  float Rho_tot=0.0f;
  float u_tot[3]={0.0f,0.0f,0.0f};
  
  #pragma unroll
  for(int ii=0;ii<LB_COMPONENTS;++ii)
  { 
      /** re-construct the real density
      * remember that the populations are stored as differences to their
      * equilibrium value */

      d_v[index].rho[ii] = mode[0 + ii * LBQ] + para.rho[ii]*para.agrid*para.agrid*para.agrid;
      Rho_tot  += mode[0 + ii * LBQ] + para.rho[ii]*para.agrid*para.agrid*para.agrid;
      u_tot[0] += mode[1 + ii * LBQ];
      u_tot[1] += mode[2 + ii * LBQ];
      u_tot[2] += mode[3 + ii * LBQ];

      /** if forces are present, the momentum density is redefined to
      * inlcude one half-step of the force action.  See the
      * Chapman-Enskog expansion in [Ladd & Verberg]. */

      u_tot[0] += 0.5f*node_f.force[(0+ii*3)*para.number_of_nodes + index];
      u_tot[1] += 0.5f*node_f.force[(1+ii*3)*para.number_of_nodes + index];
      u_tot[2] += 0.5f*node_f.force[(2+ii*3)*para.number_of_nodes + index];
  }

  u_tot[0]/=Rho_tot;
  u_tot[1]/=Rho_tot;
  u_tot[2]/=Rho_tot;

  d_v[index].v[0]=u_tot[0]; 
  d_v[index].v[1]=u_tot[1]; 
  d_v[index].v[2]=u_tot[2]; 
}

/**lb_relax_modes, means collision update of the modes
 * @param index   node index / thread index (Input)
 * @param mode    Pointer to the local register values mode (Input/Output)
 * @param node_f  Pointer to local node force (Input)
 * @param *d_v    Pointer to local device values
*/
__device__ void relax_modes(float *mode, unsigned int index, LB_node_force_gpu node_f, LB_rho_v_gpu *d_v){
  float u_tot[3]={0.0f,0.0f,0.0f};

  update_rho_v(mode, index, node_f, d_v);

  u_tot[0]=d_v[index].v[0];  
  u_tot[1]=d_v[index].v[1];  
  u_tot[2]=d_v[index].v[2];  
 
  #pragma unroll
  for(int ii=0;ii<LB_COMPONENTS;++ii)
  { 
      float Rho; float j[3]; float modes_from_pi_eq[6];

      Rho = mode[0 + ii * LBQ] + para.rho[ii]*para.agrid*para.agrid*para.agrid ;
      j[0] = Rho * u_tot[0];
      j[1] = Rho * u_tot[1];
      j[2] = Rho * u_tot[2];

      /** equilibrium part of the stress modes (eq13 schiller)*/

      modes_from_pi_eq[0] = ((j[0]*j[0])+(j[1]*j[1])+(j[2]*j[2]))/Rho;
      modes_from_pi_eq[1] = ((j[0]*j[0])-(j[1]*j[1]))/Rho;
      modes_from_pi_eq[2] = (((j[0]*j[0])+(j[1]*j[1])+(j[2]*j[2])) - 3.0f*(j[2]*j[2]))/Rho;
      modes_from_pi_eq[3] = j[0]*j[1]/Rho;
      modes_from_pi_eq[4] = j[0]*j[2]/Rho;
      modes_from_pi_eq[5] = j[1]*j[2]/Rho;
 
      /** in Shan-Chen we have to relax the momentum modes as well using the mobility, but
          the total momentum is conserved */  

#ifdef SHANCHEN
      mode[1 + ii * LBQ] = j[0] + para.gamma_mobility[0]*(mode[1 + ii * LBQ] - j[0]);
      mode[2 + ii * LBQ] = j[1] + para.gamma_mobility[0]*(mode[2 + ii * LBQ] - j[1]);
      mode[3 + ii * LBQ] = j[2] + para.gamma_mobility[0]*(mode[3 + ii * LBQ] - j[2]);
#endif
 
      /** relax the stress modes (eq14 schiller)*/

      mode[4 + ii * LBQ] = modes_from_pi_eq[0] +  para.gamma_bulk[ii]*(mode[4 + ii * LBQ] - modes_from_pi_eq[0]);
      mode[5 + ii * LBQ] = modes_from_pi_eq[1] + para.gamma_shear[ii]*(mode[5 + ii * LBQ] - modes_from_pi_eq[1]);
      mode[6 + ii * LBQ] = modes_from_pi_eq[2] + para.gamma_shear[ii]*(mode[6 + ii * LBQ] - modes_from_pi_eq[2]);
      mode[7 + ii * LBQ] = modes_from_pi_eq[3] + para.gamma_shear[ii]*(mode[7 + ii * LBQ] - modes_from_pi_eq[3]);
      mode[8 + ii * LBQ] = modes_from_pi_eq[4] + para.gamma_shear[ii]*(mode[8 + ii * LBQ] - modes_from_pi_eq[4]);
      mode[9 + ii * LBQ] = modes_from_pi_eq[5] + para.gamma_shear[ii]*(mode[9 + ii * LBQ] - modes_from_pi_eq[5]);
    
      /** relax the ghost modes (project them out) */
      /** ghost modes have no equilibrium part due to orthogonality */

      mode[10 + ii * LBQ] =  para.gamma_odd[ii]*mode[10 + ii * LBQ];
      mode[11 + ii * LBQ] =  para.gamma_odd[ii]*mode[11 + ii * LBQ];
      mode[12 + ii * LBQ] =  para.gamma_odd[ii]*mode[12 + ii * LBQ];
      mode[13 + ii * LBQ] =  para.gamma_odd[ii]*mode[13 + ii * LBQ];
      mode[14 + ii * LBQ] =  para.gamma_odd[ii]*mode[14 + ii * LBQ];
      mode[15 + ii * LBQ] =  para.gamma_odd[ii]*mode[15 + ii * LBQ];
      mode[16 + ii * LBQ] = para.gamma_even[ii]*mode[16 + ii * LBQ];
      mode[17 + ii * LBQ] = para.gamma_even[ii]*mode[17 + ii * LBQ];
      mode[18 + ii * LBQ] = para.gamma_even[ii]*mode[18 + ii * LBQ];
  }
}


/**thermalization of the modes with gaussian random numbers
 * @param index   node index / thread index (Input)
 * @param mode    Pointer to the local register values mode (Input/Output)
 * @param *rn     Pointer to randomnumber array of the local node
*/
__device__ void thermalize_modes(float *mode, unsigned int index, LB_randomnr_gpu *rn){
  float Rho;
#ifdef SHANCHEN
  float Rho_tot=0.0,c;
  #pragma unroll
  for(int ii=0;ii<LB_COMPONENTS;++ii) { 
      Rho_tot  += mode[0 + ii * LBQ]+ para.rho[ii]*para.agrid*para.agrid*para.agrid;
  }
  c = (mode[0 + 0 * LBQ]+ para.rho[0]*para.agrid*para.agrid*para.agrid ) / Rho_tot;
  random_wrapper(rn);
  for(int ii=0;ii<LB_COMPONENTS;++ii) { 
      mode[1 + ii * LBQ] +=  sqrtf(c*(1-c)*Rho_tot*(para.mu[ii]*(2.0f/3.0f)*(1.0f-(para.gamma_mobility[0]*para.gamma_mobility[0])))) * (2*ii-1) * rn->randomnr[0];
      mode[2 + ii * LBQ] +=  sqrtf(c*(1-c)*Rho_tot*(para.mu[ii]*(2.0f/3.0f)*(1.0f-(para.gamma_mobility[0]*para.gamma_mobility[0])))) * (2*ii-1) * rn->randomnr[1];
  }                                      
  random_wrapper(rn);                    
  for(int ii=0;ii<LB_COMPONENTS;++ii)    
      mode[3 + ii * LBQ] +=  sqrtf(c*(1-c)*Rho_tot*(para.mu[ii]*(2.0f/3.0f)*(1.0f-( para.gamma_mobility[0]*para.gamma_mobility[0])))) * (2*ii-1) * rn->randomnr[0];
#endif
  
  
  for(int ii=0;ii<LB_COMPONENTS;++ii)
  {  

    /** mass mode */  
    Rho = mode[0 + ii * LBQ] + para.rho[ii]*para.agrid*para.agrid*para.agrid;

    /** momentum modes */

    /** stress modes */
    random_wrapper(rn);
    mode[4 + ii * LBQ] += sqrtf(Rho*(para.mu[ii]*(2.0f/3.0f)*(1.0f-( para.gamma_bulk[ii]* para.gamma_bulk[ii])))) * rn->randomnr[0];
    mode[5 + ii * LBQ] += sqrtf(Rho*(para.mu[ii]*(4.0f/9.0f)*(1.0f-(para.gamma_shear[ii]*para.gamma_shear[ii])))) * rn->randomnr[1];

    random_wrapper(rn);
    mode[6 + ii * LBQ] += sqrtf(Rho*(para.mu[ii]*(4.0f/3.0f)*(1.0f-(para.gamma_shear[ii]*para.gamma_shear[ii])))) * rn->randomnr[0];
    mode[7 + ii * LBQ] += sqrtf(Rho*(para.mu[ii]*(1.0f/9.0f)*(1.0f-(para.gamma_shear[ii]*para.gamma_shear[ii])))) * rn->randomnr[1];

    random_wrapper(rn);
    mode[8 + ii * LBQ] += sqrtf(Rho*(para.mu[ii]*(1.0f/9.0f)*(1.0f-(para.gamma_shear[ii]*para.gamma_shear[ii])))) * rn->randomnr[0];
    mode[9 + ii * LBQ] += sqrtf(Rho*(para.mu[ii]*(1.0f/9.0f)*(1.0f-(para.gamma_shear[ii]*para.gamma_shear[ii])))) * rn->randomnr[1];

    /** ghost modes */
    random_wrapper(rn);
    mode[10 + ii * LBQ] += sqrtf(Rho*(para.mu[ii]*(2.0f/3.0f)*(1.0f-(para.gamma_odd[ii]*para.gamma_odd[ii])))) * rn->randomnr[0];
    mode[11 + ii * LBQ] += sqrtf(Rho*(para.mu[ii]*(2.0f/3.0f)*(1.0f-(para.gamma_odd[ii]*para.gamma_odd[ii])))) * rn->randomnr[1];

    random_wrapper(rn);
    mode[12 + ii * LBQ] += sqrtf(Rho*(para.mu[ii]*(2.0f/3.0f)*(1.0f-(para.gamma_odd[ii]*para.gamma_odd[ii])))) * rn->randomnr[0];
    mode[13 + ii * LBQ] += sqrtf(Rho*(para.mu[ii]*(2.0f/9.0f)*(1.0f-(para.gamma_odd[ii]*para.gamma_odd[ii])))) * rn->randomnr[1];

    random_wrapper(rn);
    mode[14 + ii * LBQ] += sqrtf(Rho*(para.mu[ii]*(2.0f/9.0f)*(1.0f-(para.gamma_odd[ii]*para.gamma_odd[ii])))) * rn->randomnr[0];
    mode[15 + ii * LBQ] += sqrtf(Rho*(para.mu[ii]*(2.0f/9.0f)*(1.0f-(para.gamma_odd[ii]*para.gamma_odd[ii])))) * rn->randomnr[1];

    random_wrapper(rn);
    mode[16 + ii * LBQ] += sqrtf(Rho*(para.mu[ii]*(2.0f)*(1.0f-(para.gamma_even[ii]*para.gamma_even[ii]))))     * rn->randomnr[0];
    mode[17 + ii * LBQ] += sqrtf(Rho*(para.mu[ii]*(4.0f/9.0f)*(1.0f-(para.gamma_even[ii]*para.gamma_even[ii])))) * rn->randomnr[1];

    random_wrapper(rn);
    mode[18 + ii * LBQ] += sqrtf(Rho*(para.mu[ii]*(4.0f/3.0f)*(1.0f-(para.gamma_even[ii]*para.gamma_even[ii])))) * rn->randomnr[0];
  }
}


/*-------------------------------------------------------*/
/**normalization of the modes need befor backtransformation into velocity space
 * @param mode    Pointer to the local register values mode (Input/Output)
*/
__device__ void normalize_modes(float* mode){
  #pragma unroll
  for(int ii=0;ii<LB_COMPONENTS;++ii)
  { 
    /** normalization factors enter in the back transformation */
    mode[ 0 + ii * LBQ] *= 1.0f;
    mode[ 1 + ii * LBQ] *= 3.0f;
    mode[ 2 + ii * LBQ] *= 3.0f;
    mode[ 3 + ii * LBQ] *= 3.0f;
    mode[ 4 + ii * LBQ] *= 3.0f/2.0f;
    mode[ 5 + ii * LBQ] *= 9.0f/4.0f;
    mode[ 6 + ii * LBQ] *= 3.0f/4.0f;
    mode[ 7 + ii * LBQ] *= 9.0f;
    mode[ 8 + ii * LBQ] *= 9.0f;
    mode[ 9 + ii * LBQ] *= 9.0f;
    mode[10 + ii * LBQ] *= 3.0f/2.0f;
    mode[11 + ii * LBQ] *= 3.0f/2.0f;
    mode[12 + ii * LBQ] *= 3.0f/2.0f;
    mode[13 + ii * LBQ] *= 9.0f/2.0f;
    mode[14 + ii * LBQ] *= 9.0f/2.0f;
    mode[15 + ii * LBQ] *= 9.0f/2.0f;
    mode[16 + ii * LBQ] *= 1.0f/2.0f;
    mode[17 + ii * LBQ] *= 9.0f/4.0f;
    mode[18 + ii * LBQ] *= 3.0f/4.0f;
  }
}



/*-------------------------------------------------------*/
/**backtransformation from modespace to desityspace and streaming with the push method using pbc
 * @param index   node index / thread index (Input)
 * @param mode    Pointer to the local register values mode (Input)
 * @param *n_b    Pointer to local node residing in array b (Output)
*/
__device__ void calc_n_from_modes_push(LB_nodes_gpu n_b, float *mode, unsigned int index){

  unsigned int xyz[3];
  index_to_xyz(index, xyz);
  unsigned int x = xyz[0];
  unsigned int y = xyz[1];
  unsigned int z = xyz[2];

  #pragma unroll
  for(int ii=0;ii<LB_COMPONENTS;++ii)
  {
 
    n_b.vd[(0 + ii*LBQ ) * para.number_of_nodes + x 
                                                + para.dim_x*y
                                                + para.dim_x*para.dim_y*z] = 
      1.0f/3.0f * (mode[0 + ii * LBQ] - mode[4 + ii * LBQ] + mode[16 + ii * LBQ]);

    n_b.vd[(1 + ii*LBQ ) * para.number_of_nodes + (x+1)%para.dim_x
                                                + para.dim_x*y 
                                                + para.dim_x*para.dim_y*z] = 
      1.0f/18.0f * (
                       mode[ 0 + ii * LBQ] + mode[ 1 + ii * LBQ]
                     + mode[ 5 + ii * LBQ] + mode[ 6 + ii * LBQ]
                     - mode[17 + ii * LBQ] - mode[18 + ii * LBQ]
                     - 2.0f*(mode[10 + ii * LBQ] + mode[16 + ii * LBQ])
                   );

    n_b.vd[(2 + ii*LBQ ) * para.number_of_nodes + (para.dim_x+x-1)%para.dim_x
                                                + para.dim_x*y
                                                + para.dim_x*para.dim_y*z] =
      1.0f/18.0f * (
                       mode[ 0 + ii * LBQ] - mode[ 1 + ii * LBQ]
                     + mode[ 5 + ii * LBQ] + mode[ 6 + ii * LBQ]
                     - mode[17 + ii * LBQ] - mode[18 + ii * LBQ]
                     + 2.0f*(mode[10 + ii * LBQ] - mode[16 + ii * LBQ])
                   );

    n_b.vd[(3 + ii*LBQ ) * para.number_of_nodes + x
                                                + para.dim_x*((y+1)%para.dim_y)
                                                + para.dim_x*para.dim_y*z] =
      1.0f/18.0f * (
                       mode[ 0 + ii * LBQ] + mode[ 2 + ii * LBQ]
                     - mode[ 5 + ii * LBQ] + mode[ 6 + ii * LBQ]
                     + mode[17 + ii * LBQ] - mode[18 + ii * LBQ]
                     - 2.0f*(mode[11 + ii * LBQ] + mode[16 + ii * LBQ])
                   );

    n_b.vd[(4 + ii*LBQ ) * para.number_of_nodes + x
                                                + para.dim_x*((para.dim_y+y-1)%para.dim_y)
                                                + para.dim_x*para.dim_y*z] =
      1.0f/18.0f * (
                       mode[ 0 + ii * LBQ] - mode[ 2 + ii * LBQ]
                     - mode[ 5 + ii * LBQ] + mode[ 6 + ii * LBQ]
                     + mode[17 + ii * LBQ] - mode[18 + ii * LBQ]
                     + 2.0f*(mode[11 + ii * LBQ] - mode[16 + ii * LBQ])
                   );

    n_b.vd[(5 + ii*LBQ ) * para.number_of_nodes + x
                                                + para.dim_x*y
                                                + para.dim_x*para.dim_y*((z+1)%para.dim_z)] =
      1.0f/18.0f * (
                       mode[0 + ii * LBQ] + mode[3 + ii * LBQ]
                     - 2.0f*(   mode[ 6 + ii * LBQ] + mode[12 + ii * LBQ]
                              + mode[16 + ii * LBQ] - mode[18 + ii * LBQ])
                   );

    n_b.vd[(6 + ii*LBQ ) * para.number_of_nodes + x
                                                + para.dim_x*y
                                                + para.dim_x*para.dim_y*((para.dim_z+z-1)%para.dim_z)] =
      1.0f/18.0f * (
                       mode[0 + ii * LBQ] - mode[3 + ii * LBQ]
                     - 2.0f*(   mode[6 + ii * LBQ] - mode[12 + ii * LBQ]
                              + mode[16 + ii * LBQ] - mode[18 + ii * LBQ])
                   );

    n_b.vd[(7 + ii*LBQ ) * para.number_of_nodes + (x+1)%para.dim_x
                                                + para.dim_x*((y+1)%para.dim_y)
                                                + para.dim_x*para.dim_y*z] =
      1.0f/36.0f * (
                       mode[ 0 + ii * LBQ] + mode[ 1 + ii * LBQ]
                     + mode[ 2 + ii * LBQ] + mode[ 4 + ii * LBQ]
                     + 2.0f*mode[ 6 + ii * LBQ] + mode[ 7 + ii * LBQ]
                     + mode[10 + ii * LBQ] + mode[11 + ii * LBQ]
                     + mode[13 + ii * LBQ] + mode[14 + ii * LBQ]
                     + mode[16 + ii * LBQ] + 2.0f*mode[18 + ii * LBQ]
                   );

    n_b.vd[(8 + ii*LBQ ) * para.number_of_nodes + (para.dim_x+x-1)%para.dim_x
                                                + para.dim_x*((para.dim_y+y-1)%para.dim_y)
                                                + para.dim_x*para.dim_y*z] =
      1.0f/36.0f * (
                       mode[ 0 + ii * LBQ] - mode[ 1 + ii * LBQ]
                     - mode[ 2 + ii * LBQ] + mode[ 4 + ii * LBQ]
                     + 2.0f*mode[ 6 + ii * LBQ] + mode[ 7 + ii * LBQ]
                     - mode[10 + ii * LBQ] - mode[11 + ii * LBQ]
                     - mode[13 + ii * LBQ] - mode[14 + ii * LBQ]
                     + mode[16 + ii * LBQ] + 2.0f*mode[18 + ii * LBQ]
                   );

    n_b.vd[(9 + ii*LBQ ) * para.number_of_nodes + (x+1)%para.dim_x
                                                + para.dim_x*((para.dim_y+y-1)%para.dim_y)
                                                + para.dim_x*para.dim_y*z] =
      1.0f/36.0f * (
                       mode[ 0 + ii * LBQ] + mode[ 1 + ii * LBQ]
                     - mode[ 2 + ii * LBQ] + mode[ 4 + ii * LBQ]
                     + 2.0f*mode[ 6 + ii * LBQ] - mode[ 7 + ii * LBQ]
                     + mode[10 + ii * LBQ] - mode[11 + ii * LBQ]
                     + mode[13 + ii * LBQ] - mode[14 + ii * LBQ]
                     + mode[16 + ii * LBQ] + 2.0f*mode[18 + ii * LBQ]
                   );

    n_b.vd[(10 + ii*LBQ ) * para.number_of_nodes + (para.dim_x+x-1)%para.dim_x
                                                 + para.dim_x*((y+1)%para.dim_y)
                                                 + para.dim_x*para.dim_y*z] = 
      1.0f/36.0f * (
                       mode[ 0 + ii * LBQ] - mode[ 1 + ii * LBQ]
                     + mode[ 2 + ii * LBQ] + mode[ 4 + ii * LBQ]
                     + 2.0f*mode[ 6 + ii * LBQ] - mode[ 7 + ii * LBQ]
                     - mode[10 + ii * LBQ] + mode[11 + ii * LBQ]
                     - mode[13 + ii * LBQ] + mode[14 + ii * LBQ]
                     + mode[16 + ii * LBQ] + 2.0f*mode[18 + ii * LBQ]
                   );

    n_b.vd[(11 + ii*LBQ ) * para.number_of_nodes + (x+1)%para.dim_x
                                                 + para.dim_x*y
                                                 + para.dim_x*para.dim_y*((z+1)%para.dim_z)] =
      1.0f/36.0f * (
                       mode[ 0 + ii * LBQ] + mode[ 1 + ii * LBQ]
                     + mode[ 3 + ii * LBQ] + mode[ 4 + ii * LBQ]
                     + mode[ 5 + ii * LBQ] - mode[ 6 + ii * LBQ]
                     + mode[ 8 + ii * LBQ] + mode[10 + ii * LBQ]
                     + mode[12 + ii * LBQ] - mode[13 + ii * LBQ]
                     + mode[15 + ii * LBQ] + mode[16 + ii * LBQ]
                     + mode[17 + ii * LBQ] - mode[18 + ii * LBQ]
                   );

    n_b.vd[(12 + ii*LBQ ) * para.number_of_nodes + (para.dim_x+x-1)%para.dim_x
                                                 + para.dim_x*y
                                                 + para.dim_x*para.dim_y*((para.dim_z+z-1)%para.dim_z)] =
      1.0f/36.0f * (
                       mode[ 0 + ii * LBQ] - mode[ 1 + ii * LBQ]
                     - mode[ 3 + ii * LBQ] + mode[ 4 + ii * LBQ]
                     + mode[ 5 + ii * LBQ] - mode[ 6 + ii * LBQ]
                     + mode[ 8 + ii * LBQ] - mode[10 + ii * LBQ]
                     - mode[12 + ii * LBQ] + mode[13 + ii * LBQ]
                     - mode[15 + ii * LBQ] + mode[16 + ii * LBQ]
                     + mode[17 + ii * LBQ] - mode[18 + ii * LBQ]
                   );

    n_b.vd[(13 + ii*LBQ ) * para.number_of_nodes + (x+1)%para.dim_x
                                                 + para.dim_x*y
                                                 + para.dim_x*para.dim_y*((para.dim_z+z-1)%para.dim_z)] =
      1.0f/36.0f * (
                       mode[ 0 + ii * LBQ] + mode[ 1 + ii * LBQ]
                     - mode[ 3 + ii * LBQ] + mode[ 4 + ii * LBQ]
                     + mode[ 5 + ii * LBQ] - mode[ 6 + ii * LBQ]
                     - mode[ 8 + ii * LBQ] + mode[10 + ii * LBQ]
                     - mode[12 + ii * LBQ] - mode[13 + ii * LBQ]
                     - mode[15 + ii * LBQ] + mode[16 + ii * LBQ]
                     + mode[17 + ii * LBQ] - mode[18 + ii * LBQ]
                   );

    n_b.vd[(14 + ii*LBQ ) * para.number_of_nodes + (para.dim_x+x-1)%para.dim_x
                                                 + para.dim_x*y
                                                 + para.dim_x*para.dim_y*((z+1)%para.dim_z)] =
      1.0f/36.0f * (
                       mode[ 0 + ii * LBQ] - mode[ 1 + ii * LBQ]
                     + mode[ 3 + ii * LBQ] + mode[ 4 + ii * LBQ]
                     + mode[ 5 + ii * LBQ] - mode[ 6 + ii * LBQ]
                     - mode[ 8 + ii * LBQ] - mode[10 + ii * LBQ]
                     + mode[12 + ii * LBQ] + mode[13 + ii * LBQ]
                     + mode[15 + ii * LBQ] + mode[16 + ii * LBQ]
                     + mode[17 + ii * LBQ] - mode[18 + ii * LBQ]
                   );

    n_b.vd[(15 + ii*LBQ ) * para.number_of_nodes + x
                                                 + para.dim_x*((y+1)%para.dim_y)
                                                 + para.dim_x*para.dim_y*((z+1)%para.dim_z)] =
      1.0f/36.0f * (
                       mode[ 0 + ii * LBQ] + mode[ 2 + ii * LBQ]
                     + mode[ 3 + ii * LBQ] + mode[ 4 + ii * LBQ]
                     - mode[ 5 + ii * LBQ] - mode[ 6 + ii * LBQ]
                     + mode[ 9 + ii * LBQ] + mode[11 + ii * LBQ]
                     + mode[12 + ii * LBQ] - mode[14 + ii * LBQ]
                     - mode[15 + ii * LBQ] + mode[16 + ii * LBQ]
                     - mode[17 + ii * LBQ] - mode[18 + ii * LBQ]
                   );

    n_b.vd[(16 + ii*LBQ ) * para.number_of_nodes + x
                                                 + para.dim_x*((para.dim_y+y-1)%para.dim_y)
                                                 + para.dim_x*para.dim_y*((para.dim_z+z-1)%para.dim_z)] =
      1.0f/36.0f * (
                       mode[ 0 + ii * LBQ] - mode[ 2 + ii * LBQ]
                     - mode[ 3 + ii * LBQ] + mode[ 4 + ii * LBQ]
                     - mode[ 5 + ii * LBQ] - mode[ 6 + ii * LBQ]
                     + mode[ 9 + ii * LBQ] - mode[11 + ii * LBQ]
                     - mode[12 + ii * LBQ] + mode[14 + ii * LBQ]
                     + mode[15 + ii * LBQ] + mode[16 + ii * LBQ]
                     - mode[17 + ii * LBQ] - mode[18 + ii * LBQ]
                   );

    n_b.vd[(17 + ii*LBQ ) * para.number_of_nodes + x
                                                 + para.dim_x*((y+1)%para.dim_y)
                                                 + para.dim_x*para.dim_y*((para.dim_z+z-1)%para.dim_z)] =
      1.0f/36.0f * (
                       mode[ 0 + ii * LBQ] + mode[ 2 + ii * LBQ]
                     - mode[ 3 + ii * LBQ] + mode[ 4 + ii * LBQ]
                     - mode[ 5 + ii * LBQ] - mode[ 6 + ii * LBQ]
                     - mode[ 9 + ii * LBQ] + mode[11 + ii * LBQ]
                     - mode[12 + ii * LBQ] - mode[14 + ii * LBQ]
                     + mode[15 + ii * LBQ] + mode[16 + ii * LBQ]
                     - mode[17 + ii * LBQ] - mode[18 + ii * LBQ]
                   );

    n_b.vd[(18 + ii*LBQ ) * para.number_of_nodes + x
                                                 + para.dim_x*((para.dim_y+y-1)%para.dim_y)
                                                 + para.dim_x*para.dim_y*((z+1)%para.dim_z)] =
      1.0f/36.0f * (
                       mode[ 0 + ii * LBQ] - mode[ 2 + ii * LBQ]
                     + mode[ 3 + ii * LBQ] + mode[ 4 + ii * LBQ]
                     - mode[ 5 + ii * LBQ] - mode[ 6 + ii * LBQ]
                     - mode[ 9 + ii * LBQ] - mode[11 + ii * LBQ]
                     + mode[12 + ii * LBQ] + mode[14 + ii * LBQ]
                     - mode[15 + ii * LBQ] + mode[16 + ii * LBQ]
                     - mode[17 + ii * LBQ] - mode[18 + ii * LBQ]
                   );
  }
}


/* bounce_back function used in bounce_back_boundaries
** called 18times, once for each velocity
*/
__device__ __inline__ void bounce_back(LB_nodes_gpu n_curr, LB_rho_v_gpu *d_v, float *boundary_force, float* torque, unsigned int index, float weight,  int *c, float *v, float* omega, int population, int inverse, short moving_flag, float *to_center, unsigned int x, unsigned int y, unsigned int z){
  
  size_t to_index, to_index_x, to_index_y, to_index_z;
  float r[3];
  float loc_v[3];
  float f[3];
  float shift, pop_to_bounce_back;
  float vc;
  
#ifndef SHANCHEN
  //+ 4.5f*uc*uc + u_sq_neg
  pop_to_bounce_back =  n_curr.vd[population*para.number_of_nodes + index ];
  to_index_x = (x+c[0]+para.dim_x)%para.dim_x;
  to_index_y = (y+c[1]+para.dim_y)%para.dim_y;
  to_index_z = (z+c[2]+para.dim_z)%para.dim_z;
  to_index = to_index_x + para.dim_x*to_index_y + para.dim_x*para.dim_y*to_index_z;
  if ( n_curr.boundary[to_index] == 0)
  {
    r[0] = to_center[0] + 0.5f * c[0];//the bounceback happens at half-way and is treated as such
    r[1] = to_center[1] + 0.5f * c[1];
    r[2] = to_center[2] + 0.5f * c[2];
    loc_v[0] = v[0] - omega[1]*r[2] + omega[2]*r[1];//change of sign since r points to the center.
    loc_v[1] = v[1] - omega[2]*r[0] + omega[0]*r[2];
    loc_v[2] = v[2] - omega[0]*r[1] + omega[1]*r[0];
    
    vc = loc_v[0]*c[0] + loc_v[1]*c[1] + loc_v[2]*c[2];
    shift = weight * d_v[to_index].rho[0] * 6.0f*vc;
    n_curr.vd[inverse*para.number_of_nodes + to_index ] = pop_to_bounce_back + shift;
    pop_to_bounce_back *= 2;
    f[0] = (pop_to_bounce_back + shift) * c[0];
    f[1] = (pop_to_bounce_back + shift) * c[1];
    f[2] = (pop_to_bounce_back + shift) * c[2];
    torque[0] -= (r[1]*f[2] - r[2]*f[1]);//change of sign because r to center
    torque[1] -= (r[2]*f[0] - r[0]*f[2]); 
    torque[2] -= (r[0]*f[1] - r[1]*f[0]); 
    boundary_force[0] += f[0];
    boundary_force[1] += f[1];
    boundary_force[2] += f[2];
  }

 
#else 
  for(int component=0; component<LB_COMPONENTS;component++){
     shift = 2.0f*para.agrid*para.agrid*para.rho[component]*3.0f*weight*para.tau*(v[0]*c[0] + v[1]*c[1] + v[2]*c[2]); 
     pop_to_bounce_back =  n_curr.vd[(population+component*LBQ)*para.number_of_nodes + index ]; 
     to_index_x = (x+c[0]+para.dim_x)%para.dim_x; 
     to_index_y = (y+c[1]+para.dim_y)%para.dim_y; 
     to_index_z = (z+c[2]+para.dim_z)%para.dim_z; 
     to_index = to_index_x + para.dim_x*to_index_y + para.dim_x*para.dim_y*to_index_z; 
     if ( n_curr.boundary[to_index] == 0) 
     { 
       boundary_force[0] += (2.0f*pop_to_bounce_back+shift)*c[0]/para.tau/para.tau/para.agrid; 
       boundary_force[1] += (2.0f*pop_to_bounce_back+shift)*c[1]/para.tau/para.tau/para.agrid; 
       boundary_force[2] += (2.0f*pop_to_bounce_back+shift)*c[2]/para.tau/para.tau/para.agrid; 
       n_curr.vd[(inverse+component*LBQ)*para.number_of_nodes + to_index ] = pop_to_bounce_back + shift; 
     } 
  }
#endif // not SHANCHEN
}

/** Bounce back boundary conditions.
 * The populations that have propagated into a boundary node
 * are bounced back to the node they came from. This results
 * in no slip boundary conditions.
 *
 * [cf. Ladd and Verberg, J. Stat. Phys. 104(5/6):1191-1251, 2001]
 * @param index   node index / thread index (Input)
 * @param n_curr  Pointer to local node which receives the current node field (Input)
 * @param lb_boundary_velocity    The constant velocity at the boundary, set by the user (Input)
 * @param lb_boundary_force       The force on the boundary nodes (Output)
*/
__device__ void bounce_back_boundaries(LB_nodes_gpu n_curr, LB_rho_v_gpu *d_v, LB_moving_boundary *lb_moving_boundary, unsigned int index, float* lb_boundary_velocity, float* lb_boundary_force){
    
  int c[3];
  float v[3];
  float omega[3];
  float weight;
  float boundary_force[3] = {0.0f,0.0f,0.0f};
  float torque[3] = {0.0f, 0.0f, 0.0f};
  int population, inverse;
  int boundary_index;
  short moving_flag;

  boundary_index= n_curr.boundary[index];
  if(boundary_index)
  {
    unsigned int xyz[3];
    float delta_xyz[3];
    index_to_xyz(index, xyz);
    
    unsigned int x = xyz[0];
    unsigned int y = xyz[1];
    unsigned int z = xyz[2];
    if(boundary_index > 0){
      moving_flag = 0;
      v[0]=lb_boundary_velocity[0 + 3*(boundary_index-1)];
      v[1]=lb_boundary_velocity[1 + 3*(boundary_index-1)];
      v[2]=lb_boundary_velocity[2 + 3*(boundary_index-1)];
      omega[0]=0;
      omega[1]=0;
      omega[2]=0;
    } 
    else {
      boundary_index = -(boundary_index + 1);
      moving_flag = 1;
      v[0] = lb_moving_boundary[boundary_index].velocity[0];
      v[1] = lb_moving_boundary[boundary_index].velocity[1];
      v[2] = lb_moving_boundary[boundary_index].velocity[2];
      omega[0] = lb_moving_boundary[boundary_index].omega[0];
      omega[1] = lb_moving_boundary[boundary_index].omega[1];
      omega[2] = lb_moving_boundary[boundary_index].omega[2];
      xyz_center_to_delta(xyz, lb_moving_boundary[boundary_index].center, delta_xyz);
    }


    /* CPU analog of shift:
       lbpar.agrid*lbpar.agrid*lbpar.agrid*lbpar.rho*2*lbmodel.c[i][l]*lb_boundaries[lbfields[k].boundary-1].velocity[l] */
  
    /** store vd temporary in second lattice to avoid race conditions */

// TODO : PUT IN EQUILIBRIUM CONTRIBUTION TO THE BOUNCE-BACK DENSITY FOR THE BOUNDARY FORCE
// TODO : INITIALIZE BOUNDARY FORCE PROPERLY, HAS NONZERO ELEMENTS IN FIRST STEP

#define W1 0.0555555555f
#define W2 0.0277777777f

    // the resting population does nothing, i.e., population 0.
    c[0]= 1;c[1]= 0;c[2]= 0; weight=W1; population= 2; inverse= 1; 
    bounce_back(n_curr, d_v, boundary_force, torque, index, weight, c, v, omega, population, inverse, moving_flag, delta_xyz, x, y, z);
    
    c[0]=-1;c[1]= 0;c[2]= 0; weight=W1; population= 1; inverse= 2; 
    bounce_back(n_curr, d_v, boundary_force, torque, index, weight, c, v, omega, population, inverse, moving_flag, delta_xyz, x, y, z);
    
    
    c[0]= 0;c[1]= 1;c[2]= 0; weight=W1; population= 4; inverse= 3; 
    bounce_back(n_curr, d_v, boundary_force, torque, index, weight, c, v, omega, population, inverse, moving_flag, delta_xyz, x, y, z);

    c[0]= 0;c[1]=-1;c[2]= 0; weight=W1; population= 3; inverse= 4; 
    bounce_back(n_curr, d_v, boundary_force, torque, index, weight, c, v, omega, population, inverse, moving_flag, delta_xyz, x, y, z);
    
    
    c[0]= 0;c[1]= 0;c[2]= 1; weight=W1; population= 6; inverse= 5; 
    bounce_back(n_curr, d_v, boundary_force, torque, index, weight, c, v, omega, population, inverse, moving_flag, delta_xyz, x, y, z);

    c[0]= 0;c[1]= 0;c[2]=-1; weight=W1; population= 5; inverse= 6; 
    bounce_back(n_curr, d_v, boundary_force, torque, index, weight, c, v, omega, population, inverse, moving_flag, delta_xyz, x, y, z); 
    
    
    c[0]= 1;c[1]= 1;c[2]= 0; weight=W2; population= 8; inverse= 7; 
    bounce_back(n_curr, d_v, boundary_force, torque, index, weight, c, v, omega, population, inverse, moving_flag, delta_xyz, x, y, z);
    
    c[0]=-1;c[1]=-1;c[2]= 0; weight=W2; population= 7; inverse= 8; 
    bounce_back(n_curr, d_v, boundary_force, torque, index, weight, c, v, omega, population, inverse, moving_flag, delta_xyz, x, y, z);
    
    
    c[0]= 1;c[1]=-1;c[2]= 0; weight=W2; population=10; inverse= 9; 
    bounce_back(n_curr, d_v, boundary_force, torque, index, weight, c, v, omega, population, inverse, moving_flag, delta_xyz, x, y, z);

    c[0]=-1;c[1]= 1;c[2]= 0; weight=W2; population= 9; inverse=10; 
    bounce_back(n_curr, d_v, boundary_force, torque, index, weight, c, v, omega, population, inverse, moving_flag, delta_xyz, x, y, z);
    
    
    c[0]= 1;c[1]= 0;c[2]= 1; weight=W2; population=12; inverse=11; 
    bounce_back(n_curr, d_v, boundary_force, torque, index, weight, c, v, omega, population, inverse, moving_flag, delta_xyz, x, y, z);
    
    c[0]=-1;c[1]= 0;c[2]=-1; weight=W2; population=11; inverse=12; 
    bounce_back(n_curr, d_v, boundary_force, torque, index, weight, c, v, omega, population, inverse, moving_flag, delta_xyz, x, y, z);
    

    c[0]= 1;c[1]= 0;c[2]=-1; weight=W2; population=14; inverse=13; 
    bounce_back(n_curr, d_v, boundary_force, torque, index, weight, c, v, omega, population, inverse, moving_flag, delta_xyz, x, y, z);
    
    c[0]=-1;c[1]= 0;c[2]= 1; weight=W2; population=13; inverse=14; 
    bounce_back(n_curr, d_v, boundary_force, torque, index, weight, c, v, omega, population, inverse, moving_flag, delta_xyz, x, y, z);
    

    c[0]= 0;c[1]= 1;c[2]= 1; weight=W2; population=16; inverse=15; 
    bounce_back(n_curr, d_v, boundary_force, torque, index, weight, c, v, omega, population, inverse, moving_flag, delta_xyz, x, y, z);
    
    c[0]= 0;c[1]=-1;c[2]=-1; weight=W2; population=15; inverse=16; 
    bounce_back(n_curr, d_v, boundary_force, torque, index, weight, c, v, omega, population, inverse, moving_flag, delta_xyz, x, y, z);
    
    
    c[0]= 0;c[1]= 1;c[2]=-1; weight=W2; population=18; inverse=17; 
    bounce_back(n_curr, d_v, boundary_force, torque, index, weight, c, v, omega, population, inverse, moving_flag, delta_xyz, x, y, z);
    
    c[0]= 0;c[1]=-1;c[2]= 1; weight=W2; population=17; inverse=18; 
    bounce_back(n_curr, d_v, boundary_force, torque, index, weight, c, v, omega, population, inverse, moving_flag, delta_xyz, x, y, z); 
    
    
#undef W1
#undef W2
    if(moving_flag){
      atomicAdd(lb_moving_boundary[boundary_index].force_hyd + 0, -boundary_force[0] );
      atomicAdd(lb_moving_boundary[boundary_index].force_hyd + 1, -boundary_force[1] );
      atomicAdd(lb_moving_boundary[boundary_index].force_hyd + 2, -boundary_force[2] );
      atomicAdd(lb_moving_boundary[boundary_index].torque + 0, -torque[0] );
      atomicAdd(lb_moving_boundary[boundary_index].torque + 1, -torque[1] );
      atomicAdd(lb_moving_boundary[boundary_index].torque + 2, -torque[2] );
    }
  }
}

/** add of (external) forces within the modespace, needed for particle-interaction
 * @param index   node index / thread index (Input)
 * @param mode    Pointer to the local register values mode (Input/Output)
 * @param node_f  Pointer to local node force (Input)
 * @param *d_v    Pointer to local device values
*/
__device__ void apply_forces(unsigned int index, float *mode, LB_node_force_gpu node_f, LB_rho_v_gpu *d_v) {
  
  float u[3]={0.0f,0.0f,0.0f},
        C[6]={0.0f,0.0f,0.0f,0.0f,0.0f,0.0f};
  /* Note: the values d_v were calculated in relax_modes() */

  u[0]=d_v[index].v[0]; 
  u[1]=d_v[index].v[1]; 
  u[2]=d_v[index].v[2]; 

  #pragma unroll
  for(int ii=0;ii<LB_COMPONENTS;++ii)
  {  
       C[0] += (1.0f + para.gamma_bulk[ii])*u[0]*node_f.force[(0 + ii*3 ) * para.number_of_nodes + index] + 
                1.0f/3.0f*(para.gamma_bulk[ii]-para.gamma_shear[ii])*(
                                                                         u[0]*node_f.force[(0 + ii*3 ) * para.number_of_nodes + index]
                                                                       + u[1]*node_f.force[(1 + ii*3 ) * para.number_of_nodes + index]
                                                                       + u[2]*node_f.force[(2 + ii*3 ) * para.number_of_nodes + index]
                                                                     );

       C[2] += (1.0f + para.gamma_bulk[ii])*u[1]*node_f.force[(1 + ii*3 ) * para.number_of_nodes + index] + 
                1.0f/3.0f*(para.gamma_bulk[ii]-para.gamma_shear[ii])*(
                                                                         u[0]*node_f.force[(0 + ii*3 ) * para.number_of_nodes + index]
                                                                       + u[1]*node_f.force[(1 + ii*3 ) * para.number_of_nodes + index]
                                                                       + u[2]*node_f.force[(2 + ii*3 ) * para.number_of_nodes + index]
                                                                     );

       C[5] += (1.0f + para.gamma_bulk[ii])*u[2]*node_f.force[(2 + ii*3 ) * para.number_of_nodes + index] + 
                1.0f/3.0f*(para.gamma_bulk[ii]-para.gamma_shear[ii])*(
                                                                         u[0]*node_f.force[(0 + ii*3 ) * para.number_of_nodes + index] 
                                                                       + u[1]*node_f.force[(1 + ii*3 ) * para.number_of_nodes + index]
                                                                       + u[2]*node_f.force[(2 + ii*3 ) * para.number_of_nodes + index]
                                                                     );

       C[1] += 1.0f/2.0f*(1.0f+para.gamma_shear[ii])*(
                                                         u[0]*node_f.force[(1 + ii*3 ) * para.number_of_nodes + index]
                                                       + u[1]*node_f.force[(0 + ii*3 ) * para.number_of_nodes + index]
                                                     );

       C[3] += 1.0f/2.0f*(1.0f+para.gamma_shear[ii])*(
                                                         u[0]*node_f.force[(2 + ii*3 ) * para.number_of_nodes + index]
                                                       + u[2]*node_f.force[(0 + ii*3 ) * para.number_of_nodes + index]
                                                     );

       C[4] += 1.0f/2.0f*(1.0f+para.gamma_shear[ii])*(
                                                         u[1]*node_f.force[(2 + ii*3 ) * para.number_of_nodes + index]
                                                       + u[2]*node_f.force[(1 + ii*3 ) * para.number_of_nodes + index]
                                                     );
  }

  #pragma unroll
  for(int ii=0;ii<LB_COMPONENTS;++ii)
  {  
      /** update momentum modes */
#ifdef SHANCHEN
      float mobility_factor=1.0f/2.0f*(1.0f+para.gamma_mobility[0]);
#else
      float mobility_factor=1.0f;
#endif 
 /** update momentum modes */
      mode[1 + ii * LBQ] += mobility_factor * node_f.force[(0 + ii*3 ) * para.number_of_nodes + index];
      mode[2 + ii * LBQ] += mobility_factor * node_f.force[(1 + ii*3 ) * para.number_of_nodes + index];
      mode[3 + ii * LBQ] += mobility_factor * node_f.force[(2 + ii*3 ) * para.number_of_nodes + index];

      /** update stress modes */
      mode[4 + ii * LBQ] += C[0] + C[2] + C[5];
      mode[5 + ii * LBQ] += C[0] - C[2];
      mode[6 + ii * LBQ] += C[0] + C[2] - 2.0f*C[5];
      mode[7 + ii * LBQ] += C[1];
      mode[8 + ii * LBQ] += C[3];
      mode[9 + ii * LBQ] += C[4];
    
  }

//#if !defined(IMMERSED_BOUNDARY)
  // This must not be done here since we need the forces after LB update for the velocity interpolation
  // It is done by calling IBM_ResetLBForces_GPU from integrate_vv
  reset_LB_forces(index, node_f);
//#endif

#ifdef SHANCHEN
  for(int ii=0;ii<LB_COMPONENTS;++ii)
  {  
     node_f.force[(0 + ii*3 ) * para.number_of_nodes + index] +=node_f.scforce[(0+ii*3)*para.number_of_nodes + index];
     node_f.force[(1 + ii*3 ) * para.number_of_nodes + index] +=node_f.scforce[(1+ii*3)*para.number_of_nodes + index];
     node_f.force[(2 + ii*3 ) * para.number_of_nodes + index] +=node_f.scforce[(2+ii*3)*para.number_of_nodes + index];
  }
#endif
}

/**function used to calculate hydrodynamic fields in MD units.
 * @param n_a     Pointer to local node residing in array a for boundary flag(Input)
 * @param mode    Pointer to the local register values mode (Input)
 * @param d_p_v   Pointer to local print values (Output)
 * @param d_v     Pointer to local device values (Input)
 * @param node_f  Pointer to local node force (Input)
 * @param index   node index / thread index (Input)
 * @param print_index   node index / thread index (Output)
*/
__device__ void calc_values_in_MD_units(LB_nodes_gpu n_a, float *mode, LB_rho_v_pi_gpu *d_p_v, LB_rho_v_gpu *d_v, LB_node_force_gpu node_f, unsigned int index, unsigned int print_index) {
  
  float j[3]; 
  float modes_from_pi_eq[6]; 
  float pi[6]={0.0f,0.0f,0.0f,0.0f,0.0f,0.0f};

  if(n_a.boundary[index] == 0)
  {
    /* Ensure we are working with the current values of d_v */

    update_rho_v(mode, index, node_f, d_v);

    for(int ii= 0; ii < LB_COMPONENTS; ii++)
    {
      d_p_v[print_index].rho[ii] = d_v[index].rho[ii] / para.agrid / para.agrid / para.agrid;
    }
      
    d_p_v[print_index].v[0] = d_v[index].v[0] * para.agrid / para.tau;
    d_p_v[print_index].v[1] = d_v[index].v[1] * para.agrid / para.tau;
    d_p_v[print_index].v[2] = d_v[index].v[2] * para.agrid / para.tau;

    /* stress calculation */ 
    for(int ii = 0; ii < LB_COMPONENTS; ii++)
    {
      float Rho = d_v[index].rho[ii];
      
      /* note that d_v[index].v[] already includes the 1/2 f term, accounting for the pre- and post-collisional average */

      j[0] = Rho * d_v[index].v[0];
      j[1] = Rho * d_v[index].v[1];
      j[2] = Rho * d_v[index].v[2];

      // equilibrium part of the stress modes, which comes from 
      // the equality between modes and stress tensor components

      /* m4 = trace(pi) - rho
         m5 = pi_xx - pi_yy
         m6 = trace(pi) - 3 pi_zz
         m7 = pi_xy
         m8 = pi_xz
         m9 = pi_yz */

      // and pluggin in the Euler stress for the equilibrium:
      // pi_eq = rho_0*c_s^2*I3 + (j \otimes j)/rho
      // with I3 the 3D identity matrix and
      // rho = \trace(rho_0*c_s^2*I3), which yields

      /* m4_from_pi_eq = j.j
         m5_from_pi_eq = j_x*j_x - j_y*j_y
         m6_from_pi_eq = j.j - 3*j_z*j_z
         m7_from_pi_eq = j_x*j_y
         m8_from_pi_eq = j_x*j_z
         m9_from_pi_eq = j_y*j_z */

      // where the / Rho term has been dropped. We thus obtain: 

      modes_from_pi_eq[0] = ( j[0]*j[0] + j[1]*j[1] + j[2]*j[2] ) / Rho;
      modes_from_pi_eq[1] = ( j[0]*j[0] - j[1]*j[1] ) / Rho;
      modes_from_pi_eq[2] = ( j[0]*j[0] + j[1]*j[1] + j[2]*j[2] - 3.0f*j[2]*j[2] ) / Rho;
      modes_from_pi_eq[3] = j[0]*j[1] / Rho;
      modes_from_pi_eq[4] = j[0]*j[2] / Rho;
      modes_from_pi_eq[5] = j[1]*j[2] / Rho;
     
      /* Now we must predict the outcome of the next collision */
      /* We immediately average pre- and post-collision.  */
      /* TODO: need a reference for this.   */

      mode[4 + ii * LBQ ] = modes_from_pi_eq[0] + (0.5f + 0.5f*para.gamma_bulk[ii]) * (mode[4 + ii * LBQ] - modes_from_pi_eq[0]);
      mode[5 + ii * LBQ ] = modes_from_pi_eq[1] + (0.5f + 0.5f*para.gamma_shear[ii]) * (mode[5 + ii * LBQ] - modes_from_pi_eq[1]);
      mode[6 + ii * LBQ ] = modes_from_pi_eq[2] + (0.5f + 0.5f*para.gamma_shear[ii]) * (mode[6 + ii * LBQ] - modes_from_pi_eq[2]);
      mode[7 + ii * LBQ ] = modes_from_pi_eq[3] + (0.5f + 0.5f*para.gamma_shear[ii]) * (mode[7 + ii * LBQ] - modes_from_pi_eq[3]);
      mode[8 + ii * LBQ ] = modes_from_pi_eq[4] + (0.5f + 0.5f*para.gamma_shear[ii]) * (mode[8 + ii * LBQ] - modes_from_pi_eq[4]);
      mode[9 + ii * LBQ ] = modes_from_pi_eq[5] + (0.5f + 0.5f*para.gamma_shear[ii]) * (mode[9 + ii * LBQ] - modes_from_pi_eq[5]);

      // Transform the stress tensor components according to the modes that
      // correspond to those used by U. Schiller. In terms of populations this
      // expression then corresponds exactly to those in Eqs. 116 - 121 in the
      // Duenweg and Ladd paper, when these are written out in populations.
      // But to ensure this, the expression in Schiller's modes has to be different!

      pi[0] += (   2.0f*(mode[0 + ii * LBQ] + mode[4 + ii * LBQ])
                + mode[6 + ii * LBQ] + 3.0f*mode[5 + ii * LBQ] )/6.0f;  // xx
      pi[1] += mode[7 + ii * LBQ];                                      // xy
      pi[2] += (   2.0f*(mode[0 + ii * LBQ] + mode[4 + ii * LBQ])
                + mode[6 + ii * LBQ] - 3.0f*mode[5 + ii * LBQ] )/6.0f;  // yy
      pi[3] += mode[8 + ii * LBQ];                                      // xz
      pi[4] += mode[9 + ii * LBQ];                                      // yz
      pi[5] += (   mode[0 + ii * LBQ] + mode[4 + ii * LBQ]
                - mode[6 + ii * LBQ] )/3.0f;                            // zz

    }
     
    for(int i = 0; i < 6; i++)
    {
      d_p_v[print_index].pi[i] = pi[i] / para.tau
                                       / para.tau
                                       / para.agrid;
    }
  }
  else
  {
    for(int ii = 0; ii < LB_COMPONENTS; ii++)
      d_p_v[print_index].rho[ii] = 0.0f;
     
    for(int i = 0; i < 3; i++)
      d_p_v[print_index].v[i] = 0.0f;

    for(int i = 0; i < 6; i++)
      d_p_v[print_index].pi[i] = 0.0f;
  }
}

/**function used to calculate hydrodynamic fields in MD units.
 * @param mode_single   Pointer to the local register values mode (Input)
 * @param d_v_single    Pointer to local device values (Input)
 * @param rho_out       Pointer to density (Output)
 * @param j_out         Pointer to momentum (Output)
 * @param pi_out        Pointer to pressure tensor (Output)
*/
__device__ void calc_values_from_m_in_LB_units(float *mode_single, LB_rho_v_gpu *d_v_single, float* rho_out, float* j_out, float* pi_out) {

  float modes_from_pi_eq[6];
  float j[6];
  float Rho; 

  // stress calculation

  for(int ii = 0; ii < LB_COMPONENTS; ii++)
  {
    // Set the rho ouput value

    Rho = d_v_single->rho[ii];
    rho_out[ii] = Rho;
    
    // note that d_v_single->v[] already includes the 1/2 f term, 
    // accounting for the pre- and post-collisional average

    j[0] = Rho * d_v_single->v[0];
    j[1] = Rho * d_v_single->v[1];
    j[2] = Rho * d_v_single->v[2];

    j_out[3*ii + 0] = j[0];
    j_out[3*ii + 1] = j[1];
    j_out[3*ii + 2] = j[2];    

    // equilibrium part of the stress modes, which comes from 
    // the equality between modes and stress tensor components

    modes_from_pi_eq[0] = ( j[0]*j[0] + j[1]*j[1] + j[2]*j[2] ) / Rho;
    modes_from_pi_eq[1] = ( j[0]*j[0] - j[1]*j[1] ) / Rho;
    modes_from_pi_eq[2] = ( j[0]*j[0] + j[1]*j[1] + j[2]*j[2] - 3.0f*j[2]*j[2] ) / Rho;
    modes_from_pi_eq[3] = j[0]*j[1] / Rho;
    modes_from_pi_eq[4] = j[0]*j[2] / Rho;
    modes_from_pi_eq[5] = j[1]*j[2] / Rho;
   
    // Now we must predict the outcome of the next collision
    // We immediately average pre- and post-collision.

    mode_single[4 + ii * LBQ ] = modes_from_pi_eq[0] + (0.5f + 0.5f* para.gamma_bulk[ii]) * (mode_single[4 + ii * LBQ] - modes_from_pi_eq[0]);
    mode_single[5 + ii * LBQ ] = modes_from_pi_eq[1] + (0.5f + 0.5f*para.gamma_shear[ii]) * (mode_single[5 + ii * LBQ] - modes_from_pi_eq[1]);
    mode_single[6 + ii * LBQ ] = modes_from_pi_eq[2] + (0.5f + 0.5f*para.gamma_shear[ii]) * (mode_single[6 + ii * LBQ] - modes_from_pi_eq[2]);
    mode_single[7 + ii * LBQ ] = modes_from_pi_eq[3] + (0.5f + 0.5f*para.gamma_shear[ii]) * (mode_single[7 + ii * LBQ] - modes_from_pi_eq[3]);
    mode_single[8 + ii * LBQ ] = modes_from_pi_eq[4] + (0.5f + 0.5f*para.gamma_shear[ii]) * (mode_single[8 + ii * LBQ] - modes_from_pi_eq[4]);
    mode_single[9 + ii * LBQ ] = modes_from_pi_eq[5] + (0.5f + 0.5f*para.gamma_shear[ii]) * (mode_single[9 + ii * LBQ] - modes_from_pi_eq[5]);

    // Transform the stress tensor components according to the mode_singles.

    pi_out[6*ii + 0] = (   2.0f*(mode_single[0 + ii * LBQ] + mode_single[4 + ii * LBQ])
                         + mode_single[6 + ii * LBQ] + 3.0f*mode_single[5 + ii * LBQ] )/6.0f;   // xx
    pi_out[6*ii + 1] = mode_single[7 + ii * LBQ];                                               // xy
    pi_out[6*ii + 2] = (   2.0f*(mode_single[0 + ii * LBQ] + mode_single[4 + ii * LBQ])
                         + mode_single[6 + ii * LBQ] - 3.0f*mode_single[5 + ii * LBQ] )/6.0f;   // yy
    pi_out[6*ii + 3] = mode_single[8 + ii * LBQ];                                               // xz
    pi_out[6*ii + 4] = mode_single[9 + ii * LBQ];                                               // yz
    pi_out[6*ii + 5] = (   mode_single[0 + ii * LBQ] + mode_single[4 + ii * LBQ]
                         - mode_single[6 + ii * LBQ] )/3.0f;                                    // zz
  }
}

/**function used to calc physical values of every node
 * @param n_a     Pointer to local node residing in array a for boundary flag(Input)
 * @param mode    Pointer to the local register values mode (Input)
 * @param d_v     Pointer to local device values (Input/Output)
 * @param node_f  Pointer to local node force (Input)
 * @param index   node index / thread index (Input)
*/

/* FIXME this function is basically un-used, think about removing/replacing it */
__device__ void calc_values(LB_nodes_gpu n_a, float *mode, LB_rho_v_gpu *d_v, LB_node_force_gpu node_f, unsigned int index){ 

  float Rho_tot=0.0f;
  float u_tot[3]={0.0f,0.0f,0.0f};

  if(n_a.boundary[index] != 1){
      #pragma unroll
      for(int ii=0;ii<LB_COMPONENTS;++ii) { 
          /** re-construct the real density
          * remember that the populations are stored as differences to their
          * equilibrium value */
          d_v[index].rho[ii]= mode[0 + ii * 4]+ para.rho[ii]*para.agrid*para.agrid*para.agrid;
          Rho_tot  += mode[0 + ii * 4]+ para.rho[ii]*para.agrid*para.agrid*para.agrid;
          u_tot[0] += mode[1 + ii * 4];
          u_tot[1] += mode[2 + ii * 4];
          u_tot[2] += mode[3 + ii * 4];
    
          /** if forces are present, the momentum density is redefined to
          * inlcude one half-step of the force action.  See the
          * Chapman-Enskog expansion in [Ladd & Verberg]. */
    
          u_tot[0] += 0.5f*node_f.force[(0+ii*3)*para.number_of_nodes + index];
          u_tot[1] += 0.5f*node_f.force[(1+ii*3)*para.number_of_nodes + index];
          u_tot[2] += 0.5f*node_f.force[(2+ii*3)*para.number_of_nodes + index];
      }
      u_tot[0]/=Rho_tot;
      u_tot[1]/=Rho_tot;
      u_tot[2]/=Rho_tot;
    
      d_v[index].v[0]=u_tot[0]; 
      d_v[index].v[1]=u_tot[1]; 
      d_v[index].v[2]=u_tot[2]; 
  } else { 
    #pragma unroll
    for(int ii=0;ii<LB_COMPONENTS;++ii) { 
       d_v[index].rho[ii]   = 1.;
    }
    d_v[index].v[0] = 0.0f;
    d_v[index].v[1] = 0.0f; 
    d_v[index].v[2] = 0.0f; 
  }   
}


/** 
 * @param node_index  node index around (8) particle (Input)
 * @param *mode       Pointer to the local register values mode (Output)
 * @param n_a         Pointer to local node residing in array a(Input)
 * @param component_index   Shanchen component index        (Input)
*/
__device__ void calc_mode(float *mode, LB_nodes_gpu n_a, unsigned int node_index, int component_index){

  /** mass mode */
  mode[0] =   n_a.vd[( 0 + component_index*LBQ ) * para.number_of_nodes + node_index]
            + n_a.vd[( 1 + component_index*LBQ ) * para.number_of_nodes + node_index] + n_a.vd[( 2 + component_index*LBQ ) * para.number_of_nodes + node_index] 
            + n_a.vd[( 3 + component_index*LBQ ) * para.number_of_nodes + node_index] + n_a.vd[( 4 + component_index*LBQ ) * para.number_of_nodes + node_index]
            + n_a.vd[( 5 + component_index*LBQ ) * para.number_of_nodes + node_index] + n_a.vd[( 6 + component_index*LBQ ) * para.number_of_nodes + node_index]
            + n_a.vd[( 7 + component_index*LBQ ) * para.number_of_nodes + node_index] + n_a.vd[( 8 + component_index*LBQ ) * para.number_of_nodes + node_index]
            + n_a.vd[( 9 + component_index*LBQ ) * para.number_of_nodes + node_index] + n_a.vd[(10 + component_index*LBQ ) * para.number_of_nodes + node_index]
            + n_a.vd[(11 + component_index*LBQ ) * para.number_of_nodes + node_index] + n_a.vd[(12 + component_index*LBQ ) * para.number_of_nodes + node_index]
            + n_a.vd[(13 + component_index*LBQ ) * para.number_of_nodes + node_index] + n_a.vd[(14 + component_index*LBQ ) * para.number_of_nodes + node_index]
            + n_a.vd[(15 + component_index*LBQ ) * para.number_of_nodes + node_index] + n_a.vd[(16 + component_index*LBQ ) * para.number_of_nodes + node_index]
            + n_a.vd[(17 + component_index*LBQ ) * para.number_of_nodes + node_index] + n_a.vd[(18 + component_index*LBQ ) * para.number_of_nodes + node_index];

  /** momentum modes */
  mode[1] =   (n_a.vd[( 1 + component_index*LBQ ) * para.number_of_nodes + node_index] - n_a.vd[( 2 + component_index*LBQ ) * para.number_of_nodes + node_index])
            + (n_a.vd[( 7 + component_index*LBQ ) * para.number_of_nodes + node_index] - n_a.vd[( 8 + component_index*LBQ ) * para.number_of_nodes + node_index])
            + (n_a.vd[( 9 + component_index*LBQ ) * para.number_of_nodes + node_index] - n_a.vd[(10 + component_index*LBQ ) * para.number_of_nodes + node_index])
            + (n_a.vd[(11 + component_index*LBQ ) * para.number_of_nodes + node_index] - n_a.vd[(12 + component_index*LBQ ) * para.number_of_nodes + node_index])
            + (n_a.vd[(13 + component_index*LBQ ) * para.number_of_nodes + node_index] - n_a.vd[(14 + component_index*LBQ ) * para.number_of_nodes + node_index]);

  mode[2] =   (n_a.vd[( 3 + component_index*LBQ ) * para.number_of_nodes + node_index] - n_a.vd[( 4 + component_index*LBQ ) * para.number_of_nodes + node_index])
            + (n_a.vd[( 7 + component_index*LBQ ) * para.number_of_nodes + node_index] - n_a.vd[( 8 + component_index*LBQ ) * para.number_of_nodes + node_index])
            - (n_a.vd[( 9 + component_index*LBQ ) * para.number_of_nodes + node_index] - n_a.vd[(10 + component_index*LBQ ) * para.number_of_nodes + node_index])
            + (n_a.vd[(15 + component_index*LBQ ) * para.number_of_nodes + node_index] - n_a.vd[(16 + component_index*LBQ ) * para.number_of_nodes + node_index])
            + (n_a.vd[(17 + component_index*LBQ ) * para.number_of_nodes + node_index] - n_a.vd[(18 + component_index*LBQ ) * para.number_of_nodes + node_index]);

  mode[3] =   (n_a.vd[( 5 + component_index*LBQ ) * para.number_of_nodes + node_index] - n_a.vd[( 6 + component_index*LBQ ) * para.number_of_nodes + node_index])
            + (n_a.vd[(11 + component_index*LBQ ) * para.number_of_nodes + node_index] - n_a.vd[(12 + component_index*LBQ ) * para.number_of_nodes + node_index])
            - (n_a.vd[(13 + component_index*LBQ ) * para.number_of_nodes + node_index] - n_a.vd[(14 + component_index*LBQ ) * para.number_of_nodes + node_index])
            + (n_a.vd[(15 + component_index*LBQ ) * para.number_of_nodes + node_index] - n_a.vd[(16 + component_index*LBQ ) * para.number_of_nodes + node_index])
            - (n_a.vd[(17 + component_index*LBQ ) * para.number_of_nodes + node_index] - n_a.vd[(18 + component_index*LBQ ) * para.number_of_nodes + node_index]);
}

/*********************************************************/
/** \name interpolation_three_point_coupling */
/*********************************************************/
/**(Eq. (12) Ahlrichs and Duenweg, JCP 111(17):8225 (1999))
 * @param n_a             Pointer to local node residing in array a (Input)
 * @param *delta          Pointer for the weighting of particle position (Output)
 * @param *particle_data  Pointer to the particle position and velocity (Input)
 * @param *particle_force Pointer to the particle force (Input)
 * @param part_index      particle id / thread id (Input)
 * @param node_index      node index around (8) particle (Output)
 * @param *d_v            Pointer to local device values
 * @param *interpolated_u Pointer to the interpolated velocity
*/
__device__ __inline__ void interpolation_three_point_coupling( LB_nodes_gpu n_a, float* particle_position, unsigned int *node_index, LB_rho_v_gpu *d_v, float *delta, float *interpolated_u ) {

  int my_center[3];
  float temp_delta[27];
  float mode[19*LB_COMPONENTS];

  /** see Duenweg and Ladd http://arxiv.org/abs/0803.2826 eqn. 301 */
  /** the i index is left node, nearest node, right node */
  for(int i=0; i<3; ++i){
    /** note the -0.5f is to account for the shift of the LB grid relative to the MD */
    float scaledpos = particle_position[i]/para.agrid-0.5f;
    /** the +0.5 is to turn the floorf into a round function */
    my_center[i] = (int)(floorf(scaledpos+0.5f));
    scaledpos = scaledpos-1.0f*my_center[i];
    temp_delta[0+3*i] = ( 5.0f - 3.0f*abs(scaledpos+1.0f) - sqrtf( -2.0f + 6.0f*abs(scaledpos+1.0f) - 3.0f*pow(scaledpos+1.0f,2) ) )/6.0f;
    temp_delta[1+3*i] = ( 1.0f + sqrtf( 1.0f - 3.0f*powf(scaledpos,2) ) )/3.0f;
    temp_delta[2+3*i] = ( 5.0f - 3.0f*abs(scaledpos-1.0f) - sqrtf( -2.0f + 6.0f*abs(scaledpos-1.0f) - 3.0f*pow(scaledpos-1.0f,2) ) )/6.0f;
  }

  for (int i=-1; i<=1; i++) {
    for (int j=-1; j<=1; j++) {
      for (int k=-1; k<=1; k++) {
        delta[i+3*j+9*k+13] = temp_delta[i+1] * temp_delta[3+j+1] * temp_delta[6+k+1];
      }
    }
  }

  // modulo for negative numbers is strange at best, shift to make sure we are positive
  int x = my_center[0] + para.dim_x;
  int y = my_center[1] + para.dim_y;
  int z = my_center[2] + para.dim_z;
  /** Here we collect the nodes for the three point coupling scheme (27 nodes in 3d) with the analogous numbering scheme of the two point coupling scheme */
  for (int i=-1; i<=1; i++) {
    for (int j=-1; j<=1; j++) {
      for (int k=-1; k<=1; k++) {
        node_index[i+3*j+9*k+13] = (x+i+para.dim_x)%para.dim_x + para.dim_x*((y+j+para.dim_y)%para.dim_y) + para.dim_x*para.dim_y*((z+k+para.dim_z)%para.dim_z);
      }
    }
  }

  interpolated_u[0] = 0.0f;
  interpolated_u[1] = 0.0f;
  interpolated_u[2] = 0.0f;
#pragma unroll
  for(int i=0; i<27; ++i){
    float totmass=0.0f;
    calc_m_from_n(n_a,node_index[i],mode);
#pragma unroll
    for(int ii=0;ii<LB_COMPONENTS;ii++){
      totmass+=mode[0]+para.rho[ii]*para.agrid*para.agrid*para.agrid;
    } 
    /* The boolean expression (n_a.boundary[node_index[i]] == 0) causes boundary nodes
       to couple with velocity 0 to particles. This is necessary, since boundary nodes
       undergo the same LB dynamics as fluid nodes do. The flow within the boundaries
       does not interact with the physical fluid, since these populations are overwritten
       by the bounce back kernel. Particles close to walls can couple to this unphysical
       flow, though.
    */
    interpolated_u[0] += (mode[1]/totmass)*delta[i] * (n_a.boundary[node_index[i]] == 0);
    interpolated_u[1] += (mode[2]/totmass)*delta[i] * (n_a.boundary[node_index[i]] == 0);
    interpolated_u[2] += (mode[3]/totmass)*delta[i] * (n_a.boundary[node_index[i]] == 0);
  }
}


/*********************************************************/
/** \name calc_viscous_force_three_point_couple */
/*********************************************************/
/**(Eq. (12) Ahlrichs and Duenweg, JCP 111(17):8225 (1999))
 * @param n_a                Pointer to local node residing in array a (Input)
 * @param *delta             Pointer for the weighting of particle position (Output)
 * @param *delta_j           Pointer for the weighting of particle momentum (Output)
 * @param *particle_position Pointer to the particle position (Input)
 * @param *rn_part           Pointer to randomnumber array of the particle
 * @param node_index         node index around (8) particle (Output)
 * @param *d_v               Pointer to local device values
 * @param flag_cs            Determine if we are at the centre (0, typical) or at the source (1, swimmer only)
*/
__device__ void calc_viscous_force_three_point_couple(LB_nodes_gpu n_a, float *delta, CUDA_particle_data *particle_data, float *particle_force, unsigned int part_index, LB_randomnr_gpu *rn_part, float *delta_j, unsigned int *node_index, LB_rho_v_gpu *d_v, int flag_cs){

  float interpolated_u[3];
  float interpolated_rho[LB_COMPONENTS];
  float viscforce[3*LB_COMPONENTS];

  // Zero out workspace
#pragma unroll
  for(int ii=0; ii<LB_COMPONENTS; ++ii){ 
#pragma unroll
    for(int jj=0; jj<3; ++jj){ 
      viscforce[jj+ii*3]=0.0f;
      delta_j[jj+ii*3]  =0.0f;
    }
  }
  // Zero out only if we are at the centre of the particle <=> flag_cs = 0
  particle_force[3*part_index+0] = flag_cs * particle_force[3*part_index+0];
  particle_force[3*part_index+1] = flag_cs * particle_force[3*part_index+1];
  particle_force[3*part_index+2] = flag_cs * particle_force[3*part_index+2];

  float position[3];
  position[0] = particle_data[part_index].p[0];
  position[1] = particle_data[part_index].p[1];
  position[2] = particle_data[part_index].p[2];

  float velocity[3];
  velocity[0] = particle_data[part_index].v[0];
  velocity[1] = particle_data[part_index].v[1];
  velocity[2] = particle_data[part_index].v[2];

#ifdef ENGINE
  // First calculate interpolated velocity for dipole source,
  // such that we don't overwrite mode, d_v, etc. for the rest of the function
  float direction = float(particle_data[part_index].swim.push_pull) * particle_data[part_index].swim.dipole_length;
  // Extrapolate position by dipole length if we are at the centre of the particle
  position[0] += flag_cs * direction * particle_data[part_index].swim.quatu[0];
  position[1] += flag_cs * direction * particle_data[part_index].swim.quatu[1];
  position[2] += flag_cs * direction * particle_data[part_index].swim.quatu[2];
#endif

  // Do the velocity interpolation
  interpolation_three_point_coupling(n_a, position, node_index, d_v, delta, interpolated_u);

#ifdef ENGINE
  velocity[0] -= (particle_data[part_index].swim.v_swim*para.time_step)*particle_data[part_index].swim.quatu[0];
  velocity[1] -= (particle_data[part_index].swim.v_swim*para.time_step)*particle_data[part_index].swim.quatu[1];
  velocity[2] -= (particle_data[part_index].swim.v_swim*para.time_step)*particle_data[part_index].swim.quatu[2];

  // The first three components are v_center, the last three v_source
  // Do not use within LB, because these have already been converted back to MD units
  particle_data[part_index].swim.v_cs[0+3*flag_cs] = interpolated_u[0] * para.agrid / para.tau;
  particle_data[part_index].swim.v_cs[1+3*flag_cs] = interpolated_u[1] * para.agrid / para.tau;
  particle_data[part_index].swim.v_cs[2+3*flag_cs] = interpolated_u[2] * para.agrid / para.tau;
#endif

  /* for LB we do not reweight the friction force */
  for(int ii=0; ii<LB_COMPONENTS; ++ii){ 
    interpolated_rho[ii]=1.0;
  }

  /** calculate viscous force
   * take care to rescale velocities with time_step and transform to MD units
   * (Eq. (9) Ahlrichs and Duenweg, JCP 111(17):8225 (1999)) */
  float rhotot=0;

#pragma unroll
  for(int ii=0; ii<LB_COMPONENTS; ++ii){ 
    rhotot+=interpolated_rho[ii];
  }


  /* Viscous force */
  for(int ii=0; ii<LB_COMPONENTS; ++ii){ 
    viscforce[0+ii*3] -= interpolated_rho[ii]*para.friction[ii]*(velocity[0]/para.time_step - interpolated_u[0]*para.agrid/para.tau)/rhotot;
    viscforce[1+ii*3] -= interpolated_rho[ii]*para.friction[ii]*(velocity[1]/para.time_step - interpolated_u[1]*para.agrid/para.tau)/rhotot;
    viscforce[2+ii*3] -= interpolated_rho[ii]*para.friction[ii]*(velocity[2]/para.time_step - interpolated_u[2]*para.agrid/para.tau)/rhotot;

#ifdef LB_ELECTROHYDRODYNAMICS
    viscforce[0+ii*3] += interpolated_rho[ii]*para.friction[ii] * particle_data[part_index].mu_E[0]/rhotot;
    viscforce[1+ii*3] += interpolated_rho[ii]*para.friction[ii] * particle_data[part_index].mu_E[1]/rhotot;
    viscforce[2+ii*3] += interpolated_rho[ii]*para.friction[ii] * particle_data[part_index].mu_E[2]/rhotot;
#endif

    /** add stochastic force of zero mean (Ahlrichs, Duenweg equ. 15)*/
#ifdef FLATNOISE
    random_01(rn_part);
    viscforce[0+ii*3] += para.lb_coupl_pref[ii]*(rn_part->randomnr[0]-0.5f);
    viscforce[1+ii*3] += para.lb_coupl_pref[ii]*(rn_part->randomnr[1]-0.5f);
    random_01(rn_part);
    viscforce[2+ii*3] += para.lb_coupl_pref[ii]*(rn_part->randomnr[0]-0.5f);
#elif defined(GAUSSRANDOMCUT)
    gaussian_random_cut(rn_part);
    viscforce[0+ii*3] += para.lb_coupl_pref2[ii]*rn_part->randomnr[0];
    viscforce[1+ii*3] += para.lb_coupl_pref2[ii]*rn_part->randomnr[1];
    gaussian_random_cut(rn_part);
    viscforce[2+ii*3] += para.lb_coupl_pref2[ii]*rn_part->randomnr[0];
#elif defined(GAUSSRANDOM)
    gaussian_random(rn_part);
    viscforce[0+ii*3] += para.lb_coupl_pref2[ii]*rn_part->randomnr[0];
    viscforce[1+ii*3] += para.lb_coupl_pref2[ii]*rn_part->randomnr[1];
    gaussian_random(rn_part);
    viscforce[2+ii*3] += para.lb_coupl_pref2[ii]*rn_part->randomnr[0];
#else
#error No noise type defined for the GPU LB
#endif    
    /** delta_j for transform momentum transfer to lattice units which is done in calc_node_force
      (Eq. (12) Ahlrichs and Duenweg, JCP 111(17):8225 (1999)) */
    // only add to particle_force for particle centre <=> (1-flag_cs) = 1
    particle_force[3*part_index+0] += (1-flag_cs) * viscforce[0+ii*3];
    particle_force[3*part_index+1] += (1-flag_cs) * viscforce[1+ii*3];
    particle_force[3*part_index+2] += (1-flag_cs) * viscforce[2+ii*3];

    // only add to particle_force for particle centre <=> (1-flag_cs) = 1
    delta_j[0+3*ii] -= (1-flag_cs)*viscforce[0+ii*3]*para.time_step*para.tau/para.agrid;
    delta_j[1+3*ii] -= (1-flag_cs)*viscforce[1+ii*3]*para.time_step*para.tau/para.agrid;
    delta_j[2+3*ii] -= (1-flag_cs)*viscforce[2+ii*3]*para.time_step*para.tau/para.agrid;
#ifdef ENGINE
    // add swimming force to source position
    delta_j[0+3*ii] -= flag_cs*particle_data[part_index].swim.f_swim*particle_data[part_index].swim.quatu[0]*para.time_step*para.tau/para.agrid;
    delta_j[1+3*ii] -= flag_cs*particle_data[part_index].swim.f_swim*particle_data[part_index].swim.quatu[1]*para.time_step*para.tau/para.agrid;
    delta_j[2+3*ii] -= flag_cs*particle_data[part_index].swim.f_swim*particle_data[part_index].swim.quatu[2]*para.time_step*para.tau/para.agrid;
#endif
  }
}

/**calcutlation of the node force caused by the particles, with atomicadd due to avoiding race conditions 
  (Eq. (14) Ahlrichs and Duenweg, JCP 111(17):8225 (1999))
 * @param *delta    Pointer for the weighting of particle position (Input)
 * @param *delta_j    Pointer for the weighting of particle momentum (Input)
 * @param node_index    node index around (8) particle (Input)
 * @param node_f        Pointer to the node force (Output).
*/
__device__ void calc_node_force_three_point_couple(float *delta, float *delta_j, unsigned int *node_index, LB_node_force_gpu node_f){
/* TODO: should the drag depend on the density?? */

  for (int i=-1; i<=1; i++) {
    for (int j=-1; j<=1; j++) {
      for (int k=-1; k<=1; k++) {
        atomicadd(&(node_f.force[0*para.number_of_nodes + node_index[i+3*j+9*k+13]]), (delta[i+3*j+9*k+13]*delta_j[0]));
        atomicadd(&(node_f.force[1*para.number_of_nodes + node_index[i+3*j+9*k+13]]), (delta[i+3*j+9*k+13]*delta_j[1]));
        atomicadd(&(node_f.force[2*para.number_of_nodes + node_index[i+3*j+9*k+13]]), (delta[i+3*j+9*k+13]*delta_j[2]));
      }
    }
  }
}


/**calculate temperature of the fluid kernel
 * @param *cpu_jsquared   Pointer to result storage value (Output)
 * @param n_a             Pointer to local node residing in array a (Input)
*/
__global__ void temperature(LB_nodes_gpu n_a, float *cpu_jsquared, int *number_of_non_boundary_nodes ) {
  float mode[4];
  float jsquared = 0.0f;
  unsigned int index = blockIdx.y * gridDim.x * blockDim.x + blockDim.x * blockIdx.x + threadIdx.x;

  if(index<para.number_of_nodes)
  {
    if(!n_a.boundary[index])
    {
      for(int ii=0;ii<LB_COMPONENTS;++ii)
      {  
         calc_mode(mode, n_a, index,ii);
         jsquared = mode[1]*mode[1]+mode[2]*mode[2]+mode[3]*mode[3];
         atomicadd(cpu_jsquared, jsquared);
         atomicAdd(number_of_non_boundary_nodes, 1);
      }
    }
  }
}

/*********************************************************/
/** \name interpolation_two_point_coupling */
/*********************************************************/
/**(Eq. (12) Ahlrichs and Duenweg, JCP 111(17):8225 (1999))
 * @param n_a                   Pointer to local node residing in array a (Input)
 * @param node_index            node index around (8) particle (Output)
 * @param *particle_position    Pointer to the particle position (Input)
 * @param *mode                 Pointer to the 19 modes for current lattice point
 * @param *d_v                  Pointer to local device values
 * @param *interpolated_u       Pointer to the interpolated velocity
 * @param *delta                Pointer for the weighting of particle position (Output)
*/
__device__ __inline__ void interpolation_two_point_coupling( LB_nodes_gpu n_a, float *particle_position, unsigned int* node_index, float* mode, LB_rho_v_gpu *d_v, float* delta, float *interpolated_u ) {
  int   left_node_index[3];
  float temp_delta[6];
  float temp_delta_half[6];

  // see ahlrichs + duenweg page 8227 equ (10) and (11)
#pragma unroll
  for(int i=0; i<3; ++i)
  {
    float scaledpos = particle_position[i]/para.agrid - 0.5f;
    left_node_index[i] = (int)(floorf(scaledpos));
    temp_delta[3+i] = scaledpos - left_node_index[i];
    temp_delta[i] = 1.0f - temp_delta[3+i];
    // further value used for interpolation of fluid velocity at part pos near boundaries
    temp_delta_half[3+i] = (scaledpos - left_node_index[i])*2.0f;
    temp_delta_half[i] = 2.0f - temp_delta_half[3+i];
  }

  delta[0] = temp_delta[0] * temp_delta[1] * temp_delta[2];
  delta[1] = temp_delta[3] * temp_delta[1] * temp_delta[2];
  delta[2] = temp_delta[0] * temp_delta[4] * temp_delta[2];
  delta[3] = temp_delta[3] * temp_delta[4] * temp_delta[2];
  delta[4] = temp_delta[0] * temp_delta[1] * temp_delta[5];
  delta[5] = temp_delta[3] * temp_delta[1] * temp_delta[5];
  delta[6] = temp_delta[0] * temp_delta[4] * temp_delta[5];
  delta[7] = temp_delta[3] * temp_delta[4] * temp_delta[5];

  // modulo for negative numbers is strange at best, shift to make sure we are positive
  int x = left_node_index[0] + para.dim_x;
  int y = left_node_index[1] + para.dim_y;
  int z = left_node_index[2] + para.dim_z;

  node_index[0] = x%para.dim_x     + para.dim_x*(y%para.dim_y)     + para.dim_x*para.dim_y*(z%para.dim_z);
  node_index[1] = (x+1)%para.dim_x + para.dim_x*(y%para.dim_y)     + para.dim_x*para.dim_y*(z%para.dim_z);
  node_index[2] = x%para.dim_x     + para.dim_x*((y+1)%para.dim_y) + para.dim_x*para.dim_y*(z%para.dim_z);
  node_index[3] = (x+1)%para.dim_x + para.dim_x*((y+1)%para.dim_y) + para.dim_x*para.dim_y*(z%para.dim_z);
  node_index[4] = x%para.dim_x     + para.dim_x*(y%para.dim_y)     + para.dim_x*para.dim_y*((z+1)%para.dim_z);
  node_index[5] = (x+1)%para.dim_x + para.dim_x*(y%para.dim_y)     + para.dim_x*para.dim_y*((z+1)%para.dim_z);
  node_index[6] = x%para.dim_x     + para.dim_x*((y+1)%para.dim_y) + para.dim_x*para.dim_y*((z+1)%para.dim_z);
  node_index[7] = (x+1)%para.dim_x + para.dim_x*((y+1)%para.dim_y) + para.dim_x*para.dim_y*((z+1)%para.dim_z);


  interpolated_u[0] = 0.0f;
  interpolated_u[1] = 0.0f;
  interpolated_u[2] = 0.0f;
#pragma unroll
  for(int i=0; i<8; ++i)
  {
    float totmass=0.0f;

    calc_m_from_n(n_a,node_index[i],mode);

#pragma unroll
    for(int ii=0;ii<LB_COMPONENTS;ii++)
    {
      totmass+=mode[0]+para.rho[ii]*para.agrid*para.agrid*para.agrid;
    } 

#ifdef SHANCHEN
    interpolated_u[0] += d_v[node_index[i]].v[0]/8.0f * (n_a.boundary[node_index[i]] == 0);  
    interpolated_u[1] += d_v[node_index[i]].v[1]/8.0f * (n_a.boundary[node_index[i]] == 0);
    interpolated_u[2] += d_v[node_index[i]].v[2]/8.0f * (n_a.boundary[node_index[i]] == 0);
#else
    /* The boolean expression (n_a.boundary[node_index[i]] == 0) causes boundary nodes
       to couple with velocity 0 to particles. This is necessary, since boundary nodes
       undergo the same LB dynamics as fluid nodes do. The flow within the boundaries
       does not interact with the physical fluid, since these populations are overwritten
       by the bounce back kernel. Particles close to walls can couple to this unphysical
       flow, though.
    */
    interpolated_u[0] += (mode[1]/totmass)*delta[i] * (n_a.boundary[node_index[i]] == 0);
    interpolated_u[1] += (mode[2]/totmass)*delta[i] * (n_a.boundary[node_index[i]] == 0);
    interpolated_u[2] += (mode[3]/totmass)*delta[i] * (n_a.boundary[node_index[i]] == 0);
#endif
  }
}

/*********************************************************/
/** \name calc_viscous_force */
/*********************************************************/
/**(Eq. (12) Ahlrichs and Duenweg, JCP 111(17):8225 (1999))
 * @param n_a                   Pointer to local node residing in array a (Input)
 * @param partgrad1             particle gradient for the Shan-Chen
 * @param partgrad2             particle gradient for the Shan-Chen
 * @param partgrad3             particle gradient for the Shan-Chen
 * @param *delta                Pointer for the weighting of particle position (Output)
 * @param *delta_j              Pointer for the weighting of particle momentum (Output)
 * @param *particle_data        Pointer to the particle position and velocity (Input)
 * @param *particle_force       Pointer to the particle force (Input)
 * @param *fluid_composition    Pointer to the fluid composition (Input)
 * @param part_index            particle id / thread id (Input)
 * @param *rn_part              Pointer to randomnumber array of the particle
 * @param node_index            node index around (8) particle (Output)
 * @param *d_v                  Pointer to local device values
 * @param flag_cs               Determine if we are at the centre (0, typical) or at the source (1, swimmer only)
*/
__device__ void calc_viscous_force(LB_nodes_gpu n_a, float *delta, float * partgrad1, float * partgrad2, float * partgrad3, CUDA_particle_data *particle_data, float *particle_force, CUDA_fluid_composition * fluid_composition, unsigned int part_index, LB_randomnr_gpu *rn_part, float *delta_j, unsigned int *node_index, LB_rho_v_gpu *d_v, int flag_cs){

  float interpolated_u[3];
  float interpolated_rho[LB_COMPONENTS];
  float viscforce[3*LB_COMPONENTS];
  float scforce[3*LB_COMPONENTS];
  float mode[19*LB_COMPONENTS];
#ifdef SHANCHEN
  float gradrho1, gradrho2, gradrho3;
  float Rho;
#endif
  // Zero out workspace
  #pragma unroll
  for(int ii=0; ii<LB_COMPONENTS; ++ii)
  { 
    #pragma unroll
    for(int jj=0; jj<3; ++jj)
    { 
      scforce[jj+ii*3]  =0.0f;
      viscforce[jj+ii*3]=0.0f;
      delta_j[jj+ii*3]  =0.0f;
    }
    
    #pragma unroll
    for(int jj=0; jj<8; ++jj)
    { 
      partgrad1[jj+ii*8]=0.0f;
      partgrad2[jj+ii*8]=0.0f;
      partgrad3[jj+ii*8]=0.0f;
    }
  }
  // Zero out only if we are at the centre of the particle <=> flag_cs = 0
  particle_force[3*part_index+0] = flag_cs * particle_force[3*part_index+0];
  particle_force[3*part_index+1] = flag_cs * particle_force[3*part_index+1];
  particle_force[3*part_index+2] = flag_cs * particle_force[3*part_index+2];

  float position[3];
  position[0] = particle_data[part_index].p[0];
  position[1] = particle_data[part_index].p[1];
  position[2] = particle_data[part_index].p[2];

  float velocity[3];
  velocity[0] = particle_data[part_index].v[0];
  velocity[1] = particle_data[part_index].v[1];
  velocity[2] = particle_data[part_index].v[2];

#ifdef ENGINE
  // First calculate interpolated velocity for dipole source,
  // such that we don't overwrite mode, d_v, etc. for the rest of the function
  float direction = float(particle_data[part_index].swim.push_pull) * particle_data[part_index].swim.dipole_length;
  // Extrapolate position by dipole length if we are at the centre of the particle
  position[0] += flag_cs * direction * particle_data[part_index].swim.quatu[0];
  position[1] += flag_cs * direction * particle_data[part_index].swim.quatu[1];
  position[2] += flag_cs * direction * particle_data[part_index].swim.quatu[2];
#endif

  // Do the velocity interpolation
  interpolation_two_point_coupling(n_a, position, node_index, mode, d_v, delta, interpolated_u);

#ifdef ENGINE
  velocity[0] -= (particle_data[part_index].swim.v_swim*para.time_step)*particle_data[part_index].swim.quatu[0];
  velocity[1] -= (particle_data[part_index].swim.v_swim*para.time_step)*particle_data[part_index].swim.quatu[1];
  velocity[2] -= (particle_data[part_index].swim.v_swim*para.time_step)*particle_data[part_index].swim.quatu[2];

  // The first three components are v_center, the last three v_source
  // Do not use within LB, because these have already been converted back to MD units
  particle_data[part_index].swim.v_cs[0+3*flag_cs] = interpolated_u[0] * para.agrid / para.tau;
  particle_data[part_index].swim.v_cs[1+3*flag_cs] = interpolated_u[1] * para.agrid / para.tau;
  particle_data[part_index].swim.v_cs[2+3*flag_cs] = interpolated_u[2] * para.agrid / para.tau;
#endif

#ifdef SHANCHEN

 #pragma unroll
  for(int ii=0; ii<LB_COMPONENTS; ++ii)
  { 
    float solvation2 = particle_data[part_index].solvation[2*ii + 1];
   
    interpolated_rho[ii]  = 0.0f;
    gradrho1 = gradrho2 = gradrho3 = 0.0f;
  
    // TODO: should one introduce a density-dependent friction ?
    calc_mode(mode, n_a, node_index[0],ii);
    Rho = mode[0] + para.rho[ii]*para.agrid*para.agrid*para.agrid;
    interpolated_rho[ii] += delta[0] * Rho; 
    partgrad1[ii*8 + 0] += Rho * solvation2;
    partgrad2[ii*8 + 0] += Rho * solvation2;
    partgrad3[ii*8 + 0] += Rho * solvation2;
    gradrho1 -=(delta[0] + delta[1]) * Rho; 
    gradrho2 -=(delta[0] + delta[2]) * Rho; 
    gradrho3 -=(delta[0] + delta[4]) * Rho; 

    calc_mode(mode, n_a, node_index[1],ii); 
    Rho = mode[0] +  para.rho[ii]*para.agrid*para.agrid*para.agrid; 
    interpolated_rho[ii] += delta[1] * Rho; 
    partgrad1[ii*8 + 1] -= Rho * solvation2;
    partgrad2[ii*8 + 1] += Rho * solvation2;
    partgrad3[ii*8 + 1] += Rho * solvation2;
    gradrho1 +=(delta[1] + delta[0]) * Rho; 
    gradrho2 -=(delta[1] + delta[3]) * Rho; 
    gradrho3 -=(delta[1] + delta[5]) * Rho; 

    calc_mode(mode, n_a, node_index[2],ii);
    Rho = mode[0] + para.rho[ii]*para.agrid*para.agrid*para.agrid;
    interpolated_rho[ii] += delta[2] * Rho; 
    partgrad1[ii*8 + 2] += Rho * solvation2;
    partgrad2[ii*8 + 2] -= Rho * solvation2;
    partgrad3[ii*8 + 2] += Rho * solvation2;
    gradrho1 -=(delta[2] + delta[3]) * Rho; 
    gradrho2 +=(delta[2] + delta[0]) * Rho; 
    gradrho3 -=(delta[2] + delta[6]) * Rho; 

    calc_mode(mode, n_a, node_index[3],ii);
    Rho = mode[0] + para.rho[ii]*para.agrid*para.agrid*para.agrid;
    interpolated_rho[ii] += delta[3] * Rho; 
    partgrad1[ii*8 + 3] -= Rho * solvation2;
    partgrad2[ii*8 + 3] -= Rho * solvation2;
    partgrad3[ii*8 + 3] += Rho * solvation2;
    gradrho1 +=(delta[3] + delta[2]) * Rho; 
    gradrho2 +=(delta[3] + delta[1]) * Rho; 
    gradrho3 -=(delta[3] + delta[7]) * Rho; 

    calc_mode(mode, n_a, node_index[4],ii);
    Rho = mode[0] + para.rho[ii]*para.agrid*para.agrid*para.agrid;
    interpolated_rho[ii] += delta[4] * Rho; 
    partgrad1[ii*8 + 4] += Rho * solvation2;
    partgrad2[ii*8 + 4] += Rho * solvation2;
    partgrad3[ii*8 + 4] -= Rho * solvation2;
    gradrho1 -=(delta[4] + delta[5]) * Rho; 
    gradrho2 -=(delta[4] + delta[6]) * Rho; 
    gradrho3 +=(delta[4] + delta[0]) * Rho; 

    calc_mode(mode, n_a, node_index[5],ii);
    Rho = mode[0] + para.rho[ii]*para.agrid*para.agrid*para.agrid;
    interpolated_rho[ii] += delta[5] * Rho; 
    partgrad1[ii*8 + 5] -= Rho * solvation2;
    partgrad2[ii*8 + 5] += Rho * solvation2;
    partgrad3[ii*8 + 5] -= Rho * solvation2;
    gradrho1 +=(delta[5] + delta[4]) * Rho; 
    gradrho2 -=(delta[5] + delta[7]) * Rho; 
    gradrho3 +=(delta[5] + delta[1]) * Rho; 

    calc_mode(mode, n_a, node_index[6],ii);
    Rho = mode[0] + para.rho[ii]*para.agrid*para.agrid*para.agrid;
    interpolated_rho[ii] += delta[6] * Rho; 
    partgrad1[ii*8 + 6] += Rho * solvation2;
    partgrad2[ii*8 + 6] -= Rho * solvation2;
    partgrad3[ii*8 + 6] -= Rho * solvation2;
    gradrho1 -=(delta[6] + delta[7]) * Rho; 
    gradrho2 +=(delta[6] + delta[4]) * Rho; 
    gradrho3 +=(delta[6] + delta[2]) * Rho; 

    calc_mode(mode, n_a, node_index[7],ii);
    Rho = mode[0] + para.rho[ii]*para.agrid*para.agrid*para.agrid;
    interpolated_rho[ii] += delta[7] * Rho; 
    partgrad1[ii*8 + 7] -= Rho * solvation2;
    partgrad2[ii*8 + 7] -= Rho * solvation2;
    partgrad3[ii*8 + 7] -= Rho * solvation2;
    gradrho1 +=(delta[7] + delta[6]) * Rho; 
    gradrho2 +=(delta[7] + delta[5]) * Rho; 
    gradrho3 +=(delta[7] + delta[3]) * Rho; 

    /* normalize the gradient to md units TODO: is that correct?*/
    gradrho1 *= para.agrid; 
    gradrho2 *= para.agrid; 
    gradrho3 *= para.agrid; 

    // scforce is 0 at the interpolated point where the swimming force gets put back on the fluid
    // only add sc_force up for particle centre <=> (1-flag_cs) = 1
    scforce[0+ii*3] += (1-flag_cs) * particle_data[part_index].solvation[2*ii] * gradrho1 ; 
    scforce[1+ii*3] += (1-flag_cs) * particle_data[part_index].solvation[2*ii] * gradrho2 ;
    scforce[2+ii*3] += (1-flag_cs) * particle_data[part_index].solvation[2*ii] * gradrho3 ;

    /* scforce is used also later...*/
    particle_force[3*part_index+0] += scforce[0+ii*3];
    particle_force[3*part_index+1] += scforce[1+ii*3];
    particle_force[3*part_index+2] += scforce[2+ii*3];
    // only set fluid_composition for particle centre <=> (1-flag_cs) = 1
    fluid_composition[part_index].weight[ii] = (1-flag_cs) * interpolated_rho[ii];
 }

#else // SHANCHEN is not defined

  /* for LB we do not reweight the friction force */
  for(int ii=0; ii<LB_COMPONENTS; ++ii){
    interpolated_rho[ii]=1.0;
  }

#endif // SHANCHEN

  /** calculate viscous force
   * take care to rescale velocities with time_step and transform to MD units
   * (Eq. (9) Ahlrichs and Duenweg, JCP 111(17):8225 (1999)) */
  float rhotot=0;

  #pragma unroll
  for(int ii=0; ii<LB_COMPONENTS; ++ii){
    rhotot+=interpolated_rho[ii];
  }

  /* Viscous force */
  for(int ii=0; ii<LB_COMPONENTS; ++ii)
  { 
    viscforce[0+ii*3] -= interpolated_rho[ii]*para.friction[ii]*(velocity[0]/para.time_step - interpolated_u[0]*para.agrid/para.tau)/rhotot;
    viscforce[1+ii*3] -= interpolated_rho[ii]*para.friction[ii]*(velocity[1]/para.time_step - interpolated_u[1]*para.agrid/para.tau)/rhotot;
    viscforce[2+ii*3] -= interpolated_rho[ii]*para.friction[ii]*(velocity[2]/para.time_step - interpolated_u[2]*para.agrid/para.tau)/rhotot;

#ifdef LB_ELECTROHYDRODYNAMICS
    viscforce[0+ii*3] += interpolated_rho[ii]*para.friction[ii] * particle_data[part_index].mu_E[0]/rhotot;
    viscforce[1+ii*3] += interpolated_rho[ii]*para.friction[ii] * particle_data[part_index].mu_E[1]/rhotot;
    viscforce[2+ii*3] += interpolated_rho[ii]*para.friction[ii] * particle_data[part_index].mu_E[2]/rhotot;
#endif

    /** add stochastic force of zero mean (Ahlrichs, Duenweg equ. 15)*/
#ifdef FLATNOISE
    random_01(rn_part);
    viscforce[0+ii*3] += para.lb_coupl_pref[ii]*(rn_part->randomnr[0]-0.5f);
    viscforce[1+ii*3] += para.lb_coupl_pref[ii]*(rn_part->randomnr[1]-0.5f);
    random_01(rn_part);
    viscforce[2+ii*3] += para.lb_coupl_pref[ii]*(rn_part->randomnr[0]-0.5f);
#elif defined(GAUSSRANDOMCUT)
    gaussian_random_cut(rn_part);
    viscforce[0+ii*3] += para.lb_coupl_pref2[ii]*rn_part->randomnr[0];
    viscforce[1+ii*3] += para.lb_coupl_pref2[ii]*rn_part->randomnr[1];
    gaussian_random_cut(rn_part);
    viscforce[2+ii*3] += para.lb_coupl_pref2[ii]*rn_part->randomnr[0];
#elif defined(GAUSSRANDOM)
    gaussian_random(rn_part);
    viscforce[0+ii*3] += para.lb_coupl_pref2[ii]*rn_part->randomnr[0];
    viscforce[1+ii*3] += para.lb_coupl_pref2[ii]*rn_part->randomnr[1];
    gaussian_random(rn_part);
    viscforce[2+ii*3] += para.lb_coupl_pref2[ii]*rn_part->randomnr[0];
#else
#error No noise type defined for the GPU LB
#endif 

    /** delta_j for transform momentum transfer to lattice units which is done in calc_node_force
      (Eq. (12) Ahlrichs and Duenweg, JCP 111(17):8225 (1999)) */

    // only add to particle_force for particle centre <=> (1-flag_cs) = 1
    particle_force[3*part_index+0] += (1-flag_cs) * viscforce[0+ii*3];
    particle_force[3*part_index+1] += (1-flag_cs) * viscforce[1+ii*3];
    particle_force[3*part_index+2] += (1-flag_cs) * viscforce[2+ii*3];

    /* the average force from the particle to surrounding nodes is transmitted back to preserve momentum */
    for(int node=0 ; node < 8 ; node++ ) { 
      particle_force[3*part_index+0] -= (1-flag_cs) * partgrad1[node+ii*8]/8.0f;
      particle_force[3*part_index+1] -= (1-flag_cs) * partgrad2[node+ii*8]/8.0f;
      particle_force[3*part_index+2] -= (1-flag_cs) * partgrad3[node+ii*8]/8.0f;
    }

    /* note that scforce is zero if SHANCHEN is not #defined */
    // only add to particle_force for particle centre <=> (1-flag_cs) = 1
    delta_j[0+3*ii] -= (scforce[0+ii*3]+(1-flag_cs)*viscforce[0+ii*3])*para.time_step*para.tau/para.agrid;
    delta_j[1+3*ii] -= (scforce[1+ii*3]+(1-flag_cs)*viscforce[1+ii*3])*para.time_step*para.tau/para.agrid;
    delta_j[2+3*ii] -= (scforce[2+ii*3]+(1-flag_cs)*viscforce[2+ii*3])*para.time_step*para.tau/para.agrid;

#ifdef ENGINE
    // add swimming force to source position
    delta_j[0+3*ii] -= flag_cs*particle_data[part_index].swim.f_swim*particle_data[part_index].swim.quatu[0]*para.time_step*para.tau/para.agrid;
    delta_j[1+3*ii] -= flag_cs*particle_data[part_index].swim.f_swim*particle_data[part_index].swim.quatu[1]*para.time_step*para.tau/para.agrid;
    delta_j[2+3*ii] -= flag_cs*particle_data[part_index].swim.f_swim*particle_data[part_index].swim.quatu[2]*para.time_step*para.tau/para.agrid;
#endif

  }

#ifdef SHANCHEN
  for(int node=0 ; node < 8 ; node++ ) { 
    for(int ii=0 ; ii < LB_COMPONENTS ; ii++ ) { 
      partgrad1[node+ii*8]*=(para.time_step*para.tau/para.agrid);
      partgrad2[node+ii*8]*=(para.time_step*para.tau/para.agrid);
      partgrad3[node+ii*8]*=(para.time_step*para.tau/para.agrid);
    }
  }
#endif

}

/**calculation of the node force caused by the particles, with atomicadd due to avoiding race conditions 
  (Eq. (14) Ahlrichs and Duenweg, JCP 111(17):8225 (1999))
 * @param *delta        Pointer for the weighting of particle position (Input)
 * @param partgrad1             particle gradient for the Shan-Chen
 * @param partgrad2             particle gradient for the Shan-Chen
 * @param partgrad3             particle gradient for the Shan-Chen
 * @param *delta_j      Pointer for the weighting of particle momentum (Input)
 * @param node_index    node index around (8) particle (Input)
 * @param node_f        Pointer to the node force (Output).
*/
__device__ void calc_node_force(float *delta, float *delta_j, float * partgrad1, float * partgrad2, float * partgrad3,  unsigned int *node_index, LB_node_force_gpu node_f){
/* TODO: should the drag depend on the density?? */
/* NOTE: partgrad is not zero only if SHANCHEN is defined. It is initialized in calc_node_force. Alternatively one could 
         specialize this function to the single component LB */ 
  for(int ii=0; ii < LB_COMPONENTS; ++ii)
  { 
    atomicadd(&(node_f.force[(0+ii*3)*para.number_of_nodes + node_index[0]]), (delta[0]*delta_j[0+ii*3] + partgrad1[ii*8+0]));
    atomicadd(&(node_f.force[(1+ii*3)*para.number_of_nodes + node_index[0]]), (delta[0]*delta_j[1+ii*3] + partgrad2[ii*8+0]));
    atomicadd(&(node_f.force[(2+ii*3)*para.number_of_nodes + node_index[0]]), (delta[0]*delta_j[2+ii*3] + partgrad3[ii*8+0]));
                                                                                                      
    atomicadd(&(node_f.force[(0+ii*3)*para.number_of_nodes + node_index[1]]), (delta[1]*delta_j[0+ii*3] + partgrad1[ii*8+1]));
    atomicadd(&(node_f.force[(1+ii*3)*para.number_of_nodes + node_index[1]]), (delta[1]*delta_j[1+ii*3] + partgrad2[ii*8+1]));
    atomicadd(&(node_f.force[(2+ii*3)*para.number_of_nodes + node_index[1]]), (delta[1]*delta_j[2+ii*3] + partgrad3[ii*8+1]));
                                                                                                      
    atomicadd(&(node_f.force[(0+ii*3)*para.number_of_nodes + node_index[2]]), (delta[2]*delta_j[0+ii*3] + partgrad1[ii*8+2]));
    atomicadd(&(node_f.force[(1+ii*3)*para.number_of_nodes + node_index[2]]), (delta[2]*delta_j[1+ii*3] + partgrad2[ii*8+2]));
    atomicadd(&(node_f.force[(2+ii*3)*para.number_of_nodes + node_index[2]]), (delta[2]*delta_j[2+ii*3] + partgrad3[ii*8+2]));
                                                                                                      
    atomicadd(&(node_f.force[(0+ii*3)*para.number_of_nodes + node_index[3]]), (delta[3]*delta_j[0+ii*3] + partgrad1[ii*8+3]));
    atomicadd(&(node_f.force[(1+ii*3)*para.number_of_nodes + node_index[3]]), (delta[3]*delta_j[1+ii*3] + partgrad2[ii*8+3]));
    atomicadd(&(node_f.force[(2+ii*3)*para.number_of_nodes + node_index[3]]), (delta[3]*delta_j[2+ii*3] + partgrad3[ii*8+3]));
                                                                                                      
    atomicadd(&(node_f.force[(0+ii*3)*para.number_of_nodes + node_index[4]]), (delta[4]*delta_j[0+ii*3] + partgrad1[ii*8+4]));
    atomicadd(&(node_f.force[(1+ii*3)*para.number_of_nodes + node_index[4]]), (delta[4]*delta_j[1+ii*3] + partgrad2[ii*8+4]));
    atomicadd(&(node_f.force[(2+ii*3)*para.number_of_nodes + node_index[4]]), (delta[4]*delta_j[2+ii*3] + partgrad3[ii*8+4]));
                                                                                                      
    atomicadd(&(node_f.force[(0+ii*3)*para.number_of_nodes + node_index[5]]), (delta[5]*delta_j[0+ii*3] + partgrad1[ii*8+5]));
    atomicadd(&(node_f.force[(1+ii*3)*para.number_of_nodes + node_index[5]]), (delta[5]*delta_j[1+ii*3] + partgrad2[ii*8+5]));
    atomicadd(&(node_f.force[(2+ii*3)*para.number_of_nodes + node_index[5]]), (delta[5]*delta_j[2+ii*3] + partgrad3[ii*8+5]));
                                                                                                      
    atomicadd(&(node_f.force[(0+ii*3)*para.number_of_nodes + node_index[6]]), (delta[6]*delta_j[0+ii*3] + partgrad1[ii*8+6]));
    atomicadd(&(node_f.force[(1+ii*3)*para.number_of_nodes + node_index[6]]), (delta[6]*delta_j[1+ii*3] + partgrad2[ii*8+6]));
    atomicadd(&(node_f.force[(2+ii*3)*para.number_of_nodes + node_index[6]]), (delta[6]*delta_j[2+ii*3] + partgrad3[ii*8+6]));
                                                                                                      
    atomicadd(&(node_f.force[(0+ii*3)*para.number_of_nodes + node_index[7]]), (delta[7]*delta_j[0+ii*3] + partgrad1[ii*8+7]));
    atomicadd(&(node_f.force[(1+ii*3)*para.number_of_nodes + node_index[7]]), (delta[7]*delta_j[1+ii*3] + partgrad2[ii*8+7]));
    atomicadd(&(node_f.force[(2+ii*3)*para.number_of_nodes + node_index[7]]), (delta[7]*delta_j[2+ii*3] + partgrad3[ii*8+7]));
  }
}

/*********************************************************/
/** \name System setup and Kernel functions */
/*********************************************************/

/**kernel to calculate local populations from hydrodynamic fields given by the tcl values.
 * The mapping is given in terms of the equilibrium distribution.
 *
 * Eq. (2.15) Ladd, J. Fluid Mech. 271, 295-309 (1994)
 * Eq. (4) in Berk Usta, Ladd and Butler, JCP 122, 094902 (2005)
 *
 * @param n_a   Pointer to the lattice site (Input).
 * @param *gpu_check additional check if gpu kernel are executed(Input).
 * @param *d_v    Pointer to local device values
 * @param *node_f          Pointer to node forces
*/
__global__ void calc_n_from_rho_j_pi(LB_nodes_gpu n_a, LB_rho_v_gpu *d_v, LB_node_force_gpu node_f, int *gpu_check) {
   /* TODO: this can handle only a uniform density, something similar, but local, 
            has to be called every time the fields are set by the user ! */ 
  unsigned int index = blockIdx.y * gridDim.x * blockDim.x + blockDim.x * blockIdx.x + threadIdx.x;
  if(index<para.number_of_nodes)
  {
    float mode[19*LB_COMPONENTS];

    #pragma unroll
    for(int ii=0;ii<LB_COMPONENTS;++ii)
    { 
      /** default values for fields in lattice units */
      gpu_check[0] = 1;
     
      float Rho = para.rho[ii]*para.agrid*para.agrid*para.agrid;
      float v[3] = { 0.0f, 0.0f, 0.0f };
      float pi[6] = { Rho*c_sound_sq, 0.0f, Rho*c_sound_sq, 0.0f, 0.0f, Rho*c_sound_sq };
     
      float rhoc_sq = Rho*c_sound_sq;
      float avg_rho = para.rho[ii]*para.agrid*para.agrid*para.agrid;
      float local_rho, local_j[3], *local_pi, trace;
     
      local_rho  = Rho;
     
      local_j[0] = Rho * v[0];
      local_j[1] = Rho * v[1];
      local_j[2] = Rho * v[2];
     
      local_pi = pi;
     
      /** reduce the pressure tensor to the part needed here. 
          NOTE: this not true anymore for SHANCHEN 
          if the densities are not uniform. FIXME*/

      local_pi[0] -= rhoc_sq;
      local_pi[2] -= rhoc_sq;
      local_pi[5] -= rhoc_sq;
     
      trace = local_pi[0] + local_pi[2] + local_pi[5];
     
      float rho_times_coeff;
      float tmp1,tmp2;
     
      /** update the q=0 sublattice */
      n_a.vd[(0 + ii*LBQ ) * para.number_of_nodes + index] = 1.0f/3.0f * (local_rho-avg_rho) - 1.0f/2.0f*trace;
     
      /** update the q=1 sublattice */
      rho_times_coeff = 1.0f/18.0f * (local_rho-avg_rho);
     
      n_a.vd[(1 + ii*LBQ ) * para.number_of_nodes + index] = rho_times_coeff + 1.0f/6.0f*local_j[0] + 1.0f/4.0f*local_pi[0] - 1.0f/12.0f*trace;
      n_a.vd[(2 + ii*LBQ ) * para.number_of_nodes + index] = rho_times_coeff - 1.0f/6.0f*local_j[0] + 1.0f/4.0f*local_pi[0] - 1.0f/12.0f*trace;
      n_a.vd[(3 + ii*LBQ ) * para.number_of_nodes + index] = rho_times_coeff + 1.0f/6.0f*local_j[1] + 1.0f/4.0f*local_pi[2] - 1.0f/12.0f*trace;
      n_a.vd[(4 + ii*LBQ ) * para.number_of_nodes + index] = rho_times_coeff - 1.0f/6.0f*local_j[1] + 1.0f/4.0f*local_pi[2] - 1.0f/12.0f*trace;
      n_a.vd[(5 + ii*LBQ ) * para.number_of_nodes + index] = rho_times_coeff + 1.0f/6.0f*local_j[2] + 1.0f/4.0f*local_pi[5] - 1.0f/12.0f*trace;
      n_a.vd[(6 + ii*LBQ ) * para.number_of_nodes + index] = rho_times_coeff - 1.0f/6.0f*local_j[2] + 1.0f/4.0f*local_pi[5] - 1.0f/12.0f*trace;
     
      /** update the q=2 sublattice */
      rho_times_coeff = 1.0f/36.0f * (local_rho-avg_rho);
     
      tmp1 = local_pi[0] + local_pi[2];
      tmp2 = 2.0f*local_pi[1];
      n_a.vd[( 7 + ii*LBQ ) * para.number_of_nodes + index] = rho_times_coeff + 1.0f/12.0f*(local_j[0]+local_j[1]) + 1.0f/8.0f*(tmp1+tmp2) - 1.0f/24.0f*trace;
      n_a.vd[( 8 + ii*LBQ ) * para.number_of_nodes + index] = rho_times_coeff - 1.0f/12.0f*(local_j[0]+local_j[1]) + 1.0f/8.0f*(tmp1+tmp2) - 1.0f/24.0f*trace;
      n_a.vd[( 9 + ii*LBQ ) * para.number_of_nodes + index] = rho_times_coeff + 1.0f/12.0f*(local_j[0]-local_j[1]) + 1.0f/8.0f*(tmp1-tmp2) - 1.0f/24.0f*trace;
      n_a.vd[(10 + ii*LBQ ) * para.number_of_nodes + index] = rho_times_coeff - 1.0f/12.0f*(local_j[0]-local_j[1]) + 1.0f/8.0f*(tmp1-tmp2) - 1.0f/24.0f*trace;
     
      tmp1 = local_pi[0] + local_pi[5];
      tmp2 = 2.0f*local_pi[3];
     
      n_a.vd[(11 + ii*LBQ ) * para.number_of_nodes + index] = rho_times_coeff + 1.0f/12.0f*(local_j[0]+local_j[2]) + 1.0f/8.0f*(tmp1+tmp2) - 1.0f/24.0f*trace;
      n_a.vd[(12 + ii*LBQ ) * para.number_of_nodes + index] = rho_times_coeff - 1.0f/12.0f*(local_j[0]+local_j[2]) + 1.0f/8.0f*(tmp1+tmp2) - 1.0f/24.0f*trace;
      n_a.vd[(13 + ii*LBQ ) * para.number_of_nodes + index] = rho_times_coeff + 1.0f/12.0f*(local_j[0]-local_j[2]) + 1.0f/8.0f*(tmp1-tmp2) - 1.0f/24.0f*trace;
      n_a.vd[(14 + ii*LBQ ) * para.number_of_nodes + index] = rho_times_coeff - 1.0f/12.0f*(local_j[0]-local_j[2]) + 1.0f/8.0f*(tmp1-tmp2) - 1.0f/24.0f*trace;
     
      tmp1 = local_pi[2] + local_pi[5];
      tmp2 = 2.0f*local_pi[4];
     
      n_a.vd[(15 + ii*LBQ ) * para.number_of_nodes + index] = rho_times_coeff + 1.0f/12.0f*(local_j[1]+local_j[2]) + 1.0f/8.0f*(tmp1+tmp2) - 1.0f/24.0f*trace;
      n_a.vd[(16 + ii*LBQ ) * para.number_of_nodes + index] = rho_times_coeff - 1.0f/12.0f*(local_j[1]+local_j[2]) + 1.0f/8.0f*(tmp1+tmp2) - 1.0f/24.0f*trace;
      n_a.vd[(17 + ii*LBQ ) * para.number_of_nodes + index] = rho_times_coeff + 1.0f/12.0f*(local_j[1]-local_j[2]) + 1.0f/8.0f*(tmp1-tmp2) - 1.0f/24.0f*trace;
      n_a.vd[(18 + ii*LBQ ) * para.number_of_nodes + index] = rho_times_coeff - 1.0f/12.0f*(local_j[1]-local_j[2]) + 1.0f/8.0f*(tmp1-tmp2) - 1.0f/24.0f*trace;
     
      /**set different seed for randomgen on every node */
      n_a.seed[index] = para.your_seed + index;
    }

    calc_m_from_n(n_a,index,mode);
    update_rho_v(mode,index,node_f,d_v);
  }
}

/** kernel to calculate local populations from hydrodynamic fields
 * from given flow field velocities. The mapping is given in terms of
 * the equilibrium distribution.
 *
 * Eq. (2.15) Ladd, J. Fluid Mech. 271, 295-309 (1994)
 * Eq. (4) in Berk Usta, Ladd and Butler, JCP 122, 094902 (2005)
 *
 * @param n_a              the current nodes array (double buffering!)
 * @param single_nodeindex the node to set the velocity for
 * @param velocity         the velocity to set
 * @param *d_v             Pointer to local device values (Input)
 * @param *node_f          Pointer to node forces (Input)
 */ 
__global__ void set_u_from_rho_v_pi( LB_nodes_gpu n_a, int single_nodeindex, float *velocity, LB_rho_v_gpu *d_v, LB_node_force_gpu node_f ) {

  unsigned int index = blockIdx.y * gridDim.x * blockDim.x + blockDim.x * blockIdx.x + threadIdx.x;

  if(index == 0)
  {
    float local_rho;
    float local_j[3];
    float local_pi[6];
    float trace, avg_rho;
    float rho_times_coeff;
    float tmp1, tmp2; 

    float mode_for_pi[19*LB_COMPONENTS];
    float rho_from_m[1*LB_COMPONENTS];
    float j_from_m[3*LB_COMPONENTS];
    float pi_from_m[6*LB_COMPONENTS];

    // Calculate the modes for this node

    calc_m_from_n(n_a, single_nodeindex, mode_for_pi);

    // Reset the d_v

    update_rho_v(mode_for_pi, single_nodeindex, node_f, d_v);

    // Calculate the density, velocity, and pressure tensor
    // in LB unit for this node

    calc_values_from_m_in_LB_units( mode_for_pi, &d_v[single_nodeindex], rho_from_m, j_from_m, pi_from_m);

    #pragma unroll
    for(int ii=0;ii<LB_COMPONENTS;++ii)
    { 
      // Take LB component density and calculate the equilibrium part

      local_rho = rho_from_m[ii];
      avg_rho = para.rho[ii]*para.agrid*para.agrid*para.agrid;

      // Take LB component velocity and make it a momentum

      local_j[0] = local_rho * velocity[0];
      local_j[1] = local_rho * velocity[1];
      local_j[2] = local_rho * velocity[2];

      // Take LB component pressure tensor and put in equilibrium

      local_pi[0] = pi_from_m[6*ii + 0];
      local_pi[1] = pi_from_m[6*ii + 1];
      local_pi[2] = pi_from_m[6*ii + 2];
      local_pi[3] = pi_from_m[6*ii + 3];
      local_pi[4] = pi_from_m[6*ii + 4];
      local_pi[5] = pi_from_m[6*ii + 5];

      trace = local_pi[0] + local_pi[2] + local_pi[5];

      // update the q=0 sublattice

      n_a.vd[(0 + ii*LBQ ) * para.number_of_nodes + single_nodeindex] = 1.0f/3.0f * (local_rho - avg_rho) - 1.0f/2.0f*trace;

      // update the q=1 sublattice

      rho_times_coeff = 1.0f/18.0f * (local_rho - avg_rho);

      n_a.vd[(1 + ii*LBQ ) * para.number_of_nodes + single_nodeindex] =   rho_times_coeff + 1.0f/6.0f*local_j[0]
                                                                        + 1.0f/4.0f*local_pi[0] - 1.0f/12.0f*trace;
      n_a.vd[(2 + ii*LBQ ) * para.number_of_nodes + single_nodeindex] =   rho_times_coeff - 1.0f/6.0f*local_j[0]
                                                                        + 1.0f/4.0f*local_pi[0] - 1.0f/12.0f*trace;
      n_a.vd[(3 + ii*LBQ ) * para.number_of_nodes + single_nodeindex] =   rho_times_coeff + 1.0f/6.0f*local_j[1]
                                                                        + 1.0f/4.0f*local_pi[2] - 1.0f/12.0f*trace;
      n_a.vd[(4 + ii*LBQ ) * para.number_of_nodes + single_nodeindex] =   rho_times_coeff - 1.0f/6.0f*local_j[1]
                                                                        + 1.0f/4.0f*local_pi[2] - 1.0f/12.0f*trace;
      n_a.vd[(5 + ii*LBQ ) * para.number_of_nodes + single_nodeindex] =   rho_times_coeff + 1.0f/6.0f*local_j[2]
                                                                        + 1.0f/4.0f*local_pi[5] - 1.0f/12.0f*trace;
      n_a.vd[(6 + ii*LBQ ) * para.number_of_nodes + single_nodeindex] =   rho_times_coeff - 1.0f/6.0f*local_j[2]
                                                                        + 1.0f/4.0f*local_pi[5] - 1.0f/12.0f*trace;

      // update the q=2 sublattice

      rho_times_coeff = 1.0f/36.0f * (local_rho - avg_rho);

      tmp1 = local_pi[0] + local_pi[2];
      tmp2 = 2.0f*local_pi[1];

      n_a.vd[( 7 + ii*LBQ ) * para.number_of_nodes + single_nodeindex] =   rho_times_coeff + 1.0f/12.0f*(local_j[0]+local_j[1])
                                                                         + 1.0f/8.0f*(tmp1+tmp2) - 1.0f/24.0f*trace;
      n_a.vd[( 8 + ii*LBQ ) * para.number_of_nodes + single_nodeindex] =   rho_times_coeff - 1.0f/12.0f*(local_j[0]+local_j[1])
                                                                         + 1.0f/8.0f*(tmp1+tmp2) - 1.0f/24.0f*trace;
      n_a.vd[( 9 + ii*LBQ ) * para.number_of_nodes + single_nodeindex] =   rho_times_coeff + 1.0f/12.0f*(local_j[0]-local_j[1])
                                                                         + 1.0f/8.0f*(tmp1-tmp2) - 1.0f/24.0f*trace;
      n_a.vd[(10 + ii*LBQ ) * para.number_of_nodes + single_nodeindex] =   rho_times_coeff - 1.0f/12.0f*(local_j[0]-local_j[1])
                                                                         + 1.0f/8.0f*(tmp1-tmp2) - 1.0f/24.0f*trace;

      tmp1 = local_pi[0] + local_pi[5];
      tmp2 = 2.0f*local_pi[3];

      n_a.vd[(11 + ii*LBQ ) * para.number_of_nodes + single_nodeindex] =   rho_times_coeff + 1.0f/12.0f*(local_j[0]+local_j[2])
                                                                         + 1.0f/8.0f*(tmp1+tmp2) - 1.0f/24.0f*trace;
      n_a.vd[(12 + ii*LBQ ) * para.number_of_nodes + single_nodeindex] =  rho_times_coeff - 1.0f/12.0f*(local_j[0]+local_j[2])
                                                                         + 1.0f/8.0f*(tmp1+tmp2) - 1.0f/24.0f*trace;
      n_a.vd[(13 + ii*LBQ ) * para.number_of_nodes + single_nodeindex] =  rho_times_coeff + 1.0f/12.0f*(local_j[0]-local_j[2])
                                                                         + 1.0f/8.0f*(tmp1-tmp2) - 1.0f/24.0f*trace;
      n_a.vd[(14 + ii*LBQ ) * para.number_of_nodes + single_nodeindex] =  rho_times_coeff - 1.0f/12.0f*(local_j[0]-local_j[2])
                                                                         + 1.0f/8.0f*(tmp1-tmp2) - 1.0f/24.0f*trace;

      tmp1 = local_pi[2] + local_pi[5];
      tmp2 = 2.0f*local_pi[4];

      n_a.vd[(15 + ii*LBQ ) * para.number_of_nodes + single_nodeindex] =   rho_times_coeff + 1.0f/12.0f*(local_j[1]+local_j[2])
                                                                         + 1.0f/8.0f*(tmp1+tmp2) - 1.0f/24.0f*trace;
      n_a.vd[(16 + ii*LBQ ) * para.number_of_nodes + single_nodeindex] =   rho_times_coeff - 1.0f/12.0f*(local_j[1]+local_j[2])
                                                                         + 1.0f/8.0f*(tmp1+tmp2) - 1.0f/24.0f*trace;
      n_a.vd[(17 + ii*LBQ ) * para.number_of_nodes + single_nodeindex] =   rho_times_coeff + 1.0f/12.0f*(local_j[1]-local_j[2])
                                                                         + 1.0f/8.0f*(tmp1-tmp2) - 1.0f/24.0f*trace;
      n_a.vd[(18 + ii*LBQ ) * para.number_of_nodes + single_nodeindex] =   rho_times_coeff - 1.0f/12.0f*(local_j[1]-local_j[2])
                                                                         + 1.0f/8.0f*(tmp1-tmp2) - 1.0f/24.0f*trace;
    }

    // Calculate the modes for this node

    calc_m_from_n(n_a, single_nodeindex, mode_for_pi);

    // Update the density and velocity field for this mode

    update_rho_v(mode_for_pi, single_nodeindex, node_f, d_v);
  }
}



/**calculate mass of the whole fluid kernel
 * @param *sum    Pointer to result storage value (Output)
 * @param n_a     Pointer to local node residing in array a (Input)
*/
__global__ void calc_mass(LB_nodes_gpu n_a, float *sum) {
  float mode[4];

  unsigned int index = blockIdx.y * gridDim.x * blockDim.x + blockDim.x * blockIdx.x + threadIdx.x;

  if(index<para.number_of_nodes)
  {
    for(int ii=0;ii<LB_COMPONENTS;++ii)
    { 
      calc_mode(mode, n_a, index,ii);
      float Rho = mode[0] + para.rho[ii]*para.agrid*para.agrid*para.agrid;
      atomicadd(&(sum[0]), Rho);
    }
  }
}

/** (re-)initialization of the node force / set up of external force in lb units
 * @param node_f  Pointer to local node force (Input)
*/
__global__ void reinit_node_force(LB_node_force_gpu node_f){

  unsigned int index = blockIdx.y * gridDim.x * blockDim.x + blockDim.x * blockIdx.x + threadIdx.x;

  if(index<para.number_of_nodes)
  {
    #pragma unroll
    for(int ii=0;ii<LB_COMPONENTS;++ii)
    {
#ifdef EXTERNAL_FORCES
      if(para.external_force)
      {
        node_f.force[(0+ii*3)*para.number_of_nodes + index] = para.ext_force[0+ii*3]*para.agrid*para.agrid*para.tau*para.tau;
        node_f.force[(1+ii*3)*para.number_of_nodes + index] = para.ext_force[1+ii*3]*para.agrid*para.agrid*para.tau*para.tau;
        node_f.force[(2+ii*3)*para.number_of_nodes + index] = para.ext_force[2+ii*3]*para.agrid*para.agrid*para.tau*para.tau;
      }
      else
      {
        node_f.force[(0+ii*3)*para.number_of_nodes + index] = 0.0f;
        node_f.force[(1+ii*3)*para.number_of_nodes + index] = 0.0f;
        node_f.force[(2+ii*3)*para.number_of_nodes + index] = 0.0f;
      }
#else
      node_f.force[(0+ii*3)*para.number_of_nodes + index] = 0.0f;
      node_f.force[(1+ii*3)*para.number_of_nodes + index] = 0.0f;
      node_f.force[(2+ii*3)*para.number_of_nodes + index] = 0.0f;
#endif
    }
  }
}


/**set extern force on single nodes kernel
 * @param n_extern_nodeforces   number of nodes (Input)
 * @param *extern_nodeforces    Pointer to extern node force array (Input)
 * @param node_f                node force struct (Output)
*/
__global__ void init_extern_nodeforces(int n_extern_nodeforces, LB_extern_nodeforce_gpu *extern_nodeforces, LB_node_force_gpu node_f){

  unsigned int index = blockIdx.y * gridDim.x * blockDim.x + blockDim.x * blockIdx.x + threadIdx.x;
  float factor=powf(para.agrid,2)*para.tau*para.tau;
  if(index<n_extern_nodeforces)
  {
    #pragma unroll
    for(int ii=0;ii<LB_COMPONENTS;++ii)
    {
      node_f.force[(0+ii*3)*para.number_of_nodes + extern_nodeforces[index].index] = extern_nodeforces[index].force[0] * factor;
      node_f.force[(1+ii*3)*para.number_of_nodes + extern_nodeforces[index].index] = extern_nodeforces[index].force[1] * factor;
      node_f.force[(2+ii*3)*para.number_of_nodes + extern_nodeforces[index].index] = extern_nodeforces[index].force[2] * factor;
    }
  }
}

#ifdef SHANCHEN

/** 
 * @param single_nodeindex  Single node index        (Input)
 * @param component_index   Shanchen component index        (Input)
 * @param n_a               Pointer to local node residing in array a(Input)
*/
__device__ __inline__ float calc_massmode(LB_nodes_gpu n_a, int single_nodeindex, int component_index){

  /** mass mode */
  float mode;
  mode =   n_a.vd[( 0 + component_index*LBQ ) * para.number_of_nodes + single_nodeindex]
         + n_a.vd[( 1 + component_index*LBQ ) * para.number_of_nodes + single_nodeindex] + n_a.vd[( 2 + component_index*LBQ ) * para.number_of_nodes + single_nodeindex] 
         + n_a.vd[( 3 + component_index*LBQ ) * para.number_of_nodes + single_nodeindex] + n_a.vd[( 4 + component_index*LBQ ) * para.number_of_nodes + single_nodeindex]
         + n_a.vd[( 5 + component_index*LBQ ) * para.number_of_nodes + single_nodeindex] + n_a.vd[( 6 + component_index*LBQ ) * para.number_of_nodes + single_nodeindex]
         + n_a.vd[( 7 + component_index*LBQ ) * para.number_of_nodes + single_nodeindex] + n_a.vd[( 8 + component_index*LBQ ) * para.number_of_nodes + single_nodeindex]
         + n_a.vd[( 9 + component_index*LBQ ) * para.number_of_nodes + single_nodeindex] + n_a.vd[(10 + component_index*LBQ ) * para.number_of_nodes + single_nodeindex]
         + n_a.vd[(11 + component_index*LBQ ) * para.number_of_nodes + single_nodeindex] + n_a.vd[(12 + component_index*LBQ ) * para.number_of_nodes + single_nodeindex]
         + n_a.vd[(13 + component_index*LBQ ) * para.number_of_nodes + single_nodeindex] + n_a.vd[(14 + component_index*LBQ ) * para.number_of_nodes + single_nodeindex]
         + n_a.vd[(15 + component_index*LBQ ) * para.number_of_nodes + single_nodeindex] + n_a.vd[(16 + component_index*LBQ ) * para.number_of_nodes + single_nodeindex]
         + n_a.vd[(17 + component_index*LBQ ) * para.number_of_nodes + single_nodeindex] + n_a.vd[(18 + component_index*LBQ ) * para.number_of_nodes + single_nodeindex];

 mode += para.rho[component_index]*para.agrid*para.agrid*para.agrid;

 return mode;
}

__device__ __inline__ void calc_shanchen_contribution(LB_nodes_gpu n_a,int component_index, int x, int y, int z, float *p){ 

  float tmp_p[3]={0.0f,0.0f,0.0f};
  float pseudo;
  int index;

  index  = (x+1)%para.dim_x + para.dim_x*y + para.dim_x*para.dim_y*z;
  pseudo =  calc_massmode(n_a,index,component_index);
  tmp_p[0]+=pseudo/18.0f;

  index  = (para.dim_x+x-1)%para.dim_x + para.dim_x*y + para.dim_x*para.dim_y*z;
  pseudo =  calc_massmode(n_a,index,component_index);
  tmp_p[0]-=pseudo/18.0f;

  index  = x + para.dim_x*((y+1)%para.dim_y) + para.dim_x*para.dim_y*z;
  pseudo =  calc_massmode(n_a,index,component_index);
  tmp_p[1]+=pseudo/18.0f;

  index  = x + para.dim_x*((para.dim_y+y-1)%para.dim_y) + para.dim_x*para.dim_y*z;
  pseudo =  calc_massmode(n_a,index,component_index);
  tmp_p[1]-=pseudo/18.0f;

  index  = x + para.dim_x*y + para.dim_x*para.dim_y*((z+1)%para.dim_z);
  pseudo =  calc_massmode(n_a,index,component_index);
  tmp_p[2]+=pseudo/18.0f;

  index  = x + para.dim_x*y + para.dim_x*para.dim_y*((para.dim_z+z-1)%para.dim_z);
  pseudo =  calc_massmode(n_a,index,component_index);
  tmp_p[2]-=pseudo/18.0f;

  index  = (x+1)%para.dim_x + para.dim_x*((y+1)%para.dim_y) + para.dim_x*para.dim_y*z;
  pseudo =  calc_massmode(n_a,index,component_index);
  tmp_p[0]+=pseudo/36.0f;
  tmp_p[1]+=pseudo/36.0f;

  index  = (para.dim_x+x-1)%para.dim_x + para.dim_x*((para.dim_y+y-1)%para.dim_y) + para.dim_x*para.dim_y*z;
  pseudo =  calc_massmode(n_a,index,component_index);
  tmp_p[0]-=pseudo/36.0f;
  tmp_p[1]-=pseudo/36.0f;

  index  = (x+1)%para.dim_x + para.dim_x*((para.dim_y+y-1)%para.dim_y) + para.dim_x*para.dim_y*z;
  pseudo =  calc_massmode(n_a,index,component_index);
  tmp_p[0]+=pseudo/36.0f;
  tmp_p[1]-=pseudo/36.0f;

  index  = (para.dim_x+x-1)%para.dim_x + para.dim_x*((y+1)%para.dim_y) + para.dim_x*para.dim_y*z;
  pseudo =  calc_massmode(n_a,index,component_index);
  tmp_p[0]-=pseudo/36.0f;
  tmp_p[1]+=pseudo/36.0f;

  index  = (x+1)%para.dim_x + para.dim_x*y + para.dim_x*para.dim_y*((z+1)%para.dim_z);
  pseudo =  calc_massmode(n_a,index,component_index);
  tmp_p[0]+=pseudo/36.0f;
  tmp_p[2]+=pseudo/36.0f;

  index  = (para.dim_x+x-1)%para.dim_x + para.dim_x*y + para.dim_x*para.dim_y*((para.dim_z+z-1)%para.dim_z);
  pseudo =  calc_massmode(n_a,index,component_index);
  tmp_p[0]-=pseudo/36.0f;
  tmp_p[2]-=pseudo/36.0f;

  index  = (x+1)%para.dim_x + para.dim_x*y + para.dim_x*para.dim_y*((para.dim_z+z-1)%para.dim_z);
  pseudo =  calc_massmode(n_a,index,component_index);
  tmp_p[0]+=pseudo/36.0f;
  tmp_p[2]-=pseudo/36.0f;

  index  = (para.dim_x+x-1)%para.dim_x + para.dim_x*y + para.dim_x*para.dim_y*((z+1)%para.dim_z);
  pseudo =  calc_massmode(n_a,index,component_index);
  tmp_p[0]-=pseudo/36.0f;
  tmp_p[2]+=pseudo/36.0f;

  index  = x + para.dim_x*((y+1)%para.dim_y) + para.dim_x*para.dim_y*((z+1)%para.dim_z);
  pseudo =  calc_massmode(n_a,index,component_index);
  tmp_p[1]+=pseudo/36.0f;
  tmp_p[2]+=pseudo/36.0f;

  index  = x + para.dim_x*((para.dim_y+y-1)%para.dim_y) + para.dim_x*para.dim_y*((para.dim_z+z-1)%para.dim_z);
  pseudo =  calc_massmode(n_a,index,component_index);
  tmp_p[1]-=pseudo/36.0f;
  tmp_p[2]-=pseudo/36.0f;

  index  = x + para.dim_x*((y+1)%para.dim_y) + para.dim_x*para.dim_y*((para.dim_z+z-1)%para.dim_z);
  pseudo =  calc_massmode(n_a,index,component_index);
  tmp_p[1]+=pseudo/36.0f;
  tmp_p[2]-=pseudo/36.0f;
  
  index  = x + para.dim_x*((para.dim_y+y-1)%para.dim_y) + para.dim_x*para.dim_y*((z+1)%para.dim_z);
  pseudo =  calc_massmode(n_a,index,component_index);
  tmp_p[1]-=pseudo/36.0f;
  tmp_p[2]+=pseudo/36.0f;
 
  p[0]=tmp_p[0];
  p[1]=tmp_p[1];
  p[2]=tmp_p[2];
}

/** function to calc shanchen forces 
 * @param n_a     Pointer to local node residing in array a(Input)
 * @param node_f  Pointer to local node force (Input)
*/
__global__ void lb_shanchen_GPU(LB_nodes_gpu n_a,LB_node_force_gpu node_f){
#ifndef D3Q19
#error Lattices other than D3Q19 not supported
#endif
#if ( LB_COMPONENTS == 1  ) 
  #warning shanchen forces not implemented 
#else  

  unsigned int index = blockIdx.y * gridDim.x * blockDim.x + blockDim.x * blockIdx.x + threadIdx.x;
  unsigned int xyz[3];
  float pseudo;

  if(index<para.number_of_nodes)
  if( n_a.boundary[index] == 0 )
  {

    /* ShanChen forces are not reset at the end of the integration cycle, 
       in order to compute properly the hydrodynamic fields, so we have
       to reset them here. For the standard LB this is not needed */
     reset_LB_forces(index, node_f) ;
     /*Let's first identify the neighboring nodes */
     index_to_xyz(index, xyz);
     int x = xyz[0];
     int y = xyz[1];
     int z = xyz[2];
     
     #pragma unroll
     for(int ii=0;ii<LB_COMPONENTS;ii++)
     { 
       float p[3]={0.0f,0.0f,0.0f};
       pseudo =  calc_massmode(n_a,index,ii);

       #pragma unroll
       for(int jj=0;jj<LB_COMPONENTS;jj++)
       { 
             float tmpp[3]={0.0f,0.0f,0.0f};
             calc_shanchen_contribution(n_a, jj, x,y,z, tmpp);

             // FIXME  coupling HAS to be rescaled with agrid....
             p[0] += - para.coupling[(LB_COMPONENTS)*ii+jj]  * pseudo  * tmpp[0];
             p[1] += - para.coupling[(LB_COMPONENTS)*ii+jj]  * pseudo  * tmpp[1];
             p[2] += - para.coupling[(LB_COMPONENTS)*ii+jj]  * pseudo  * tmpp[2];
       }

       node_f.force[(0+ii*3)*para.number_of_nodes + index]+=p[0];
       node_f.force[(1+ii*3)*para.number_of_nodes + index]+=p[1];
       node_f.force[(2+ii*3)*para.number_of_nodes + index]+=p[2];
/* copy to be used when resetting forces */
       node_f.scforce[(0+ii*3)*para.number_of_nodes + index]=p[0];
       node_f.scforce[(1+ii*3)*para.number_of_nodes + index]=p[1];
       node_f.scforce[(2+ii*3)*para.number_of_nodes + index]=p[2];
    }
  }
#endif 
  return; 
}

#endif //SHANCHEN

/** kernel to set the local density
 *
 * @param n_a               the current nodes array (double buffering!)
 * @param single_nodeindex  the node to set the velocity for
 * @param rho               the density to set
 * @param d_v                Pointer to the local modes
*/
__global__ void set_rho(LB_nodes_gpu n_a,  LB_rho_v_gpu *d_v, int single_nodeindex,float *rho) {

  unsigned int index = blockIdx.y * gridDim.x * blockDim.x + blockDim.x * blockIdx.x + threadIdx.x;
  /*Note: this sets the velocities to zero */
  if(index == 0)
  {
    float local_rho;

    #pragma unroll
    for(int ii=0;ii<LB_COMPONENTS;++ii)
    { 
       /** default values for fields in lattice units */
       local_rho = (rho[ii]-para.rho[ii])*para.agrid*para.agrid*para.agrid;
       d_v[single_nodeindex].rho[ii]=rho[ii];

       n_a.vd[(0  + ii*LBQ ) * para.number_of_nodes + single_nodeindex] = 1.0f/ 3.0f * local_rho ;
       n_a.vd[(1  + ii*LBQ ) * para.number_of_nodes + single_nodeindex] = 1.0f/18.0f * local_rho ;
       n_a.vd[(2  + ii*LBQ ) * para.number_of_nodes + single_nodeindex] = 1.0f/18.0f * local_rho ;
       n_a.vd[(3  + ii*LBQ ) * para.number_of_nodes + single_nodeindex] = 1.0f/18.0f * local_rho ;
       n_a.vd[(4  + ii*LBQ ) * para.number_of_nodes + single_nodeindex] = 1.0f/18.0f * local_rho ;
       n_a.vd[(5  + ii*LBQ ) * para.number_of_nodes + single_nodeindex] = 1.0f/18.0f * local_rho ;
       n_a.vd[(6  + ii*LBQ ) * para.number_of_nodes + single_nodeindex] = 1.0f/18.0f * local_rho ;
       n_a.vd[(7  + ii*LBQ ) * para.number_of_nodes + single_nodeindex] = 1.0f/36.0f * local_rho ;
       n_a.vd[(8  + ii*LBQ ) * para.number_of_nodes + single_nodeindex] = 1.0f/36.0f * local_rho ;
       n_a.vd[(9  + ii*LBQ ) * para.number_of_nodes + single_nodeindex] = 1.0f/36.0f * local_rho ;
       n_a.vd[(10 + ii*LBQ ) * para.number_of_nodes + single_nodeindex] = 1.0f/36.0f * local_rho ;
       n_a.vd[(11 + ii*LBQ ) * para.number_of_nodes + single_nodeindex] = 1.0f/36.0f * local_rho ;
       n_a.vd[(12 + ii*LBQ ) * para.number_of_nodes + single_nodeindex] = 1.0f/36.0f * local_rho ;
       n_a.vd[(13 + ii*LBQ ) * para.number_of_nodes + single_nodeindex] = 1.0f/36.0f * local_rho ;
       n_a.vd[(14 + ii*LBQ ) * para.number_of_nodes + single_nodeindex] = 1.0f/36.0f * local_rho ;
       n_a.vd[(15 + ii*LBQ ) * para.number_of_nodes + single_nodeindex] = 1.0f/36.0f * local_rho ;
       n_a.vd[(16 + ii*LBQ ) * para.number_of_nodes + single_nodeindex] = 1.0f/36.0f * local_rho ;
       n_a.vd[(17 + ii*LBQ ) * para.number_of_nodes + single_nodeindex] = 1.0f/36.0f * local_rho ;
       n_a.vd[(18 + ii*LBQ ) * para.number_of_nodes + single_nodeindex] = 1.0f/36.0f * local_rho ;
    }
  }
}



/** adds the soon to be deleted momentum to the deleting boundaries momentum. necessary for momentum conservation
 * @param index             node index 
 * @param sign              1, -1, 0. add, subtract, don't do anything
 * @param current_boundary  pointer to the boundary struct which is about to occupy the node
 * @param n_in              node data which is used as input
*/
__device__ __inline__ void add_n_to_b(unsigned int index, int sign, LB_moving_boundary *current_boundary, LB_nodes_gpu n_in) {

  float v_temp[3] = {0, 0, 0};
  int c[3];
  float delta_n;
  short ii;
  if(sign){
    ii = 1;
    c[0] = 1; c[1] = 0; c[2] = 0;
    delta_n = (n_in.vd[ii * para.number_of_nodes + index] - n_in.vd[(ii+1) * para.number_of_nodes + index]);
      v_temp[0] += delta_n * c[0];
      v_temp[1] += delta_n * c[1];
      v_temp[2] += delta_n * c[2]; 
      
    ii += 2;//3
    c[0] = 0; c[1] = 1; c[2] = 0;
    delta_n = (n_in.vd[ii * para.number_of_nodes + index] - n_in.vd[(ii+1) * para.number_of_nodes + index]);
      v_temp[0] += delta_n * c[0];
      v_temp[1] += delta_n * c[1];
      v_temp[2] += delta_n * c[2]; 

    ii += 2;//5
    c[0] = 0; c[1] = 0; c[2] = 1;
    delta_n = (n_in.vd[ii * para.number_of_nodes + index] - n_in.vd[(ii+1) * para.number_of_nodes + index]);
      v_temp[0] += delta_n * c[0];
      v_temp[1] += delta_n * c[1];
      v_temp[2] += delta_n * c[2]; 

    ii += 2;//7
    c[0] = 1; c[1] = 1; c[2] = 0;
    delta_n = (n_in.vd[ii * para.number_of_nodes + index] - n_in.vd[(ii+1) * para.number_of_nodes + index]);
      v_temp[0] += delta_n * c[0];
      v_temp[1] += delta_n * c[1];
      v_temp[2] += delta_n * c[2]; 

    ii += 2;//9
    c[0] = 1; c[1] =-1; c[2] = 0;
    delta_n = (n_in.vd[ii * para.number_of_nodes + index] - n_in.vd[(ii+1) * para.number_of_nodes + index]);
      v_temp[0] += delta_n * c[0];
      v_temp[1] += delta_n * c[1];
      v_temp[2] += delta_n * c[2]; 

    ii += 2;//11
    c[0] = 1; c[1] = 0; c[2] = 1;
    delta_n = (n_in.vd[ii * para.number_of_nodes + index] - n_in.vd[(ii+1) * para.number_of_nodes + index]);
      v_temp[0] += delta_n * c[0];
      v_temp[1] += delta_n * c[1];
      v_temp[2] += delta_n * c[2]; 
      
    ii += 2;//13
    c[0] = 1; c[1] = 0; c[2] =-1;
    delta_n = (n_in.vd[ii * para.number_of_nodes + index] - n_in.vd[(ii+1) * para.number_of_nodes + index]);
      v_temp[0] += delta_n * c[0];
      v_temp[1] += delta_n * c[1];
      v_temp[2] += delta_n * c[2]; 
      
    ii += 2;//15
    c[0] = 0; c[1] = 1; c[2] = 1;
    delta_n = (n_in.vd[ii * para.number_of_nodes + index] - n_in.vd[(ii+1) * para.number_of_nodes + index]);
      v_temp[0] += delta_n * c[0];
      v_temp[1] += delta_n * c[1];
      v_temp[2] += delta_n * c[2]; 
      
    ii += 2;//17
    c[0] = 0; c[1] = 1; c[2] =-1;
    delta_n = (n_in.vd[ii * para.number_of_nodes + index] - n_in.vd[(ii+1) * para.number_of_nodes + index]);
      v_temp[0] += delta_n * c[0];
      v_temp[1] += delta_n * c[1];
      v_temp[2] += delta_n * c[2]; 
      
  atomicAdd(current_boundary->force_add + 0, sign * v_temp[0] );
  atomicAdd(current_boundary->force_add + 1, sign * v_temp[1] );
  atomicAdd(current_boundary->force_add + 2, sign * v_temp[2] );
  }
}
/**Determines local density by averaging over pre-push densities of neighboring fluid nodes
 **This is used to determine population of soon-to-be-fluid boundary node
 **[Reference goes here]
 * @param index             node index 
 * @param xyz               local node xyz computed by global fct
 * @param current_boundary  pointer to the boundary struct which is about to occupy the node
 * @param d_v               Pointer to local density array
 * @param n_in              node data which is used as input
*/
__device__ __inline__ void set_equilibrium_from_u_rho(unsigned int index, unsigned int *xyz, float* delta_xyz, LB_moving_boundary *current_boundary, LB_rho_v_gpu *d_v, LB_nodes_gpu n_in){
  
  short empty_neighbor_count = 0;
  size_t to_index, to_index_x, to_index_y, to_index_z;
  float u[3];
  float u_sq_neg;
  int c[3];
  float uc;
  float rho_avg = 0.0f;
  
  u[0] = current_boundary->velocity[0] - current_boundary->omega[1]*delta_xyz[2] + current_boundary->omega[2]*delta_xyz[1];//change of sign since r points to the center.
  u[1] = current_boundary->velocity[1] - current_boundary->omega[2]*delta_xyz[0] + current_boundary->omega[0]*delta_xyz[2];
  u[2] = current_boundary->velocity[2] - current_boundary->omega[0]*delta_xyz[1] + current_boundary->omega[1]*delta_xyz[0];
 
  u_sq_neg = -1.5f * (u[0]*u[0] + u[1]*u[1] + u[2]*u[2]);
  n_in.vd[0 * para.number_of_nodes + index] = u_sq_neg;
  
  c[0] = 1; c[1] = 0; c[2] = 0;
  
  uc = u[0]*c[0] + u[1]*c[1] + u[2]*c[2];
  n_in.vd[1 * para.number_of_nodes + index] = 3.0f*uc + 4.5f*uc*uc + u_sq_neg;
  
  to_index_x = (xyz[0] + c[0] + para.dim_x) % para.dim_x;
  to_index_y = (xyz[1] + c[1] + para.dim_y) % para.dim_y;
  to_index_z = (xyz[2] + c[2] + para.dim_z) % para.dim_z;
  to_index = to_index_x + to_index_y * para.dim_x + to_index_z * para.dim_x * para.dim_y;
  
  if(n_in.boundary[to_index]) empty_neighbor_count++;
  else rho_avg += d_v[to_index].rho[0];

  //this code repeats 18 times, the average is computed at the end of the function.
  
  c[0] =-1; c[1] = 0; c[2] = 0;
  
  uc = u[0]*c[0] + u[1]*c[1] + u[2]*c[2];
  n_in.vd[2 * para.number_of_nodes + index] = 3.0f*uc + 4.5f*uc*uc + u_sq_neg;
  
  to_index_x = (xyz[0] + c[0] + para.dim_x) % para.dim_x;
  to_index_y = (xyz[1] + c[1] + para.dim_y) % para.dim_y;
  to_index_z = (xyz[2] + c[2] + para.dim_z) % para.dim_z;
  to_index = to_index_x + to_index_y * para.dim_x + to_index_z * para.dim_x * para.dim_y;
  
  if(n_in.boundary[to_index]) empty_neighbor_count++;
  else rho_avg += d_v[to_index].rho[0];
  
  
  
  c[0] = 0; c[1] = 1; c[2] = 0;
  
  uc = u[0]*c[0] + u[1]*c[1] + u[2]*c[2];
  n_in.vd[3 * para.number_of_nodes + index] = 3.0f*uc + 4.5f*uc*uc + u_sq_neg;
  
  to_index_x = (xyz[0] + c[0] + para.dim_x) % para.dim_x;
  to_index_y = (xyz[1] + c[1] + para.dim_y) % para.dim_y;
  to_index_z = (xyz[2] + c[2] + para.dim_z) % para.dim_z;
  to_index = to_index_x + to_index_y * para.dim_x + to_index_z * para.dim_x * para.dim_y;
  
  if(n_in.boundary[to_index]) empty_neighbor_count++;
  else rho_avg += d_v[to_index].rho[0];
  
  
  
  c[0] = 0; c[1] =-1; c[2] = 0;
  
  uc = u[0]*c[0] + u[1]*c[1] + u[2]*c[2];
  n_in.vd[4 * para.number_of_nodes + index] = 3.0f*uc + 4.5f*uc*uc + u_sq_neg;
  
  to_index_x = (xyz[0] + c[0] + para.dim_x) % para.dim_x;
  to_index_y = (xyz[1] + c[1] + para.dim_y) % para.dim_y;
  to_index_z = (xyz[2] + c[2] + para.dim_z) % para.dim_z;
  to_index = to_index_x + to_index_y * para.dim_x + to_index_z * para.dim_x * para.dim_y;
  
  if(n_in.boundary[to_index]) empty_neighbor_count++;
  else rho_avg += d_v[to_index].rho[0];
  
  
  
  c[0] = 0; c[1] = 0; c[2] = 1;
   
  uc = u[0]*c[0] + u[1]*c[1] + u[2]*c[2];
  n_in.vd[5 * para.number_of_nodes + index] = 3.0f*uc + 4.5f*uc*uc + u_sq_neg;
  
  to_index_x = (xyz[0] + c[0] + para.dim_x) % para.dim_x;
  to_index_y = (xyz[1] + c[1] + para.dim_y) % para.dim_y;
  to_index_z = (xyz[2] + c[2] + para.dim_z) % para.dim_z;
  to_index = to_index_x + to_index_y * para.dim_x + to_index_z * para.dim_x * para.dim_y;
  
  if(n_in.boundary[to_index]) empty_neighbor_count++;
  else rho_avg += d_v[to_index].rho[0];
  
  
  
  c[0] = 0; c[1] = 0; c[2] =-1;
  
  uc = u[0]*c[0] + u[1]*c[1] + u[2]*c[2];
  n_in.vd[6 * para.number_of_nodes + index] = 3.0f*uc + 4.5f*uc*uc + u_sq_neg;
  
  to_index_x = (xyz[0] + c[0] + para.dim_x) % para.dim_x;
  to_index_y = (xyz[1] + c[1] + para.dim_y) % para.dim_y;
  to_index_z = (xyz[2] + c[2] + para.dim_z) % para.dim_z;
  to_index = to_index_x + to_index_y * para.dim_x + to_index_z * para.dim_x * para.dim_y;
  
  if(n_in.boundary[to_index]) empty_neighbor_count++;
  else rho_avg += d_v[to_index].rho[0];
  
  
  
  c[0] = 1; c[1] = 1; c[2] = 0;
  
  uc = u[0]*c[0] + u[1]*c[1] + u[2]*c[2];
  n_in.vd[7 * para.number_of_nodes + index] = 3.0f*uc + 4.5f*uc*uc + u_sq_neg;
  
  to_index_x = (xyz[0] + c[0] + para.dim_x) % para.dim_x;
  to_index_y = (xyz[1] + c[1] + para.dim_y) % para.dim_y;
  to_index_z = (xyz[2] + c[2] + para.dim_z) % para.dim_z;
  to_index = to_index_x + to_index_y * para.dim_x + to_index_z * para.dim_x * para.dim_y;
  
  if(n_in.boundary[to_index]) empty_neighbor_count++;
  else rho_avg += d_v[to_index].rho[0];
  
  
  
  c[0] =-1; c[1] =-1; c[2] = 0;
  
  uc = u[0]*c[0] + u[1]*c[1] + u[2]*c[2];
  n_in.vd[8 * para.number_of_nodes + index] = 3.0f*uc + 4.5f*uc*uc + u_sq_neg;
  
  to_index_x = (xyz[0] + c[0] + para.dim_x) % para.dim_x;
  to_index_y = (xyz[1] + c[1] + para.dim_y) % para.dim_y;
  to_index_z = (xyz[2] + c[2] + para.dim_z) % para.dim_z;
  to_index = to_index_x + to_index_y * para.dim_x + to_index_z * para.dim_x * para.dim_y;
  
  if(n_in.boundary[to_index]) empty_neighbor_count++;
  else rho_avg += d_v[to_index].rho[0];
  
  
  
  c[0] = 1; c[1] =-1; c[2] = 0;
  
  uc = u[0]*c[0] + u[1]*c[1] + u[2]*c[2];
  n_in.vd[9 * para.number_of_nodes + index] = 3.0f*uc + 4.5f*uc*uc + u_sq_neg;
  
  to_index_x = (xyz[0] + c[0] + para.dim_x) % para.dim_x;
  to_index_y = (xyz[1] + c[1] + para.dim_y) % para.dim_y;
  to_index_z = (xyz[2] + c[2] + para.dim_z) % para.dim_z;
  to_index = to_index_x + to_index_y * para.dim_x + to_index_z * para.dim_x * para.dim_y;
  
  if(n_in.boundary[to_index]) empty_neighbor_count++;
  else rho_avg += d_v[to_index].rho[0];
  
  
  
  c[0] =-1; c[1] = 1; c[2] = 0;
  
  uc = u[0]*c[0] + u[1]*c[1] + u[2]*c[2];
  n_in.vd[10* para.number_of_nodes + index] = 3.0f*uc + 4.5f*uc*uc + u_sq_neg;
  
  to_index_x = (xyz[0] + c[0] + para.dim_x) % para.dim_x;
  to_index_y = (xyz[1] + c[1] + para.dim_y) % para.dim_y;
  to_index_z = (xyz[2] + c[2] + para.dim_z) % para.dim_z;
  to_index = to_index_x + to_index_y * para.dim_x + to_index_z * para.dim_x * para.dim_y;
  
  if(n_in.boundary[to_index]) empty_neighbor_count++;
  else rho_avg += d_v[to_index].rho[0];
  
  
  
  c[0] = 1; c[1] = 0; c[2] = 1;
  
  uc = u[0]*c[0] + u[1]*c[1] + u[2]*c[2];
  n_in.vd[11* para.number_of_nodes + index] = 3.0f*uc + 4.5f*uc*uc + u_sq_neg;
  
  to_index_x = (xyz[0] + c[0] + para.dim_x) % para.dim_x;
  to_index_y = (xyz[1] + c[1] + para.dim_y) % para.dim_y;
  to_index_z = (xyz[2] + c[2] + para.dim_z) % para.dim_z;
  to_index = to_index_x + to_index_y * para.dim_x + to_index_z * para.dim_x * para.dim_y;
  
  if(n_in.boundary[to_index]) empty_neighbor_count++;
  else rho_avg += d_v[to_index].rho[0];
  
  
  
  c[0] =-1; c[1] = 0; c[2] =-1;
  
  uc = u[0]*c[0] + u[1]*c[1] + u[2]*c[2];
  n_in.vd[12* para.number_of_nodes + index] = 3.0f*uc + 4.5f*uc*uc + u_sq_neg;
  
  to_index_x = (xyz[0] + c[0] + para.dim_x) % para.dim_x;
  to_index_y = (xyz[1] + c[1] + para.dim_y) % para.dim_y;
  to_index_z = (xyz[2] + c[2] + para.dim_z) % para.dim_z;
  to_index = to_index_x + to_index_y * para.dim_x + to_index_z * para.dim_x * para.dim_y;
  
  if(n_in.boundary[to_index]) empty_neighbor_count++;
  else rho_avg += d_v[to_index].rho[0];
  
  
  
  c[0] = 1; c[1] = 0; c[2] =-1;
  
  uc = u[0]*c[0] + u[1]*c[1] + u[2]*c[2];
  n_in.vd[13* para.number_of_nodes + index] = 3.0f*uc + 4.5f*uc*uc + u_sq_neg;
  
  to_index_x = (xyz[0] + c[0] + para.dim_x) % para.dim_x;
  to_index_y = (xyz[1] + c[1] + para.dim_y) % para.dim_y;
  to_index_z = (xyz[2] + c[2] + para.dim_z) % para.dim_z;
  to_index = to_index_x + to_index_y * para.dim_x + to_index_z * para.dim_x * para.dim_y;
  
  if(n_in.boundary[to_index]) empty_neighbor_count++;
  else rho_avg += d_v[to_index].rho[0];
  
  
  
  c[0] =-1; c[1] = 0; c[2] = 1;
  
  uc = u[0]*c[0] + u[1]*c[1] + u[2]*c[2];
  n_in.vd[14* para.number_of_nodes + index] = 3.0f*uc + 4.5f*uc*uc + u_sq_neg;
  
  to_index_x = (xyz[0] + c[0] + para.dim_x) % para.dim_x;
  to_index_y = (xyz[1] + c[1] + para.dim_y) % para.dim_y;
  to_index_z = (xyz[2] + c[2] + para.dim_z) % para.dim_z;
  to_index = to_index_x + to_index_y * para.dim_x + to_index_z * para.dim_x * para.dim_y;
  
  if(n_in.boundary[to_index]) empty_neighbor_count++;
  else rho_avg += d_v[to_index].rho[0];
  
  
  
  c[0] = 0; c[1] = 1; c[2] = 1;
  
  uc = u[0]*c[0] + u[1]*c[1] + u[2]*c[2];
  n_in.vd[15* para.number_of_nodes + index] = 3.0f*uc + 4.5f*uc*uc + u_sq_neg;
  
  to_index_x = (xyz[0] + c[0] + para.dim_x) % para.dim_x;
  to_index_y = (xyz[1] + c[1] + para.dim_y) % para.dim_y;
  to_index_z = (xyz[2] + c[2] + para.dim_z) % para.dim_z;
  to_index = to_index_x + to_index_y * para.dim_x + to_index_z * para.dim_x * para.dim_y;
  
  if(n_in.boundary[to_index]) empty_neighbor_count++;
  else rho_avg += d_v[to_index].rho[0];
  
  
  
  c[0] = 0; c[1] =-1; c[2] =-1;
  
  uc = u[0]*c[0] + u[1]*c[1] + u[2]*c[2];
  n_in.vd[16* para.number_of_nodes + index] = 3.0f*uc + 4.5f*uc*uc + u_sq_neg;
  
  to_index_x = (xyz[0] + c[0] + para.dim_x) % para.dim_x;
  to_index_y = (xyz[1] + c[1] + para.dim_y) % para.dim_y;
  to_index_z = (xyz[2] + c[2] + para.dim_z) % para.dim_z;
  to_index = to_index_x + to_index_y * para.dim_x + to_index_z * para.dim_x * para.dim_y;
  
  if(n_in.boundary[to_index]) empty_neighbor_count++;
  else rho_avg += d_v[to_index].rho[0];
  
  
  
  c[0] = 0; c[1] = 1; c[2] =-1;
  
  uc = u[0]*c[0] + u[1]*c[1] + u[2]*c[2];
  n_in.vd[17* para.number_of_nodes + index] = 3.0f*uc + 4.5f*uc*uc + u_sq_neg;
  
  to_index_x = (xyz[0] + c[0] + para.dim_x) % para.dim_x;
  to_index_y = (xyz[1] + c[1] + para.dim_y) % para.dim_y;
  to_index_z = (xyz[2] + c[2] + para.dim_z) % para.dim_z;
  to_index = to_index_x + to_index_y * para.dim_x + to_index_z * para.dim_x * para.dim_y;
  
  if(n_in.boundary[to_index]) empty_neighbor_count++;
  else rho_avg += d_v[to_index].rho[0];
  
  
  
  c[0] = 0; c[1] =-1; c[2] = 1;
  
  uc = u[0]*c[0] + u[1]*c[1] + u[2]*c[2];
  n_in.vd[18* para.number_of_nodes + index] = 3.0f*uc + 4.5f*uc*uc + u_sq_neg;
  
  to_index_x = (xyz[0] + c[0] + para.dim_x) % para.dim_x;
  to_index_y = (xyz[1] + c[1] + para.dim_y) % para.dim_y;
  to_index_z = (xyz[2] + c[2] + para.dim_z) % para.dim_z;
  to_index = to_index_x + to_index_y * para.dim_x + to_index_z * para.dim_x * para.dim_y;
  
  if(n_in.boundary[to_index]) empty_neighbor_count++;
  else rho_avg += d_v[to_index].rho[0];
  
  //repetition done.
  rho_avg /= 18 - empty_neighbor_count;
#define W0 0.333333333f 
#define W1 0.0555555555f
#define W2 0.0277777777f
  n_in.vd[ 0 * para.number_of_nodes + index] *= rho_avg * W0;
  
  n_in.vd[ 1 * para.number_of_nodes + index] *= rho_avg * W1;
  n_in.vd[ 2 * para.number_of_nodes + index] *= rho_avg * W1;
  n_in.vd[ 3 * para.number_of_nodes + index] *= rho_avg * W1;
  n_in.vd[ 4 * para.number_of_nodes + index] *= rho_avg * W1;
  n_in.vd[ 5 * para.number_of_nodes + index] *= rho_avg * W1;
  n_in.vd[ 6 * para.number_of_nodes + index] *= rho_avg * W1;
  
  n_in.vd[ 7 * para.number_of_nodes + index] *= rho_avg * W2;
  n_in.vd[ 8 * para.number_of_nodes + index] *= rho_avg * W2;
  n_in.vd[ 9 * para.number_of_nodes + index] *= rho_avg * W2;
  n_in.vd[10 * para.number_of_nodes + index] *= rho_avg * W2;
  n_in.vd[11 * para.number_of_nodes + index] *= rho_avg * W2;
  n_in.vd[12 * para.number_of_nodes + index] *= rho_avg * W2;
  n_in.vd[13 * para.number_of_nodes + index] *= rho_avg * W2;
  n_in.vd[14 * para.number_of_nodes + index] *= rho_avg * W2;
  n_in.vd[15 * para.number_of_nodes + index] *= rho_avg * W2;
  n_in.vd[16 * para.number_of_nodes + index] *= rho_avg * W2;
  n_in.vd[17 * para.number_of_nodes + index] *= rho_avg * W2;
  n_in.vd[18 * para.number_of_nodes + index] *= rho_avg * W2;
#undef W0
#undef W1
#undef W2

}
__global__ void debug_set_eq(LB_moving_boundary *current_boundary, LB_rho_v_gpu *d_v, LB_nodes_gpu n_in){
  
  unsigned int index = blockIdx.y * gridDim.x * blockDim.x + blockDim.x * blockIdx.x + threadIdx.x;
  if(!index) printf("\n init velocity %f ", current_boundary->velocity[0]);
  unsigned int xyz[3];
  float dummy_delta[3] = {0, 0, 0};
  index_to_xyz(index, xyz);
  if(!n_in.boundary[index])
    set_equilibrium_from_u_rho(index, xyz, dummy_delta, current_boundary, d_v, n_in);
}
  
/**inits moving boundaries
 * @param index                 node index
 * @param lb_moving_boundary    Pointer to all structs holding boundary info
 * @param n_curr                pointer to LB node data which is overwritten
*/
__device__ void init_moving_boundaries(unsigned int index, LB_moving_boundary *lb_moving_boundary, LB_nodes_gpu n_curr, float *d_lab_anchors){

  unsigned int xyz[3];
  index_to_xyz(index, xyz);
  int ii, jj; //iterable
  float dist;
  float delta_xyz[3];
  float anchor_pos[3];
  int sign = 0;
    float *migrating_anchors = d_lab_anchors;
    if(n_curr.boundary[index]<=0){
      for(ii=0; ii<n_lb_moving_boundaries_gpu && !sign; ii++){
        xyz_center_to_delta(xyz, lb_moving_boundary[ii].center, delta_xyz);
        dist = sqrtf((delta_xyz[0])*(delta_xyz[0])
                   + (delta_xyz[1])*(delta_xyz[1])
                   + (delta_xyz[2])*(delta_xyz[2]));
        if(dist<=lb_moving_boundary[ii].radius){
          sign++;
          break;
        }
        else {
          for(jj=0; jj<lb_moving_boundary[ii].n_anchors; jj++){
      
            anchor_pos[0] = lb_moving_boundary[ii].center[0] + migrating_anchors[0];
            anchor_pos[1] = lb_moving_boundary[ii].center[1] + migrating_anchors[1];
            anchor_pos[2] = lb_moving_boundary[ii].center[2] + migrating_anchors[2];
            xyz_center_to_delta(xyz, anchor_pos, delta_xyz);
            dist = sqrtf((delta_xyz[0])*(delta_xyz[0])
                       + (delta_xyz[1])*(delta_xyz[1])
                       + (delta_xyz[2])*(delta_xyz[2]));
            if(dist<=migrating_anchors[3]){
              sign++;
              ii-=1;
	      break;
            }
          migrating_anchors += 4;
          }
        }  
      }
    if(sign) {
      n_curr.boundary_buffer[index] = -ii -1; 
      n_curr.boundary[index] = -ii -1;
    }
    else {
      n_curr.boundary_buffer[index] = 0;
      n_curr.boundary[index] = 0;
    }
  }
}
  
/**set the boundary flag for all boundary nodes
 * @param boundary_node_list    The indices of the boundary nodes
 * @param boundary_index_list   The flag representing the corresponding boundary
 * @param number_of_boundnodes  The number of boundary nodes
 * @param n_a                   Pointer to local node residing in array a (Input)
 * @param n_b                   Pointer to local node residing in array b (Input)
*/
__global__ void init_boundaries(int *boundary_node_list, int *boundary_index_list, int number_of_boundnodes, LB_moving_boundary *lb_moving_boundary, float *d_lab_anchors, LB_nodes_gpu n_a, LB_nodes_gpu n_b){

  unsigned int index = blockIdx.y * gridDim.x * blockDim.x + blockDim.x * blockIdx.x + threadIdx.x;
  
  
  if(index<number_of_boundnodes)
  {
    n_a.boundary[boundary_node_list[index]] = boundary_index_list[index];
    n_a.boundary_buffer[boundary_node_list[index]] = boundary_index_list[index];
    n_b.boundary[boundary_node_list[index]] = boundary_index_list[index];
    n_b.boundary_buffer[boundary_node_list[index]] = boundary_index_list[index];
  }
  if(index<para.number_of_nodes) 
  {
    init_moving_boundaries(index, lb_moving_boundary, n_a, d_lab_anchors);
    init_moving_boundaries(index, lb_moving_boundary, n_b, d_lab_anchors);
  }

}

/**reset the boundary flag of every node
 * @param n_a   Pointer to local node residing in array a (Input)
 * @param n_b   Pointer to local node residing in array b (Input)
*/
__global__ void reset_boundaries(LB_nodes_gpu n_a, LB_nodes_gpu n_b){

  size_t index = blockIdx.y * gridDim.x * blockDim.x + blockDim.x * blockIdx.x + threadIdx.x;

  if(index<para.number_of_nodes)
    n_a.boundary[index] = n_b.boundary[index] = 0;
}

/**sets boundary vd to zero
 * @param n_a   Pointer to local node residing in array a (Input)
 * @param n_b   Pointer to local node residing in array b (Input)
*/
__global__ void reset_boundary_vd(LB_nodes_gpu n_a, LB_nodes_gpu n_b){
  
  unsigned int index = blockIdx.y * gridDim.x * blockDim.x + blockDim.x * blockIdx.x + threadIdx.x;
  short ii; //iterable
  for(ii = 0; ii<19 && n_a.boundary[index]; ii++)
    n_a.vd[ii * para.number_of_nodes + index] = 0;
  for(ii = 0; ii<19 && n_b.boundary[index]; ii++)
    n_b.vd[ii * para.number_of_nodes + index] = 0;

}
/** reinits boundaries from new positions and interpolates new fluid nodes from their neighbors
 ** also adds absorbed nodes' momentum to boundaries && 
 ** subtracts created nodes' momentum from boundaries for conservation purposes.
 * @param n_curr                pointer to current node array after lb integration
 * @param *lb_moving_boundary   pointer to structs holding boundary information 
                                after lb integration && boundary integration
*/
__global__ void propagate_boundaries(LB_nodes_gpu n_curr, LB_rho_v_gpu *d_v, LB_moving_boundary *lb_moving_boundary, float* d_lab_anchors){

  unsigned int xyz[3], index;
  index = blockIdx.y * gridDim.x * blockDim.x + blockDim.x * blockIdx.x + threadIdx.x;
  if(index > para.dim_x*para.dim_y*para.dim_z) return;
  index_to_xyz(index, xyz);
  int ii, jj; //iterable
  float dist;
  float delta_xyz[3];
  float delta_xyz_center[3];
  float anchor_pos[3];
  int sign = 0;//for add_n_to_b
  if(n_curr.boundary[index] <= 0){ //when colliding with a stationary boundary, don't delete it
    float *migrating_anchors = d_lab_anchors;
    for(ii=0; ii<n_lb_moving_boundaries_gpu && !sign; ii++){
      xyz_center_to_delta(xyz, lb_moving_boundary[ii].center, delta_xyz_center);
      dist = sqrtf((delta_xyz_center[0])*(delta_xyz_center[0])
                 + (delta_xyz_center[1])*(delta_xyz_center[1])
                 + (delta_xyz_center[2])*(delta_xyz_center[2]));
      if(dist<=lb_moving_boundary[ii].radius){
        sign++;
        break;
      }
      else {
        for(jj=0; jj<lb_moving_boundary[ii].n_anchors; jj++){
          anchor_pos[0] = lb_moving_boundary[ii].center[0] + migrating_anchors[0];
          anchor_pos[1] = lb_moving_boundary[ii].center[1] + migrating_anchors[1];
          anchor_pos[2] = lb_moving_boundary[ii].center[2] + migrating_anchors[2];
          xyz_center_to_delta(xyz, anchor_pos, delta_xyz);
          dist = sqrtf((delta_xyz[0])*(delta_xyz[0])
                     + (delta_xyz[1])*(delta_xyz[1])
                     + (delta_xyz[2])*(delta_xyz[2]));
          if(dist<=migrating_anchors[3]){
            sign++;
            ii -= 1;
            break;
          }
        migrating_anchors += 4;
        }
      }  
    }
    if(sign) {
      if(n_curr.boundary[index]) sign = 0;//add momentum
      n_curr.boundary_buffer[index] = -ii-1; 
    }
    else {
      if(n_curr.boundary[index]) {
        sign = -1;//subtract momentum
        ii = -n_curr.boundary[index]-1;
        set_equilibrium_from_u_rho(index, xyz, delta_xyz_center, lb_moving_boundary +ii, d_v, n_curr);
      }
      n_curr.boundary_buffer[index] = 0;
    }
#ifdef LB_GPU_MOVING_MOMENTUM    
  add_n_to_b(index, sign, lb_moving_boundary +ii , n_curr);
#endif
  }
}
/** sets boundary flag from buffer and sets corresponding vd to zero. 
 ** necessary to avoid blockwise race-condition in set_equilibrium_from_u_rho(...) fkt
 ** which is called in propagate_boundaries(...)
 * @param n_curr  pointer to current node array after boundary propagation
*/
__global__ void update_boundaries_from_buffer(LB_nodes_gpu n_curr, LB_nodes_gpu n_next, LB_moving_boundary *lb_moving_boundary){//at first this function is extreme at zeroing out everything. this is inefficient but safe. afterwards I will test what exactly is necessary.
  
  unsigned int index = blockIdx.y * gridDim.x * blockDim.x + blockDim.x * blockIdx.x + threadIdx.x;
  short ii; //iterable
  
  if(index < para.number_of_nodes){
    
    n_curr.boundary[index] = n_curr.boundary_buffer[index];
    n_next.boundary[index] = n_curr.boundary_buffer[index];
    if (n_next.boundary[index]){
#pragma unroll
      for(ii = 0; ii<19; ii++){
        n_curr.vd[ii * para.number_of_nodes + index] = 0;
        n_next.vd[ii * para.number_of_nodes + index] = 0;
      }
    }
  }
}
/** integrationstep of the lb-fluid-solver
 * @param n_a     Pointer to local node residing in array a (Input)
 * @param n_b     Pointer to local node residing in array b (Input)
 * @param *d_v    Pointer to local device values (Input)
 * @param node_f  Pointer to local node force (Input)
 * @param ek_parameters_gpu  Pointer to the parameters for the electrokinetics (Input)
*/


__global__ void integrate(LB_nodes_gpu n_a, LB_nodes_gpu n_b, LB_rho_v_gpu *d_v, LB_node_force_gpu node_f, EK_parameters* ek_parameters_gpu) {
  /**every node is connected to a thread via the index*/
  unsigned int index = blockIdx.y * gridDim.x * blockDim.x + blockDim.x * blockIdx.x + threadIdx.x;
  /**the 19 moments (modes) are only temporary register values */
  float mode[19*LB_COMPONENTS];
  LB_randomnr_gpu rng;
  
  if( index < para.number_of_nodes )
  {
    /** storing the seed into a register value*/
    rng.seed = n_a.seed[index];
    /**calc_m_from_n*/
    calc_m_from_n(n_a, index, mode);
    /**lb_relax_modes*/
    relax_modes(mode, index, node_f,d_v);
    /**lb_thermalize_modes */
    if (para.fluct)
    {
      thermalize_modes(mode, index, &rng);
    }
#if  defined(EXTERNAL_FORCES)  ||   defined (SHANCHEN)  
    /**if external force is used apply node force */
    apply_forces(index, mode, node_f,d_v);
#else
    /**if particles are used apply node forces*/
    if (para.number_of_particles) apply_forces(index, mode, node_f,d_v); 
#endif
    /**lb_calc_n_from_modes_push*/
    normalize_modes(mode);
    /**calc of velocity densities and streaming with pbc*/
    calc_n_from_modes_push(n_b, mode, index);
    /** rewriting the seed back to the global memory*/
    n_b.seed[index] = rng.seed;
  }  
}

/** part interaction kernel
 * @param n_a                Pointer to local node residing in array a (Input)
 * @param *particle_data     Pointer to the particle position and velocity (Input)
 * @param *particle_force    Pointer to the particle force (Input)
 * @param *part              Pointer to the rn array of the particles (Input)
 * @param node_f             Pointer to local node force (Input)
 * @param *fluid_composition Pointer to the local fluid composition for the Shanchen
 * @param *d_v               Pointer to local device values
*/
__global__ void calc_fluid_particle_ia(LB_nodes_gpu n_a, CUDA_particle_data *particle_data, float *particle_force, CUDA_fluid_composition * fluid_composition, LB_node_force_gpu node_f, CUDA_particle_seed *part, LB_rho_v_gpu *d_v){

  unsigned int part_index = blockIdx.y * gridDim.x * blockDim.x + blockDim.x * blockIdx.x + threadIdx.x;
  unsigned int node_index[8];
  float delta[8];
  float delta_j[3*LB_COMPONENTS]; 
  float partgrad1[8*LB_COMPONENTS]; 
  float partgrad2[8*LB_COMPONENTS]; 
  float partgrad3[8*LB_COMPONENTS]; 
  LB_randomnr_gpu rng_part;
  if(part_index<para.number_of_particles)
  {
#if defined(IMMERSED_BOUNDARY) || defined(VIRTUAL_SITES_COM)
    if ( !particle_data[part_index].isVirtual )
#endif
    {
      rng_part.seed = part[part_index].seed;

      /**force acting on the particle. delta_j will be used later to compute the force that acts back onto the fluid. */
      calc_viscous_force(n_a, delta, partgrad1, partgrad2, partgrad3, particle_data, particle_force, fluid_composition,part_index, &rng_part, delta_j, node_index, d_v, 0);
      calc_node_force(delta, delta_j, partgrad1, partgrad2, partgrad3, node_index, node_f); 

#ifdef ENGINE
      if ( particle_data[part_index].swim.swimming ) {
        calc_viscous_force(n_a, delta, partgrad1, partgrad2, partgrad3, particle_data, particle_force, fluid_composition,part_index, &rng_part, delta_j, node_index, d_v, 1);
        calc_node_force(delta, delta_j, partgrad1, partgrad2, partgrad3, node_index, node_f);
      }
#endif

      /**force which acts back to the fluid node */
      part[part_index].seed = rng_part.seed;
    }
  }
}

/** part interaction kernel
 * @param n_a       Pointer to local node residing in array a (Input)
 * @param *particle_data    Pointer to the particle position and velocity (Input)
 * @param *particle_force   Pointer to the particle force (Input)
 * @param *part       Pointer to the rn array of the particles (Input)
 * @param node_f      Pointer to local node force (Input)
 * @param *d_v    Pointer to local device values
*/
__global__ void calc_fluid_particle_ia_three_point_couple(LB_nodes_gpu n_a, CUDA_particle_data *particle_data, float *particle_force, 
                                                                LB_node_force_gpu node_f, CUDA_particle_seed *part, LB_rho_v_gpu *d_v){

  unsigned int part_index = blockIdx.y * gridDim.x * blockDim.x + blockDim.x * blockIdx.x + threadIdx.x;
  unsigned int node_index[27];
  float delta[27];
  float delta_j[3*LB_COMPONENTS]; 
  LB_randomnr_gpu rng_part;
  if(part_index<para.number_of_particles){

    rng_part.seed = part[part_index].seed;
    /**force acting on the particle. delta_j will be used later to compute the force that acts back onto the fluid. */
    calc_viscous_force_three_point_couple(n_a, delta, particle_data, particle_force, part_index, &rng_part, delta_j, node_index,d_v,0);
    calc_node_force_three_point_couple(delta, delta_j, node_index, node_f);

#ifdef ENGINE
    if ( particle_data[part_index].swim.swimming ) {
      calc_viscous_force_three_point_couple(n_a, delta, particle_data, particle_force, part_index, &rng_part, delta_j, node_index,d_v,1);
      calc_node_force_three_point_couple(delta, delta_j, node_index, node_f);
    }
#endif

    /**force which acts back to the fluid node */
    part[part_index].seed = rng_part.seed;    
  }
}


#ifdef LB_BOUNDARIES_GPU
/**Bounce back boundary kernel
 * @param n_a         Pointer to local node residing in array a (Input)
 * @param n_b         Pointer to local node residing in array b (Input)
 * @param lb_boundary_velocity    The constant velocity at the boundary, set by the user (Input)
 * @param lb_boundary_force       The force on the boundary nodes (Output)
*/
__global__ void apply_boundaries(LB_nodes_gpu n_curr, LB_rho_v_gpu *d_v, LB_moving_boundary *lb_moving_boundary, float* lb_boundary_velocity, float* lb_boundary_force){

  unsigned int index = blockIdx.y * gridDim.x * blockDim.x + blockDim.x * blockIdx.x + threadIdx.x;

  if(index<para.number_of_nodes)
    bounce_back_boundaries(n_curr, d_v, lb_moving_boundary, index, lb_boundary_velocity, lb_boundary_force);
}

#ifdef SHANCHEN
__global__ void lb_shanchen_set_boundaries(LB_nodes_gpu n_curr){
/* This implements neutral boundary conditions for the shanchen fluid (i.e., 90 deg contact angle) */

  unsigned int index = blockIdx.y * gridDim.x * blockDim.x + blockDim.x * blockIdx.x + threadIdx.x;
  unsigned int xyz[3];
  if(index<para.number_of_nodes){
   if(n_curr.boundary[index] != 0 ) { 
    index_to_xyz(index, xyz);
    unsigned int x = xyz[0];
    unsigned int y = xyz[1];
    unsigned int z = xyz[2];
    unsigned int to_index_x,to_index_y,to_index_z,to_index;
    int c[3],count=0;

    for(int ii=0;ii<LB_COMPONENTS;ii++)
       for(int comp=0;comp<19;comp++)
          n_curr.vd[(comp + ii*LBQ ) * para.number_of_nodes + index]  =  0.0 ; 
    for(c[0]=-1;c[0]<=1;c[0]++){
       for(c[1]=-1;c[1]<=1;c[1]++){
          for(c[2]=-1;c[2]<=1;c[2]++){
             to_index_x = (x+c[0]+para.dim_x)%para.dim_x; 
             to_index_y = (y+c[1]+para.dim_y)%para.dim_y; 
             to_index_z = (z+c[2]+para.dim_z)%para.dim_z; 
             to_index = to_index_x + para.dim_x*to_index_y + para.dim_x*para.dim_y*to_index_z;  
	     if(n_curr.boundary[to_index] == 0 ) { 
                  for(int ii=0;ii<LB_COMPONENTS;ii++){
                     for(int comp=0;comp<19;comp++){ /* We copy all velocities: at the end we will need 
                                                        only the density mode, but this introduces no overhead anyway */
	                 n_curr.vd[(comp + ii*LBQ ) * para.number_of_nodes + index] += 
	                    n_curr.vd[(comp + ii*LBQ ) * para.number_of_nodes + to_index] ;
                         count++;
                     }
                  }
             }	      	
          }
       }
    }
    if(count>0)
      for(int ii=0;ii<LB_COMPONENTS;ii++)
        for(int comp=0;comp<19;comp++)
           n_curr.vd[(comp + ii*LBQ ) * para.number_of_nodes + index]  /= count ; 
  }
 }
}
#endif /* SHANCHEN */


#endif

/** get physical values of the nodes (density, velocity, ...)
 * @param n_a     Pointer to local node residing in array a (Input)
 * @param *p_v    Pointer to local print values (Output)
 * @param *d_v    Pointer to local device values (Input)
 * @param node_f  The forces on the LB nodes
*/
__global__ void get_mesoscopic_values_in_MD_units(LB_nodes_gpu n_a, LB_rho_v_pi_gpu *p_v,LB_rho_v_gpu *d_v, LB_node_force_gpu node_f) {
  unsigned int index = blockIdx.y * gridDim.x * blockDim.x + blockDim.x * blockIdx.x + threadIdx.x;

  if(index < para.number_of_nodes)
  {
    float mode[19*LB_COMPONENTS];
    calc_m_from_n(n_a, index, mode);
    calc_values_in_MD_units(n_a, mode, p_v, d_v, node_f, index, index);
  }
}

/** get boundary flags
 *  @param n_a                Pointer to local node residing in array a (Input)
 *  @param device_bound_array Pointer to local device values (Input)
 */
__global__ void lb_get_boundaries(LB_nodes_gpu n_a, int *device_bound_array){

  unsigned int index = blockIdx.y * gridDim.x * blockDim.x + blockDim.x * blockIdx.x + threadIdx.x;

  if(index<para.number_of_nodes)
   device_bound_array[index] = n_a.boundary[index];
}



/**print single node values kernel
 * @param single_nodeindex  index of the node (Input)
 * @param *d_p_v            Pointer to result storage array (Input)
 * @param n_a               Pointer to local node residing in array a (Input)
 * @param *d_v    Pointer to local device values
 * @param node_f  Pointer to local node force
*/
__global__ void lb_print_node(int single_nodeindex, LB_rho_v_pi_gpu *d_p_v, LB_nodes_gpu n_a, LB_rho_v_gpu * d_v, LB_node_force_gpu node_f){

  float mode[19*LB_COMPONENTS];
  unsigned int index = blockIdx.y * gridDim.x * blockDim.x + blockDim.x * blockIdx.x + threadIdx.x;

  if(index == 0)
  {
    calc_m_from_n(n_a, single_nodeindex, mode);
     
    /* the following actually copies rho and v from d_v, and calculates pi */
    calc_values_in_MD_units(n_a, mode, d_p_v, d_v, node_f, single_nodeindex, 0);
  }
}
__global__ void momentum(LB_nodes_gpu n_a, LB_rho_v_gpu * d_v, LB_node_force_gpu node_f, float *sum) {

  unsigned int index = blockIdx.y * gridDim.x * blockDim.x + blockDim.x * blockIdx.x + threadIdx.x;

  if(index<para.number_of_nodes)
  {
    float j[3]={0.0f,0.0f,0.0f};
    float mode[4];

    for(int ii=0 ; ii < LB_COMPONENTS ; ii++ )
    { 
      calc_mode(mode, n_a, index,ii);

      j[0] += mode[1]+node_f.force[(0+ii*3)*para.number_of_nodes + index];
      j[1] += mode[2]+node_f.force[(1+ii*3)*para.number_of_nodes + index];
      j[2] += mode[3]+node_f.force[(2+ii*3)*para.number_of_nodes + index];
    }

#ifdef LB_BOUNDARIES_GPU
    if(n_a.boundary[index])
      j[0]=j[1]=j[2]=0.0f;
#endif

    atomicadd(&(sum[0]), j[0]); 
    atomicadd(&(sum[1]), j[1]); 
    atomicadd(&(sum[2]), j[2]); 
  }
}
__global__ void remove_momentum(LB_nodes_gpu n_a, LB_rho_v_gpu * d_v, LB_node_force_gpu node_f, float *sum) {

  unsigned int index = blockIdx.y * gridDim.x * blockDim.x + blockDim.x * blockIdx.x + threadIdx.x;
  if(index<para.number_of_nodes){
    for(int ii=0 ; ii < LB_COMPONENTS ; ii++ ) { 
        node_f.force[(0+ii*3)*para.number_of_nodes + index]-=sum[0]/para.number_of_nodes;
        node_f.force[(1+ii*3)*para.number_of_nodes + index]-=sum[1]/para.number_of_nodes;
        node_f.force[(2+ii*3)*para.number_of_nodes + index]-=sum[2]/para.number_of_nodes;
    }
  }
}

/**print single node boundary flag
 * @param single_nodeindex  index of the node (Input)
 * @param *device_flag      Pointer to result storage array (Input)
 * @param n_a               Pointer to local node residing in array a (Input)
*/
__global__ void lb_get_boundary_flag(int single_nodeindex, int *device_flag, LB_nodes_gpu n_a){

  unsigned int index = blockIdx.y * gridDim.x * blockDim.x + blockDim.x * blockIdx.x + threadIdx.x;

  if(index == 0)
    device_flag[0] = n_a.boundary[single_nodeindex];
}

/**********************************************************************/
/* Host functions to setup and call kernels*/
/**********************************************************************/

void lb_get_para_pointer(LB_parameters_gpu** pointeradress) {
  if(cudaGetSymbolAddress((void**) pointeradress, para) != cudaSuccess)
  {
    fprintf(stderr, "Trouble getting address of LB parameters.\n"); //TODO give proper error message
    errexit();
  }
}

void lb_get_lbpar_pointer(LB_parameters_gpu** pointeradress) {
  *pointeradress = &lbpar_gpu;
}


void lb_get_boundary_force_pointer(float** pointeradress) {
#ifdef LB_BOUNDARIES_GPU
  *pointeradress = lb_boundary_force;
#endif
}

void lb_get_device_values_pointer(LB_rho_v_gpu** pointeradress) {
  *pointeradress = device_rho_v;
}

/**initialization for the lb gpu fluid called from host
 * @param *lbpar_gpu  Pointer to parameters to setup the lb field
*/
void lb_init_GPU(LB_parameters_gpu *lbpar_gpu){
#define free_and_realloc(var,size) \
  { if( (var) != NULL ) cudaFree((var)); cuda_safe_mem(cudaMalloc((void**)&var, size)); } //this isn't realloc, it's free and malloc.

  size_of_rho_v     = lbpar_gpu->number_of_nodes * sizeof(LB_rho_v_gpu);
  size_of_rho_v_pi  = lbpar_gpu->number_of_nodes * sizeof(LB_rho_v_pi_gpu);


  /** Allocate structs in device memory*/
  /* see the notes to the stucture device_rho_v_pi above...*/
  if(extended_values_flag==0) 
  {
    free_and_realloc(device_rho_v, size_of_rho_v);
  }
  else 
  {
    free_and_realloc(device_rho_v_pi, size_of_rho_v_pi);
  }

  /* TODO: this is a almost a copy copy of  device_rho_v think about eliminating it, and maybe pi can be added to device_rho_v in this case*/
  free_and_realloc(print_rho_v_pi  , size_of_rho_v_pi);
  free_and_realloc(nodes_a.vd      , lbpar_gpu->number_of_nodes * 19 * LB_COMPONENTS * sizeof(float));
  free_and_realloc(nodes_b.vd      , lbpar_gpu->number_of_nodes * 19 * LB_COMPONENTS * sizeof(float)); 
  free_and_realloc(node_f.force    , lbpar_gpu->number_of_nodes *  3 * LB_COMPONENTS * sizeof(lbForceFloat));
#if defined(IMMERSED_BOUNDARY) || defined(EK_DEBUG)
  free_and_realloc(node_f.force_buf    , lbpar_gpu->number_of_nodes *  3 * LB_COMPONENTS * sizeof(lbForceFloat));
#endif
#ifdef SHANCHEN
  free_and_realloc(node_f.scforce  , lbpar_gpu->number_of_nodes *  3 * LB_COMPONENTS * sizeof(float));
#endif

  free_and_realloc(nodes_a.seed           , lbpar_gpu->number_of_nodes * sizeof( unsigned int));
  free_and_realloc(nodes_a.boundary       , lbpar_gpu->number_of_nodes * sizeof( int));
  free_and_realloc(nodes_a.boundary_buffer, lbpar_gpu->number_of_nodes * sizeof( int));
  free_and_realloc(nodes_b.seed           , lbpar_gpu->number_of_nodes * sizeof( unsigned int));
  free_and_realloc(nodes_b.boundary       , lbpar_gpu->number_of_nodes * sizeof( int));
  free_and_realloc(nodes_b.boundary_buffer, lbpar_gpu->number_of_nodes * sizeof( int));
  /**write parameters in const memory*/
  cuda_safe_mem(cudaMemcpyToSymbol(para, lbpar_gpu, sizeof(LB_parameters_gpu)));

  /**check flag if lb gpu init works*/
  free_and_realloc(gpu_check, sizeof(int));

  if(h_gpu_check!=NULL)
    free(h_gpu_check);  

  h_gpu_check = (int*)Utils::malloc(sizeof(int));

  /** values for the kernel call */
  int threads_per_block = 64;
  int blocks_per_grid_y = 4;
  int blocks_per_grid_x = (lbpar_gpu->number_of_nodes + threads_per_block * blocks_per_grid_y - 1) /(threads_per_block * blocks_per_grid_y);
  dim3 dim_grid = make_uint3(blocks_per_grid_x, blocks_per_grid_y, 1);

  KERNELCALL(reset_boundaries, dim_grid, threads_per_block, (nodes_a, nodes_b));
  //TODO: xyz are constant and don't need double buffer. Consider storing them in params or somewhere else
  //          --Re: Watch out for moving boundaries, they may need an internal buffer  
  #ifdef SHANCHEN
  // TODO FIXME: 
  /* We must add shan-chen forces, which are zero only if the densities are uniform*/
  #endif

#if defined(ELECTROKINETICS)
  // We need to know if the electrokinetics is being used or not
  cuda_safe_mem(cudaMemcpyToSymbol(ek_initialized_gpu, &ek_initialized, sizeof(int)));
#endif


  /** calc of velocitydensities from given parameters and initialize the Node_Force array with zero */
  KERNELCALL(reinit_node_force, dim_grid, threads_per_block, (node_f));
  KERNELCALL(calc_n_from_rho_j_pi, dim_grid, threads_per_block, (nodes_a, device_rho_v, node_f, gpu_check));
 
  intflag = 1;
  current_nodes = &nodes_a;
  h_gpu_check[0] = 0;
  cuda_safe_mem(cudaMemcpy(h_gpu_check, gpu_check, sizeof(int), cudaMemcpyDeviceToHost));
//fprintf(stderr, "initialization of lb gpu code %i\n", lbpar_gpu->number_of_nodes);
  cudaThreadSynchronize();

#if __CUDA_ARCH__ >= 200
  if(!h_gpu_check[0])
  {
    fprintf(stderr, "initialization of lb gpu code failed! \n");
    errexit();
  }
#endif
}

/** reinitialization for the lb gpu fluid called from host
 * @param *lbpar_gpu  Pointer to parameters to setup the lb field
*/
void lb_reinit_GPU(LB_parameters_gpu *lbpar_gpu){

  /**write parameters in const memory*/
  cuda_safe_mem(cudaMemcpyToSymbol(para, lbpar_gpu, sizeof(LB_parameters_gpu)));
  
  /** values for the kernel call */
  int threads_per_block = 64;
  int blocks_per_grid_y = 4;
  int blocks_per_grid_x = (lbpar_gpu->number_of_nodes + threads_per_block * blocks_per_grid_y - 1) /(threads_per_block * blocks_per_grid_y);
  dim3 dim_grid = make_uint3(blocks_per_grid_x, blocks_per_grid_y, 1);

  /** calc of velocity densities from given parameters and initialize the Node_Force array with zero */
  KERNELCALL(calc_n_from_rho_j_pi, dim_grid, threads_per_block, (nodes_a, device_rho_v, node_f, gpu_check));
}

void lb_realloc_particles_GPU_leftovers(LB_parameters_gpu *lbpar_gpu){

  //copy parameters, especially number of parts to gpu mem
  cuda_safe_mem(cudaMemcpyToSymbol(para, lbpar_gpu, sizeof(LB_parameters_gpu)));
}

#ifdef LB_BOUNDARIES_GPU
/**
*/

void set_virtual_particles_and_lab_space_moving(LB_moving_boundary *lbb_in, Particle *part_out, int n_parts){
  int ii;
  float dummy[3];
  unsigned long int anchor_count = 0;
  unsigned long int old_anchor_count = 0;
  h_lab_anchors = (float*) realloc(h_lab_anchors, 1*sizeof(float));
  for(ii = 0; ii < n_parts; ii++){
    part_out[ii].r.p[0] = lbb_in[ii].center[0];
    part_out[ii].r.p[1] = lbb_in[ii].center[1];
    part_out[ii].r.p[2] = lbb_in[ii].center[2];
    
    part_out[ii].r.quat[0] = lbb_in[ii].quat[0];
    part_out[ii].r.quat[1] = lbb_in[ii].quat[1];
    part_out[ii].r.quat[2] = lbb_in[ii].quat[2];
    part_out[ii].r.quat[3] = lbb_in[ii].quat[3];
      lb_convert_quat_to_quatu(part_out[ii].r.quat, part_out[ii].r.quatu);
    
    part_out[ii].m.v[0] = lbb_in[ii].velocity[0];
    part_out[ii].m.v[1] = lbb_in[ii].velocity[1];
    part_out[ii].m.v[2] = lbb_in[ii].velocity[2];
    convert_omega_space_to_body(lbb_in + ii, part_out + ii);
      part_out[ii].m.somega[0] /= lbpar_gpu.tau;
      part_out[ii].m.somega[1] /= lbpar_gpu.tau;
      part_out[ii].m.somega[2] /= lbpar_gpu.tau;
    
    part_out[ii].f.sf[0] = 0;
    part_out[ii].f.sf[1] = 0;
    part_out[ii].f.sf[2] = 0;
    part_out[ii].f.storque[0] = 0;
    part_out[ii].f.storque[1] = 0;
    part_out[ii].f.storque[2] = 0;

    part_out[ii].p.mass = lbb_in[ii].scaled_mass;//this is mass 
      part_out[ii].p.rinertia[0] = lbb_in[ii].rinertia[0];
      part_out[ii].p.rinertia[1] = lbb_in[ii].rinertia[1];
      part_out[ii].p.rinertia[2] = lbb_in[ii].rinertia[2];
#ifdef EXTERNAL_FORCES
    part_out[ii].p.ext_force[0] = lbb_in[ii].ext_force[0]; 
    part_out[ii].p.ext_force[1] = lbb_in[ii].ext_force[1]; 
    part_out[ii].p.ext_force[2] = lbb_in[ii].ext_force[2]; 
    part_out[ii].p.ext_torque[0] = lbb_in[ii].ext_torque[0]; 
    part_out[ii].p.ext_torque[1] = lbb_in[ii].ext_torque[1]; 
    part_out[ii].p.ext_torque[2] = lbb_in[ii].ext_torque[2]; 
    part_out[ii].p.lab_force[0] = 0;
    part_out[ii].p.lab_force[1] = 0;
    part_out[ii].p.lab_force[2] = 0;
    part_out[ii].p.body_torque[0] = lbb_in[ii].body_torque[0];
    part_out[ii].p.body_torque[1] = lbb_in[ii].body_torque[1];
    part_out[ii].p.body_torque[2] = lbb_in[ii].body_torque[2];
    part_out[ii].p.body_force[0] = lbb_in[ii].body_force[0];
    part_out[ii].p.body_force[1] = lbb_in[ii].body_force[1];
    part_out[ii].p.body_force[2] = lbb_in[ii].body_force[2];
    
#endif    
    part_out[ii].p.n_anchors = lbb_in[ii].n_anchors;
    part_out[ii].p.anchors = lbb_in[ii].anchors;//!free: lbb_in anchors pointer will be deleted here and a copy kept in part_out 
    lbb_in[ii].anchors = NULL;//So that a deep free implementation does not delete p.p.anchors
      old_anchor_count = anchor_count;
      anchor_count += lbb_in[ii].n_anchors;
    h_lab_anchors = (float*) realloc(h_lab_anchors, (4*anchor_count+1)*sizeof(float)); // +1 for good luck ...and because cuda doesn't like nullpointers & memcpy(0)
    convert_omega_anchors_body_to_space(part_out + ii, dummy, h_lab_anchors+old_anchor_count);


  }
  h_lab_anchors[4*anchor_count] = anchor_count;//might come in handy
  free_and_realloc(d_lab_anchors, (4*anchor_count+5)*sizeof(float));//this isn't a realloc, it's actually free and malloc
  //need redundant size because cuda doesn't like to have a pointer increment in an if statement, even when it doesn't activate.
  cuda_safe_mem(cudaMemcpy(d_lab_anchors, h_lab_anchors, (4*anchor_count+1)*sizeof(float), cudaMemcpyHostToDevice));
  lbpar_gpu.number_of_anchors = anchor_count;
}
    
    
    
/** setup and call boundaries from the host
 * @param host_n_lb_boundaries        number of LB boundaries
 * @param h_n_lb_moving_boundaries number of LB moving boundaries
 * @param number_of_boundnodes        number of boundnodes
 * @param host_boundary_node_list     The indices of the boundary nodes
 * @param host_boundary_index_list    The flag representing the corresponding boundary
 * @param host_lb_boundary_velocity   The constant velocity at the boundary, set by the user (Input)
*/
void lb_init_boundaries_GPU(int host_n_lb_boundaries, int h_n_lb_moving_boundaries, int number_of_boundnodes, int *host_boundary_node_list, int* host_boundary_index_list, float* host_lb_boundary_velocity, LB_moving_boundary *host_moving_boundary){
  int temp = host_n_lb_boundaries;
  int temp2 = h_n_lb_moving_boundaries;
  //set some time vars for rotation
  tau_half = lbpar_gpu.tau / 2;
  tau_sq = lbpar_gpu.tau * lbpar_gpu.tau;
  
  //perpare host structs for moving boundaries for integration and cudamalloc array which holds lab space anchor data
  host_lb_moving_boundary = (Particle*) realloc(host_lb_moving_boundary, h_n_lb_moving_boundaries * sizeof(Particle));
  host_n_lb_moving_boundaries = h_n_lb_moving_boundaries;
  set_virtual_particles_and_lab_space_moving(host_moving_boundary, host_lb_moving_boundary, h_n_lb_moving_boundaries);
  
  size_of_boundindex = number_of_boundnodes*sizeof(int);
  cuda_safe_mem(cudaMalloc((void**)&boundary_node_list, size_of_boundindex));
  cuda_safe_mem(cudaMalloc((void**)&boundary_index_list, size_of_boundindex));
  cuda_safe_mem(cudaMemcpy(boundary_index_list, host_boundary_index_list, size_of_boundindex, cudaMemcpyHostToDevice));
  cuda_safe_mem(cudaMemcpy(boundary_node_list, host_boundary_node_list, size_of_boundindex, cudaMemcpyHostToDevice));
  cuda_safe_mem(cudaMalloc((void**)&lb_boundary_force   , 3*host_n_lb_boundaries*sizeof(float)));
  cuda_safe_mem(cudaMalloc((void**)&lb_boundary_velocity, 3*host_n_lb_boundaries*sizeof(float)));
  cuda_safe_mem(cudaMemcpy(lb_boundary_velocity, host_lb_boundary_velocity, 3*n_lb_boundaries*sizeof(float), cudaMemcpyHostToDevice));
  cuda_safe_mem(cudaMalloc((void**)&lb_moving_boundary  , h_n_lb_moving_boundaries * sizeof(LB_moving_boundary)));
  cuda_safe_mem(cudaMemcpy(lb_moving_boundary, host_moving_boundary, h_n_lb_moving_boundaries * sizeof(LB_moving_boundary), cudaMemcpyHostToDevice));
  cuda_safe_mem(cudaMemcpyToSymbol(n_lb_boundaries_gpu, &temp, sizeof(int)));
  cuda_safe_mem(cudaMemcpyToSymbol(n_lb_moving_boundaries_gpu, &temp2, sizeof(int)));
  host_lb_moving_boundary = (Particle*) realloc(host_lb_moving_boundary, h_n_lb_moving_boundaries * sizeof(Particle));
  host_n_lb_moving_boundaries = h_n_lb_moving_boundaries;

  /** values for the kernel call */
  int threads_per_block = 64;
  int blocks_per_grid_y = 4;
  int blocks_per_grid_x = (lbpar_gpu.number_of_nodes + threads_per_block * blocks_per_grid_y - 1) /(threads_per_block * blocks_per_grid_y);
  dim3 dim_grid = make_uint3(blocks_per_grid_x, blocks_per_grid_y, 1);
  KERNELCALL(reset_boundaries, dim_grid, threads_per_block, (nodes_a, nodes_b));

  if (n_lb_boundaries == 0 && !pdb_boundary_lattice)
  {
    cudaThreadSynchronize();
    return;
  }

  if(number_of_boundnodes == 0 && n_lb_moving_boundaries == 0)
  {
    fprintf(stderr, "WARNING: boundary cmd executed but no boundary node found!\n");
  }
  else
  {
    int threads_per_block_bound = 64;
    int blocks_per_grid_bound_y = 4;
    int blocks_per_grid_bound_x = (lbpar_gpu.number_of_nodes + threads_per_block * blocks_per_grid_y - 1) /(threads_per_block * blocks_per_grid_y);
    dim3 dim_grid_bound = make_uint3(blocks_per_grid_bound_x, blocks_per_grid_bound_y, 1);

    KERNELCALL(init_boundaries, dim_grid_bound, threads_per_block_bound, (boundary_node_list, boundary_index_list, number_of_boundnodes, lb_moving_boundary, d_lab_anchors, nodes_a, nodes_b));
    KERNELCALL(reset_boundary_vd, dim_grid_bound, threads_per_block_bound, (nodes_a, nodes_b));

  }

  cudaThreadSynchronize();
}
#endif
/**setup and call extern single node force initialization from the host
 * @param *lbpar_gpu    Pointer to host parameter struct
*/
void lb_reinit_extern_nodeforce_GPU(LB_parameters_gpu *lbpar_gpu){

  cuda_safe_mem(cudaMemcpyToSymbol(para, lbpar_gpu, sizeof(LB_parameters_gpu))); 

  /** values for the kernel call */
  int threads_per_block = 64;
  int blocks_per_grid_y = 4;
  int blocks_per_grid_x = (lbpar_gpu->number_of_nodes + threads_per_block * blocks_per_grid_y - 1) /(threads_per_block * blocks_per_grid_y);
  dim3 dim_grid = make_uint3(blocks_per_grid_x, blocks_per_grid_y, 1);

  KERNELCALL(reinit_node_force, dim_grid, threads_per_block, (node_f));

}
/**setup and call extern single node force initialization from the host
 * @param n_extern_nodeforces       number of nodes on which the external force has to be applied
 * @param *host_extern_nodeforces   Pointer to the host extern node forces
 * @param *lbpar_gpu                Pointer to host parameter struct
*/
void lb_init_extern_nodeforces_GPU(int n_extern_nodeforces, LB_extern_nodeforce_gpu *host_extern_nodeforces, LB_parameters_gpu *lbpar_gpu){

  size_of_extern_nodeforces = n_extern_nodeforces*sizeof(LB_extern_nodeforce_gpu);
  cuda_safe_mem(cudaMalloc((void**)&extern_nodeforces, size_of_extern_nodeforces));
  cuda_safe_mem(cudaMemcpy(extern_nodeforces, host_extern_nodeforces, size_of_extern_nodeforces, cudaMemcpyHostToDevice));

  if(lbpar_gpu->external_force == 0)
    cuda_safe_mem(cudaMemcpyToSymbol(para, lbpar_gpu, sizeof(LB_parameters_gpu))); 

  int threads_per_block_exf = 64;
  int blocks_per_grid_exf_y = 4;
  int blocks_per_grid_exf_x = (n_extern_nodeforces + threads_per_block_exf * blocks_per_grid_exf_y - 1) / 
                              (threads_per_block_exf * blocks_per_grid_exf_y);
  dim3 dim_grid_exf = make_uint3(blocks_per_grid_exf_x, blocks_per_grid_exf_y, 1);

  KERNELCALL(init_extern_nodeforces, dim_grid_exf, threads_per_block_exf, (n_extern_nodeforces, extern_nodeforces, node_f));
  cudaFree(extern_nodeforces);
}

/**setup and call particle kernel from the host
*/
void lb_calc_particle_lattice_ia_gpu(){
  if (lbpar_gpu.number_of_particles) 
  {
    /** call of the particle kernel */
    /** values for the particle kernel */
    int threads_per_block_particles = 64;
    int blocks_per_grid_particles_y = 4;
    int blocks_per_grid_particles_x = (lbpar_gpu.number_of_particles + threads_per_block_particles * blocks_per_grid_particles_y - 1) / 
                                      (threads_per_block_particles * blocks_per_grid_particles_y);
    dim3 dim_grid_particles = make_uint3(blocks_per_grid_particles_x, blocks_per_grid_particles_y, 1);

    if ( lbpar_gpu.lb_couple_switch & LB_COUPLE_TWO_POINT )
    {
      KERNELCALL( calc_fluid_particle_ia, dim_grid_particles, threads_per_block_particles, 
                  ( *current_nodes, gpu_get_particle_pointer(), 
                    gpu_get_particle_force_pointer(), gpu_get_fluid_composition_pointer(),
                    node_f, gpu_get_particle_seed_pointer(), device_rho_v )
                );
    }
    else { /** only other option is the three point coupling scheme */
#ifdef SHANCHEN
#if __CUDA_ARCH__ >= 200
      fprintf (stderr, "The three point particle coupling is not currently compatible with the Shan-Chen implementation of the LB\n");
      errexit(); 
#endif
#endif
      KERNELCALL( calc_fluid_particle_ia_three_point_couple, dim_grid_particles, threads_per_block_particles,
                   ( *current_nodes, gpu_get_particle_pointer(),
                     gpu_get_particle_force_pointer(), node_f,
                     gpu_get_particle_seed_pointer(), device_rho_v )
                );
    }
  }
}

/** setup and call kernel for getting macroscopic fluid values of all nodes
 * @param *host_values struct to save the gpu values
*/
void lb_get_values_GPU(LB_rho_v_pi_gpu *host_values){

  /** values for the kernel call */
  int threads_per_block = 64;
  int blocks_per_grid_y = 4;
  int blocks_per_grid_x = (lbpar_gpu.number_of_nodes + threads_per_block * blocks_per_grid_y - 1) / 
                          (threads_per_block * blocks_per_grid_y);
  dim3 dim_grid = make_uint3(blocks_per_grid_x, blocks_per_grid_y, 1);

  KERNELCALL( get_mesoscopic_values_in_MD_units, dim_grid, threads_per_block,
              ( *current_nodes, print_rho_v_pi, device_rho_v, node_f ) );
  cuda_safe_mem( cudaMemcpy( host_values, print_rho_v_pi, size_of_rho_v_pi, cudaMemcpyDeviceToHost ) );

}

/** get all the boundary flags for all nodes
 *  @param host_bound_array here go the values of the boundary flag
 */
void lb_get_boundary_flags_GPU(int* host_bound_array){
   
  int* device_bound_array;
  cuda_safe_mem(cudaMalloc((void**)&device_bound_array, lbpar_gpu.number_of_nodes*sizeof(int)));
  /** values for the kernel call */
  int threads_per_block = 64;
  int blocks_per_grid_y = 4;
  int blocks_per_grid_x = (lbpar_gpu.number_of_nodes + threads_per_block * blocks_per_grid_y - 1) / (threads_per_block * blocks_per_grid_y);
  dim3 dim_grid = make_uint3(blocks_per_grid_x, blocks_per_grid_y, 1);

  KERNELCALL(lb_get_boundaries, dim_grid, threads_per_block, (*current_nodes, device_bound_array));

  cuda_safe_mem(cudaMemcpy(host_bound_array, device_bound_array, lbpar_gpu.number_of_nodes*sizeof(int), cudaMemcpyDeviceToHost));

  cudaFree(device_bound_array);

}

/** setup and call kernel for getting macroscopic fluid values of a single node*/
void lb_print_node_GPU(int single_nodeindex, LB_rho_v_pi_gpu *host_print_values){ 
      
  LB_rho_v_pi_gpu *device_print_values;
  cuda_safe_mem(cudaMalloc((void**)&device_print_values, sizeof(LB_rho_v_pi_gpu)));
  int threads_per_block_print = 1;
  int blocks_per_grid_print_y = 1;
  int blocks_per_grid_print_x = 1;
  dim3 dim_grid_print = make_uint3(blocks_per_grid_print_x, blocks_per_grid_print_y, 1);

  KERNELCALL(lb_print_node, dim_grid_print, threads_per_block_print, (single_nodeindex, device_print_values, *current_nodes, device_rho_v, node_f));

  cuda_safe_mem(cudaMemcpy(host_print_values, device_print_values, sizeof(LB_rho_v_pi_gpu), cudaMemcpyDeviceToHost));
  cudaFree(device_print_values);

}
int lb_lbfluid_print_moving_pos(int part_num, double *print_val){
  if(part_num < 0 || part_num >= host_n_lb_moving_boundaries) return 1;
  print_val[0]=host_lb_moving_boundary[part_num].r.p[0]*lbpar_gpu.agrid;
  print_val[1]=host_lb_moving_boundary[part_num].r.p[1]*lbpar_gpu.agrid;
  print_val[2]=host_lb_moving_boundary[part_num].r.p[2]*lbpar_gpu.agrid;

  return 0;
}

int lb_lbfluid_print_moving_vel(int part_num, double *print_val){
  if(part_num < 0 || part_num >= host_n_lb_moving_boundaries) return 1;   
  print_val[0]=host_lb_moving_boundary[part_num].m.v[0] / lbpar_gpu.tau * lbpar_gpu.agrid;
  print_val[1]=host_lb_moving_boundary[part_num].m.v[1] / lbpar_gpu.tau * lbpar_gpu.agrid;
  print_val[2]=host_lb_moving_boundary[part_num].m.v[2] / lbpar_gpu.tau * lbpar_gpu.agrid;

  return 0;
}

int lb_lbfluid_print_moving_omega_body(int part_num, double *print_val){
  if(part_num < 0 || part_num >= host_n_lb_moving_boundaries) return 1;  
  print_val[0]=host_lb_moving_boundary[part_num].m.somega[0];
  print_val[1]=host_lb_moving_boundary[part_num].m.somega[1];
  print_val[2]=host_lb_moving_boundary[part_num].m.somega[2];

  return 0;
}

int lb_lbfluid_print_moving_omega_lab(int part_num, double *print_val){
  if(part_num < 0 || part_num >= host_n_lb_moving_boundaries) return 1;
  lb_convert_omega_to_space_for_print(host_lb_moving_boundary + part_num, print_val);
  return 0;
}

int lb_lbfluid_print_moving_torque_body(int part_num, double *print_val){
  if(part_num < 0 || part_num >= host_n_lb_moving_boundaries) return 1;
  print_val[0]=host_lb_moving_boundary[part_num].f.storque[0];
  print_val[1]=host_lb_moving_boundary[part_num].f.storque[1];
  print_val[2]=host_lb_moving_boundary[part_num].f.storque[2];

  return 0;
}

int lb_lbfluid_print_moving_torque_lab(int part_num, double *print_val){
  if(part_num < 0 || part_num >= host_n_lb_moving_boundaries) return 1;
  lb_convert_torque_to_space_for_print(host_lb_moving_boundary + part_num, print_val);
  return 0;
}

int lb_lbfluid_print_moving_force_body(int part_num, double *print_val){
  if(part_num < 0 || part_num >= host_n_lb_moving_boundaries) return 1;  
  lb_convert_force_to_body_for_print(host_lb_moving_boundary + part_num, print_val);
  return 0;
}

int lb_lbfluid_print_moving_force_lab(int part_num, double *print_val){
  if(part_num < 0 || part_num >= host_n_lb_moving_boundaries) return 1;  
  print_val[0]=host_lb_moving_boundary[part_num].f.sf[0]/tau_sq;
  print_val[1]=host_lb_moving_boundary[part_num].f.sf[1]/tau_sq;
  print_val[2]=host_lb_moving_boundary[part_num].f.sf[2]/tau_sq;
  return 0;
}
/** setup and call kernel to calculate the total momentum of the hole fluid
 * @param *mass value of the mass calcutated on the GPU
*/
void lb_calc_fluid_mass_GPU(double* mass){

  float* tot_mass;
  float cpu_mass =  0.0f ;
  cuda_safe_mem(cudaMalloc((void**)&tot_mass, sizeof(float)));
  cuda_safe_mem(cudaMemcpy(tot_mass, &cpu_mass, sizeof(float), cudaMemcpyHostToDevice));

  /** values for the kernel call */
  int threads_per_block = 64;
  int blocks_per_grid_y = 4;
  int blocks_per_grid_x = (lbpar_gpu.number_of_nodes + threads_per_block * blocks_per_grid_y - 1) /(threads_per_block * blocks_per_grid_y);
  dim3 dim_grid = make_uint3(blocks_per_grid_x, blocks_per_grid_y, 1);

  KERNELCALL(calc_mass, dim_grid, threads_per_block,(*current_nodes, tot_mass));

  cuda_safe_mem(cudaMemcpy(&cpu_mass, tot_mass, sizeof(float), cudaMemcpyDeviceToHost));
  
  cudaFree(tot_mass);
  mass[0] = (double)(cpu_mass);
}

/** setup and call kernel to calculate the total momentum of the whole fluid
 *  @param host_mom value of the momentum calcutated on the GPU
 */
void lb_calc_fluid_momentum_GPU(double* host_mom){

  float* tot_momentum;
  float host_momentum[3] = { 0.0f, 0.0f, 0.0f};
  cuda_safe_mem(cudaMalloc((void**)&tot_momentum, 3*sizeof(float)));
  cuda_safe_mem(cudaMemcpy(tot_momentum, host_momentum, 3*sizeof(float), cudaMemcpyHostToDevice));

  /** values for the kernel call */
  int threads_per_block = 64;
  int blocks_per_grid_y = 4;
  int blocks_per_grid_x = (lbpar_gpu.number_of_nodes + threads_per_block * blocks_per_grid_y - 1) /(threads_per_block * blocks_per_grid_y);
  dim3 dim_grid = make_uint3(blocks_per_grid_x, blocks_per_grid_y, 1);

  KERNELCALL(momentum, dim_grid, threads_per_block,(*current_nodes, device_rho_v, node_f, tot_momentum));
  
  cuda_safe_mem(cudaMemcpy(host_momentum, tot_momentum, 3*sizeof(float), cudaMemcpyDeviceToHost));
  
  cudaFree(tot_momentum);
  host_mom[0] = (double)(host_momentum[0]* lbpar_gpu.agrid/lbpar_gpu.tau);
  host_mom[1] = (double)(host_momentum[1]* lbpar_gpu.agrid/lbpar_gpu.tau);
  host_mom[2] = (double)(host_momentum[2]* lbpar_gpu.agrid/lbpar_gpu.tau);
}

/** setup and call kernel to remove the net momentum of the whole fluid
 */
void lb_remove_fluid_momentum_GPU(void){
  float* tot_momentum;
  float host_momentum[3] = { 0.0f, 0.0f, 0.0f};
  cuda_safe_mem(cudaMalloc((void**)&tot_momentum, 3*sizeof(float)));
  cuda_safe_mem(cudaMemcpy(tot_momentum, host_momentum, 3*sizeof(float), cudaMemcpyHostToDevice));

  /** values for the kernel call */
  int threads_per_block = 64;
  int blocks_per_grid_y = 4;
  int blocks_per_grid_x = (lbpar_gpu.number_of_nodes + threads_per_block * blocks_per_grid_y - 1) /(threads_per_block * blocks_per_grid_y);
  dim3 dim_grid = make_uint3(blocks_per_grid_x, blocks_per_grid_y, 1);

  KERNELCALL(momentum, dim_grid, threads_per_block,(*current_nodes, device_rho_v, node_f, tot_momentum));
  
  cuda_safe_mem(cudaMemcpy(host_momentum, tot_momentum, 3*sizeof(float), cudaMemcpyDeviceToHost));

  KERNELCALL(remove_momentum, dim_grid, threads_per_block,(*current_nodes, device_rho_v, node_f, tot_momentum));
  
  cudaFree(tot_momentum);
}


/** setup and call kernel to calculate the temperature of the hole fluid
 *  @param host_temp value of the temperatur calcutated on the GPU
*/
void lb_calc_fluid_temperature_GPU(double* host_temp){

  int host_number_of_non_boundary_nodes = 0;
  int *device_number_of_non_boundary_nodes;
  cuda_safe_mem(cudaMalloc((void**)&device_number_of_non_boundary_nodes, sizeof(int)));
  cuda_safe_mem(cudaMemcpy(device_number_of_non_boundary_nodes, &host_number_of_non_boundary_nodes, sizeof(int), cudaMemcpyHostToDevice));

  float host_jsquared = 0.0f;
  float* device_jsquared;
  cuda_safe_mem(cudaMalloc((void**)&device_jsquared, sizeof(float)));
  cuda_safe_mem(cudaMemcpy(device_jsquared, &host_jsquared, sizeof(float), cudaMemcpyHostToDevice));

  /** values for the kernel call */
  int threads_per_block = 64;
  int blocks_per_grid_y = 4;
  int blocks_per_grid_x = (lbpar_gpu.number_of_nodes + threads_per_block * blocks_per_grid_y - 1) /(threads_per_block * blocks_per_grid_y);
  dim3 dim_grid = make_uint3(blocks_per_grid_x, blocks_per_grid_y, 1);

  KERNELCALL(temperature, dim_grid, threads_per_block,(*current_nodes, device_jsquared, device_number_of_non_boundary_nodes));

  cuda_safe_mem(cudaMemcpy(&host_number_of_non_boundary_nodes, device_number_of_non_boundary_nodes, sizeof(int), cudaMemcpyDeviceToHost));
  cuda_safe_mem(cudaMemcpy(&host_jsquared, device_jsquared, sizeof(float), cudaMemcpyDeviceToHost));

  // TODO: check that temperature calculation is properly implemented for shanchen
  *host_temp=0;

  #pragma unroll
  for(int ii=0;ii<LB_COMPONENTS;++ii)
  { 
      *host_temp += (double)(host_jsquared*1./(3.0f*lbpar_gpu.rho[ii]*host_number_of_non_boundary_nodes*lbpar_gpu.tau*lbpar_gpu.tau*lbpar_gpu.agrid));
  }
}


#ifdef SHANCHEN
void lb_calc_shanchen_GPU(){
  /** values for the kernel call */
  int threads_per_block = 64;
  int blocks_per_grid_y = 4;
  int blocks_per_grid_x = (lbpar_gpu.number_of_nodes + threads_per_block * blocks_per_grid_y - 1) /(threads_per_block * blocks_per_grid_y);
  dim3 dim_grid = make_uint3(blocks_per_grid_x, blocks_per_grid_y, 1);

#ifdef LB_BOUNDARIES_GPU
  if (n_lb_boundaries != 0)
  {
    KERNELCALL(lb_shanchen_set_boundaries, dim_grid, threads_per_block,(*current_nodes));
    cudaThreadSynchronize();
  }
#endif
  KERNELCALL(lb_shanchen_GPU, dim_grid, threads_per_block,(*current_nodes, node_f));
}

#endif // SHANCHEN

/** setup and call kernel for getting macroscopic fluid values of all nodes
 * @param *host_checkpoint_vd struct to save the gpu populations
 * @param *host_checkpoint_seed struct to save the nodes' seeds for the lb on the gpu
 * @param *host_checkpoint_boundary struct to save the boundary nodes
 * @param *host_checkpoint_force struct to save the forces on the nodes
 */
void lb_save_checkpoint_GPU(float *host_checkpoint_vd, unsigned int *host_checkpoint_seed, int *host_checkpoint_boundary, lbForceFloat *host_checkpoint_force){

  cuda_safe_mem(cudaMemcpy(host_checkpoint_vd, current_nodes->vd, lbpar_gpu.number_of_nodes * 19 * sizeof(float), cudaMemcpyDeviceToHost));
  cuda_safe_mem(cudaMemcpy(host_checkpoint_seed, current_nodes->seed, lbpar_gpu.number_of_nodes * sizeof(unsigned int), cudaMemcpyDeviceToHost));
  cuda_safe_mem(cudaMemcpy(host_checkpoint_boundary, current_nodes->boundary, lbpar_gpu.number_of_nodes * sizeof(int), cudaMemcpyDeviceToHost));
  cuda_safe_mem(cudaMemcpy(host_checkpoint_force, node_f.force, lbpar_gpu.number_of_nodes * 3 * sizeof(lbForceFloat), cudaMemcpyDeviceToHost));

}

/** setup and call kernel for setting macroscopic fluid values of all nodes
 * @param *host_checkpoint_vd struct to save the gpu populations
 * @param *host_checkpoint_seed struct to save the nodes' seeds for the lb on the gpu
 * @param *host_checkpoint_boundary struct to save the boundary nodes
 * @param *host_checkpoint_force struct to save the forces on the nodes
*/
void lb_load_checkpoint_GPU(float *host_checkpoint_vd, unsigned int *host_checkpoint_seed, int *host_checkpoint_boundary, lbForceFloat *host_checkpoint_force){

  current_nodes = &nodes_a;
  intflag = 1;

  cuda_safe_mem(cudaMemcpy(current_nodes->vd, host_checkpoint_vd, lbpar_gpu.number_of_nodes * 19 * sizeof(float), cudaMemcpyHostToDevice));

  cuda_safe_mem(cudaMemcpy(current_nodes->seed, host_checkpoint_seed, lbpar_gpu.number_of_nodes * sizeof(unsigned int), cudaMemcpyHostToDevice));
  cuda_safe_mem(cudaMemcpy(current_nodes->boundary, host_checkpoint_boundary, lbpar_gpu.number_of_nodes * sizeof(int), cudaMemcpyHostToDevice));
  cuda_safe_mem(cudaMemcpy(node_f.force, host_checkpoint_force, lbpar_gpu.number_of_nodes * 3 * sizeof(lbForceFloat), cudaMemcpyHostToDevice));
}

/** setup and call kernel to get the boundary flag of a single node
 *  @param single_nodeindex number of the node to get the flag for
 *  @param host_flag her goes the value of the boundary flag
 */
void lb_get_boundary_flag_GPU(int single_nodeindex, int* host_flag){
   
  int* device_flag;
  cuda_safe_mem(cudaMalloc((void**)&device_flag, sizeof(int)));
  int threads_per_block_flag = 1;
  int blocks_per_grid_flag_y = 1;
  int blocks_per_grid_flag_x = 1;
  dim3 dim_grid_flag = make_uint3(blocks_per_grid_flag_x, blocks_per_grid_flag_y, 1);

  KERNELCALL(lb_get_boundary_flag, dim_grid_flag, threads_per_block_flag, (single_nodeindex, device_flag, *current_nodes));

  cuda_safe_mem(cudaMemcpy(host_flag, device_flag, sizeof(int), cudaMemcpyDeviceToHost));

  cudaFree(device_flag);
}

/** set the density at a single node
 *  @param single_nodeindex the node to set the velocity for 
 *  @param *host_rho the density to set
 */
void lb_set_node_rho_GPU(int single_nodeindex, float* host_rho){
   
  float* device_rho;
  cuda_safe_mem(cudaMalloc((void**)&device_rho, LB_COMPONENTS*sizeof(float)));
  cuda_safe_mem(cudaMemcpy(device_rho, host_rho, LB_COMPONENTS*sizeof(float), cudaMemcpyHostToDevice));
  int threads_per_block_flag = 1;
  int blocks_per_grid_flag_y = 1;
  int blocks_per_grid_flag_x = 1;
  dim3 dim_grid_flag = make_uint3(blocks_per_grid_flag_x, blocks_per_grid_flag_y, 1);
  KERNELCALL(set_rho, dim_grid_flag, threads_per_block_flag, (*current_nodes, device_rho_v, single_nodeindex, device_rho)); 
  cudaFree(device_rho);
}

/** set the net velocity at a single node
 *  @param single_nodeindex the node to set the velocity for 
 *  @param host_velocity the velocity to set
 */
void lb_set_node_velocity_GPU(int single_nodeindex, float* host_velocity){
   
  float* device_velocity;
  cuda_safe_mem(cudaMalloc((void**)&device_velocity, 3*sizeof(float)));
  cuda_safe_mem(cudaMemcpy(device_velocity, host_velocity, 3*sizeof(float), cudaMemcpyHostToDevice));
  int threads_per_block_flag = 1;
  int blocks_per_grid_flag_y = 1;
  int blocks_per_grid_flag_x = 1;
  dim3 dim_grid_flag = make_uint3(blocks_per_grid_flag_x, blocks_per_grid_flag_y, 1);

  KERNELCALL(set_u_from_rho_v_pi, dim_grid_flag, threads_per_block_flag, (*current_nodes, single_nodeindex, device_velocity, device_rho_v, node_f));

  cudaFree(device_velocity);
}

/** reinit of params 
 * @param *lbpar_gpu struct containing the paramters of the fluid
*/
void reinit_parameters_GPU(LB_parameters_gpu *lbpar_gpu){
  /**write parameters in const memory*/
  cuda_safe_mem(cudaMemcpyToSymbol(para, lbpar_gpu, sizeof(LB_parameters_gpu)));
}


/** function to determine parameters for rotational integration rewritten for floats.
    See \ref rotation.cpp for original implementation for doubles
    no ifndef rotation necessary since these values are float instead of double */
    
void define_Qdd(Particle *p, float Qd[4], float Qdd[4], float S[3], float Wd[3])
{
  float S1;
  
    /* calculate the first derivative of the quaternion */
  Qd[0] = 0.5 * ( -p->r.quat[1] * p->m.somega[0] -
                   p->r.quat[2] * p->m.somega[1] -
		   p->r.quat[3] * p->m.somega[2] );

  Qd[1] = 0.5 * (  p->r.quat[0] * p->m.somega[0] -
                   p->r.quat[3] * p->m.somega[1] +
		   p->r.quat[2] * p->m.somega[2] );

  Qd[2] = 0.5 * (  p->r.quat[3] * p->m.somega[0] +
                   p->r.quat[0] * p->m.somega[1] -
		   p->r.quat[1] * p->m.somega[2] );

  Qd[3] = 0.5 * ( -p->r.quat[2] * p->m.somega[0] +
                   p->r.quat[1] * p->m.somega[1] +
		   p->r.quat[0] * p->m.somega[2] );
  
  /* calculate the second derivative of the quaternion */
  Wd[0] =  (p->f.storque[0] + p->m.somega[1]*p->m.somega[2]*(p->p.rinertia[1]-p->p.rinertia[2]))/p->p.rinertia[0];
  Wd[1] =  (p->f.storque[1] + p->m.somega[2]*p->m.somega[0]*(p->p.rinertia[2]-p->p.rinertia[0]))/p->p.rinertia[1];
  Wd[2] =  (p->f.storque[2] + p->m.somega[0]*p->m.somega[1]*(p->p.rinertia[0]-p->p.rinertia[1]))/p->p.rinertia[2];
  
  S1 = Qd[0]*Qd[0] + Qd[1]*Qd[1] + Qd[2]*Qd[2] + Qd[3]*Qd[3];

  Qdd[0] = 0.5 * ( -p->r.quat[1] * Wd[0] -
		    p->r.quat[2] * Wd[1] -
	  	    p->r.quat[3] * Wd[2] ) - p->r.quat[0] * S1;

  Qdd[1] = 0.5 * (  p->r.quat[0] * Wd[0] -
		    p->r.quat[3] * Wd[1] +
		    p->r.quat[2] * Wd[2] ) - p->r.quat[1] * S1;
  
  Qdd[2] = 0.5 * (  p->r.quat[3] * Wd[0] +
		    p->r.quat[0] * Wd[1] -
		    p->r.quat[1] * Wd[2] ) - p->r.quat[2] * S1;
  
  Qdd[3] = 0.5 * ( -p->r.quat[2] * Wd[0] +
                    p->r.quat[1] * Wd[1] +
		    p->r.quat[0] * Wd[2] ) - p->r.quat[3] * S1;

  S[0] = S1;
  S[1] = Qd[0]*Qdd[0]  + Qd[1]*Qdd[1]  + Qd[2]*Qdd[2]  + Qd[3]*Qdd[3];
  S[2] = Qdd[0]*Qdd[0] + Qdd[1]*Qdd[1] + Qdd[2]*Qdd[2] + Qdd[3]*Qdd[3];
}


/** Here we use quaternions to calculate the rotation matrix which
    will be used then to transform torques from the laboratory to
    the body-fixed frames */ 
void lb_define_rotation_matrix(Particle *p, double A[9])
{
  double q0q0 = p->r.quat[0];
  q0q0 *= q0q0;
  
  double q1q1 = p->r.quat[1];
  q1q1 *= q1q1;
  
  double q2q2 = p->r.quat[2];
  q2q2 *= q2q2;
  
  double q3q3 = p->r.quat[3];
  q3q3 *= q3q3;
  
  A[0 + 3*0] = q0q0 + q1q1 - q2q2 - q3q3;
  A[1 + 3*1] = q0q0 - q1q1 + q2q2 - q3q3;
  A[2 + 3*2] = q0q0 - q1q1 - q2q2 + q3q3;

  A[0 + 3*1] = 2*(p->r.quat[1]*p->r.quat[2] + p->r.quat[0]*p->r.quat[3]);
  A[0 + 3*2] = 2*(p->r.quat[1]*p->r.quat[3] - p->r.quat[0]*p->r.quat[2]);
  A[1 + 3*0] = 2*(p->r.quat[1]*p->r.quat[2] - p->r.quat[0]*p->r.quat[3]);

  A[1 + 3*2] = 2*(p->r.quat[2]*p->r.quat[3] + p->r.quat[0]*p->r.quat[1]);
  A[2 + 3*0] = 2*(p->r.quat[1]*p->r.quat[3] + p->r.quat[0]*p->r.quat[2]);
  A[2 + 3*1] = 2*(p->r.quat[2]*p->r.quat[3] - p->r.quat[0]*p->r.quat[1]);
}
void convert_omega_space_to_body(LB_moving_boundary *lbb, Particle *p)//used for init
{
  double A[9];
  lb_define_rotation_matrix(p, A);
  
  
  p->m.somega[0] = A[0 + 3*0]*lbb->omega[0] + A[0 + 3*1]*lbb->omega[1] + A[0 + 3*2]*lbb->omega[2]; 
  p->m.somega[1] = A[1 + 3*0]*lbb->omega[0] + A[1 + 3*1]*lbb->omega[1] + A[1 + 3*2]*lbb->omega[2]; 
  p->m.somega[2] = A[2 + 3*0]*lbb->omega[0] + A[2 + 3*1]*lbb->omega[1] + A[2 + 3*2]*lbb->omega[2]; 
}
void convert_omega_anchors_body_to_space(Particle *p, float *omega, float *h_anchors)
{
  double A[9];
  lb_define_rotation_matrix(p, A);
  omega[0] = A[0 + 3*0]*p->m.somega[0] + A[1 + 3*0]*p->m.somega[1] + A[2 + 3*0]*p->m.somega[2];
  omega[1] = A[0 + 3*1]*p->m.somega[0] + A[1 + 3*1]*p->m.somega[1] + A[2 + 3*1]*p->m.somega[2];
  omega[2] = A[0 + 3*2]*p->m.somega[0] + A[1 + 3*2]*p->m.somega[1] + A[2 + 3*2]*p->m.somega[2];
  int ii;
  for(ii=0; ii<p->p.n_anchors; ii++){
    h_anchors[4*ii + 0] = A[0 + 3*0]*p->p.anchors[4*ii + 0] + A[1 + 3*0]*p->p.anchors[4*ii + 1] + A[2 + 3*0]*p->p.anchors[4*ii + 2];
    h_anchors[4*ii + 1] = A[0 + 3*1]*p->p.anchors[4*ii + 0] + A[1 + 3*1]*p->p.anchors[4*ii + 1] + A[2 + 3*1]*p->p.anchors[4*ii + 2];
    h_anchors[4*ii + 2] = A[0 + 3*2]*p->p.anchors[4*ii + 0] + A[1 + 3*2]*p->p.anchors[4*ii + 1] + A[2 + 3*2]*p->p.anchors[4*ii + 2];
    h_anchors[4*ii + 3] = p->p.anchors[4*ii + 3];
  }
}

#ifdef EXTERNAL_FORCES
void convert_force_body_to_space(Particle *p)
{
  double A[9];
  lb_define_rotation_matrix(p, A);
  p->p.lab_force[0] = A[0 + 3*0]*p->p.body_force[0] + A[1 + 3*0]*p->p.body_force[1] + A[2 + 3*0]*p->p.body_force[2];
  p->p.lab_force[1] = A[0 + 3*1]*p->p.body_force[0] + A[1 + 3*1]*p->p.body_force[1] + A[2 + 3*1]*p->p.body_force[2];
  p->p.lab_force[2] = A[0 + 3*2]*p->p.body_force[0] + A[1 + 3*2]*p->p.body_force[1] + A[2 + 3*2]*p->p.body_force[2];
}
#endif

inline void lb_convert_quat_to_quatu(double quat[4], double quatu[3])
{
  /* director */
  quatu[0] = 2*(quat[1]*quat[3] + quat[0]*quat[2]);
  quatu[1] = 2*(quat[2]*quat[3] - quat[0]*quat[1]);
  quatu[2] =   (quat[0]*quat[0] - quat[1]*quat[1] -
		quat[2]*quat[2] + quat[3]*quat[3]); 
}

void lb_propagate_omega_quat_particle(Particle* p)
{
  double lambda;

  float Qd[4], Qdd[4], S[3], Wd[3];
 
  define_Qdd(p, Qd, Qdd, S, Wd);

  lambda = 1 - S[0]*tau_sq / 2.0 - sqrt(1 - tau_sq*(S[0] + lbpar_gpu.tau*(S[1] + tau_half/2.0*(S[2]-S[0]*S[0]))));
  
  for(int j=0; j < 3; j++){
    p->m.somega[j]+= tau_half *Wd[j];
  }
  ONEPART_TRACE(if(p->p.identity==check_id) fprintf(stderr,"%d: OPT: PV_1 v_new = (%.3e,%.3e,%.3e)\n",this_node,p->m.v[0],p->m.v[1],p->m.v[2]));
  
  p->r.quat[0]+= lbpar_gpu.tau*(Qd[0] + tau_half*Qdd[0]) - lambda*p->r.quat[0];
  p->r.quat[1]+= lbpar_gpu.tau*(Qd[1] + tau_half*Qdd[1]) - lambda*p->r.quat[1];
  p->r.quat[2]+= lbpar_gpu.tau*(Qd[2] + tau_half*Qdd[2]) - lambda*p->r.quat[2];
  p->r.quat[3]+= lbpar_gpu.tau*(Qd[3] + tau_half*Qdd[3]) - lambda*p->r.quat[3];
  // Update the director
  lb_convert_quat_to_quatu(p->r.quat, p->r.quatu);
}

void lb_convert_torques_propagate_omega(Particle *p)
{


  double tx, ty, tz;
  double A[9];
  lb_define_rotation_matrix(p, A);

  tx = A[0 + 3*0]*p->f.storque[0] + A[0 + 3*1]*p->f.storque[1] + A[0 + 3*2]*p->f.storque[2];
  ty = A[1 + 3*0]*p->f.storque[0] + A[1 + 3*1]*p->f.storque[1] + A[1 + 3*2]*p->f.storque[2];
  tz = A[2 + 3*0]*p->f.storque[0] + A[2 + 3*1]*p->f.storque[1] + A[2 + 3*2]*p->f.storque[2];
#ifdef EXTERNAL_FORCES
  p->f.storque[0] = tx + p->p.body_torque[0];
  p->f.storque[1] = ty + p->p.body_torque[1];
  p->f.storque[2] = tz + p->p.body_torque[2];
#else
  p->f.storque[0] = tx;
  p->f.storque[1] = ty;
  p->f.storque[2] = tz;
#endif


  p->m.somega[0]+= tau_half*p->f.storque[0]/p->p.rinertia[0];
  p->m.somega[1]+= tau_half*p->f.storque[1]/p->p.rinertia[1];
  p->m.somega[2]+= tau_half*p->f.storque[2]/p->p.rinertia[2];

      /* if the tensor of inertia is isotropic, the following refinement is not needed.
         Otherwise repeat this loop 2-3 times depending on the required accuracy */
  for(int times=0; times <= 3; times++)
  { 
    double Wd[3];


    Wd[0] = (p->m.somega[1]*p->m.somega[2]*(p->p.rinertia[1]-p->p.rinertia[2]))/p->p.rinertia[0];
    Wd[1] = (p->m.somega[2]*p->m.somega[0]*(p->p.rinertia[2]-p->p.rinertia[0]))/p->p.rinertia[1]; 
    Wd[2] = (p->m.somega[0]*p->m.somega[1]*(p->p.rinertia[0]-p->p.rinertia[1]))/p->p.rinertia[2];


    p->m.somega[0]+= tau_half*Wd[0];
    p->m.somega[1]+= tau_half*Wd[1];
    p->m.somega[2]+= tau_half*Wd[2];
  }

}

void lb_convert_torque_to_space_for_print(Particle *p, double *print_val)
{
  double A[9];
  lb_define_rotation_matrix(p, A);
  
  print_val[0] = A[0 + 3*0]*p->f.storque[0] + A[1 + 3*0]*p->f.storque[1] + A[2 + 3*0]*p->f.storque[2];
  print_val[1] = A[0 + 3*1]*p->f.storque[0] + A[1 + 3*1]*p->f.storque[1] + A[2 + 3*1]*p->f.storque[2];
  print_val[2] = A[0 + 3*2]*p->f.storque[0] + A[1 + 3*2]*p->f.storque[1] + A[2 + 3*2]*p->f.storque[2];
}

void lb_convert_omega_to_space_for_print(Particle *p, double *print_val)
{
  double A[9];
  lb_define_rotation_matrix(p, A);
  
  print_val[0] = A[0 + 3*0]*p->m.somega[0] + A[1 + 3*0]*p->m.somega[1] + A[2 + 3*0]*p->m.somega[2];
  print_val[1] = A[0 + 3*1]*p->m.somega[0] + A[1 + 3*1]*p->m.somega[1] + A[2 + 3*1]*p->m.somega[2];
  print_val[2] = A[0 + 3*2]*p->m.somega[0] + A[1 + 3*2]*p->m.somega[1] + A[2 + 3*2]*p->m.somega[2];
}

void lb_convert_force_to_body_for_print(Particle *p, double *print_val)
{
  double A[9];
  lb_define_rotation_matrix(p,A);
  
  print_val[0] = A[0 + 3*0]*p->f.sf[0] + A[0 + 3*1]*p->f.sf[1] + A[0 + 3*2]*p->f.sf[2];
  print_val[1] = A[1 + 3*0]*p->f.sf[0] + A[1 + 3*1]*p->f.sf[1] + A[1 + 3*2]*p->f.sf[2];
  print_val[2] = A[2 + 3*0]*p->f.sf[0] + A[2 + 3*1]*p->f.sf[1] + A[2 + 3*2]*p->f.sf[2];
}


//this should run in parallel to integrate()
//this needs external_forces feature
//this cannot (and should not) do torques&rotation!
#ifdef EXTERNAL_FORCES
void calc_MD_force_and_set_ia_vel(){
  // TODO this needs to be extended for ia forces between particles
  int ii;
  float cast_v[3];
  for(ii = 0; ii<host_n_lb_moving_boundaries; ii++){
    
    //FORCECALC()
    convert_force_body_to_space(host_lb_moving_boundary + ii);
    
    host_lb_moving_boundary[ii].p.lab_force[0] += host_lb_moving_boundary[ii].p.ext_force[0];
    host_lb_moving_boundary[ii].m.v[0] += host_lb_moving_boundary[ii].p.lab_force[0];
    
    host_lb_moving_boundary[ii].p.lab_force[1] += host_lb_moving_boundary[ii].p.ext_force[1];
    host_lb_moving_boundary[ii].m.v[1] += host_lb_moving_boundary[ii].p.lab_force[1];
    
    host_lb_moving_boundary[ii].p.lab_force[2] += host_lb_moving_boundary[ii].p.ext_force[2];
    host_lb_moving_boundary[ii].m.v[2] += host_lb_moving_boundary[ii].p.lab_force[2];
    
    cast_v[0] = host_lb_moving_boundary[ii].m.v[0];
    cast_v[1] = host_lb_moving_boundary[ii].m.v[1];
    cast_v[2] = host_lb_moving_boundary[ii].m.v[2];
    
    cudaMemcpy(lb_moving_boundary[ii].velocity, cast_v, 3*sizeof(float), cudaMemcpyHostToDevice);
    
    
  }
}
#endif
 /* Integrator for the boundaries equations of motion */
void lb_integrate_moving_boundaries(){
  float fold_p[3];
  float cast_v[3];
  float omega[3];
  float force_add[3];
  int ii;
  float *migrating_anchors = h_lab_anchors;
  for(ii = 0; ii<host_n_lb_moving_boundaries; ii++){

    lb_propagate_omega_quat_particle(host_lb_moving_boundary + ii); //use old torque in body space (vv)
  //no safe mem, since there are no mallocs and initialisations happening. also these communnications happen a lot so it needs to be efficient.
    cudaMemcpy(host_lb_moving_boundary[ii].f.sf, lb_moving_boundary[ii].force_hyd, 3*sizeof(float), cudaMemcpyDeviceToHost);
    cudaMemset(lb_moving_boundary[ii].force_hyd, 0, 3*sizeof(float));
    cudaMemcpy(force_add, lb_moving_boundary[ii].force_add, 3*sizeof(float), cudaMemcpyDeviceToHost);
    cudaMemset(lb_moving_boundary[ii].force_add, 0, 3*sizeof(float));
    
    
    cudaMemcpy(host_lb_moving_boundary[ii].f.storque, lb_moving_boundary[ii].torque, 3*sizeof(float), cudaMemcpyDeviceToHost);
    cudaMemset(lb_moving_boundary[ii].torque, 0, 3*sizeof(float));
    //getting forces/torques
      host_lb_moving_boundary[ii].f.storque[0] /= tau_sq;//transform from kg a^2/tau^2 to kg a^2/s^2
      #ifdef EXTERNAL_FORCES
      host_lb_moving_boundary[ii].f.storque[0] += host_lb_moving_boundary[ii].p.ext_torque[0];//add ext_torques labspace.
      #endif
      host_lb_moving_boundary[ii].f.storque[1] /= tau_sq;
      #ifdef EXTERNAL_FORCES
      host_lb_moving_boundary[ii].f.storque[1] += host_lb_moving_boundary[ii].p.ext_torque[1];
      #endif
      host_lb_moving_boundary[ii].f.storque[2] /= tau_sq;
      #ifdef EXTERNAL_FORCES
      host_lb_moving_boundary[ii].f.storque[2] += host_lb_moving_boundary[ii].p.ext_torque[2];
      #endif
      

    lb_convert_torques_propagate_omega(host_lb_moving_boundary + ii); //calc new torques and do the rest of the rotation
    
    convert_omega_anchors_body_to_space(host_lb_moving_boundary + ii, omega, migrating_anchors);
    migrating_anchors += 4*host_lb_moving_boundary[ii].p.n_anchors;
      omega[0] *= lbpar_gpu.tau;//back to lb time units
      omega[1] *= lbpar_gpu.tau;
      omega[2] *= lbpar_gpu.tau;
    
    
    #ifdef LB_GPU_MOVING_MOMENTUM //when using momentum conservation
    host_lb_moving_boundary[ii].m.v[0] += (force_add[0] + host_lb_moving_boundary[ii].f.sf[0])/host_lb_moving_boundary[ii].p.mass; 
    #else
    host_lb_moving_boundary[ii].m.v[0] += host_lb_moving_boundary[ii].f.sf[0]/host_lb_moving_boundary[ii].p.mass;
    #endif
    #ifdef EXTERNAL_FORCES
    host_lb_moving_boundary[ii].m.v[0] += host_lb_moving_boundary[ii].p.lab_force[0];
    #endif
    cast_v[0] = host_lb_moving_boundary[ii].m.v[0];
    host_lb_moving_boundary[ii].r.p[0] += host_lb_moving_boundary[ii].m.v[0];//do the propagation w/ absolute coords
    fold_p[0] = host_lb_moving_boundary[ii].r.p[0] - lbpar_gpu.dim_x * floor(host_lb_moving_boundary[ii].r.p[0]/lbpar_gpu.dim_x);//fold absolute coords to avoid single-prec shenanigans
    
    
    #ifdef LB_GPU_MOVING_MOMENTUM
    host_lb_moving_boundary[ii].m.v[1] += (force_add[1] + host_lb_moving_boundary[ii].f.sf[1])/host_lb_moving_boundary[ii].p.mass; 
    #else
    host_lb_moving_boundary[ii].m.v[1] += host_lb_moving_boundary[ii].f.sf[1]/host_lb_moving_boundary[ii].p.mass;
    #endif
    #ifdef EXTERNAL_FORCES
    host_lb_moving_boundary[ii].m.v[1] += host_lb_moving_boundary[ii].p.lab_force[1];
    #endif
    cast_v[1] = host_lb_moving_boundary[ii].m.v[1];
    host_lb_moving_boundary[ii].r.p[1] += host_lb_moving_boundary[ii].m.v[1];
    fold_p[1] = host_lb_moving_boundary[ii].r.p[1] - lbpar_gpu.dim_y * floor(host_lb_moving_boundary[ii].r.p[1]/lbpar_gpu.dim_y);
    
    
    #ifdef LB_GPU_MOVING_MOMENTUM
    host_lb_moving_boundary[ii].m.v[2] += (force_add[2] + host_lb_moving_boundary[ii].f.sf[2])/host_lb_moving_boundary[ii].p.mass; 
    #else
    host_lb_moving_boundary[ii].m.v[2] += host_lb_moving_boundary[ii].f.sf[2]/host_lb_moving_boundary[ii].p.mass;
    #endif
    #ifdef EXTERNAL_FORCES
    host_lb_moving_boundary[ii].m.v[2] += host_lb_moving_boundary[ii].p.lab_force[2];
    #endif
    cast_v[2] = host_lb_moving_boundary[ii].m.v[2];
    host_lb_moving_boundary[ii].r.p[2] += host_lb_moving_boundary[ii].m.v[2];
    fold_p[2] = host_lb_moving_boundary[ii].r.p[2] - lbpar_gpu.dim_z * floor(host_lb_moving_boundary[ii].r.p[2]/lbpar_gpu.dim_z);


    
    cudaMemcpy(lb_moving_boundary[ii].omega, omega, 3*sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(lb_moving_boundary[ii].center, fold_p, 3*sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(lb_moving_boundary[ii].velocity, cast_v, 3*sizeof(float), cudaMemcpyHostToDevice);


   //DEBUG
   n_integration_steps++;
   if(n_integration_steps > 1){
               add_track += force_add[0];
	       hyd_track += host_lb_moving_boundary[ii].f.sf[0]/host_lb_moving_boundary[ii].p.mass;
	       torque_track[0] += host_lb_moving_boundary[ii].f.storque[0];
	       torque_track[1] += host_lb_moving_boundary[ii].f.storque[1];
	       torque_track[2] += host_lb_moving_boundary[ii].f.storque[2];
	       
   }
   //if ( !n_integration_steps % 5000)
       // printf("\n %d: boundary velocity %f %f %f; boundary omega(lab) %f %f %f,force_hyd %e, track %e ,force_add %f, track %f, diff %f torque_track %e %e %e\n",n_integration_steps, host_lb_moving_boundary[ii].m.v[0], host_lb_moving_boundary[ii].m.v[1], host_lb_moving_boundary[ii].m.v[2],omega[0],omega[1],omega[2], host_lb_moving_boundary[ii].f.sf[0], hyd_track, force_add[0], add_track, -hyd_track+add_track,torque_track[0] ,torque_track[1] ,torque_track[2] );
   }
  cudaMemcpy(d_lab_anchors, h_lab_anchors, (4* *migrating_anchors +1)*sizeof(float), cudaMemcpyHostToDevice);
}

/**integration kernel for the lb gpu fluid update called from host */
void lb_integrate_GPU() {

  /** values for the kernel call */
  int threads_per_block = 64;
  int blocks_per_grid_y = 4;
  int blocks_per_grid_x = (lbpar_gpu.number_of_nodes + threads_per_block * blocks_per_grid_y - 1) /(threads_per_block * blocks_per_grid_y);
  dim3 dim_grid = make_uint3(blocks_per_grid_x, blocks_per_grid_y, 1);
#ifdef LB_BOUNDARIES_GPU
  if (n_lb_boundaries > 0)
  {
    cuda_safe_mem(cudaMemset( lb_boundary_force, 0, 3*n_lb_boundaries*sizeof(float)));
  }
#endif

  /**call of fluid step*/
  /* NOTE: if pi is needed at every integration step, one should call an extended version 
           of the integrate kernel, or pass also device_rho_v_pi and make sure that either 
           it or device_rho_v are NULL depending on extended_values_flag */ 
  if (intflag == 1)
  {
    KERNELCALL(integrate, dim_grid, threads_per_block, (nodes_a, nodes_b, device_rho_v, node_f, lb_ek_parameters_gpu));
    current_nodes = &nodes_b;
    next_nodes = &nodes_a;
    intflag = 0;
  }
  else
  {
    KERNELCALL(integrate, dim_grid, threads_per_block, (nodes_b, nodes_a, device_rho_v, node_f, lb_ek_parameters_gpu));
    current_nodes = &nodes_a;
    next_nodes = &nodes_b;
    intflag = 1;
  }

#ifdef LB_BOUNDARIES_GPU
    if (n_lb_boundaries > 0)
    {
      #ifdef EXTERNAL_FORCES
      calc_MD_force_and_set_ia_vel();
      #endif
      KERNELCALL(apply_boundaries, dim_grid, threads_per_block, (*current_nodes, device_rho_v, lb_moving_boundary, lb_boundary_velocity, lb_boundary_force));
      lb_integrate_moving_boundaries();
      KERNELCALL(propagate_boundaries, dim_grid, threads_per_block, (*current_nodes, device_rho_v, lb_moving_boundary, d_lab_anchors));
      KERNELCALL(update_boundaries_from_buffer, dim_grid, threads_per_block,(*current_nodes, *next_nodes, lb_moving_boundary));
    }
#endif
}

void lb_gpu_get_boundary_forces(double* forces) {
#ifdef LB_BOUNDARIES_GPU
  float* temp = (float*) Utils::malloc(3*n_lb_boundaries*sizeof(float));
  cuda_safe_mem(cudaMemcpy(temp, lb_boundary_force, 3*n_lb_boundaries*sizeof(float), cudaMemcpyDeviceToHost));

  for (int i =0; i<3*n_lb_boundaries; i++)
  {
    forces[i]=(double)temp[i];
  }
  free(temp);
#endif
}

__device__ void get_interpolated_velocity(LB_nodes_gpu n_a, float* r, float* u, LB_node_force_gpu node_f, int asdf) {

  /** see ahlrichs + duenweg page 8227 equ (10) and (11) */
  float temp_delta[6];
  float delta[8];
  int my_left[3];
  int node_index[8];
  float mode[4];
  u[0]=u[1]=u[2]=0;
  #pragma unroll
  for(int i=0; i<3; ++i){
    float scaledpos = r[i]/para.agrid - 0.5f;
    my_left[i] = (int)(floorf(scaledpos));
    temp_delta[3+i] = scaledpos - my_left[i];
    temp_delta[i] = 1.f - temp_delta[3+i];
  }

  delta[0] = temp_delta[0] * temp_delta[1] * temp_delta[2];
  delta[1] = temp_delta[3] * temp_delta[1] * temp_delta[2];
  delta[2] = temp_delta[0] * temp_delta[4] * temp_delta[2];
  delta[3] = temp_delta[3] * temp_delta[4] * temp_delta[2];
  delta[4] = temp_delta[0] * temp_delta[1] * temp_delta[5];
  delta[5] = temp_delta[3] * temp_delta[1] * temp_delta[5];
  delta[6] = temp_delta[0] * temp_delta[4] * temp_delta[5];
  delta[7] = temp_delta[3] * temp_delta[4] * temp_delta[5];

  // modulo for negative numbers is strange at best, shift to make sure we are positive
  int x = my_left[0] + para.dim_x;
  int y = my_left[1] + para.dim_y;
  int z = my_left[2] + para.dim_z;

  node_index[0] = x%para.dim_x     + para.dim_x*(y%para.dim_y)     + para.dim_x*para.dim_y*(z%para.dim_z);
  node_index[1] = (x+1)%para.dim_x + para.dim_x*(y%para.dim_y)     + para.dim_x*para.dim_y*(z%para.dim_z);
  node_index[2] = x%para.dim_x     + para.dim_x*((y+1)%para.dim_y) + para.dim_x*para.dim_y*(z%para.dim_z);
  node_index[3] = (x+1)%para.dim_x + para.dim_x*((y+1)%para.dim_y) + para.dim_x*para.dim_y*(z%para.dim_z);
  node_index[4] = x%para.dim_x     + para.dim_x*(y%para.dim_y)     + para.dim_x*para.dim_y*((z+1)%para.dim_z);
  node_index[5] = (x+1)%para.dim_x + para.dim_x*(y%para.dim_y)     + para.dim_x*para.dim_y*((z+1)%para.dim_z);
  node_index[6] = x%para.dim_x     + para.dim_x*((y+1)%para.dim_y) + para.dim_x*para.dim_y*((z+1)%para.dim_z);
  node_index[7] = (x+1)%para.dim_x + para.dim_x*((y+1)%para.dim_y) + para.dim_x*para.dim_y*((z+1)%para.dim_z);

  for(int i=0; i<8; ++i)
  {
      float totmass=0.0f;

      if(n_a.boundary[node_index[i]])
        continue;


      calc_m_from_n(n_a,node_index[i],mode);

      #pragma unroll
      for(int ii=0;ii<LB_COMPONENTS;ii++)
      {
        totmass+=mode[0]+para.rho[ii]*para.agrid*para.agrid*para.agrid;
      } 

#ifndef SHANCHEN
      u[0] += (mode[1]/totmass)*delta[i];
      u[1] += (mode[2]/totmass)*delta[i];
      u[2] += (mode[3]/totmass)*delta[i];
#else //SHANCHEN
//      u[0] += d_v[node_index[i]].v[0]/8.0f;  
//      u[1] += d_v[node_index[i]].v[1]/8.0f;
//      u[2] += d_v[node_index[i]].v[2]/8.0f;
#warning "lb_radial_velocity_profile does not work with SHANCHEN yet/"
        u[0] = 0;
        u[1] = 0;
        u[2] = 0;
#endif

//      mode[1]+=0.5f*node_f.force[0*para.number_of_nodes + node_index[i]];
//      mode[2]+=0.5f*node_f.force[1*para.number_of_nodes + node_index[i]];
//      mode[3]+=0.5f*node_f.force[2*para.number_of_nodes + node_index[i]];
//
    }

  #pragma unroll
  for(int i=0; i<3; ++i){
    u[i]*=para.agrid/para.tau;
  }
}
__global__ void fill_lb_radial_velocity_profile(LB_nodes_gpu n_a, radial_profile_data* pdata, float* data, LB_node_force_gpu node_f){

  unsigned int rbin=threadIdx.x;
  unsigned int phibin=blockIdx.x;
  unsigned int zbin=blockIdx.y;

  float roffset=pdata->minr;
  float r_incr=(pdata->maxr-pdata->minr)/(pdata->rbins-1);

  float r = roffset + rbin*r_incr;

  unsigned int maxj;
  float phioffset, phi_incr;
  if ( pdata->phibins == 1 ) {
    maxj = (int)floorf( 2*3.1415f*pdata->maxr/para.agrid ) ; 
    phioffset=0;
    phi_incr=2*3.1415f/maxj;
  } else {
    maxj = pdata->phibins;
    phioffset=pdata->minphi;
    phi_incr=(pdata->maxphi-pdata->minphi)/(pdata->phibins);
  }
  float phi = phioffset + phibin*phi_incr;

  unsigned int maxk;
  float zoffset, z_incr;
  if ( pdata->zbins == 1 ) {
    maxk = (int) para.dim_z;
    zoffset=-pdata->center[2];
    z_incr=para.agrid;
  } else {
    maxk = (int) pdata->zbins;
    zoffset=pdata->minz;
    z_incr=(pdata->maxz-pdata->minz)/(pdata->zbins-1);
  }

  float z = zoffset + zbin*z_incr;

  float p[3];
  p[0]=r*cos(phi)+pdata->center[0];
  p[1]=r*sin(phi)+pdata->center[1];
  p[2]=z+pdata->center[2];

  float v[3];
  get_interpolated_velocity(n_a, p, v, node_f, 0);
  unsigned int linear_index = rbin*maxj*maxk + phibin*maxk + zbin;

 float v_r,v_phi;

  if (r==0) {
    v_r=0;
    v_phi=0;
  } else {
    v_r = 1/r*((p[0]-pdata->center[0])*v[0] + (p[1]-pdata->center[1])*v[1]); 
    v_phi = 1/r/r*((p[0]-pdata->center[0])*v[1]-(p[1]-pdata->center[1])*v[0]);
  }
  data[3*linear_index+0]=v_r;
  data[3*linear_index+1]=v_phi;
  data[3*linear_index+2]=v[2];

}


__global__ void fill_lb_velocity_profile(LB_nodes_gpu n_a, profile_data* pdata, float* data, LB_node_force_gpu node_f){

  unsigned int xbin=threadIdx.x;
  unsigned int ybin=blockIdx.x;
  unsigned int zbin=blockIdx.y;

  float xoffset, yoffset, zoffset;
  float x_incr, y_incr, z_incr;
  unsigned int maxj, maxk;



  if ( pdata->xbins == 1 ) {
    /* maxi = (int) floor(gridDim.x/para.agrid); */
    xoffset=0;
    x_incr=para.agrid;
  } else {
    /* maxi = pdata->xbins; */
    xoffset=pdata->minx;
    x_incr=(pdata->maxx-pdata->minx)/(pdata->xbins-1);
  }
  float x = xoffset + xbin*x_incr;
  if ( pdata->ybins == 1 ) {
    maxj = (int) floorf(para.dim_y/para.agrid);
    yoffset=0;
    y_incr=para.agrid;
  } else {
    maxj = pdata->ybins;
    yoffset=pdata->miny;
    y_incr=(pdata->maxy-pdata->miny)/(pdata->ybins-1);
  }
  float y = yoffset + ybin*y_incr;
  if ( pdata->zbins == 1 ) {
    maxk = (int) floorf(para.dim_z/para.agrid);
    zoffset=0;
    z_incr=para.agrid;
  } else {
    maxk = (int) pdata->zbins;
    zoffset=pdata->minz;
    z_incr=(pdata->maxz-pdata->minz)/(pdata->zbins-1);
  }
  float z = zoffset + zbin*z_incr;

  float p[3];
  p[0]=x;
  p[1]=y;
  p[2]=z;

  float v[3];
  get_interpolated_velocity(n_a, p, v, node_f, 0);
  unsigned int linear_index = xbin*maxj*maxk + ybin*maxk + zbin;

  data[3*linear_index+0]=v[0];
  data[3*linear_index+1]=v[1];
  data[3*linear_index+2]=v[2];

}

int statistics_observable_lbgpu_radial_velocity_profile(radial_profile_data* pdata, double* A, unsigned int n_A){

  unsigned int maxj, maxk;
  float normalization_factor=1;
  
  if ( pdata->rbins == 1 ) {
    return 1;
  }

  unsigned int maxi=pdata->rbins;
  
  if ( pdata->phibins == 1 ) {
    maxj = (int)floorf( 2*3.1415f*pdata->maxr/lbpar_gpu.agrid ) ; 
    normalization_factor/=maxj;
  } else {
    maxj = pdata->phibins;
  }
  if ( pdata->zbins == 1 ) {
    maxk = (int) lbpar_gpu.dim_z;
    normalization_factor/=maxk;
  } else {
    maxk = pdata->zbins;
  }

  for (int i = 0; i<n_A; i++) {
    A[i]=0;
  }

  
  // copy radial profile to device
  radial_profile_data* pdata_device;
  cuda_safe_mem(cudaMalloc((void**)&pdata_device, sizeof(radial_profile_data)));
  cuda_safe_mem(cudaMemcpy(pdata_device, pdata,  sizeof(radial_profile_data), cudaMemcpyHostToDevice));

  // allocate data on device
  float* data_device;
  cuda_safe_mem(cudaMalloc((void**)&data_device, sizeof(float)*3*maxi*maxj*maxk));
  // kernellcall
  int threads_per_block = maxi;
  int blocks_per_grid_x = maxj;
  int blocks_per_grid_y = maxk;
  dim3 dim_grid = make_uint3(blocks_per_grid_x, blocks_per_grid_y, 1);
  KERNELCALL(fill_lb_radial_velocity_profile, dim_grid, threads_per_block, (nodes_a, pdata_device, data_device, node_f ));

  // allocate data on host
  float* host_data;
  host_data = (float*) Utils::malloc(sizeof(float)*3*maxi*maxj*maxk);

  // copy data back
  cuda_safe_mem(cudaMemcpy(host_data, data_device,  sizeof(float)*3*maxi*maxj*maxk, cudaMemcpyDeviceToHost));

  // average (or store)
  unsigned int linear_index;
  for (int i =0; i<maxi; i++)
    for (int j =0; j<maxj; j++)
      for (int k =0; k<maxk; k++) {
        linear_index = 0;
        if (pdata->rbins > 1)
          linear_index += i*pdata->phibins*pdata->zbins;
        if (pdata->phibins > 1)
          linear_index += j*pdata->zbins;
        if (pdata->zbins > 1)
          linear_index +=k;
        A[3*linear_index+0]+=host_data[3*(i*maxj*maxk + j*maxk + k)+0]*normalization_factor*lbpar_gpu.tau/lbpar_gpu.agrid;
        A[3*linear_index+1]+=host_data[3*(i*maxj*maxk + j*maxk + k)+1]*normalization_factor*lbpar_gpu.tau/lbpar_gpu.agrid;
        A[3*linear_index+2]+=host_data[3*(i*maxj*maxk + j*maxk + k)+2]*normalization_factor*lbpar_gpu.tau/lbpar_gpu.agrid;
      }

  // free device data
  cudaFree(pdata_device);
  cudaFree(data_device);

  // free host data
  free(host_data);

  return 0;
}

int statistics_observable_lbgpu_velocity_profile(profile_data* pdata, double* A, unsigned int n_A){
  unsigned int maxi, maxj, maxk;
  int linear_index;
  float normalization_factor=1;


  if ( pdata->xbins == 1 ) {
    maxi = (int) floor(lbpar_gpu.dim_x/lbpar_gpu.agrid);
    normalization_factor/=maxi;
  } else {
    maxi = pdata->xbins;
  }
  if ( pdata->ybins == 1 ) {
    maxj = (int) floor(lbpar_gpu.dim_y/lbpar_gpu.agrid);
    normalization_factor/=maxj;
  } else {
    maxj = pdata->ybins;
  }
  if ( pdata->zbins == 1 ) {
    maxk = (int) floor(lbpar_gpu.dim_z/lbpar_gpu.agrid);
    normalization_factor/=maxk;
  } else {
    maxk = pdata->zbins;
  }

  for (int i = 0; i<n_A; i++) {
    A[i]=0;
  }

  
  // copy  profile to device
  profile_data* pdata_device;
  cuda_safe_mem(cudaMalloc((void**)&pdata_device, sizeof(profile_data)));
  cuda_safe_mem(cudaMemcpy(pdata_device, pdata,  sizeof(profile_data), cudaMemcpyHostToDevice));

  // allocate data on device
  float* data_device;
  cuda_safe_mem(cudaMalloc((void**)&data_device, sizeof(float)*3*maxi*maxj*maxk));
  // kernellcall
  int threads_per_block = maxi;
  int blocks_per_grid_x = maxj;
  int blocks_per_grid_y = maxk;
  dim3 dim_grid = make_uint3(blocks_per_grid_x, blocks_per_grid_y, 1);

  KERNELCALL(fill_lb_velocity_profile, dim_grid, threads_per_block, (nodes_a, pdata_device, data_device, node_f));
  

  // allocate data on host
  float* host_data;
  host_data = (float*) Utils::malloc(sizeof(float)*3*maxi*maxj*maxk);

  // copy data back
  cuda_safe_mem(cudaMemcpy(host_data, data_device,  sizeof(float)*3*maxi*maxj*maxk, cudaMemcpyDeviceToHost));

  // average (or store)
  unsigned int i, j, k;
  for ( i = 0; i < maxi; i++ ) {
    for ( j = 0; j < maxj; j++ ) {
      for ( k = 0; k < maxk; k++ ) {
        linear_index = 0;
        if (pdata->xbins > 1)
          linear_index += i*pdata->ybins*pdata->zbins;
        if (pdata->ybins > 1)
          linear_index += j*pdata->zbins;
        if (pdata->zbins > 1)
          linear_index +=k;

        A[3*linear_index+0]+=host_data[3*(i*maxj*maxk + j*maxk + k)+0]*normalization_factor*lbpar_gpu.tau/lbpar_gpu.agrid;
        A[3*linear_index+1]+=host_data[3*(i*maxj*maxk + j*maxk + k)+1]*normalization_factor*lbpar_gpu.tau/lbpar_gpu.agrid;
        A[3*linear_index+2]+=host_data[3*(i*maxj*maxk + j*maxk + k)+2]*normalization_factor*lbpar_gpu.tau/lbpar_gpu.agrid;
      }
    }
  }

  // free device data
  cudaFree(pdata_device);
  cudaFree(data_device);

  // free host data
  free(host_data);

  return 0;
}

#endif /* LB_GPU */
