#ifndef __CUDA_KERNEL_LY__
#define __CUDA_KERNEL_LY__

#include "mapparticle.h"
#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#define MAX_PIX_PER_PAR 200 //might be incorrect
//#define MAX_ANGULAR_RADIUS 0.001

struct mapMap{
	int pix;
	Real signal;
	Real factor;
};


/*cudaError_t doWithCuda(const long MAX_Num_Paritcle, const long Npix_map,
		 const long Nside, 
		 const Real theta0, const Real fluxfactor, const long nmax,
		 Real *allskymap, MapParticle * particles,
		 Real * rotmatrix, Real * opos);*/
/*cudaError_t doWithCuda(const long MAX_Num_Paritcle, const long Nside, 
		 const Real theta0, const Real fluxfactor, const long nmax,
		 Real *allskymap, MapParticle * particles,
		 Real * rotmatrix, Real * opos, mapMap *dev_maplist, 
		 mapMap *host_maplist);*/

cudaError_t doWithCuda_Par(const long MAX_Num_Paritcle, const long Nside, 
		 const Real theta0, const Real fluxfactor, const long nmax,
		 Real *allskymap, MapParticle * dev_par, MapParticle * host_par,
		 Real * dev_rotm, Real * dev_opos);

/*__global__
void generate_map(const long Nside, 
		 const Real theta0, const Real fluxfactor, const long nmax,
		 Real *allskymap, MapParticle * particles,
		 Real * rotmatrix, Real * opos);*/

#endif

