#ifndef __CUDACC__
#define __CUDACC__
#endif
#include "cuda.h"
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <stdio.h>
#include "kernel.h"
#include <iostream>
#include "structures.h"
#include <ctime>

#define MAX_PIXES 500

///////////
#ifndef M_PI
#define M_PI           3.14159265358979323846
#endif

#ifndef PI
#define PI           M_PI
#endif

const dim3 dimBlock(256, 1, 1);


//prototype:
__device__ void ang2pix_ring( const long nside, Real theta, Real phi, long *ipix);
__device__ 
void calc_angles( Real xpos, Real ypos, Real zpos, Real &distances, 
					 Real * opos, Real *rotmatrix, Real & costheta, Real &phi);

__device__
void pix2ang_ring( long nside, long ipix, Real *theta, Real *phi);

__device__ __host__
void query_disc_getnorm(long nside_, const Real theta, const Real phi, 
	const Real radius, const Real vecx,  const Real vecy, const Real vecz, 
	Real &weight_norm, int &Npix);


__global__
void generate_map_step1(const long Nside, 
		 const Real theta0, const long nmax,
		 MapParticle * particles,
		 Real * rotmatrix, Real * opos);


__global__
void generate_map_step2(const long Nside, 
		 const Real theta0, const long nmax,
		 MapParticle * particles, mapMap *skymap,
		 Real * rotmatrix, Real * opos);


__global__
void generate_map_step0(const long Nside, 
		 const Real theta0, const long nmax,
		 MapParticle * particles,
		 Real * rotmatrix, Real * opos,
		 int * pd_key, int * pd_val, int stpoint);

__global__
void generate_map_step3(const long nmax, 
	mapMap *maplist);



__device__ 
const Real inv_twopi = 1.0 / (2.0 * PI);
__device__ 
const Real pi = PI;
__device__ 
const Real twothird = 2.0 / 3.0;


#include "ang2pix.cpp"
#include "calcangles.cpp"
#include "pix2ang.cpp"
#include "inring.cpp"
#include "ringabove.cpp"
#include "querydisc.cpp"
#include "generatemap.cpp"



cudaError_t doWithCuda_pre(const long MAX_Num_Paritcle, const long Nside, 
		 const Real theta0, const Real fluxfactor, const long nmax,
		 Real *allskymap, MapParticle * dev_par, MapParticle * host_par,
		 Real * dev_rotm, Real * dev_opos, int * pd_key, int * pd_val, int stpoint)
{
	cudaError_t cudaStatus;
	dim3 dimGrid(MAX_Num_Paritcle/dimBlock.x + 1, 1, 1);
	//std::cout << 0 << " sec -> step1 start" << std::endl;
	clock_t starttime = clock();
	generate_map_step0<<<dimGrid, dimBlock>>>(Nside, theta0, nmax,
		dev_par, dev_rotm, dev_opos, pd_key, pd_val, stpoint);
	cudaStatus = cudaDeviceSynchronize();
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaDeviceSynchronize returned error code %d after launching generate_map_step0!\n", cudaStatus);
		return cudaStatus;
	}
	//std::cout <<"step0 cost: " << (clock() - starttime) / 1000.0 << " secs. "<< std::endl;
	cudaStatus = cudaMemcpy(host_par, dev_par, sizeof(MapParticle) * nmax, cudaMemcpyDeviceToHost);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMemcpy failed!");
		exit(0);
	}
	return cudaStatus;
}


cudaError_t run_kernel(const long MAX_Num_Paritcle, const long Nside, 
		 const Real theta0, const Real fluxfactor, const long nmax,
		 Real *allskymap, MapParticle * dev_par, MapParticle * host_par,
		 Real * dev_rotm, Real * dev_opos, int init){
	//std::cout << "Add map length cost" << (clock() -starttime) / 1000.0 << " sec"<< std::endl;
	//copy the host_par back to device_par
	//Real fluxfactor = master->codeunits.annihilation_flux_to_cgs;
	dim3 dimGrid(MAX_Num_Paritcle/dimBlock.x + 1, 1, 1);
	mapMap * maplist;
	mapMap * dev_maplist;

	cudaError_t cudaStatus;
	cudaStatus = cudaMemcpy(dev_par, host_par, sizeof(MapParticle) * nmax, cudaMemcpyHostToDevice);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMemcpy failed!");
		return cudaStatus;
	}
	//relocate the memory for skymap
	maplist = (mapMap *)calloc(init,sizeof(mapMap));
	cudaStatus = cudaMalloc((void**)&dev_maplist, sizeof(mapMap) * init);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMalloc failed!");
		return cudaStatus;
	}
	//initiate maplist
	cudaStatus = cudaMemcpy(dev_maplist, maplist, sizeof(mapMap) * init, cudaMemcpyHostToDevice);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMemcpy failed!");
		return cudaStatus;
	}

	//std::cout << (clock() -starttime) / 1000.0 << " sec -> step2 start" << std::endl;
	generate_map_step2<<<dimGrid, dimBlock>>>(Nside, theta0, nmax,
		dev_par, dev_maplist, dev_rotm, dev_opos);
	cudaStatus = cudaDeviceSynchronize();
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaDeviceSynchronize returned error code %d after launching generate_map_step2!\n", cudaStatus);
		return cudaStatus;
	}

	cudaStatus = cudaMemcpy(maplist, dev_maplist, sizeof(mapMap) * init, cudaMemcpyDeviceToHost);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMemcpy failed!");
		return cudaStatus;
	}
	//add to allskymap
	double norm = 1;
	for(int i = 0; i < init; i++){
		if(maplist[i].pix < 12*512L*512L && maplist[i].pix >= 0){
			if(maplist[i].factor > 0 )
				{norm = maplist[i].factor;}
			double ff = maplist[i].signal/norm;
			allskymap[maplist[i].pix] += ff;
		}
		else{
			//int k = maplist[i].pix;
		}
	}
	//std::cout << " maplist: " << init << std::endl;
	free(maplist);
	cudaFree(dev_maplist);
	//std::cout << (clock() -starttime) / 1000.0 << " sec -> kernel end" << std::endl;
    return cudaStatus;
}


// Helper function for using CUDA to add vectors in parallel.
cudaError_t doWithCuda_Par(const long MAX_Num_Paritcle, const long Nside, 
		 const Real theta0, const Real fluxfactor, const long nmax,
		 Real *allskymap, MapParticle * dev_par, MapParticle * host_par,
		 Real * dev_rotm, Real * dev_opos)
{
	cudaError_t cudaStatus;
	dim3 dimGrid(MAX_Num_Paritcle/dimBlock.x + 1, 1, 1);

	//std::cout << 0 << " sec -> step1 start" << std::endl;
	clock_t starttime = clock();
	cudaStatus = cudaMemcpy(dev_par, host_par, sizeof(MapParticle) * nmax, cudaMemcpyHostToDevice);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMemcpy failed!");
		return cudaStatus;
	}
    generate_map_step1<<<dimGrid, dimBlock>>>( Nside, 
		 theta0, nmax,
		 dev_par,
		dev_rotm, dev_opos);
	cudaStatus = cudaDeviceSynchronize();
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaDeviceSynchronize returned error code %d after launching generate_map_step1!\n", cudaStatus);
		return cudaStatus;
	}
	cudaStatus = cudaMemcpy(host_par, dev_par, sizeof(MapParticle) * nmax, cudaMemcpyDeviceToHost);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMemcpy failed!");
		return cudaStatus;
	}
	//calculting the index for maplist
	int init = 0;
	int _tenmax = 0;
	int maxPixx = MAX_Num_Paritcle * MAX_PIXES;
	for(int i=0; i < nmax; i++){
		//accumulate the particles and deal with the memery overflow
		int init1 = init;
		init += (int)(host_par[i].xpos);
		host_par[i].xpos = init1; 
		if(init1 > maxPixx){
			//run the kernel onec
			int st_p = _tenmax;
			_tenmax = i - _tenmax + 1;
			MapParticle * __dev_par = dev_par + st_p;
			cudaStatus = run_kernel(MAX_Num_Paritcle, Nside, theta0, fluxfactor, _tenmax, allskymap,
				__dev_par, host_par + st_p, dev_rotm, dev_opos, init);
			init = 0;
			_tenmax = i + 1;
		}

	}
	int _nmax = nmax - _tenmax;
	if( _nmax > 0 ){
		cudaStatus = run_kernel(MAX_Num_Paritcle,Nside,theta0,fluxfactor,_nmax,allskymap,
			dev_par + _tenmax, host_par + _tenmax, dev_rotm, dev_opos, init);
	}
	return cudaStatus;
}




