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
#include "querydisc.cpp"



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
			int k = maplist[i].pix;
		}
	}
	//std::cout << " maplist: " << init << std::endl;
	free(maplist);
	cudaFree(dev_maplist);
	//std::cout << (clock() -starttime) / 1000.0 << " sec -> kernel end" << std::endl;
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





__global__
void generate_map_step0(const long Nside, 
		 const Real theta0, const long nmax,
		 MapParticle * particles,
		 Real * rotmatrix, Real * opos,
		 int * pd_key, int * pd_val, int stpoint){
		//fluxfactor = master->codeunits.annihilation_flux_to_cgs 
	//int i = threadIdx.x;
	int i = blockIdx.x * blockDim.x + threadIdx.x;
	if( i >= nmax){
		return;
	}
	MapParticle mp = particles[i];
	Real density = mp.density;
	Real xpos = mp.xpos;
	Real ypos = mp.ypos;
	Real zpos = mp.zpos;
	Real mass = mp.mass;
	Real hsmooth = mp.hsmooth;
	Real distances;
	Real fluxes;
	Real costheta;
	Real phi;
	Real theta;
	Real angular_radius;

#ifdef POS_FLOAT
	float3 _opos = {opos[0], opos[1], opos[2]};
#else
	double3 _opos = {opos[0], opos[1], opos[2]};
#endif

	/*since if density <0, the threads becomes null*/
	/*please make sure all the density is greater or equal to 0*/
	/*this is very important*/
	if( density < 0.0){
		particles[i].xpos = 0;
		//xpos holds pixel num
		particles[i].hsmooth = 0;//;angular_radius;
		particles[i].ypos = 0;//theta;
		particles[i].zpos = 0;//phi;
		particles[i].mass = 0;//flux;
		//particles[i].density = weight_norm;
		pd_val[i] = stpoint + i;
		pd_key[i] = 0;
	}else if( density >= 0.0){
		distances = sqrt( (xpos-_opos.x) * (xpos-_opos.x) + 
			(ypos-_opos.y) * (ypos-_opos.y) +
			(zpos-_opos.z) *(zpos-_opos.z) );
		fluxes = density * mass / (4.0 * PI * distances * distances);

		calc_angles( xpos-_opos.x, ypos-_opos.y, zpos-_opos.z, distances, 
			opos, rotmatrix, costheta, phi);
		theta = acos(costheta);
		if(distances != 0.0){
			angular_radius = hsmooth / distances;}
		else
			{angular_radius = PI;}
		//hsmooth holds angular_radius
		particles[i].hsmooth = angular_radius;
		particles[i].ypos = theta;
		particles[i].zpos = phi;
		particles[i].mass = fluxes;
		//delete step1
		Real ad =  2 * angular_radius * (1 + 0.1 );
		int pixss = floor(ad * ad * ((double)(12L * Nside * Nside)) / 4) + 1; //sort by anglar_radius
		particles[i].xpos = pixss;
		particles[i].density = 1;//2.0 * angular_radius * angular_radius * 12 * 512 * 512 / 4 / PI;
		pd_key[i] = pixss;
		pd_val[i] = stpoint + i;
		//density holds norm 1.62707
	}
}




__global__
void generate_map_step1(const long Nside, 
		 const Real theta0, const long nmax,
		 MapParticle * particles,
		 Real * rotmatrix, Real * opos){
		//fluxfactor = master->codeunits.annihilation_flux_to_cgs 
	//int i = threadIdx.x;
	int i = blockIdx.x * blockDim.x + threadIdx.x;
	if( i >= nmax){
		return;
	}
	MapParticle mp = particles[i];
	Real density = mp.density;
	Real hsmooth = mp.hsmooth;
	Real fluxes = mp.mass;
	Real phi = mp.zpos;
	Real theta = mp.ypos;
	Real angular_radius = mp.hsmooth;

#ifdef POS_FLOAT
	float3 _opos = {opos[0], opos[1], opos[2]};
#else
	double3 _opos = {opos[0], opos[1], opos[2]};
#endif

	if( density >= 0.0){
		//density holds norm
	/*------------------------------------------add particles------------------------------------------------------*/
		{
			//pointing p (theta,phi);
			Real vecx = sin(theta) * cos(phi);
			Real vecy = sin(theta) * sin(phi);
			Real vecz = cos(theta);
			if( 2.0*angular_radius < theta0 ) {
				/*should be atomic plus */
				//allskymap[pix] += fluxes;
				particles[i].xpos = 1;
				particles[i].density = 1;
			}else{
				Real weight_norm = 0.0;
				int Npix = 0;
				query_disc_getnorm(Nside, theta, phi, 2.0*angular_radius, vecx, vecy, vecz, weight_norm, Npix);
				particles[i].xpos = Npix;
				particles[i].density = weight_norm;
			}
		}
	}
}




__global__
void generate_map_step2(const long Nside, 
		 const Real theta0, const long nmax,
		  MapParticle * particles, mapMap *maplist,
		 Real * rotmatrix, Real * opos){

	int i = blockIdx.x * blockDim.x + threadIdx.x;

	if( i >= nmax){
		return;
	}
	MapParticle mp = particles[i];
	Real fluxes = mp.mass;
	int startp = (int)mp.xpos;
	Real theta = mp.ypos;
	Real phi = mp.zpos;
	Real angular_radius = mp.hsmooth;
	Real weight_norm = mp.density;
	Real density = mp.density;

	/*since if density <0, the threads becomes null*/
	/*please make sure all the density is greater or equal to 0*/
	/*this is very important*/
	if( density >= 0.0){

		Real costheta = cos(theta);
	/*------------------------------------------add particles------------------------------------------------------*/
		{
			//pointing p (theta,phi);
			Real vecx = sin(theta) * cos(phi);
			Real vecy = sin(theta) * sin(phi);
			Real vecz = costheta;
			//ang2vec(theta,phi,&vec);

			long pix;
			ang2pix_ring(Nside,theta,phi,&pix);
			if( 2.0*angular_radius < theta0 ) {
				/*should be atomic plus */
				//allskymap[pix] += fluxes;
				maplist[startp].pix = pix;
				maplist[startp].signal = fluxes;
				maplist[startp].factor = 1;
			}else{
				long index = startp;
				double getnorm = 0;
				query_disc_calcflux(Nside, theta, phi, 2.0*angular_radius, vecx, vecy, vecz, 
					weight_norm, fluxes, maplist, index, getnorm);
				maplist[startp].factor = getnorm;
				if(getnorm == 0){ // consider Npix <2 --> keep the flux conservative
					maplist[startp].pix = pix;
					maplist[startp].signal = fluxes;
					maplist[startp].factor = 1;
				}
			}
		}
	}
}
