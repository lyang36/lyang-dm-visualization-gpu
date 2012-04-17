/**************************************************************/
/*                    Class: skymap                           */
/*	              Generate the skymap.                        */ 
/*           Author: Lin Yang      03/07/2012                 */
/**************************************************************/
#include <cmath>
#include <fstream>
#include <string>
//#include "tipsydefs.h"
#include "mapparticle.h"
#include "info.h"
#include "structures.h"
#include "info.h"
#include "VL2_debug.h"
#include <iostream>
#include <vector>
#include "skymap.h"
#include <ctime>
#include "kernel.h"

#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <thrust/sort.h>
#include <thrust/copy.h>


//#include "cuda_runtime.h"
//#include "device_launch_parameters.h"
#include <cuda.h>
#include <cstdio>
//#include <b40c/radix_sort/enactor.cuh>
//#include <b40c/util/multiple_buffering.cuh>

//#include 


using namespace std;;



Skymap::Skymap(){
	_reload = NULL;
	_rotate = NULL;
}

bool Skymap::creat_map(){
	cudaError_t cudaStatus;
	if(_reload == NULL)	reload = false;
	else reload = *_reload;

	if(_rotate == NULL) rotate = false;
	else rotate = *_rotate;

	//allocate memory for map
	//const int LP = 10000;
	int Nparts = 0;

	//get rotation matrix
	if ( rotate ){
		for(int _i = 0; _i < 3; _i++){
			rotmatrix[0 + _i] = master->rotmatrix[0][_i];
			rotmatrix[3 + _i] = master->rotmatrix[1][_i];
			rotmatrix[6 + _i] = master->rotmatrix[2][_i];
		}
	}else{
		for(int _i = 0; _i < 3; _i++){
			rotmatrix[0 + _i] = 0;
			rotmatrix[3 + _i] = 0;
			rotmatrix[6 + _i] = 0;
		}
		rotmatrix[0] = 1.0;
		rotmatrix[4] = 1.0;
		rotmatrix[8] = 1.0;
	}
#ifdef _DEBUG__LY__
	//cout << "good1" <<endl;
#endif
	// Read particle_numbers
    ifstream data_input_file((*datafile).c_str(), ios::binary);
	if(data_input_file.bad()){
		cout << "Data Error!!!" << endl;
		return false;
	}
    data_input_file.read((char*)&Nparts, sizeof(Nparts));  
	cout << "Particles: " << Nparts << endl;
	Np = Nparts;
	num_p = 0;
	//particles = new MapParticle[CPU_trunk];
	//MapParticle * sorted_particles = new MapParticle[CPU_trunk];
	Real * opos = master->params.opos;
//	Real fluxes;//  = master.codeunits.annihilation_flux_to_cgs * density * mass / (4.0 * !pi * distances^2)
	long Nside = master->map.Nside;
	long Npix_in_map = master->map.Npix;
	int Npix_map = 12 * Nside*Nside;
	Real dOmega = 4.0 * PI / Npix_map;
	Real theta0 = acos( 1.0 - dOmega/(2.0*PI) );
	//cudaError_t cudaStatus;
	
	Real * allskymap = (Real *) calloc(Npix_map, sizeof(Real));
	//Real * dev_allskymap;
	Real * dev_rotm = 0; //(should be constant)/
	Real * dev_opos = 0; //(should be constant)/
	double fluxfactor = master->codeunits.annihilation_flux_to_cgs;
	MapParticle * dev_par = 0;
	//int * host_keys;
	//int * dev_keys;
	//int * dev_values;

	cudaStatus = cudaSuccess;
	cudaStatus = cudaSetDevice(0);
	if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaSetDevice failed!  Do you have a CUDA-capable GPU installed?");
        return false;
    	}
		
	particles = new MapParticle[CPU_trunk];
	MapParticle * sorted_particles = new MapParticle[CPU_trunk];

	//copy rotmatrix into GPU
	cudaStatus = cudaMalloc((void**)&dev_rotm, sizeof(Real) * 9);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMalloc failed!");
		return false;
    }
	cudaStatus = cudaMemcpy(dev_rotm, rotmatrix, sizeof(Real) * 9, cudaMemcpyHostToDevice);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMemcpy failed!");
        return false;
    }

	//copy o_pos into GPU
	cudaStatus = cudaMalloc((void**)&dev_opos, sizeof(Real) * 3);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMalloc failed!");
		return false;
    }
	cudaStatus = cudaMemcpy(dev_opos, opos, sizeof(Real) * 3, cudaMemcpyHostToDevice);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMemcpy failed!");
        return false;
    }

	//allocate particle memery into GPU
	int parsize = PRE_trunk > MAX_Num_Particle ? PRE_trunk : MAX_Num_Particle;
	cudaStatus = cudaMalloc((void**)&dev_par, sizeof(MapParticle) * parsize);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMalloc failed!");
		return false;
	}

	//<debug
	//int c = sizeof(MapParticle);
	//int d = sizeof(MapParticle) * nmax;
	//debug>
	//thrust::host_vector<int> host_key(CPU_trunk);
	//use for sorting
	thrust::device_vector<int> dev_key(CPU_trunk);
	thrust::device_vector<int> dev_val(CPU_trunk);
	thrust::host_vector<int> host_val(CPU_trunk);
	// obtain raw pointer to device vectors memory
	int * pd_key = thrust::raw_pointer_cast(&dev_key[0]);
	int * pd_val = thrust::raw_pointer_cast(&dev_val[0]);


	cout << "Creating map!!!" << endl;
	//cout << "---10---20---30---40---50---60---70---80---90--100%\n";
	//int rec = Nparts / MAX_Num_Particle / 50;
	clock_t time;
	time = clock(); 
	
#ifdef _DEBUG__LY__
	for(int _ip = 0, _jp=0 ; _ip < Nparts, _jp<1; _jp ++ ){
		//if(_jp != 1017) continue;
#else
	for(int _ip = 0, _jp=0 ; _ip < Nparts; _jp ++ ){ 
#endif
		clock_t time_a = clock();
		cout << "CPU_trunck " << _jp << "--- Particles: " << CPU_trunk + _ip << "/"<< Nparts << "..." << endl;
		clock_t time_b = clock();
		int nmax = 0;
		int tnmax = 0;

		//read to CPU trunk
		if( (Nparts - _ip) >= CPU_trunk ){//read a block of data
			data_input_file.read((char*)particles, sizeof(MapParticle) * CPU_trunk); 	
			_ip += CPU_trunk;
			tnmax = CPU_trunk;
		}else{
			tnmax = (Nparts - _ip);
			data_input_file.read((char*)particles, sizeof(MapParticle) * tnmax);
			_ip += tnmax;
		}
		//if(_jp < 6) continue;
		//step 1: pre-deal with particles
		//get the start point of pre-process data
		for(int _pt =0; _pt < CPU_trunk; ){
			if( (Nparts - _pt) >= PRE_trunk ){//read a block of data
				nmax = PRE_trunk;
			}else{
				nmax = (CPU_trunk - _pt);
			}

			cudaStatus = cudaMemcpy(dev_par, particles + _pt, sizeof(MapParticle) * nmax, cudaMemcpyHostToDevice);
			if (cudaStatus != cudaSuccess) {
				fprintf(stderr, "cudaMemcpy failed!");
				return false;
			}
			cudaStatus = doWithCuda_pre(PRE_trunk, Nside, theta0, 1, nmax, allskymap,
			 dev_par, particles + _pt, dev_rotm, dev_opos, pd_key + _pt, pd_val + _pt, _pt);
			if (cudaStatus != cudaSuccess) {
				fprintf(stderr, "Kernel failed!");
				return false;
			}
			_pt += nmax;
		}
		//cudaFree(dev_par);
		std::cout <<"step1 cost: " << (clock() - time_b) / 1000.0 << " secs. "<< std::endl;

		//step 2: sort
		time_b = clock();
		// interface to CUDA code
		thrust::sort_by_key(dev_key.begin(), dev_key.end(), dev_val.begin());
		thrust::copy(dev_val.begin(), dev_val.end(),host_val.begin());
		for(int _pkk = CPU_trunk - 1; _pkk >=0; _pkk --){
			int pg =host_val[_pkk];
			sorted_particles[_pkk] = particles[pg];		
		}
		{//swape
			MapParticle * temp;
			temp = particles;
			particles = sorted_particles;
			sorted_particles = temp;
		}
/*		{//clear
			dev_key.clear();
			dev_val.clear();
			host_val.clear();
			cudaFree(pd_key);
			cudaFree(pd_val);
			free(ph_val);
		}*/
		std::cout <<"sort cost: " << (clock() - time_b) / 1000.0 << " secs. "<< std::endl;

		//step3: calculate flux
		time_b = clock();
		/*cudaStatus = cudaMalloc((void**)&dev_par, sizeof(MapParticle) * MAX_Num_Particle);
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaMalloc failed!");
			exit(0);
		}*/

		for(int _pt =0, _ptn = 0; _pt < CPU_trunk; _ptn++ ){
			if( (Nparts - _pt) >= MAX_Num_Particle ){//read a block of data
				nmax = MAX_Num_Particle;
			}else{
				nmax = (CPU_trunk - _pt);
			}
			//if(_pt < 2031616){ _pt += nmax; continue;}
			cudaStatus = cudaMemcpy(dev_par, particles + _pt, sizeof(MapParticle) * nmax, cudaMemcpyHostToDevice);
			if (cudaStatus != cudaSuccess) {
				fprintf(stderr, "cudaMemcpy failed!");
				return false;
			}
			//decide whether or not do step 1
			//if > 5000 exceeds 
			cudaStatus = doWithCuda_Par(MAX_Num_Particle, Nside, theta0, 1, nmax, allskymap,
			dev_par, particles + _pt, dev_rotm, dev_opos);
			if (cudaStatus != cudaSuccess) {
				fprintf(stderr, "Kernel failed!");
				return false;
			}
			_pt += nmax;
		    if(_ptn % 10 ==0 ){
				std::cout << ".";// << _pt << "/" << CPU_trunk << endl;
				std::cout.flush();
		
			}
		}
		std::cout << endl;
		std::cout <<"step3 cost: " << (clock() - time_b) / 1000.0 << " secs. "<< std::endl;
		std::cout << "trunk " << _jp << ": "<< (float)_ip / Nparts *100<<"% finished, costs " << (Real)(clock() - time_a) / 1000.0 << 
			" secs, escaped: " << (Real)(clock() - time) / 1000.0 << " secs\n" << endl;
		_jp =_jp;
/*****************************************************************************************************************/
	}

	/*cudaStatus = cudaMemcpy(allskymap, dev_allskymap, sizeof(Real) * Npix_map, cudaMemcpyDeviceToHost);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMemcpy failed!");
		return false;
	}*/

	cout << endl;
	time = clock() - time;
	cout << "Time cosumed: " << (Real) time / 1000.0  << " seconds" << endl; 
#ifdef _DEBUG__LY__
	//cout << "good3" <<endl;
	print_out_master(master);
#endif

	//divide by solid angle of map's pixels
	//conversion of allsky from g^2 cm^-5 to Gev^2 cm^-6 kpc

	Real unit_factor = pow(pow((master -> natconst.c_in_cgs), 2) /
	(master->units.eV_in_cgs * 1.0e9), 2) / (master->units.Mpc_in_cgs * 1.0e-3);
	for(int i = 0; i < Npix_in_map; i++){
		float amap = (float)((double)(unit_factor) / 
			(double)(master->map.dOmega) * double (fluxfactor) * (double) allskymap[i]);
		allskymap[i] = amap;
	}
//	cblas_dscal( Npix_in_map, unit_factor / master->map.dOmega, allskymap, 1);


#ifdef _DEBUG__LY__

	Real rot1[9];
	cudaStatus = cudaMemcpy(rot1, dev_rotm, 9*sizeof(Real), cudaMemcpyDeviceToHost);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMemcpy failed!");
		return false;
	}

//	cout << "flux" << fluxes << endl;
	cout << "unit_factor: " << unit_factor << endl;
	cout << "dOmega: " << master->map.dOmega << endl;
	cout << "All Skymap:" <<endl;
	cout << "skymap[0]: " << allskymap[0] << endl;
	cout << "skymap[1]: " << allskymap[1] << endl;
	cout << "skymap[2]: " << allskymap[2] << endl;
	cout << "skymap[100]: " << allskymap[100] << endl;
	cout << "skymap[10000]: " << allskymap[10000] << endl;
	cout << "skymap[1000000]: " << allskymap[1000000] << endl;
	cout << "skymap[2977220]: " << allskymap[2977220] << endl;
	/*for(int i =0; i < 12*512*512; i++){
		if(allskymap[i] !=0 ){
			cout << i <<"-----"<< allskymap[i]<<endl;		
		}
	}*/
#endif

	//float * newskymap = (float *) calloc(Npix_map, sizeof(float));
	//for (int i=0; i < Npix_map; i++){
	//	newskymap[i] = allskymap[i];
	//}
	//write_healpix_map(newskymap,Nside,"trymap.fits", 0, "G");
	cout << "Writing to file \"" << *fits_filename << "\":" << endl;
	ofstream output_file (fits_filename -> c_str(), ios::out | ios::binary);
	if(output_file.good()){
		output_file.write ((char *)allskymap, Npix_map * sizeof(Real));
	}else{
		cout << "Writing Error!";
	}
	output_file.close();
	cout << "success!" << endl;
	free(allskymap);
	data_input_file.close();
	//free(newskymap);
   	cudaFree(dev_rotm);
    	cudaFree(dev_opos);
	cudaFree(dev_par);
	//cudaFree(dev_allskymap);
	
	cudaStatus = cudaDeviceReset();
	if (cudaStatus != cudaSuccess){
		fprintf(stderr, "cudaDeviceReset failed!");
                return false;
	}
	return true;

}




