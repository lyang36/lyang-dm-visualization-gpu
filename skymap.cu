/**************************************************************/
/*                    Class: skymap                           */
/*	              Generate the skymap.                        */ 
/*           Author: Lin Yang      03/07/2012                 */
/**************************************************************/
#include <cmath>
#include <fstream>
#include <string>
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


#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <cuda.h>
#include <cstdio>


using namespace std;;



Skymap::Skymap(){
    //if set _reload, then reload
	_reload = NULL;
    //if set _rotate, then rotate
    //if not set, not doing rotate
	_rotate = NULL;
}

bool Skymap::creat_map(){

	if(_reload == NULL)	reload = false;
	else reload = *_reload;

	if(_rotate == NULL) rotate = false;
	else rotate = *_rotate;

	//allocate memory for map
	//const int LP = 10000;
    
    //the particle numbers
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
		exit(0);
	}
    data_input_file.read((char*)&Nparts, sizeof(Nparts));  
	cout << "Particles: " << Nparts << endl;
	Np = Nparts;
    
    //setup cpu-memory for particles
	num_p = 0;
	particles = new MapParticle[CPU_chunk];
    //sorted cpu-memory for particles
	MapParticle * sorted_particles = new MapParticle[CPU_chunk];
    
    //setup observation position
	Real * opos = master->params.opos;
    
    //setup Nside of Healpix map
	long Nside = master->map.Nside;
    
    //setup Total Pix number of Healpix map
	long Npix_in_map = master->map.Npix;
    
    //this guy is the same of Npix_in_map
	//int Npix_map = 12 * Nside * Nside;
    
    //the omega element of each pix
	Real dOmega = 4.0 * PI / Npix_in_map;
	Real theta0 = acos( 1.0 - dOmega/(2.0*PI) );

	
    //alloc and initialize the allksymap array
	Real * allskymap = (Real *) calloc(Npix_in_map, sizeof(Real));
    
    //pointers to the rotmatrix in the GPU memory
	Real * dev_rotm = 0; //(should be constant)/
	Real * dev_opos = 0; //(should be constant)/
    
    //final factor for the flux
	double fluxfactor = master->codeunits.annihilation_flux_to_cgs;
    
    //pointers to the particle memory in the GPU memory
	MapParticle * dev_par = 0;
    
    //key and values for sorting
	//int * host_keys;
	//int * dev_keys;
	//int * dev_values;

    //checkout is there any GPU
	cudaError_t cudaStatus = cudaSetDevice(0);
	if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaSetDevice failed!  Do you have a CUDA-capable GPU installed?");
        return false;
    }

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
	int parsize = PRE_chunk > GPU_chunk ? PRE_chunk : GPU_chunk;
	cudaStatus = cudaMalloc((void**)&dev_par, sizeof(MapParticle) * parsize);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMalloc failed!");
		return false;
	}


	//use for sorting
	thrust::device_vector<int> dev_key(CPU_chunk);
	thrust::device_vector<int> dev_val(CPU_chunk);
	thrust::host_vector<int> host_val(CPU_chunk);
	// obtain raw pointer to device vectors memory
	int * pd_key = thrust::raw_pointer_cast(&dev_key[0]);
	int * pd_val = thrust::raw_pointer_cast(&dev_val[0]);


	cout << "Creating map!!!" << endl;
	//int rec = Nparts / GPU_chunk / 50;
    
    //recording time
	clock_t time_start;
	time_start = clock(); 
	
#ifdef _DEBUG__LY__
	for(int _ip = 0, _jp=0 ; _ip < Nparts, _jp<1; _jp ++ ){
#else
	for(int _ip = 0, _jp=0 ; _ip < Nparts; _jp ++ ){ 
#endif
/*****************************************************************************************************************/
		clock_t time_loop_start = clock();
		cout << ">>>>CPU_chunk " << _jp << "--- Particles: " << CPU_chunk + _ip << "/"<< Nparts << "..." << endl;
        
		clock_t time_step_start = clock();
                
		int nmax = 0;
        
        //tnmax is the number of particles read into the CPU memory from the hard drive
		int tnmax = 0;

		//read to CPU chunk
		if( (Nparts - _ip) >= CPU_chunk ){//read a block of data
			data_input_file.read((char*)particles, sizeof(MapParticle) * CPU_chunk); 	
			_ip += CPU_chunk;
			tnmax = CPU_chunk;
		}else{
			tnmax = (Nparts - _ip);
			data_input_file.read((char*)particles, sizeof(MapParticle) * tnmax);
			_ip += tnmax;
		}
        std::cout <<"1) read from disk cost: " << (clock() - time_step_start) / 1000.0 << " secs. "<< std::endl;
        
		//step 1: pre-deal with particles
		//get the start point of pre-process data
        time_step_start = clock();
		for(int _pt =0; _pt < CPU_chunk; ){
			if( (Nparts - _pt) >= PRE_chunk ){//read a block of data
				nmax = PRE_chunk;
			}else{
				nmax = (CPU_chunk - _pt);
			}

			cudaStatus = cudaMemcpy(dev_par, particles + _pt, sizeof(MapParticle) * nmax, cudaMemcpyHostToDevice);
			if (cudaStatus != cudaSuccess) {
				fprintf(stderr, "cudaMemcpy failed!");
				return false;
			}
			cudaStatus = doWithCuda_pre(PRE_chunk, Nside, theta0, 1, nmax, allskymap,
			 dev_par, particles + _pt, dev_rotm, dev_opos, pd_key + _pt, pd_val + _pt, _pt);
			if (cudaStatus != cudaSuccess) {
				fprintf(stderr, "Kernel failed!");
				return false;
			}
			_pt += nmax;
		}
		//cudaFree(dev_par);
		std::cout <<"2) pre-calculating cost: " << (clock() - time_step_start) / 1000.0 << " secs. "<< std::endl;

		//step 2: sort
		time_step_start = clock();
		// interface to CUDA code
		thrust::sort_by_key(dev_key.begin(), dev_key.end(), dev_val.begin());
		thrust::copy(dev_val.begin(), dev_val.end(),host_val.begin());
        //actually sort on the particles
		for(int _pkk = CPU_chunk - 1; _pkk >=0; _pkk --){
			int pg =host_val[_pkk];
			sorted_particles[_pkk] = particles[pg];		
		}
		{//swape the sorted particles with the unsorted particles
			MapParticle * temp;
			temp = particles;
			particles = sorted_particles;
			sorted_particles = temp;
		}

		std::cout <<"3) sort cost: " << (clock() - time_step_start) / 1000.0 << " secs. "<< std::endl;
        
		//step3: calculate flux
		time_step_start = clock();
		for(int _pt =0, _ptn = 0; _pt < CPU_chunk; _ptn++ ){
			if( (Nparts - _pt) >= GPU_chunk ){//read a block of data
				nmax = GPU_chunk;
			}else{
				nmax = (CPU_chunk - _pt);
			}
			//if(_pt < 2031616){ _pt += nmax; continue;}
			cudaStatus = cudaMemcpy(dev_par, particles + _pt, sizeof(MapParticle) * nmax, cudaMemcpyHostToDevice);
			if (cudaStatus != cudaSuccess) {
				fprintf(stderr, "cudaMemcpy failed!");
				return false;
			}
			//decide whether or not do step 1
			//if > 5000 exceeds 
			cudaStatus = doWithCuda_Par(GPU_chunk, Nside, theta0, 1, nmax, allskymap,
			dev_par, particles + _pt, dev_rotm, dev_opos);
			if (cudaStatus != cudaSuccess) {
				fprintf(stderr, "Kernel failed!");
				return false;
			}
			_pt += nmax;
		    if(_ptn % 10 ==0 ){
				std::cout << ".";// << _pt << "/" << CPU_chunk << endl;
				std::cout.flush();
		
			}
		}
		std::cout << endl;
		std::cout <<"4) flux calculating cost: " << (clock() - time_step_start) / 1000.0 << " secs. "<< std::endl;
		std::cout << ">>>>chunk " << _jp << ": "<< (float)_ip / Nparts *100<<"% finished, costs " << (Real)(clock() - time_loop_start) / 1000.0 << 
			" secs, escaped: " << (Real)(clock() - time_start) / 1000.0 << " secs\n" << endl;
		_jp =_jp;
/*****************************************************************************************************************/
	}


	cout << endl;
	cout << "Time cosumed: " << (Real)(clock() - time_start) / 1000.0  << " seconds" << endl; 
#ifdef _DEBUG__LY__
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


	cout << "Writing to file \"" << *fits_filename << "\":" << endl;
	ofstream output_file (fits_filename -> c_str(), ios::out | ios::binary);
	if(output_file.good()){
		output_file.write ((char *)allskymap, Npix_in_map * sizeof(Real));
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
	return true;

}




