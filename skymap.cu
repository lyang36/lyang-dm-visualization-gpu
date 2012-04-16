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
#include <algorithm>
#include <thrust/sort.h>

#include "cuda_runtime.h"
#include "device_launch_parameters.h"


using namespace std;;



Skymap::Skymap(){
	_reload = NULL;
	_rotate = NULL;
}

bool Skymap::creat_map(){

	if(_reload == NULL)	reload = false;
	else reload = *_reload;

	if(_rotate == NULL) rotate = false;
	else rotate = *_rotate;

	//allocate memory for map
	const int LP = 10000;
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
	//cout << "File name: " << datafile->c_str() << endl;
    data_input_file.read((char*)&Nparts, sizeof(Nparts));  
	cout << "Particles: " << Nparts << endl;
	Np = Nparts;
	num_p = 0;
	particles = new MapParticle[MAX_Num_Particle];
	Real * opos = master->params.opos;
//	Real fluxes;//  = master.codeunits.annihilation_flux_to_cgs * density * mass / (4.0 * !pi * distances^2)
	long Nside = master->map.Nside;
	long Npix_in_map = master->map.Npix;
	int Npix_map = 12 * Nside*Nside;
	Real dOmega = 4.0 * PI / Npix_map;
	Real theta0 = acos( 1.0 - dOmega/(2.0*PI) );

	
	Real * allskymap = (Real *) calloc(Npix_map, sizeof(Real));
	//Real * dev_allskymap;
	Real * dev_rotm = 0; //(should be constant)/
	Real * dev_opos = 0; //(should be constant)/
	double fluxfactor = master->codeunits.annihilation_flux_to_cgs;
	MapParticle * dev_par = 0;
	//mapMap * host_maplist;
	//mapMap * dev_maplist;

	cudaError_t cudaStatus = cudaSetDevice(0);
	if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaSetDevice failed!  Do you have a CUDA-capable GPU installed?");
        exit(1);
    }

	//host_maplist = (mapMap *)calloc(MAX_PIX_PER_PAR * MAX_Num_Particle,sizeof(mapMap));
	//cudaStatus = cudaMalloc((void**)&dev_maplist, sizeof(mapMap) * MAX_PIX_PER_PAR * MAX_Num_Particle);
    //if (cudaStatus != cudaSuccess) {
    //    fprintf(stderr, "cudaMalloc failed!");
	//	exit(0);
    //}


	//initiate the allskymap in GPU
/*	cudaStatus = cudaMalloc((void**)&dev_allskymap, sizeof(Real) * Npix_map);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMalloc failed!");
		exit(0);
    }
	cudaStatus = cudaMemcpy(dev_allskymap, allskymap,  sizeof(Real) * Npix_map, cudaMemcpyHostToDevice);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMemcpy failed!");
        exit(0);
    }*/

	//copy rotmatrix into GPU
	cudaStatus = cudaMalloc((void**)&dev_rotm, sizeof(Real) * 9);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMalloc failed!");
		exit(0);
    }
	cudaStatus = cudaMemcpy(dev_rotm, rotmatrix, sizeof(Real) * 9, cudaMemcpyHostToDevice);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMemcpy failed!");
        exit(0);
    }

	//copy o_pos into GPU
	cudaStatus = cudaMalloc((void**)&dev_opos, sizeof(Real) * 3);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMalloc failed!");
		exit(0);
    }
	cudaStatus = cudaMemcpy(dev_opos, opos, sizeof(Real) * 3, cudaMemcpyHostToDevice);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMemcpy failed!");
        exit(0);
    }
	
	//allocate particle memery into GPU
	cudaStatus = cudaMalloc((void**)&dev_par, sizeof(MapParticle) * MAX_Num_Particle);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMalloc failed!");
		exit(0);
    }
	//<debug
	//int c = sizeof(MapParticle);
	//int d = sizeof(MapParticle) * nmax;
	//debug>

	cout << "Creating map!!!" << endl;
	//cout << "---10---20---30---40---50---60---70---80---90--100%\n";
	int rec = Nparts / MAX_Num_Particle / 50;
	clock_t time;
	time = clock(); 
	
#ifdef _DEBUG__LY__
	for(int _ip = 0, _jp=0 ; _ip < Nparts, _jp<1; _jp ++ ){
		//if(_jp != 1017) continue;
#else
	for(int _ip = 0, _jp=0 ; _ip < Nparts; _jp ++ ){ /*use _jp < 100 to specify loops*///int _ip = 0; 
#endif
		cout << "block " << _jp << "--- Particles: " << MAX_Num_Particle + _ip << "/"<< Nparts << "..." << endl;
		clock_t time_b = clock();
		bool islast = false;
		int nmax = 0;
	
		if( (Nparts - _ip) >= MAX_Num_Particle ){//read a block of data
			data_input_file.read((char*)particles, sizeof(MapParticle) * MAX_Num_Particle); 	
			_ip += MAX_Num_Particle;
			nmax = MAX_Num_Particle;
			//sort particles, so that the parallelism is more efficient
		}else{
			islast = true;
			nmax = (Nparts - _ip);
			data_input_file.read((char*)particles, sizeof(MapParticle) * nmax);
			_ip += nmax;

		}
		//if(_jp != 168 ) continue;
		//sort particles, so that the parallelism is more efficient
		//cout << "sorting ... " << endl;
		//sort(particles, particles + nmax, (&compare));
		//thrust :: stable_sort (particles, particles + nmax, thrust::greater<MapParticle>());
		//cout << "starting ..." << endl;

		//copy particle data into GPU 
		cudaStatus = cudaMemcpy(dev_par, particles, sizeof(MapParticle) * nmax, cudaMemcpyHostToDevice);
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaMemcpy failed!");
			exit(0);
		}
		cudaStatus = doWithCuda_Par(MAX_Num_Particle,Nside,theta0,1,nmax,allskymap,
			dev_par,particles,dev_rotm,dev_opos);
		
		//cudaStatus = doWithCuda(MAX_Num_Particle,Nside,theta0,1.0,nmax,
		//	allskymap,dev_par,dev_rotm,dev_opos,dev_maplist,host_maplist);
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "Kernel failed!");
			return false;
		}

		//<debug
/*		cudaStatus = cudaMemcpy(allskymap, dev_allskymap, 12 * Nside * Nside, cudaMemcpyDeviceToHost);
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaMemcpy failed!");
			goto Err;
		}
		

		Real sk[100];
		Real a = allskymap[1000];
		int j=0;
		for(int i =0; i < 12*512*512 && j<100; i++){
			if(allskymap[i] !=0 ){
				sk[j] = allskymap[i];
				j++;
			}
		
		}
		a = a;*/
		//debug>
		cout << "block " << _jp << ": "<< (float)_ip / Nparts *100<<"% finished, costs " << (Real)(clock() - time_b) / 1000.0 << 
			" secs, escaped: " << (Real)(clock() - time) / 1000.0 << " secs" << endl;
		if(_jp % rec == 0){
			//cout << "#";
			cout.flush();
		}
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
	return true;

}




