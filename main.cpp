//#include <stdio>
#include <iostream>
#include <sstream>
#include <fstream>
#include <cmath>
#include <string>
#include "structures.h"
#include "setparams.h"
#include "skymap.h"



using namespace std;

int main(int argc, const char **argv){
	int m=1;
	string datafile = "/home/gpuclass4/data/data_new.bin";
        string info_file =  "/home/gpuclass4/data/VL2_info.txt";
        string fits_file = "/home/gpuclass4/data/skymap.fits";
	long Nside = 512;
	long Npix = 12 * Nside * Nside;
	Real dOmega = 4.0 * M_PI / Npix;
	Real 	theta0 = acos(1.0 - dOmega / (2.0 * PI));
	Real observerpos_t[] = {0.36116062, 0.21047031, 0.0069339937};
	Real align_vector_t[] = {-0.0265062, -0.285643, -0.957970};


	while (m<argc)
	{
		string arg = argv[m];
		if (arg == "-data") { datafile = argv[m+1]; m+=2;}
		if (arg == "-info") { info_file = argv[m+1]; m+=2;}
		if (arg == "-fits") { fits_file = argv[m+1]; m+=2;}
		if (arg == "-ob"){
			Real inp;
			cout << "Input observer X:" << endl;
			cin >> inp;
			observerpos_t[0] = inp;
			cout << "Input observer Y:" << endl;
			cin >> inp;
			observerpos_t[1] = inp;
			cout << "Input observer Z:" << endl;
			cin >> inp;
			observerpos_t[2] = inp;
		}
		if (arg == "-al"){
			Real the;
			Real phi;
			cout << "Input observer align theta in degree:" << endl;
			cin >> the;
			cout << "Input observer align phi in degree:" << endl;
			cin >> phi;
			the = the / 180.0 * PI;
			phi = phi / 180.0 * PI;
			align_vector_t[0] = sin(the)*cos(phi);//{-0.0265062, -0.285643, -0.957970};
			align_vector_t[1] = sin(the)*sin(phi);
			align_vector_t[2] = cos(the);
		}
	}
	cout << "BASE DIR: " << datafile << endl;
	cout << "INFO FILE: " << info_file << endl;
	cout << "OUTPUT FITS: " << fits_file << endl;

/*----------------------------------------------------*/
	Skymap * skymap = new Skymap();
/*----------------------------------------------------*/
	Master master;

/*----------prepare for map---------------------------*/
	Map map;
	map.projection = "mollweide";
	map.Nside = Nside;
	map.Npix = Npix;
	map.dOmega = dOmega;
	map.theta0 = theta0;

/*--------------set master structure-------------------------------------------------------*/
	SetParams * setparams = new SetParams();
	//master = new Master;
	setparams -> map = & map;
	
	Real shapes = 0;
	setparams -> halo_shape = & shapes;
	master = (*setparams).set_master(& info_file);
	(*setparams).set_observer_position( & master, observerpos_t, false);
	
	(*setparams).rotation_matrix( & master, align_vector_t);
	//cout <<"MASTER:   "<< master -> params.Lbox_in_Mpc << endl;
/*-----------------------------------------------------------------------------------------*/

//*--------------generating fits map------------------*/
	skymap->MAX_Num_Particle = 16 * 1024;
	skymap->CPU_trunk = 1024 / 64 * 1024 * 1024; //~1Gb 
	skymap->PRE_trunk = skymap->CPU_trunk / 64; 
	skymap->fits_filename = & fits_file;
	skymap->master = &master;
	skymap->datafile = &datafile;
	skymap->map = &(master.map);
	bool rotate = true;
	skymap->_rotate = &rotate; 
	skymap->creat_map();
	getchar();
/*----------------------------------------------------*/
	return 0;	
}
