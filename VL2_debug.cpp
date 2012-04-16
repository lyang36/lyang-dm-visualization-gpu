#include <iostream>
#include <sstream>
#include <string>
//#include "tipsydefs.h"
#include "info.h"
#include "structures.h"
#include "VL2_debug.h"

using namespace std;


string getstring_Real( double * d, long numx){
        #ifdef _DEBUG__LY__
                //stringstream ss(stringstream::in | stringstream::out);
                string oputs;
                string temps;
                for(long i = 0; i< numx; i++){
                        stringstream ss;
                        ss << d[i];
                        //cout << d[i] << endl;
                        ss >> temps;
                        oputs += temps + " ";
                }
                //ss >> oputs;
                return oputs;
        #endif
                return "";
}
string getstring_long( long * d, long numx){
        #ifdef _DEBUG__LY__
                //stringstream ss(stringstream::in | stringstream::out);
                string oputs;
                string temps;
                for(long i = 0; i< numx; i++){
                        stringstream ss;
                        ss << d[i];
                        //cout << d[i] << endl;
                        ss >> temps;
                        oputs += temps + " ";
                }
                //ss >> oputs;
                return oputs;
        #endif
                return "";
}


void print_out_natconst( Natconst * natconst){
	#ifdef _DEBUG__LY__
	cout << "Out put Natconst:-----------------------------------------------------" << endl;
		cout << "h100 " << natconst -> h100 << endl;
		cout << "G_in_cgs " << natconst -> G_in_cgs << endl;
		cout << "rho_crit_in_cgs " << natconst -> rho_crit_in_cgs << endl;
		cout << "Omega_m " << natconst -> Omega_m << endl;
		cout << "Delta_vir " << natconst -> Delta_vir << endl;
		cout << "Rvir_MW_in_Mpc " << natconst -> Rvir_MW_in_Mpc << endl;
		cout << "c_in_cgs " << natconst -> c_in_cgs << endl;
	#endif
}

void print_out_halo( Halo * halo){
	#ifdef _DEBUG__LY__
	cout << "Out put Halo:-----------------------------------------------------" << endl;
	cout << "Mvir_in_Msun " << halo -> Mvir_in_Msun << endl;
	cout << "M200_in_Msun " << halo -> M200_in_Msun << endl;
	cout << "M200crit_in_Msun " << halo -> M200crit_in_Msun << endl;
	cout << "Rvir_in_Mpc " << halo -> Rvir_in_Mpc << endl;
	cout << "R200_in_Mpc " << halo -> R200_in_Mpc << endl;
	cout << "R200crit_in_Mpc " << halo -> R200crit_in_Mpc << endl;
	cout << "Vmax_in_kms " << halo -> Vmax_in_kms << endl;
	cout << "RVmax_in_kpc " << halo -> RVmax_in_kpc << endl;
	cout << "rconverg_in_kpc " << halo -> rconverg_in_kpc << endl;
	cout << "params_NFW " << getstring_Real( halo -> params_NFW, 3) << endl;
	cout << "params_GNFW " << getstring_Real(  halo -> params_GNFW, 3) << endl;
	cout << "params_Einasto " << getstring_Real( halo -> params_Einasto, 3) << endl;
	cout << "shape_r " << halo -> shape_r << endl;
	cout << "shape " << halo -> shape << endl;
	#endif
}

void print_out_dm( Dm * dm){
	#ifdef _DEBUG__LY__
	cout << "DM:-----------------------------------------------------" << endl;
	cout << "M_in_GeV " << dm -> M_in_GeV << endl;// = 46.0;
	cout << "M_in_cgs " << dm -> M_in_cgs << endl;// = 0.0;
	cout << "sigma_v_in_cgs "<< dm -> sigma_v_in_cgs << endl;// = 5.0e-26;

	#endif
}

void print_out_codeunits( Codeunits * codeunits){
	#ifdef _DEBUG__LY__
	cout << "Codeunits:-----------------------------------------------------" << endl;
	cout << "mass_to_cgs " <<  codeunits -> mass_to_cgs << endl;
	cout << "mass_to_Msun " <<  codeunits -> mass_to_Msun << endl;
	cout << "length_to_Mpc " << codeunits -> length_to_Mpc << endl;
	cout << "length_to_cgs " << codeunits -> length_to_cgs << endl;
	cout << "time_to_cgs " << codeunits -> time_to_cgs << endl;
	cout << "velocity_to_cgs " << codeunits ->  velocity_to_cgs << endl;
	cout << "density_to_cgs " << codeunits -> density_to_cgs << endl;
	cout << "annihilation_flus_to_cgs " << codeunits -> annihilation_flux_to_cgs << endl;

	#endif
}


void print_out_params( Params * params){
	#ifdef _DEBUG__LY__
	cout << "Params:-----------------------------------------------------" << endl;
	cout << "z " <<  params -> z << endl;
//	cout << "cpos " << getstring_Real(params -> cpos, 3) << endl;
//	cout << "cvel " << getstring_Real(params -> cvel, 3) << endl;//cvel[3];
//	cout << "opos " << getstring_Real(params -> opos, 3) << endl; // opos[3];
	cout << "otheta " << params -> otheta << endl;
	cout << "ophi " << params -> ophi << endl;//;
	cout << "Lbox_in_Mpc " << params -> Lbox_in_Mpc << endl;//;
	
	//cout << "particle_numbers " <<  getstring_long( (params -> particle_numbers), 10) << endl;//[10];
	//cout << "particle_masses " <<  getstring_Real( (params -> particle_masses), 10) << endl;//[10];
	//cout << "particle_masses_in_Msun " << getstring_Real( (params -> particle_masses_in_Msun), 10) << endl;// [10];

	#endif
}


void print_out_master( Master * master){
	#ifdef _DEBUG__LY__
	cout << "Master:-----------------------------------------------------" << endl;
	print_out_natconst( &(master -> natconst) );
	print_out_dm( &(master -> dm) );
	print_out_halo( &(master -> halo) );
	print_out_params( &(master -> params) );
	//print_out_units units;
	print_out_codeunits( &(master -> codeunits));
	print_out_map( &(master -> map));
	cout << "detector " << master -> detector << endl;
	//cout << "map " << master -> map << endl;
	cout << "files " << master -> files << endl;
	cout << "analytical " << master -> analytical << endl;
	cout << "other " << master -> other << endl;
//	cout << "rotmatrix " << getstring_Real( master -> rotmatrix[0], 3) << endl;//[3][3];
//	cout << "rotmatrix " << getstring_Real( master -> rotmatrix[1], 3) << endl;
//	cout << "rotmatrix " << getstring_Real( master -> rotmatrix[2], 3) << endl;
	#endif
}
/*
void print_out_tipsyheader( Tipsyheader * tipsyheader){
	#ifdef _DEBUG__LY__
	cout << "Tipsyheader:-----------------------------------------------------" << endl;
	cout << "time " << tipsyheader -> time << endl;
	cout << "nboides " << tipsyheader -> nbodies << endl;//;
	cout << "ndim " << tipsyheader -> ndim << endl;//;
	cout << "nsph " << tipsyheader -> nsph << endl;//;
	cout << "ndark " << tipsyheader -> ndark << endl;//;
	cout << "nstar " << tipsyheader -> nstar << endl;//;
//	cout << "dummy " << tipsyheader -> dummy << endl;//;
	#endif

}
*/
void print_out_particle( Particle * particle){
	#ifdef _DEBUG__LY__
	cout << "Particle:-----------------------------------------------------" << endl;
	cout << "mass " << particle -> mass << endl;//;
	cout << "density " << particle -> density << endl;//;
	cout << "hsmooth " << particle -> hsmooth << endl;//;
	cout << "xpos " << particle -> xpos << endl;//;
	cout << "ypos " << particle -> ypos << endl;//;
	cout << "zpos " << particle -> zpos << endl;//;

	#endif
}


void print_out_map( Map * map){
	#ifdef _DEBUG__LY__
	cout << "map:-----------------------------------------------------" << endl;
	cout << "projection " << map -> projection << endl;//;
	cout << "Nside " << map -> Nside << endl;//;
	cout << "Npix " << map -> Npix << endl;//;
	cout << "dOmega " << map -> dOmega << endl;//;
	cout << "theta0 " << map -> theta0 << endl;//;

	#endif
}
/*
void print_out_allparts( Allparts * allparts){
	#ifdef _DEBUG__LY__

	Real * mass;
	Real * density;
	Real * hsmooth;
	Real * xpos;
	Real * ypos;
	Real * zpos;

	#endif
}
*/
