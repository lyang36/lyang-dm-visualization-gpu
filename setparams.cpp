/********************************************************************/
//                Class: SetParams
//				Set up all the parameters
//					Author: Lin Yang
//				     03/08/2012
/********************************************************************/


#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <stdlib.h>
#include <cmath>
//#include "tipsydefs.h"
#include "info.h"
#include "structures.h"
#include "VL2_debug.h"
#include "setparams.h"
using namespace std;

SetParams::SetParams(){
	//keywords
	centerpos = NULL;
	centvel = NULL;
	observerpos = NULL;
	particle_masses = NULL;
	particle_numbers = NULL;
	redshift = NULL;
	scaled = NULL;
	length_scale_factor = NULL;
	Mvir_in_Msun = NULL;
	M200_in_Msun = NULL;
	M200crit_in_Msun = NULL;
	Rvir_in_Mpc = NULL;
	R200_in_Mpc = NULL;
	R200crit_in_Mpc = NULL;
	Lbox_in_Mpc = NULL;
	rconverg_in_kpc = NULL;
	shape_radii = NULL;
	halo_shape = NULL;
	detector = NULL;
	dm = NULL;
	map = NULL;
	analytical = NULL;
	other = NULL;
	files = NULL;
}
Master SetParams::set_master(string * info_filename){
//	return set_master(info_filename, NULL, new Map);
//}

//Master SetParams::set_master(string * info_filename, Real * haloshape, Map * map){
	Master master;
	Halo halo;
	units.Msun_in_cgs = 1.9889212e+33;
	units.Mpc_in_cgs = 3.08568025e+24;
	units.eV_in_cgs = 1.602e-12;

	Real redshift = 0;
	
	const int MAXLENGTH = 300;
	if ((info_filename->length()) == 0){
		cout << "NO INFO FILE" << endl;
	}
	
	Info info;
	//initiate the long 
	for( int i = 0; i<10; i++){
		info.particle_masses[i] = 0;
		info.particle_numbers[i] = 0;
	}
	
	string line;
	ifstream infofile (info_filename -> c_str());
	if (infofile.is_open()){
		while (infofile.good()){
			getline (infofile, line);
			char linestr[MAXLENGTH];
			strcpy (linestr, line.c_str());
			//cout << "OK3" << endl;
			if (line.length() == 0) 
				continue;
			char * toks = strtok(linestr, " ");
			//cout << "why: " << toks << endl;
			if (toks == NULL)
				continue;
			string temp = toks;
			if (toks[0] == '#')
				continue;
			//cout << toks << endl;
			char * val = strtok(NULL, " ");
			val = strtok(NULL, " ");
			//cout << "OK4" << endl;
			if (temp == "Omega_M"){
				info.Omega_M = atof(val);
				//cout << "Omega_M " << info.Omega_M << endl;
			}else if (temp == "h100"){
				info.h100 = atof(val);
				//cout << "h100 " << info.h100<< endl;
			}else if (temp == "n_s"){
                info.n_s = atof(val);
                //cout << "n_s " << info.n_s << endl;
            }else if (temp == "sigma_8"){
                info.sigma_8 = atof(val);
                //cout << "sigma_8 " << info.sigma_8<< endl;
            }else if (temp == "Lbox_in_Mpc"){
                info.Lbox_in_Mpc = atof(val);
                //cout << "Lbox_in_Mpc " << info.Lbox_in_Mpc << endl;
            }else if (temp == "particle_masses"){
				int i = 0;
				do{
					info.particle_masses[i] = atof(val);
#ifdef _DEBUG__LY__
					cout << "VAL:->" << info.particle_masses[i]  << endl;
#endif
                    
					i++;
					val = strtok(NULL, " ");
				}while (val != NULL && i<10);
                                //cout << "particle_masses " << * info.particle_masses << endl;
            }else if (temp == "particle_numbers"){
                int i = 0;
			    do{
                 	info.particle_numbers[i] = atol(val);
					i++;
					val = strtok(NULL, " ");
				}while (val != NULL && i<10);
                                //cout << "particle_numbers" << * info.particle_numbers << endl;
            }else if (temp == "centerpos"){
                info.centerpos[0] = atof(val);
				val = strtok(NULL, " ");
				info.centerpos[1] = atof(val);
				val = strtok(NULL, " ");
				info.centerpos[2] = atof(val);
                                //cout << "centerpos " << * info.centerpos << endl;
			}else if (temp == "centervel"){
                info.centervel[0] = atof(val);
				val = strtok(NULL, " ");
				info.centervel[1] = atof(val);
				val = strtok(NULL, " ");
				info.centervel[2] = atof(val);
					//cout << "centervel " << * info.centervel << endl;
			}else if (temp == "Mvir_in_Msun"){
				info.Mvir_in_Msun = atof(val);
					//cout << "Mvir_in_Msun " << info.Mvir_in_Msun << endl;
			}else if (temp == "M200_in_Msun"){
					info.M200_in_Msun = atof(val);
					//cout << "M200_in_Msun " << info.M200_in_Msun << endl;
			}else if (temp == "M200crit_in_Msun"){
					info.M200crit_in_Msun = atof(val);
					//cout << "M200crit_in_Msun " << info.M200crit_in_Msun << endl;
			}else if (temp == "Rvir_in_Mpc"){
					info.Rvir_in_Mpc = atof(val);
					//cout << "Rvir_in_Mpc " << info.Rvir_in_Mpc << endl;
			}else if (temp == "R200_in_Mpc"){
					info.R200_in_Mpc = atof(val);
					//cout << "R200_in_Mpc " << info.R200_in_Mpc << endl;
			}else if (temp == "R200crit_in_Mpc"){
					info.R200crit_in_Mpc = atof(val);
					//cout << "R200crit_in_Mpc " << info.R200crit_in_Mpc << endl;
			}else if (temp == "Vmax_in_kms"){
					info.Vmax_in_kms = atof(val);
					//cout << "Vmax_in_kms " << info.Vmax_in_kms << endl;
			}else if (temp == "RVmax_in_kpc"){
					info.RVmax_in_kpc = atof(val);
					//cout << "RVmax_in_kpc " << info.RVmax_in_kpc << endl;
			}else if (temp == "rconverg_in_kpc"){
					info.rconverg_in_kpc = atof(val);
					//cout << "rconverg_in_kpc " << info.rconverg_in_kpc << endl;
			}else if (temp == "params_NFW"){
					info.params_NFW[0] = atof(val);
					val = strtok(NULL, " ");
					info.params_NFW[1] = atof(val);
					val = strtok(NULL, " ");
					info.params_NFW[2] = atof(val);
	#ifdef _DEBUG__LY__
					cout << "params_NFW " << info.params_NFW[ 0 ] << endl;
					cout << "params_NFW " << info.params_NFW[ 1 ] << endl;
					cout << "params_NFW " << info.params_NFW[ 2 ] << endl;
	#endif
			}else if (temp == "params_GNFW"){
					info.params_GNFW[0] = atof(val);
					val = strtok(NULL, " ");
					info.params_GNFW[1] = atof(val);
					val = strtok(NULL, " ");
					info.params_GNFW[2] = atof(val);
					//cout << "params_GNFW " << * info.params_GNFW << endl;
			}else if (temp == "params_Einasto"){
					info.params_Einasto[0] = atof(val);
					val = strtok(NULL, " ");
					info.params_Einasto[1] = atof(val);
					val = strtok(NULL, " ");
					info.params_Einasto[2] = atof(val);
					//cout << "params_Einasto " << * info.params_Einasto << endl;
			}
		}
	}
 
	//no scale

	infofile.close();
	natconst.h100 = info.h100;
	natconst.rho_crit_in_cgs = 1.8788309e-29 * pow(info.h100, 2);
	natconst.Omega_m = info.Omega_M;
	natconst.Rvir_MW_in_Mpc = 0.275;
	natconst.c_in_cgs = 2.99792458e10;
	natconst.G_in_cgs = 6.67259e-8;

	// definition of the virial overdensity (Bryan & Norman 1998)  (rho/rho_0)
  	Real Omega_m_at_z = natconst.Omega_m * pow((1.0 + redshift), 3) / (natconst.Omega_m * pow((1.0 + redshift), 3) + (1.0 - natconst.Omega_m));
  	Real x = Omega_m_at_z - 1.0;
  	natconst.Delta_vir = (18.0 * PI * PI + 82.0 * x - 39.0 * x * x) / (1.0 + x);
  	natconst.rho_crit_in_cgs = 3.0 * pow((100.0 * natconst.h100/(units.Mpc_in_cgs * 1e-5)), 2) / (8 * PI * natconst.G_in_cgs);


	//no keywords for dm
	if( dm == NULL){
		dm = new Dm;
		dm -> M_in_GeV = 46.0;
		dm -> M_in_cgs = 0.0;
		dm -> sigma_v_in_cgs = 5.0e-26;
	}

	dm -> M_in_cgs = dm -> M_in_GeV * 1e9 * units.eV_in_cgs / pow(natconst.c_in_cgs, 2);

	Real shape_radii = 0.0;
	Real halo_shape =0.0;	

	halo.Mvir_in_Msun = info.Mvir_in_Msun;
    halo.M200_in_Msun = info.M200_in_Msun;
    halo.M200crit_in_Msun = info.M200crit_in_Msun;
    halo.Rvir_in_Mpc = info.Rvir_in_Mpc;
    halo.R200_in_Mpc = info.R200_in_Mpc;
    halo.R200crit_in_Mpc = info.R200crit_in_Mpc;
    halo.Vmax_in_kms = info.Vmax_in_kms;
    halo.RVmax_in_kpc = info.RVmax_in_kpc;
    halo.rconverg_in_kpc = info.rconverg_in_kpc;
	for(int i = 0; i < 3; i++){
		halo.params_NFW[i] = info.params_NFW[i];
		halo.params_GNFW[i] = info.params_GNFW[i];
		halo.params_Einasto[i] = info.params_Einasto[i];
	}
    halo.shape_r = shape_radii;
    halo.shape = halo_shape;

    codeunits.length_to_Mpc = info.Lbox_in_Mpc;
    codeunits.length_to_cgs = info.Lbox_in_Mpc*units.Mpc_in_cgs;
    codeunits.time_to_cgs = pow((natconst.rho_crit_in_cgs * natconst.G_in_cgs), (-0.5));
    codeunits.density_to_cgs = natconst.rho_crit_in_cgs;

  	codeunits.mass_to_cgs = codeunits.density_to_cgs * codeunits.length_to_cgs * codeunits.length_to_cgs *codeunits.length_to_cgs;
  	codeunits.mass_to_Msun = codeunits.mass_to_cgs / units.Msun_in_cgs;
  	codeunits.velocity_to_cgs = codeunits.length_to_cgs /codeunits.time_to_cgs;
	codeunits.annihilation_flux_to_cgs = (codeunits.density_to_cgs * codeunits.mass_to_cgs / pow(codeunits.length_to_cgs, 2));


	//no setting for observerpos
	Real observerpos[] = {0, 0, 0};
	params.z = redshift;
	for(int i = 0; i < 3; i++){
		params.cpos[i] = info.centerpos[i];
		params.cvel[i] = info.centervel[i];
		params.opos[i] = observerpos[i];
	}
	params.otheta = 0.0;
	params.ophi = 0.0;
	params.Lbox_in_Mpc = info.Lbox_in_Mpc;
	for( int i = 0, j = 0; i < 10; i++){
		params.particle_numbers[i] = 0;
		params.particle_masses[i] = 0;
		params.particle_masses_in_Msun[i] = 0;
		if (info.particle_masses[i] > 0){
#ifdef _DEBUG__LY__
			cout << "INFO: " << i << " " << j << " "<<info.particle_masses[i] << endl;
#endif
		    params.particle_numbers[j] = info.particle_numbers[i];
    		params.particle_masses[j] = info.particle_masses[i];
    		params.particle_masses_in_Msun[j] = info.particle_masses[i] * codeunits.mass_to_cgs / units.Msun_in_cgs;
			j++;
		}
	}

	

	Real other = 0;
	Real analytical =0;
	Real files = 0;
	Real detector = 0;
	//Real map = 0;
	
	print_out_halo( &halo );
	static Real rotmatrix[][3] = {1, 0, 0, 0, 1, 0, 0, 0, 1};
	master.rotmatrix = rotmatrix;
	master.natconst = natconst;
    master.dm = *dm;
    master.halo = halo;
    master.params = params;
    master.units = units;
    master.codeunits = codeunits;
    master.detector = detector;
	if( map == NULL){
		map = new Map;
		
	}
	master.map = (* map);
    master.files = files;
    master.analytical = analytical;
    master.other = other;
	return master;
}



void SetParams::set_observer_position(Master * master, Real observerpos[3], bool randompos){
  	// calculate ophi and otheta from observerpos and centerpos

  	master -> params.opos[0] = observerpos[0];
	master -> params.opos[1] = observerpos[1];
	master -> params.opos[2] = observerpos[2];

  	Real radius = sqrt(pow(master -> params.opos[0] - master -> params.cpos[0], (Real)2.0)
			+ pow(master -> params.opos[1] - master -> params.cpos[1], (Real)2.0)
			+ pow(master -> params.opos[2] - master -> params.cpos[2], (Real)2.0));
	//cout << "R = " << radius <<endl;
	//cout << "cpos = " << master -> params.cpos[0] <<endl;
  	Real otheta = acos((master -> params.opos[2] - master -> params.cpos[2])/radius );
  	Real ophi = -PI + atan((master -> params.opos[1]-master -> params.cpos[1]) / (master -> params.opos[0]-master -> params.cpos[0]));
  	if (ophi < 0.0)
		ophi += 2.0 * PI;
	//cout << "Params Otheta:" << otheta << endl;

  	master -> params.otheta = otheta;
  	master -> params.ophi = ophi;
	
	//cout << "theta+phi" << otheta << " " << ophi << endl;	
	return;

}

void SetParams::rotation_matrix(Master * master, Real align_vector[3]){
	Real otheta = master -> params.otheta;
  	Real ophi = master -> params.ophi;

  	Real rotaxis[3];
	rotaxis[0] = 0.0;
  	rotaxis[1] = cos(otheta) / sqrt(pow((sin(otheta) * sin(ophi)), (Real)2.0) + pow(cos(otheta), (Real)2));
  	rotaxis[2] = -sin(otheta) * sin(ophi) / sqrt(pow((sin(otheta)*sin(ophi)), (Real)2) + pow(cos(otheta), (Real)2));

	Real c = -sin(otheta) * cos(ophi);
  	Real s = -sqrt(pow(sin(otheta) * sin(ophi), 2)+ pow(cos(otheta), 2));
  	Real t = 1.0 - c;
	
	//cout << "OTHETA" << otheta <<endl;
	//cout << "ophi" << ophi <<endl;
	//cout << "rotaxis "<< rotaxis[1] <<endl;

  	rotmatrix[0][0] = t * pow(rotaxis[0], 2) + c;
  	rotmatrix[1][0] = t * rotaxis[0] * rotaxis[1] - s * rotaxis[2];
  	rotmatrix[2][0] = t * rotaxis[0] * rotaxis[2] + s * rotaxis[1];

  	rotmatrix[0][1] = t * rotaxis[0] * rotaxis[1] + s * rotaxis[2];
  	rotmatrix[1][1] = t * pow(rotaxis[1], 2) + c;
  	rotmatrix[2][1] = t * rotaxis[1] * rotaxis[2] - s * rotaxis[0];

 	rotmatrix[0][2] = t * rotaxis[0] * rotaxis[2] - s * rotaxis[1];
  	rotmatrix[1][2] = t * rotaxis[1] * rotaxis[2] + s * rotaxis[0];
  	rotmatrix[2][2] = t * pow(rotaxis[2], 2) + c;

	//cout << "ROT[0][0] " << rotmatrix[0][0] << endl;

  	Real distance = sqrt(pow(master -> params.opos[0] - master -> params.cpos[0], (Real)2.0)
			+ pow(master -> params.opos[1] - master -> params.cpos[1], (Real)2.0)
			+ pow(master -> params.opos[2] - master -> params.cpos[2], (Real)2.0));
  	Real tmpvec[] = {master -> params.cpos[0] + distance * align_vector[0] - master -> params.opos[0], 
				master -> params.cpos[1] + distance * align_vector[1] - master -> params.opos[1], 
				master -> params.cpos[2] + distance * align_vector[2] - master -> params.opos[2]};

  	Real xtmp = tmpvec[0] * rotmatrix[0][0] + 
        	tmpvec[1] * rotmatrix[1][0] + 
        	tmpvec[2] * rotmatrix[2][0];

  	Real ytmp = tmpvec[0] * rotmatrix[0][1] + 
        	tmpvec[1] * rotmatrix[1][1] + 
        	tmpvec[2] * rotmatrix[2][1];

  	Real ztmp = tmpvec[0] * rotmatrix[0][2] + 
         	tmpvec[1] * rotmatrix[1][2] + 
         	tmpvec[2] * rotmatrix[2][2];

  	tmpvec[0] = xtmp;
	tmpvec[1] = ytmp;
	tmpvec[2] = ztmp;
  	Real gamma_t = fabs(atan(ztmp / ytmp));

  	if (gamma_t > PI/4.0)
		gamma_t *= -1.0;
	//cout << "GAMMA = " << gamma_t <<endl;

  	Real rotmatrix2[][3] = { 1.0, 0.0, 0.0,
				0.0, cos(gamma_t), sin(gamma_t),
				0.0, -sin(gamma_t), cos(gamma_t)};
	static Real rottmp[3][3];	
	for (int i = 0; i<3; i++){
		for (int j =0; j<3; j++){
			Real val = 0;
			for (int k=0; k<3; k++)
				val += rotmatrix[i][k] * rotmatrix2[k][j];
			
			rottmp[i][j] = val;
		}
	}
	
	master -> rotmatrix = rottmp;
 	return;                       
}

