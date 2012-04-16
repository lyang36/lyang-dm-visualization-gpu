#ifndef __STRUCTURES__
#define __STRUCTURES__

#include <string>
//#define POS_FLOAT
#ifdef POS_FLOAT
  #define Real float
#else
  #define Real double
#endif

#ifndef M_PI
#define M_PI           3.14159265358979323846
#endif

const Real PI  = M_PI; 

struct Units{//static const 
	double Msun_in_cgs;// = 1.9889212e+33;
	double Mpc_in_cgs;// = 3.08568025e+24;
	double eV_in_cgs;// = 1.602e-12;
};

typedef struct{
	double h100;
	double G_in_cgs;// = 6.67259e-8;
	double rho_crit_in_cgs;
	double Omega_m;
	double Delta_vir;
    double Rvir_MW_in_Mpc;// = 0.275;
	double c_in_cgs;// = 2.99792458e10;
} Natconst;

typedef struct{
	double Mvir_in_Msun;
	double M200_in_Msun;
	double M200crit_in_Msun;
	double Rvir_in_Mpc;
	double R200_in_Mpc;
	double R200crit_in_Mpc;
	double Vmax_in_kms;
	double RVmax_in_kpc;
	double rconverg_in_kpc;
	double params_NFW[3];
	double params_GNFW[3];
	double params_Einasto[3];
	double shape_r;
	double shape;
} Halo;

typedef struct{
	double M_in_GeV;// = 46.0;
	double M_in_cgs;// = 0.0;
	double sigma_v_in_cgs;// = 5.0e-26;
	
} Dm;

typedef struct{
	double mass_to_cgs;
	double mass_to_Msun;
	double length_to_Mpc;
	double length_to_cgs;
	double time_to_cgs;
	double velocity_to_cgs;
	double density_to_cgs;
	double annihilation_flux_to_cgs;
} Codeunits;

typedef struct{
	Real z;
	Real cpos[3];
	Real cvel[3];
	Real opos[3];
	Real otheta;
	Real ophi;
	Real Lbox_in_Mpc;
	long particle_numbers[10];
	Real particle_masses[10];
	Real particle_masses_in_Msun[10];
} Params;


//typedef tipsy_header Tipsyheader;

typedef struct{
	Real mass;
	Real density;
	Real hsmooth;
	Real xpos;
	Real ypos;
	Real zpos;
} Particle;

typedef struct{
	std::string projection;
	long Nside;
	long Npix;
	Real dOmega;
	Real theta0;
}Map;

typedef struct{
	Real * mass;
	Real * density;
	Real * hsmooth;
	Real * xpos;
	Real * ypos;
	Real * zpos;
	Real * distance;
} Allparts;

typedef struct{
	Natconst natconst;
	Dm dm;
	Halo halo;
	Params params;
	Units units;
	Codeunits codeunits;
	Real detector;
	Map map;
	Real files;
	Real analytical;
	Real other;
	Real (*rotmatrix)[3];//[3][3];	
} Master;

#endif

