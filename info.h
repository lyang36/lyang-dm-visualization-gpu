#include "structures.h"
#ifndef __INFO__
#define __INFO__

// Info from info file
typedef struct{
	Real Omega_M;
	Real h100;
	Real n_s;
	Real sigma_8;
	Real Lbox_in_Mpc;
	Real particle_masses[10];
	long particle_numbers[10];
	Real centerpos[3];
	Real centervel[3];
	Real Mvir_in_Msun;
	Real M200_in_Msun;
	Real M200crit_in_Msun;
	Real Rvir_in_Mpc;
	Real R200_in_Mpc;
	Real R200crit_in_Mpc;
	Real Vmax_in_kms;
	Real RVmax_in_kpc;
	Real rconverg_in_kpc;
	Real params_NFW[3];
	Real params_GNFW[3];
	Real params_Einasto[3];
}Info;
#endif

