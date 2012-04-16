#ifndef __SETPARAMS__
#define __SETPARAMS__
#include "info.h"
using namespace std;


class SetParams{
	Info info;
	Units units;
	Natconst natconst;
	Halo halo;
	Codeunits codeunits;
	Params params;
	Real rotmatrix[3][3];
	Master master;
	Real randomu();
	Real Dsun_in_kpc;

	/*---------------------*/
	public:

	SetParams();
	//keywords
	Real *centerpos;
	Real *centvel;
	Real *observerpos;
	Real *particle_masses;
	Real *particle_numbers;
	Real *redshift;
	Real *scaled;
	Real *length_scale_factor;
	Real *Mvir_in_Msun;
	Real *M200_in_Msun;
	Real *M200crit_in_Msun;
	Real *Rvir_in_Mpc;
	Real *R200_in_Mpc;
	Real *R200crit_in_Mpc;
	Real *Lbox_in_Mpc;
	Real *rconverg_in_kpc;
	Real *shape_radii;
	Real *halo_shape;
	Real *detector;
	Dm *dm;
	Map *map;
	Real *analytical;
	Real *other;
	Real *files;
	////////////////////////////////////////


	//functions
	Master set_master(string * info_filename);
	
	void set_observer_position(Master * master, Real observerpos[3], bool randompos);
	void rotation_matrix(Master * master, Real align_vector[3]);
};

#endif
