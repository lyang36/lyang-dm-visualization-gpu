//#include <stdio>
#include "structures.h"
#ifndef __DOITVL2__
#define __DOITVL2__
class Prepare{
	
private:
	long Npix;
	long Nside;
	bool reload;
    bool master_only;
	Real dOmega;
	Real theta0;
	Map map;
	SetParams * setparams;
	ReadFiles * readfiles;
	Real shapes;
	static Real * observerpos;//[] = {0.36116062, 0.21047031, 0.0069339937};
	static Real * align_vector;//[] = {-0.0265062, -0.285643, -0.957970};
	Tipsyheader header;
	Real redshift;
	Particle particle;
	Allparts * allp;
    //wher find the comtar with less than compare to be true, and edit those value in editar with trueval, otherwise falsval
	long findwhere(Real comtar[], Real lesscompare, Real editar[], Real falsval, int num);
	long findwhere(Real comtar[], Real lesscompare, float editar[], Real falsval, int num);
	

public:
	long Nparticle;
	/***********************/
	//Parameter
	long *_Nside;
	bool *_reload;
	bool *_master_only;
	/***********************/
	Prepare();
	string info_file;//[] = "VL2_info.txt";
	string basedir;//[] = "";
	string basename;//[] = "";
	string projection;//[] = "mollweide";
	void doit_VL2(Master & master, Allparts * allp);
	Real min(Real ar[], int num);
}; 



#endif

