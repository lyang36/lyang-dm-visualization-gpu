#include "structures.h"
#ifndef __READFILES__
#define __READFILES__
using namespace std;

class ReadFiles{
	
	public:
		//keywords
		//read tipsyheader
		string * filename;
		Tipsyheader * tipsyheader;
		bool native;
		long * Nparticles;
		Real * redshift;

		//readparticles
		Allparts *particles;
		bool mass;
		Real * value;
        bool Realpos;
		bool relative_positions;
		Real centerpos[3];


		//read_scarlar

		//function
		void read_tipsyheader();
		void read_particles();
		void read_scalar(float * &); 
};

#endif

