#ifndef __SKYMAP__
#define __SKYMAP__

#include <string>
#include "mapparticle.h"
#include <sys/time.h>
using namespace std;

class Skymap{
private:
	Real rotmatrix[9];
	bool rotate;
    bool reload;
	int num_p; //number of particles in a paritcle block
	MapParticle *particles;
    long clock_get_msec();
    void clock_ini();
    struct timeval initialtime;
    

public:
	long Np;
    long GPU_chunk;  //number of particles in the GPU memory
	//long MAX_Num_Particle; 
	long CPU_chunk;  //number of paritcles in the CPU memory 
	long PRE_chunk;  //number of particles of pre-calculation in GPU memory
                     // with the data
		             //depend on the memory, set up a maximum number of particles that could exist in the memery
	Master *master;
	//Allparts * particles;
	string *fits_filename;
    string *png_filename;
	string *datafile;
    Real *allskymap;
	Map * map;
	Skymap();
	//bool compare(MapParticle a, MapParticle b);
	bool creat_map();
	/*parameters*/
	bool *_rotate;
	bool *_reload;

};


#endif

