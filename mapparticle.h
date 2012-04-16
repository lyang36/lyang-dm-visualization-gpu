#include "structures.h"

#ifndef __MAP_PARTICLE_STRUCTURE__
#define __MAP_PARTICLE_STRUCTURE__

struct MapParticle{
public:
	Real mass;
	Real density;
	Real hsmooth;
	Real xpos;
	Real ypos;
	Real zpos;
/*	bool operator()(const MapParticle &a, const MapParticle &b) const{
		return a.zpos > b.zpos;
	}*/
	bool operator>(const MapParticle &b) const{
		return hsmooth > b.hsmooth;//hsmooth > b.hsmooth
			/*(xpos-0.36116062)*(xpos-0.36116062) + 
				(ypos-0.21047031)*(ypos-0.21047031) + 
				(zpos-0.0069339937)*(zpos-0.0069339937) > 
				(b.xpos-0.36116062)*(b.xpos-0.36116062) + 
				(b.ypos-0.21047031)*(b.ypos-0.21047031) + 
				(b.zpos-0.0069339937)*(b.zpos-0.0069339937)*/
	}
};


#endif