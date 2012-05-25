
__global__
void generate_map_step0(const long Nside, 
                        const Real theta0, const long nmax,
                        MapParticle * particles,
                        Real * rotmatrix, Real * opos,
                        int * pd_key, int * pd_val, int stpoint){
    //fluxfactor = master->codeunits.annihilation_flux_to_cgs 
	//int i = threadIdx.x;
	int i = blockIdx.x * blockDim.x + threadIdx.x;
	if( i >= nmax){
		return;
	}
	MapParticle mp = particles[i];
	Real density = mp.density;
	Real xpos = mp.xpos;
	Real ypos = mp.ypos;
	Real zpos = mp.zpos;
	Real mass = mp.mass;
	Real hsmooth = mp.hsmooth;
	Real distances;
	Real fluxes;
	Real costheta;
	Real phi;
	Real theta;
	Real angular_radius;
    
#ifdef POS_FLOAT
	float3 _opos = {opos[0], opos[1], opos[2]};
#else
	double3 _opos = {opos[0], opos[1], opos[2]};
#endif
    
	/*since if density <0, the threads becomes null*/
	/*please make sure all the density is greater or equal to 0*/
	/*this is very important*/
	if( density < 0.0){
		particles[i].xpos = 0;
		//xpos holds pixel num
		particles[i].hsmooth = 0;//;angular_radius;
		particles[i].ypos = 0;//theta;
		particles[i].zpos = 0;//phi;
		particles[i].mass = 0;//flux;
		//particles[i].density = weight_norm;
		pd_val[i] = stpoint + i;
		pd_key[i] = 0;
	}else if( density >= 0.0){
		distances = sqrt( (xpos-_opos.x) * (xpos-_opos.x) + 
                         (ypos-_opos.y) * (ypos-_opos.y) +
                         (zpos-_opos.z) *(zpos-_opos.z) );
		fluxes = density * mass / (4.0 * PI * distances * distances);
        
		calc_angles( xpos-_opos.x, ypos-_opos.y, zpos-_opos.z, distances, 
                    opos, rotmatrix, costheta, phi);
		theta = acos(costheta);
		if(distances != 0.0){
			angular_radius = hsmooth / distances;}
		else
        {angular_radius = PI;}
		//hsmooth holds angular_radius
		particles[i].hsmooth = angular_radius;
		particles[i].ypos = theta;
		particles[i].zpos = phi;
		particles[i].mass = fluxes;
		//delete step1
		Real ad =  2 * angular_radius * (1 + 0.1 );
		int pixss = floor(ad * ad * ((double)(12L * Nside * Nside)) / 4) + 1; //sort by anglar_radius
		particles[i].xpos = pixss;
		particles[i].density = 1;//2.0 * angular_radius * angular_radius * 12 * 512 * 512 / 4 / PI;
		pd_key[i] = pixss;
		pd_val[i] = stpoint + i;
		//density holds norm 1.62707
	}
}




__global__
void generate_map_step1(const long Nside, 
                        const Real theta0, const long nmax,
                        MapParticle * particles,
                        Real * rotmatrix, Real * opos){
    //fluxfactor = master->codeunits.annihilation_flux_to_cgs 
	//int i = threadIdx.x;
	int i = blockIdx.x * blockDim.x + threadIdx.x;
	if( i >= nmax){
		return;
	}
	MapParticle mp = particles[i];
	Real density = mp.density;
	Real hsmooth = mp.hsmooth;
	Real fluxes = mp.mass;
	Real phi = mp.zpos;
	Real theta = mp.ypos;
	Real angular_radius = mp.hsmooth;
    
#ifdef POS_FLOAT
	float3 _opos = {opos[0], opos[1], opos[2]};
#else
	double3 _opos = {opos[0], opos[1], opos[2]};
#endif
    
	if( density >= 0.0){
		//density holds norm
        /*------------------------------------------add particles------------------------------------------------------*/
		{
			//pointing p (theta,phi);
			Real vecx = sin(theta) * cos(phi);
			Real vecy = sin(theta) * sin(phi);
			Real vecz = cos(theta);
			if( 2.0*angular_radius < theta0 ) {
				/*should be atomic plus */
				//allskymap[pix] += fluxes;
				particles[i].xpos = 1;
				particles[i].density = 1;
			}else{
				Real weight_norm = 0.0;
				int Npix = 0;
				query_disc_getnorm(Nside, theta, phi, 2.0*angular_radius, vecx, vecy, vecz, weight_norm, Npix);
				particles[i].xpos = Npix;
				particles[i].density = weight_norm;
			}
		}
	}
}




__global__
void generate_map_step2(const long Nside, 
                        const Real theta0, const long nmax,
                        MapParticle * particles, mapMap *maplist,
                        Real * rotmatrix, Real * opos){
    
	int i = blockIdx.x * blockDim.x + threadIdx.x;
    
	if( i >= nmax){
		return;
	}
	MapParticle mp = particles[i];
	Real fluxes = mp.mass;
	int startp = (int)mp.xpos;
	Real theta = mp.ypos;
	Real phi = mp.zpos;
	Real angular_radius = mp.hsmooth;
	Real weight_norm = mp.density;
	Real density = mp.density;
    
	/*since if density <0, the threads becomes null*/
	/*please make sure all the density is greater or equal to 0*/
	/*this is very important*/
	if( density >= 0.0){
        
		Real costheta = cos(theta);
        /*------------------------------------------add particles------------------------------------------------------*/
		{
			//pointing p (theta,phi);
			Real vecx = sin(theta) * cos(phi);
			Real vecy = sin(theta) * sin(phi);
			Real vecz = costheta;
			//ang2vec(theta,phi,&vec);
            
			long pix;
			ang2pix_ring(Nside,theta,phi,&pix);
			if( 2.0*angular_radius < theta0 ) {
				/*should be atomic plus */
				//allskymap[pix] += fluxes;
				maplist[startp].pix = pix;
				maplist[startp].signal = fluxes;
				maplist[startp].factor = 1;
			}else{
				long index = startp;
				double getnorm = 0;
				query_disc_calcflux(Nside, theta, phi, 2.0*angular_radius, vecx, vecy, vecz, 
                                    weight_norm, fluxes, maplist, index, getnorm);
				maplist[startp].factor = getnorm;
				if(getnorm == 0){ // consider Npix <2 --> keep the flux conservative
					maplist[startp].pix = pix;
					maplist[startp].signal = fluxes;
					maplist[startp].factor = 1;
				}
			}
		}
	}
}
