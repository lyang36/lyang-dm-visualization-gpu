

__device__ __host__
void query_disc_getnorm(long nside_, const Real theta, const Real phi, 
                        const Real radius, const Real vecx,  const Real vecy, const Real vecz, 
                        Real &weight_norm, int &Npix){
	Npix = 0;
	weight_norm = 0;
	long npface_ = nside_ * nside_;
    //	long ncap_   = (npface_-nside_)<<1;
	long npix_   = 12 * npface_;
	Real fact2_  = 4. / (Real) npix_;
	Real fact1_  = (nside_<<1) * (Real) fact2_;
    
	Real dth1 = fact2_;
	Real dth2 = fact1_;
	Real cosang = cos(radius);
    
	Real z0 = cos(theta);
	Real xa = 1./sqrt((1-z0)*(1+z0));
    
	Real rlat1  = theta - radius;
	Real zmax = cos(rlat1);
	int irmin = ring_above (nside_, zmax)+1;
    
	if (rlat1<=0) // north pole in the disc
        for (int m=1; m<irmin; ++m) // rings completely in the disc
            in_ring_count (nside_, m, 0, pi, vecx, vecy, vecz, radius / 2.0, weight_norm, Npix);
    
	Real rlat2  = theta + radius;
	Real zmin = cos(rlat2);
	int irmax = ring_above (nside_, zmin);
    
	// ------------- loop on ring number ---------------------
	for (int iz=irmin; iz<=irmax; ++iz) // rings partially in the disc
	{
		Real z;
		if (iz<nside_) // north polar cap
            z = 1.0 - iz*iz*dth1;
		else if (iz <= (3*nside_)) // tropical band + equat.
            z = (2*nside_-iz) * dth2;
		else
            z = -1.0 + (4*nside_-iz)*(4*nside_-iz)*dth1;
        
		// --------- phi range in the disc for each z ---------
		Real x = (cosang-z*z0)*xa;
		Real ysq = 1-z*z-x*x;
		//planck_assert(ysq>=0, "error in query_disc()");
		Real dphi=atan2(sqrt(ysq),x);
		in_ring_count (nside_, iz, phi, dphi, vecx, vecy, vecz, radius / 2.0, weight_norm, Npix);
	}
    
	if (rlat2>=pi) // south pole in the disc
		for (int m=irmax+1; m<(4*nside_); ++m)  // rings completely in the disc
			in_ring_count (nside_, m, 0, pi,  vecx, vecy, vecz, radius / 2.0, weight_norm, Npix);
}

__device__
void query_disc_calcflux(long nside_, const Real theta, const Real phi, const Real radius, 
                         const Real vecx,  const Real vecy, const Real vecz, const Real weight_norm,
                         const Real aflux, mapMap *maplist, long &indexp, double &getnorm){
	long npface_ = nside_ * nside_;
    //	long ncap_   = (npface_-nside_)<<1;
	long npix_   = 12 * npface_;
	Real fact2_  = 4. / (Real) npix_;
	Real fact1_  = (nside_<<1) * (Real) fact2_;
    
	Real dth1 = fact2_;
	Real dth2 = fact1_;
	Real cosang = cos(radius);
    
	Real z0 = cos(theta);
	Real xa = 1./sqrt((1-z0)*(1+z0));
    
	Real rlat1  = theta - radius;
	Real zmax = cos(rlat1);
	int irmin = ring_above (nside_, zmax)+1;
    
	if (rlat1<=0) // north pole in the disc
        for (int m=1; m<irmin; ++m) // rings completely in the disc
            in_ring_allskymap (nside_, m, 0, pi, vecx, vecy, vecz, radius / 2.0, 
                               weight_norm, aflux, maplist, indexp, getnorm);
    
	Real rlat2  = theta + radius;
	Real zmin = cos(rlat2);
	int irmax = ring_above (nside_, zmin);
    
	// ------------- loop on ring number ---------------------
	for (int iz=irmin; iz<=irmax; ++iz) // rings partially in the disc
	{
		Real z;
		if (iz<nside_) // north polar cap
            z = 1.0 - iz*iz*dth1;
		else if (iz <= (3*nside_)) // tropical band + equat.
            z = (2*nside_-iz) * dth2;
		else
            z = -1.0 + (4*nside_-iz)*(4*nside_-iz)*dth1;
        
		// --------- phi range in the disc for each z ---------
		Real x = (cosang-z*z0)*xa;
		Real ysq = 1-z*z-x*x;
		//planck_assert(ysq>=0, "error in query_disc()");
		Real dphi=atan2(sqrt(ysq),x);
		in_ring_allskymap (nside_, iz, phi, dphi, vecx, vecy, vecz,
                           radius / 2.0, weight_norm, aflux, maplist, indexp, getnorm);
	}
    
	if (rlat2>=pi) // south pole in the disc
		for (int m=irmax+1; m<(4*nside_); ++m)  // rings completely in the disc
			in_ring_allskymap (nside_, m, 0, pi,  vecx, vecy, 
                               vecz, radius / 2.0, weight_norm, aflux, maplist, indexp, getnorm);
}