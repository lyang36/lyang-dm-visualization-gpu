

__device__ __host__
void pix2ang_ring( long nside, long ipix, Real *theta, Real *phi) {
    /*
     c=======================================================================
     c     gives theta and phi corresponding to pixel ipix (RING) 
     c     for a parameter nside
     c=======================================================================
     */
    
    int nl2, nl4, npix, ncap, iring, iphi, ip, ipix1;
    Real  fact1, fact2, fodd, hip, fihip;
    //      PARAMETER (pi     = 3.1415926535897932384626434d0)
    //      parameter (ns_max = 8192) ! 2^13 : largest nside available
    
    //  int ns_max=8192;
    
    /*  if( nside<1 || nside>ns_max ) {
     fprintf(stderr, "%s (%d): nside out of range: %ld\n", __FILE__, __LINE__, nside);
     exit(0);
     }*/
    npix = 12*nside*nside;      // ! total number of points
    /*  if( ipix<0 || ipix>npix-1 ) {
     fprintf(stderr, "%s (%d): ipix out of range: %ld\n", __FILE__, __LINE__, ipix);
     exit(0);
     }*/
    
    ipix1 = ipix + 1; // in {1, npix}
    nl2 = 2*nside;
    nl4 = 4*nside;
    ncap = 2*nside*(nside-1);// ! points in each polar cap, =0 for nside =1
    fact1 = 1.5*nside;
    fact2 = 3.0*nside*nside;
    
    if( ipix1 <= ncap ) {  //! North Polar cap -------------
        
        hip   = ipix1/2.;
        fihip = floor(hip);
        iring = (int)floor( sqrt( hip - sqrt(fihip) ) ) + 1;// ! counted from North pole
        iphi  = ipix1 - 2*iring*(iring - 1);
        
        *theta = acos( 1. - iring*iring / fact2 );
        *phi   = (1.*iphi - 0.5) * PI/(2.*iring);
    }
    else if( ipix1 <= nl2*(5*nside+1) ) {//then ! Equatorial region ------
        
        ip    = ipix1 - ncap - 1;
        iring =  ip / nl4 + nside;// ! counted from North pole
        iphi  =  ip % nl4 + 1;
        
        fodd  = 0.5 * (1 + fmod((Real)(iring+nside),2));//  ! 1 if iring+nside is odd, 1/2 otherwise
        *theta = acos( (nl2 - iring) / fact1 );
        *phi   = (1.*iphi - fodd) * PI /(2.*nside);
    }
    else {//! South Polar cap -----------------------------------
        
        ip    = npix - ipix1 + 1;
        hip   = ip/2.;
        /* bug corrige floor instead of 1.* */
        fihip = floor(hip);
        iring = (int)floor( sqrt( hip - sqrt(fihip) ) ) + 1;//     ! counted from South pole
        iphi  = (int)(4.*iring + 1 - (ip - 2.*iring*(iring-1)));
        
        *theta = acos( -1. + iring*iring / fact2 );
        *phi   = (1.*iphi - 0.5) * PI/(2.*iring);
    }
}


////////////////////////query disk//////////////////////////////

__device__ __host__
int ring_above (long nside_, Real z){
	Real az=abs(z);
	if (az>twothird) // polar caps
	{
		int iring = int(nside_*sqrt(3*(1-az)));
		return (z>0) ? iring : 4*nside_-iring-1;
	}
	else // ----- equatorial region ---------
        return int(nside_*(2-1.5*z));
}



__device__ __host__
Real calc_weight(const long nside_, const int pix, const Real vecx, 
                 const Real vecy, const Real vecz, const Real norm, const Real angular_radius){
	Real _t, _p;
	pix2ang_ring(nside_, pix, &_t, &_p);
	Real this_vecx = sin(_t) * cos(_p);
	Real this_vecy = sin(_t) * sin(_p);
	Real this_vecz = cos(_t);
	Real d2 = acos( this_vecx * vecx + this_vecy * vecy + this_vecz * vecz) / angular_radius;
	d2 = d2*d2;
	//norm += exp(-0.5 * d2 / 0.333) / norm;  
	return exp(-0.5 * d2 / 0.333) / norm;
}

__device__ __host__
void in_ring_count(long nside_, int iz, Real phi0, Real dphi,
                   const Real vecx, const Real vecy, const Real vecz, 
                   const Real angular_radius, Real &norm, int &num){
	long npface_ = nside_ * nside_;
	long ncap_   = (npface_-nside_)<<1;
	long npix_   = 12 * npface_;
    //	Real fact2_  = 4. / (Real) npix_;
    //	Real fact1_  = (nside_<<1) * (Real) fact2_;
    
	int nr, ir, ipix1;
	Real shift=0.5;
    
	if (iz<nside_) // north pole
	{
		ir = iz;
		nr = ir*4;
		ipix1 = 2*ir*(ir-1);        //    lowest pixel number in the ring
	}
	else if (iz>(3*nside_)) // south pole
	{
		ir = 4*nside_ - iz;
		nr = ir*4;
		ipix1 = npix_ - 2*ir*(ir+1); // lowest pixel number in the ring
	}
	else // equatorial region
	{
		ir = iz - nside_ + 1;           //    within {1, 2*nside + 1}
		nr = nside_*4;
		if ((ir&1)==0) shift = 0;
		ipix1 = ncap_ + (ir-1)*nr; // lowest pixel number in the ring
	}
    
	int ipix2 = ipix1 + nr - 1;       //    highest pixel number in the ring
    
	// ----------- constructs the pixel list --------------
	if (dphi > (PI-1e-7))
		for (int i=ipix1; i<=ipix2; ++i){
			num++;		
			/*if( ++pos < maxnum){
             listir[(pos)] = i;
             }*/
			norm += calc_weight(nside_, i, vecx, vecy, vecz, 1.0, angular_radius);
        }
	else
	{
		int ip_lo = floor(nr*inv_twopi*(phi0-dphi) - shift)+1;
		int ip_hi = floor(nr*inv_twopi*(phi0+dphi) - shift);
		int pixnum = ip_lo+ipix1;
		if (pixnum<ipix1) pixnum += nr;
		for (int i=ip_lo; i<=ip_hi; ++i, ++pixnum)
		{
			if (pixnum>ipix2) pixnum -= nr;
			num ++;
			norm += calc_weight(nside_, pixnum, vecx, vecy, vecz, 1.0, angular_radius);
			/*if( ++pos < maxnum){
             listir[(pos)] = pixnum;
             }*/
		}
	}
}

__device__
void in_ring_allskymap(long nside_, int iz, Real phi0, Real dphi,
                       const Real vecx, const Real vecy, const Real vecz, 
                       const Real angular_radius, const Real norm, 
                       const Real aflux, mapMap * maplist, long &indexp, double &getnorm){
	long npface_ = nside_ * nside_;
	long ncap_   = (npface_-nside_)<<1;
	long npix_   = 12 * npface_;
    
	int nr, ir, ipix1;
	Real shift=0.5;
    
	if (iz<nside_) // north pole
	{
		ir = iz;
		nr = ir*4;
		ipix1 = 2*ir*(ir-1);        //    lowest pixel number in the ring
	}
	else if (iz>(3*nside_)) // south pole
	{
		ir = 4*nside_ - iz;
		nr = ir*4;
		ipix1 = npix_ - 2*ir*(ir+1); // lowest pixel number in the ring
	}
	else // equatorial region
	{
		ir = iz - nside_ + 1;           //    within {1, 2*nside + 1}
		nr = nside_*4;
		if ((ir&1)==0) shift = 0;
		ipix1 = ncap_ + (ir-1)*nr; // lowest pixel number in the ring
	}
    
	int ipix2 = ipix1 + nr - 1;       //    highest pixel number in the ring
    
	// ----------- constructs the pixel list --------------
	if (dphi > (PI-1e-7))
		for (int i=ipix1; i<=ipix2; ++i){
			Real _flux = calc_weight(nside_, i, vecx, vecy, vecz, norm, angular_radius);
			maplist[indexp].pix = i;
			maplist[indexp].signal = _flux * aflux;
			maplist[indexp].factor = 0;
			getnorm += _flux;
			indexp = indexp +1 ;
        }
	else
	{
		int ip_lo = floor(nr*inv_twopi*(phi0-dphi) - shift)+1;
		int ip_hi = floor(nr*inv_twopi*(phi0+dphi) - shift);
		int pixnum = ip_lo+ipix1;
		if (pixnum<ipix1) pixnum += nr;
		for (int i=ip_lo; i<=ip_hi; ++i, ++pixnum)
		{
			if (pixnum>ipix2) pixnum -= nr;
			//num ++;
			//atomic plus
			//allskymap[pixnum] += 
			Real _flux = calc_weight(nside_, pixnum, vecx, vecy, vecz, norm, angular_radius);
			maplist[indexp].pix = pixnum;
			maplist[indexp].signal = _flux * aflux;
			maplist[indexp].factor = 0;
			getnorm += _flux;
			indexp = indexp + 1;
		}
	}
}


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