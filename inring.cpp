//calculate the weight for a given pixel
__device__ __host__
Real calc_weight(const long nside_, const int pix, const Real vecx, 
                 const Real vecy, const Real vecz, const Real norm, const Real angular_radius){
	Real _t, _p;
	pix2ang_ring(nside_, pix, &_t, &_p);
	Real this_vecx = sin(_t) * cos(_p);
	Real this_vecy = sin(_t) * sin(_p);
	Real this_vecz = cos(_t);
	Real d2 = acos(this_vecx * vecx + this_vecy * vecy + this_vecz * vecz) / angular_radius;
	d2 = d2*d2;
	return exp(-0.5 * d2 / 0.333) / norm;
}


//count how many pixels in this ring
__device__ __host__
void in_ring_count(long nside_, int iz, Real phi0, Real dphi,
                   const Real vecx, const Real vecy, const Real vecz, 
                   const Real angular_radius, Real &norm, int &num){
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
			num++;		
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
			if (pixnum > ipix2) pixnum -= nr;
			num ++;
			norm += calc_weight(nside_, pixnum, vecx, vecy, vecz, 1.0, angular_radius);
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
