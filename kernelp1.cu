#ifndef __CUDACC__
#define __CUDACC__
#endif
#include "cuda.h"
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <stdio.h>
#include "kernel.h"
#include <iostream>
#include "structures.h"


//prototype:
__global__ void doKernel(int *c, const int *a, const int *b);
__device__ void ang2pix_ring( const long nside, Real theta, Real phi, long *ipix);
__device__ 
void calc_angles( Real xpos, Real ypos, Real zpos, Real &distances, 
					 Real * opos, Real *rotmatrix, Real & costheta, Real &phi);
__device__
void query_disc (long nside_, const Real theta, const Real phi, 
	const Real radius, const int maxlength, int *disk, int & disknum); 
//return the disk in disk, disknum is the number of pixels
//disk is a group of shared memery in the device with the length smaller then the maximum length
__device__
void pix2ang_ring( long nside, long ipix, Real *theta, Real *phi);

__global__
void generate_map(const long Nside, 
		 const Real theta0, const Real fluxfactor, const long nmax,
		 mapMap *maplist, MapParticle * particles,
		 Real * rotmatrix, Real * opos);


///////////
#ifndef M_PI
#define M_PI           3.14159265358979323846
#endif

#ifndef PI
#define PI           M_PI
#endif

__device__ 
const Real inv_twopi = 1.0 / (2.0 * PI);
__device__ 
const Real pi = PI;
__device__ 
const Real twothird = 2.0 / 3.0;

// Helper function for using CUDA to add vectors in parallel.
cudaError_t doWithCuda(const long MAX_Num_Particle, const long Nside, 
		 const Real theta0, const Real fluxfactor, const long nmax,
		 Real *allskymap, MapParticle * particles,
		 Real * rotmatrix, Real * opos, mapMap *dev_maplist, 
		 mapMap *host_maplist){
	cudaError_t cudaStatus;
	dim3 dimBlock(512, 1, 1);
	dim3 dimGrid(MAX_Num_Particle/dimBlock.x + 1, 1, 1);
	//Real fluxfactor = master->codeunits.annihilation_flux_to_cgs;

	generate_map<<<dimGrid, dimBlock>>>(Nside, theta0, fluxfactor , nmax, dev_maplist, 
		particles, rotmatrix, opos);
		
	cudaStatus = cudaDeviceSynchronize();
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaDeviceSynchronize returned error code %d after launching addKernel!\n", cudaStatus);
		return cudaStatus;
	}
	cudaStatus = cudaMemcpy(host_maplist, dev_maplist, sizeof(mapMap) * MAX_PIX_PER_PAR * nmax, cudaMemcpyDeviceToHost);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMemcpy failed!");
		exit(0);
	}
	//add to allskymap
	for(int i = 0; i < nmax; i++){
		if(host_maplist[i*MAX_PIX_PER_PAR].pix > MAX_PIX_PER_PAR - 1){
			std::cout<<"Warning: particle " << i << " contributes to more than " << MAX_PIX_PER_PAR 
				<< "pixels. Only " <<  MAX_PIX_PER_PAR << " pixels are considerd!" << std::endl;
		}
		for(int j=i*MAX_PIX_PER_PAR+1; j<=host_maplist[i*MAX_PIX_PER_PAR].pix; j++){
			allskymap[host_maplist[j].pix] += host_maplist[j].signal;
			int t = host_maplist[j].pix;
			double kk = allskymap[t];
			t=t;
		}

	}

	return cudaStatus;
}

__device__ inline void ang2pix_ring( const long nside, Real theta, Real phi, long *ipix) {
  /*
    c=======================================================================
    c     gives the pixel number ipix (RING) 
    c     corresponding to angles theta and phi
    c=======================================================================
  */
  

  int nl2, nl4, ncap, npix, jp, jm, ipix1;
  Real  z, za, tt, tp, tmp;
  int ir, ip, kshift;
  
  Real piover2 = 0.5 * M_PI;
  Real twopi=2.0*M_PI;
  Real z0=2.0/3.0;
  
  
  z = cos(theta);
  za = fabs(z);
  if( phi >= twopi)  phi = phi - twopi;
  if (phi < 0.)     phi = phi + twopi;
  tt = phi / piover2;//  ! in [0,4)
  
  nl2 = 2*nside;
  nl4 = 4*nside;
  ncap  = nl2*(nside-1);// ! number of pixels in the north polar cap
  npix  = 12*nside*nside;
  
  if( za <= z0 ) {
    
    jp = (int)floor(nside*(0.5 + tt - z*0.75)); /*index of ascending edge line*/
    jm = (int)floor(nside*(0.5 + tt + z*0.75)); /*index of descending edge line*/
    
    ir = nside + 1 + jp - jm;// ! in {1,2n+1} (ring number counted from z=2/3)
    kshift = 0;
    if (ir % 2==0) kshift = 1;// ! kshift=1 if ir even, 0 otherwise

    ip = (( jp+jm - nside + kshift + 1 ) / 2)  + 1;// ! in {1,4n}
    if( ip>nl4 ) ip = ip - nl4;
    
    ipix1 = ncap + nl4*(ir-1) + ip ;
  }
  else {
    
    tp = tt - floor(tt);//      !MOD(tt,1.d0)
    tmp = sqrt( 3.*(1. - za) );
    
    jp = (int)floor( nside * tp * tmp );// ! increasing edge line index
    jm = (int)floor( nside * (1. - tp) * tmp );// ! decreasing edge line index
    
    ir = jp + jm + 1;//        ! ring number counted from the closest pole
    ip = (int)floor( tt * ir ) + 1;// ! in {1,4*ir}
    if( ip>4*ir ) ip = ip - 4*ir;
    
    ipix1 = 2*ir*(ir-1) + ip;
    if( z<=0. ) {
      ipix1 = npix - 2*ir*(ir+1) + ip;
    }
  }
  *ipix = ipix1 - 1;// ! in {0, npix-1}
  
}


__device__ 
void calc_angles( Real xpos, Real ypos, Real zpos, Real &distances, 
						 Real * opos, Real *rotmatrix, Real & costheta, Real &phi){
	Real vec[] = { xpos, ypos, zpos};
	//Real temp[] = {0, 0, 0};
	//cblas_dgemv(CblasColMajor, CblasNoTrans, 3, 3, 1.0, rotmatrix, 3, vec, 1, 1.0, temp, 1 );
	Real xcoord;
	Real ycoord;
	Real zcoord;
	xcoord = rotmatrix[0] * vec[0] + rotmatrix[3] * vec[1]  + rotmatrix[6] * vec[2];
	ycoord = rotmatrix[1] * vec[0] + rotmatrix[4] * vec[1]  + rotmatrix[7] * vec[2];
	zcoord = rotmatrix[2] * vec[0] + rotmatrix[5] * vec[1]  + rotmatrix[8] * vec[2];
  
	costheta = zcoord / distances;

	phi = atan2( ycoord, xcoord );
	
	if( phi < 0 ){
		phi += 2.0 * PI;
	}

	//a few adjustments...
  
	//one more rotation by 180 degrees...
	phi -= PI;
  
	//force phi to lie between 0 and 2*pi  
	if( phi < 0 ){
		phi = 2.0 * PI + phi;
	}
  
	if( phi > 2 * PI){
		phi = phi - 2.0 * PI;
	}
}

__device__
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

__device__
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

/*
__device__
void in_ring(long nside_, int iz, Real phi0, Real dphi,
  int* listir, int maxnum, int & pos){
	long npface_ = nside_ * nside_;
	long ncap_   = (npface_-nside_)<<1;
	long npix_   = 12 * npface_;
	Real fact2_  = 4. / (Real) npix_;
	Real fact1_  = (nside_<<1) * (Real) fact2_;

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
			if( ++pos < maxnum){
				listir[(pos)] = i;
			}
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
			if( ++pos < maxnum){
				listir[(pos)] = pixnum;
			}
		}
	}
}
*/


__device__
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

__device__
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
  const Real angular_radius, const Real norm, const Real aflux, 
  mapMap * maplist, int & indexp, const int pmax){
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
			/*if( ++pos < maxnum){
				listir[(pos)] = i;
			}*/
			//atomic plus
			//allskymap[i] += 
			Real _flux = calc_weight(nside_, i, vecx, vecy, vecz, norm, angular_radius);
			if(indexp < pmax){
				maplist[indexp].pix = i;
				maplist[indexp].signal = _flux * aflux;
			}
			indexp ++;
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
			if(indexp < pmax){
				maplist[indexp].pix = pixnum;
				maplist[indexp].signal = _flux * aflux;
				indexp ++;
			}
			/*if( ++pos < maxnum){
				listir[(pos)] = pixnum;
			}*/
		}
	}
}


__device__
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
	const Real vecx,  const Real vecy, const Real vecz, const Real weight_norm, const Real aflux, mapMap *maplist, 
	int &indexp, const int pmax){
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
		in_ring_allskymap (nside_, m, 0, pi, vecx, vecy, vecz, radius / 2.0, weight_norm, aflux, maplist, indexp, pmax);

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
		in_ring_allskymap (nside_, iz, phi, dphi, vecx, vecy, vecz, radius / 2.0, weight_norm, aflux, maplist, indexp, pmax);
	}

	if (rlat2>=pi) // south pole in the disc
		for (int m=irmax+1; m<(4*nside_); ++m)  // rings completely in the disc
			in_ring_allskymap (nside_, m, 0, pi,  vecx, vecy, vecz, radius / 2.0, weight_norm, aflux, maplist, indexp, pmax);
}

__global__
void generate_map(const long Nside, 
		 const Real theta0, const Real fluxfactor, const long nmax,
		 mapMap *maplist, MapParticle * particles,
		 Real * rotmatrix, Real * opos){
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
	int startp = MAX_PIX_PER_PAR * i;

	// Real3 opos...
	// opos.x, .y, .z
#ifdef POS_FLOAT
	float3 _opos = {opos[0], opos[1], opos[2]};
#else
	double3 _opos = {opos[0], opos[1], opos[2]};
#endif
	/*since if density <0, the threads becomes null*/
	/*please make sure all the density is greater or equal to 0*/
	/*this is very important*/
	if( density >= 0.0){
		distances = sqrt( (xpos-_opos.x) * (xpos-_opos.x) + 
			(ypos-_opos.y) * (ypos-_opos.y) +
			(zpos-_opos.z) *(zpos-_opos.z) );
		fluxes = fluxfactor * density * mass / (4.0 * PI * distances * distances);

		calc_angles( xpos-_opos.x, ypos-_opos.y, zpos-_opos.z, distances, 
		opos, rotmatrix, costheta, phi);
		theta = acos(costheta);
		angular_radius = hsmooth / distances;
	/*------------------------------------------add particles------------------------------------------------------*/
		{
			//pointing p (theta,phi);
			Real vecx = sin(theta) * cos(phi);
			Real vecy = sin(theta) * sin(phi);
			Real vecz = cos(theta);
			//ang2vec(theta,phi,&vec);

			long pix;
			ang2pix_ring(Nside,theta,phi,&pix);
			if( 2.0*angular_radius < theta0 ) {
				/*should be atomic plus */
				//allskymap[pix] += fluxes;
				maplist[startp + 0].pix = 1;
				maplist[startp + 1].pix = pix;
				maplist[startp + 1].signal = fluxes;
			}else{
				Real weight_norm = 0.0;
				int Npix = 0;
				query_disc_getnorm(Nside, theta, phi, 2.0*angular_radius, vecx, vecy, vecz, weight_norm, Npix);
				// if either zero or one pixel are covered by particle (this should be avoided by the above condition...)
					//<debug
					//allskymap[i] = weight_norm;
					//return;
					//debug>
				//store the number of pixels
				maplist[startp + 0].pix = Npix;
				if(Npix < 2) {
					/*should be atomic plus */
					//allskymap[pix] += fluxes;
					//atomicAdd((allskymap + pix), fluxes);
					maplist[startp + 1].pix = pix;
					maplist[startp + 1].signal = fluxes;
				}else{
					int indexp = startp + 1;
					int pmax = indexp + MAX_PIX_PER_PAR - 1;
					// get here only if the particle covers more than one pixel
					query_disc_calcflux(Nside, theta, phi, 2.0*angular_radius, vecx, vecy, vecz, weight_norm, fluxes, maplist, indexp, pmax);

				}
			}
		}
	}
}

		