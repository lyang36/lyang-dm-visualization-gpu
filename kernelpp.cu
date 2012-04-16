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
#include <ctime>

#define MAX_PIXES 500

///////////
#ifndef M_PI
#define M_PI           3.14159265358979323846
#endif

#ifndef PI
#define PI           M_PI
#endif

const dim3 dimBlock(256, 1, 1);


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

__device__ __host__
void query_disc_getnorm(long nside_, const Real theta, const Real phi, 
	const Real radius, const Real vecx,  const Real vecy, const Real vecz, 
	Real &weight_norm, int &Npix);


__global__
void generate_map_step1(const long Nside, 
		 const Real theta0, const long nmax,
		 MapParticle * particles,
		 Real * rotmatrix, Real * opos);


__global__
void generate_map_step2(const long Nside, 
		 const Real theta0, const long nmax,
		 MapParticle * particles, mapMap *skymap,
		 Real * rotmatrix, Real * opos);


__global__
void generate_map_step0(const long Nside, 
		 const Real theta0, const long nmax,
		 MapParticle * particles,
		 Real * rotmatrix, Real * opos,
		 int * pd_key, int * pd_val, int stpoint);

__global__
void generate_map_step3(const long nmax, 
	mapMap *maplist);



__device__ 
const Real inv_twopi = 1.0 / (2.0 * PI);
__device__ 
const Real pi = PI;
__device__ 
const Real twothird = 2.0 / 3.0;

cudaError_t doWithCuda_pre(const long MAX_Num_Paritcle, const long Nside, 
		 const Real theta0, const Real fluxfactor, const long nmax,
		 Real *allskymap, MapParticle * dev_par, MapParticle * host_par,
		 Real * dev_rotm, Real * dev_opos, int * pd_key, int * pd_val, int stpoint)
{
	cudaError_t cudaStatus;
	dim3 dimGrid(MAX_Num_Paritcle/dimBlock.x + 1, 1, 1);
	//std::cout << 0 << " sec -> step1 start" << std::endl;
	clock_t starttime = clock();
	generate_map_step0<<<dimGrid, dimBlock>>>(Nside, theta0, nmax,
		dev_par, dev_rotm, dev_opos, pd_key, pd_val, stpoint);
	cudaStatus = cudaDeviceSynchronize();
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaDeviceSynchronize returned error code %d after launching generate_map_step0!\n", cudaStatus);
		return cudaStatus;
	}
	//std::cout <<"step0 cost: " << (clock() - starttime) / 1000.0 << " secs. "<< std::endl;
	cudaStatus = cudaMemcpy(host_par, dev_par, sizeof(MapParticle) * nmax, cudaMemcpyDeviceToHost);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMemcpy failed!");
		exit(0);
	}
	return cudaStatus;
}


cudaError_t run_kernel(const long MAX_Num_Paritcle, const long Nside, 
		 const Real theta0, const Real fluxfactor, const long nmax,
		 Real *allskymap, MapParticle * dev_par, MapParticle * host_par,
		 Real * dev_rotm, Real * dev_opos, int init){
	//std::cout << "Add map length cost" << (clock() -starttime) / 1000.0 << " sec"<< std::endl;
	//copy the host_par back to device_par
	//Real fluxfactor = master->codeunits.annihilation_flux_to_cgs;
	dim3 dimGrid(MAX_Num_Paritcle/dimBlock.x + 1, 1, 1);
	mapMap * maplist;
	mapMap * dev_maplist;

	cudaError_t cudaStatus;
	cudaStatus = cudaMemcpy(dev_par, host_par, sizeof(MapParticle) * nmax, cudaMemcpyHostToDevice);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMemcpy failed!");
		return cudaStatus;
	}
	//relocate the memory for skymap
	maplist = (mapMap *)calloc(init,sizeof(mapMap));
	cudaStatus = cudaMalloc((void**)&dev_maplist, sizeof(mapMap) * init);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMalloc failed!");
		return cudaStatus;
	}
	//initiate maplist
	cudaStatus = cudaMemcpy(dev_maplist, maplist, sizeof(mapMap) * init, cudaMemcpyHostToDevice);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMemcpy failed!");
		return cudaStatus;
	}

	//std::cout << (clock() -starttime) / 1000.0 << " sec -> step2 start" << std::endl;
	generate_map_step2<<<dimGrid, dimBlock>>>(Nside, theta0, nmax,
		dev_par, dev_maplist, dev_rotm, dev_opos);
	cudaStatus = cudaDeviceSynchronize();
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaDeviceSynchronize returned error code %d after launching generate_map_step2!\n", cudaStatus);
		return cudaStatus;
	}

	cudaStatus = cudaMemcpy(maplist, dev_maplist, sizeof(mapMap) * init, cudaMemcpyDeviceToHost);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMemcpy failed!");
		return cudaStatus;
	}
	//add to allskymap
	double norm = 1;
	for(int i = 0; i < init; i++){
		if(maplist[i].pix < 12*512L*512L && maplist[i].pix >= 0){
			if(maplist[i].factor > 0 )
				{norm = maplist[i].factor;}
			double ff = maplist[i].signal/norm;
			allskymap[maplist[i].pix] += ff;
		}
		else{
			int k = maplist[i].pix;
		}
	}
	//std::cout << " maplist: " << init << std::endl;
	free(maplist);
	cudaFree(dev_maplist);
	//std::cout << (clock() -starttime) / 1000.0 << " sec -> kernel end" << std::endl;
}


// Helper function for using CUDA to add vectors in parallel.
cudaError_t doWithCuda_Par(const long MAX_Num_Paritcle, const long Nside, 
		 const Real theta0, const Real fluxfactor, const long nmax,
		 Real *allskymap, MapParticle * dev_par, MapParticle * host_par,
		 Real * dev_rotm, Real * dev_opos)
{
	cudaError_t cudaStatus;
	dim3 dimGrid(MAX_Num_Paritcle/dimBlock.x + 1, 1, 1);

	//std::cout << 0 << " sec -> step1 start" << std::endl;
	clock_t starttime = clock();
	cudaStatus = cudaMemcpy(dev_par, host_par, sizeof(MapParticle) * nmax, cudaMemcpyHostToDevice);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMemcpy failed!");
		return cudaStatus;
	}
    generate_map_step1<<<dimGrid, dimBlock>>>( Nside, 
		 theta0, nmax,
		 dev_par,
		dev_rotm, dev_opos);
	cudaStatus = cudaDeviceSynchronize();
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaDeviceSynchronize returned error code %d after launching generate_map_step1!\n", cudaStatus);
		return cudaStatus;
	}
	cudaStatus = cudaMemcpy(host_par, dev_par, sizeof(MapParticle) * nmax, cudaMemcpyDeviceToHost);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMemcpy failed!");
		return cudaStatus;
	}
	//calculting the index for maplist
	int init = 0;
	int _tenmax = 0;
	int maxPixx = MAX_Num_Paritcle * MAX_PIXES;
	for(int i=0; i < nmax; i++){
		//accumulate the particles and deal with the memery overflow
		int init1 = init;
		init += (int)(host_par[i].xpos);
		host_par[i].xpos = init1; 
		if(init1 > maxPixx){
			//run the kernel onec
			int st_p = _tenmax;
			_tenmax = i - _tenmax + 1;
			MapParticle * __dev_par = dev_par + st_p;
			cudaStatus = run_kernel(MAX_Num_Paritcle, Nside, theta0, fluxfactor, _tenmax, allskymap,
				__dev_par, host_par + st_p, dev_rotm, dev_opos, init);
			init = 0;
			_tenmax = i + 1;
		}

	}
	int _nmax = nmax - _tenmax;
	if( _nmax > 0 ){
		cudaStatus = run_kernel(MAX_Num_Paritcle,Nside,theta0,fluxfactor,_nmax,allskymap,
			dev_par + _tenmax, host_par + _tenmax, dev_rotm, dev_opos, init);
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
