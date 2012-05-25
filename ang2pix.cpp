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
