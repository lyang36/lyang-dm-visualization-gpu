
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

