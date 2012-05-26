//calculate the angle for given position in terms with the rotation matrix

__device__ 
void calc_angles( Real xpos, Real ypos, Real zpos, Real &distances, 
                 Real * opos, Real *rotmatrix, Real & costheta, Real &phi){
    
	Real vec[] = { xpos, ypos, zpos};

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
