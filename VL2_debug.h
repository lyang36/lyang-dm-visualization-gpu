
#ifndef __DEBUGGGG__
#define __DEBUGGGG__
//#define _DEBUG__LY__
#include <string>
using namespace std;

string getstring_Real( double * d, long numx);
string getstring_long( long * d, long numx);

void print_out_natconst( Natconst * natconst);
void print_out_halo( Halo * halo);
void print_out_dm( Dm * dm);

void print_out_codeunits( Codeunits * codeunits);

void print_out_params( Params * params);

void print_out_master( Master * master);
//void print_out_tipsyheader( Tipsyheader * tipsyheader);

void print_out_particle( Particle * particle);


void print_out_map( Map * map);/*
void print_out_allparts( Allparts * allparts){
	#ifdef _DEBUG__LY__

	Real * mass;
	Real * density;
	Real * hsmooth;
	Real * xpos;
	Real * ypos;
	Real * zpos;

	#endif
}
*/
#endif

