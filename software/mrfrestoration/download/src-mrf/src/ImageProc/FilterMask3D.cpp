#include <stdlib.h>
#include <stdio.h>

#include "FilterMask3D.h"

/************************************************************************
 * Constructor:	Empty Mask
 ************************************************************************/

FilterMask3D::FilterMask3D(int x, int y, int z){
	M =	new	double[x*y*z];
	xsize =	x;
	ysize =	y;
	zsize =	z;
}

/************************************************************************
 * Normalize the filter
 ************************************************************************/

void FilterMask3D::normalize(){
	double sum = 0;

	for	(int z=0; z<zsize; ++z)
	for	(int y=0; y<ysize; ++y)
	for	(int x=0; x<xsize; ++x)	{
		sum	+= get (x,y,z);
	}

	if (sum==0)
		return;

	for	(int z=0; z<zsize; ++z)
	for	(int y=0; y<ysize; ++y)
	for	(int x=0; x<xsize; ++x)
		set(x,y,z, get(x,y,z)/sum);

}

