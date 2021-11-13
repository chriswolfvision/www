/***************************************************************************
                          FloatMatrixClassMethods.cpp  -  description
                             -------------------
    begin                : Fri Mar 14 2003
    copyright            : (C) 2003 by Christian Wolf
    email                : wolf@rfv.insa-lyon.fr
 ***************************************************************************/

// From this module
#include "FloatMatrix.h"

// *************************************************************
// Class method:
// allocate an image plane
// *************************************************************

float ** FloatMatrix::_AllocPlane (int xsize, int ysize)	{
	float **im;
	float *big;
	typedef float * pfloat;

#ifdef CHECK_CODE
	if (xsize<=0 || ysize<=0) {
		cerr << "Invalid image sizes in FloatMatrix::_AllocPlane():\n"
			 << "xsize: " << xsize << ", ysize: " << ysize << endl;
		CORE_DUMP;
	}
#endif	

	im = new pfloat [xsize];
	big	= new float [xsize*ysize];

	for	(int i=0 ; i<xsize ; i++)
		im[i] =	big	+ i*ysize;

	return (im);
}

