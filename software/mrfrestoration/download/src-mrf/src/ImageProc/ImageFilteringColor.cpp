/***************************************************************************
                          ImageFilteringColor.cc  -  description
                             -------------------
    begin                : Mon Feb 4 2002
    copyright            : (C) 2002 by Christian Wolf
    email                : wolf@rfv.insa-lyon.fr
 ***************************************************************************/

 // C
#include <math.h>

// From the IMAGE library
#include "FloatMatrix.h"
#include "FloatMatrixColor.h"

// From the IMAGE PROCESSING library
#include "FilterMask.h"

// From our own libary
#include "BaseMath.h"
#include "ImageProc.h"

// *************************************************************
// Use Deriche's filter on a given plane of a Float Image
// to calculate the derivatives in x and y direction
// *************************************************************

void filterDericheGetDerivatives (FloatMatrixColor &fim, int plane, float alpha,
	FloatMatrix *&outDerivX, FloatMatrix *&outDerivY) {

	int	maxSize;
	float  **A;
	double cste;
	float a, b1, b2, *ym, *yp;
	float S;

	// Allocate the return images
	outDerivX = new FloatMatrix (fim.xsize, fim.ysize);
	outDerivY = new FloatMatrix (fim.xsize, fim.ysize);

	// Creation	of the tables
	cste = 255.*255. ;
	A =	MAT_FLOAT(fim.ysize, fim.xsize) ;
	maxSize	= (fim.ysize > fim.xsize ? fim.ysize : fim.xsize) ;
	yp = VECT_FLOAT(maxSize) ;
	ym = VECT_FLOAT(maxSize) ;

	// Initializations
	b1 = exp(-alpha) ;
	S =	-(1.-b1)*(1.-b1)/b1	;
	a =	S*b1;
	b2 = -b1 * b1 ;
	b1 = 2*b1 ;

	for	(int i = 0 ; i < fim.xsize ; i++)
	for	(int j = 0 ; j < fim.ysize ; j++)
		A[i][j]	= fim.get(plane, i, j);

	// ------------------------------------------------
	// Gradient	in x direction
	// ------------------------------------------------

	for	(int i = 0 ; i < fim.xsize ; i++)	{

		yp[0] =	0 ;
		yp[1] =	A[i][0]+b1*yp[0] ;
		for	(int j = 2 ; j < fim.ysize ; j++)
			yp[j] =	A[i][j-1]+b1*yp[j-1]+b2*yp[j-2]	;

		ym[fim.ysize-1] =	0 ;
		ym[fim.ysize-2] =	-A[i][fim.ysize-1]+b1*ym[fim.ysize-1] ;
		for	(int j = fim.ysize-3 ; j >= 0	; j--)
			ym[j] =	-A[i][j+1]+b1*ym[j+1]+b2*ym[j+2] ;

		for	(int j = 0 ; j < fim.ysize ; j++)
			outDerivX->set(i, j, a*(yp[j]+ym[j]));
	}

	// ------------------------------------------------
	// Gradient	in y direction
	// ------------------------------------------------

	for	(int i = 0 ; i < fim.ysize ; i++)	{
		yp[0] =	0 ;
		yp[1] =	A[0][i]+b1*yp[0] ;
		for	(int j = 2 ; j < fim.xsize ; j++)
			yp[j] =	A[j-1][i]+b1*yp[j-1]+b2*yp[j-2]	;

		ym[fim.xsize-1] =	0 ;
		ym[fim.xsize-2] =	-A[fim.xsize-1][i]+b1*ym[fim.xsize-1] ;
		for	(int j = fim.xsize-3 ; j >= 0	; j--)
			ym[j] =	-A[j+1][i]+b1*ym[j+1]+b2*ym[j+2] ;

		for	(int j = 0 ; j < fim.xsize ; j++)
			outDerivY->set(j, i, a*(yp[j]+ym[j]));
	}

	// Clean up
	FREE_MAT_FLOAT (fim.xsize, A);
	FREE_VECT_FLOAT(yp);
	FREE_VECT_FLOAT(ym)	;
}



