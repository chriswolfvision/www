/**********************************************************
 * ImageFilter.cc
 * Filtering algorithms for the Image class
 *
 * Christian Wolf
 **********************************************************/

// C
#include <iostream>
#include <errno.h>
#include <string.h>
#include <math.h>

// From the IMAGE module
#include <color.h>

// From the MATHEMATICS module
#include <Histogram.h>

// From our own module
#include "BaseMath.h"
#include "FilterMask.h"
#include "ImageProc.h"

// *************************************************************
// Filter the image with a deriche edge filter
// Version for Grayscale Images
// *************************************************************

void filterDeriche (Image &im, float alpha) {
	int	maxSize;
	unsigned char  **I ;
	float  **A,	**Ix, **Iy ;
	double cste	;
	float	g, a, b1, b2, *ym, *yp ;
	float S;

	I = im.R;

	// Creation	of the tables
	cste = 255.*255. ;
	A =	MAT_FLOAT(im.ysize,im.xsize) ;
	Ix = MAT_FLOAT(im.ysize,im.xsize)	;
	Iy = MAT_FLOAT(im.ysize,im.xsize)	;
	maxSize	= (im.ysize > im.xsize ? im.ysize : im.xsize) ;
	yp = VECT_FLOAT(maxSize) ;
	ym = VECT_FLOAT(maxSize) ;

	// Initializations
	b1 = exp(-alpha) ;
	S =	-(1.-b1)*(1.-b1)/b1	;
	a =	S*b1;
	b2 = -b1 * b1 ;
	b1 = 2*b1 ;

	// fprintf (stderr,	"Deriche - b1: %.4f	b2:	%.4f a:	%.4f\n", b1, b2, a);

	for	(int i = 0 ; i < im.xsize ; i++)
	for	(int j = 0 ; j < im.ysize ; j++)
		A[i][j]	= I[i][j] ;

	// ------------------------------------------------
	// Gradient	in x direction
	// ------------------------------------------------

	for	(int i = 0 ; i < im.xsize ; i++)	{

		yp[0] =	0 ;
		yp[1] =	A[i][0]+b1*yp[0] ;
		for	(int j = 2 ; j < im.ysize ; j++)
			yp[j] =	A[i][j-1]+b1*yp[j-1]+b2*yp[j-2]	;

		ym[im.ysize-1] =	0 ;
		ym[im.ysize-2] =	-A[i][im.ysize-1]+b1*ym[im.ysize-1] ;
		for	(int j = im.ysize-3 ; j >= 0	; j--)
			ym[j] =	-A[i][j+1]+b1*ym[j+1]+b2*ym[j+2] ;

		for	(int j = 0 ; j < im.ysize ; j++)
			Ix[i][j] = a*(yp[j]+ym[j]) ;
	}

	// ------------------------------------------------
	// Gradient	in y direction
	// ------------------------------------------------

	for	(int i = 0 ; i < im.ysize ; i++)	{
		yp[0] =	0 ;
		yp[1] =	A[0][i]+b1*yp[0] ;
		for	(int j = 2 ; j < im.xsize ; j++)
			yp[j] =	A[j-1][i]+b1*yp[j-1]+b2*yp[j-2]	;

		ym[im.xsize-1] =	0 ;
		ym[im.xsize-2] =	-A[im.xsize-1][i]+b1*ym[im.xsize-1] ;
		for	(int j = im.xsize-3 ; j >= 0	; j--)
			ym[j] =	-A[j+1][i]+b1*ym[j+1]+b2*ym[j+2] ;

		for	(int j = 0 ; j < im.xsize ; j++)
			Iy[j][i] = a*(yp[j]+ym[j]) ;
	}

	// ------------------------------------------------
	// Collect the results
	// ------------------------------------------------

	for	(int i = 0 ; i < im.xsize ; i++)
	for	(int j = 0 ; j < im.ysize ; j++)	{

		g =	sqrt(Ix[i][j]*Ix[i][j]+Iy[i][j]*Iy[i][j]) ;
		// g =	abs(Iy[i][j]);

		/*
		Ix[i][j] = abs ((int) Ix[i][j]);
		Iy[i][j] = abs ((int) Iy[i][j]);
		g =	Ix[i][j] >	Iy[i][j] ?	Ix[i][j] :	Iy[i][j];
		*/

		I[i][j]	= (byte) (g	< 0. ? 0 : (g >	255	? 255 :	g))	/1;
	}

	// Clean up
	FREE_MAT_FLOAT (im.xsize, A);
	FREE_MAT_FLOAT (im.xsize, Ix);
	FREE_MAT_FLOAT (im.xsize, Iy);
	FREE_VECT_FLOAT(yp);
	FREE_VECT_FLOAT(ym)	;
}


// *************************************************************
// Filter the image with a Sobel edge filter
// The static engine function, is not called externally
// *************************************************************

static byte ** _sobel (Image &im, bool horiz, bool	vertical) 
{
	byte **J;
	int	gx,gy,g;
	int	nbFilters;

	// The number of filters to	weight the results
	nbFilters=0;
	if (horiz) ++nbFilters;
	if (vertical) ++nbFilters;
	if (nbFilters==0)
		return NULL;

	J =	CREATE_IMAGE (im.ysize, im.xsize);

	// Process the border
	for	(int i = 0 ; i < im.xsize ; i++)
		J[i][0]	= J[i][im.ysize-1] =	0 ;
	for	(int i = 0 ; i < im.ysize ; i++)
		J[0][i]	= J[im.xsize-1][i] =	0 ;

	// Travers all pixels
	for	(int i = 1 ; i < im.xsize-1 ; i++) 
	for	(int j = 1 ; j < im.ysize-1 ; j++) 
	{
		if (horiz) 
		{
			gx =  im.R[i+1][j+1] +	im.R[i+1][j-1]
				-	im.R[i-1][j+1] -	im.R[i-1][j-1]
				+2*( im.R[i+1][j]	-	im.R[i-1][j]) ;

			gx /= 4;
		}
		else
			gx = 0;

		if (vertical) 
		{
			gy =  im.R[i+1][j]	+	im.R[i-1][j+1]
				-	im.R[i-1][j-1] -	im.R[i+1][j-1]
				+2*( im.R[i][j+1]	-	im.R[i][j-1]) ;

			gy /= 4;
		}
		else

			gy = 0;

		g =	(abs(gx) + abs(gy))	/ nbFilters;
		if (g>255) g = 255;

		J[i][j]	= g;
	}

	return J;
}

// *************************************************************
// Filter the image with a Sobel edge filter
// Apply results to the current image
// *************************************************************

void filterSobel (Image &im, bool horiz, bool vertical) {
    if (im.type==3) {
		cerr << "*** Internal ERROR in filterSobel():\n"
			"Cannot process color image!!!!\n";
		exit (1);
	}
	byte **J = _sobel (im, horiz, vertical);
	if (J!=NULL) {
		FREE_IMAGE (im.R);
		im.R =	J;
	}
}

// *************************************************************
// Filter the image with a Sobel edge filter
// *************************************************************

Image *filterSobelReturn (Image &im, bool horiz, bool vertical) 
{
	Image *ret;
	byte **J;

	if (im.type==3) 
	{
		cerr << "*** Internal ERROR in filterSobel()\n"
			"Cannot process color image!!!!\n";
		exit (1);
	}

	ret = new Image;
	J = _sobel (im, horiz, vertical);
	if (J==NULL) 
	{
		delete ret;
		return NULL;
	}

	ret->R = J;
	ret->type = 2;
	ret->xsize = im.xsize;
	ret->ysize = im.ysize;

	return ret;
}

// *************************************************************
// Filter the image	with a given filter	mask
// Border treatment: Recopying of the border pixels.
// *************************************************************

void filterMask (Image &im, FilterMask &fm) 
{
	byte **FIL[3];
	int	x, y, ix, iy, bx, fx, fy, by;
	int trimed_x, trimed_y;
	int	planes;
	double newval;
	int inewval;

	// How many	color planes?
	planes = im.nbColorPlanes();

	for	(int pl=0; pl<planes; ++pl)
		FIL[pl]	= CREATE_IMAGE (im.ysize, im.xsize);

	for	(y=0; y<im.ysize; ++y ) 
	{
		for	(x=0; x<im.xsize; ++x ) 
		{
			/* Calculate the filter	boundaries */
			bx = x-(fm.xsize/2);
			by = y-(fm.ysize/2);

			// Traverse	all	colour planes
			for	(int pl=0; pl<planes; ++pl)	
			{
				/* convolve	*/
				newval = 0;
				for	(iy=by,	fy=0; fy<fm.ysize; ++iy, ++fy) 
				{
					for	(ix=bx,	fx=0; fx<fm.xsize; ++ix, ++fx) 
					{
						trimed_x = (ix < 0 ? 0 : ix);
						trimed_y = (iy < 0 ? 0 : iy);
						trimed_x = (trimed_x >= im.xsize ? im.xsize-1 : trimed_x);
						trimed_y = (trimed_y >= im.ysize ? im.ysize-1 : trimed_y);
						newval += ((double)im.get(pl+1,trimed_x,trimed_y)*fm.get(fx,fy));
					}
				}
				inewval = (int) rint(newval);
				if (inewval<0) inewval = 0;
				if (inewval>255) inewval = 255;
				FIL[pl][x][y] =	inewval;
			}
		}
	}

 	FREE_IMAGE(im.R);
	im.R =	FIL[0];


	if (im.type==3)
	{
		FREE_IMAGE(im.G);
		FREE_IMAGE(im.B);
		im.G =	FIL[1];
		im.B =	FIL[2];
	}
}

// *************************************************************
// Filter the image	with a median Filter
// *************************************************************

void filterMedian (Image &im, int plane, int fsize) 
{
	filterGeneralizedMedian (im, plane, fsize, 0.5, false, false);
}

// *************************************************************
// Filter the image	with a kind of parametrable median filter,
// which does not peak the middle of the sorted values, but
// the one determined by a parameter.
// *************************************************************

void filterGeneralizedMedian (Image &im, int plane, int fsize, double k, bool horiz, bool vert) 
{
	byte *arr;
	byte **FIL;
	int	i, x, y, lx, ly, bx, by, ex, ey, curSort, valCount;
	byte foo;
	double fi;
	int index;

	if (horiz && vert) {
		cerr << "Internal error in filterGeneralizedMedian():\n"
			 << "The filter cannot be horizontal AND vertical!\n";
		exit (1);
	}

	arr	= new byte[fsize*fsize];
	FIL	= CREATE_IMAGE (im.ysize, im.xsize);

	/* Calculate the filter	boundaries */
	if (horiz || vert)
		valCount = fsize;
	else
		valCount = fsize*fsize;
	
	fi = k*(double)valCount;
	if ((fsize*fsize)%2==0)
		index = (int) floor(fi)-1;
	else
		index = (int) floor(fi);

	// Border treatment
	for	(y=0; y<im.ysize; ++y) 
	for	(x=0; x<=fsize/2; ++x) 	
	{
		FIL[x][y]=im.get(x,y);
		FIL[im.xsize-x-1][y]=im.get(im.xsize-x-1,y);		
	}
	for	(x=0; x<im.xsize; ++x) 
	for	(y=0; y<=fsize/2; ++y) 	
	{
		FIL[x][y]=im.get(x,y);
		FIL[x][im.ysize-y-1]=im.get(x,im.ysize-y-1);		
	}
	
	// The standard case
	for	(y=0+fsize/2; y<(im.ysize-fsize/2); ++y ) 
	for	(x=0+fsize/2; x<(im.xsize-fsize/2); ++x )
	{

		/* Calculate the filter	boundaries */
		if (horiz) 
		{
			bx = x-(fsize/2);
			ex = x+(fsize/2);
			by = y;
			ey = y;
		}

		else if (vert) 
		{
			bx = x;
			ex = x;
			by = y-(fsize/2);
			ey = y+(fsize/2);
		}

		else 
		{
			bx = x-(fsize/2);
			ex = x+(fsize/2);
			by = y-(fsize/2);    			
			ey = y+(fsize/2);
		}						

		/* Read	the	values into	the	array */
		i=0;
		for	(ly=by;	ly<=ey;	++ly) 
		for	(lx=bx;	lx<=ex;	++lx) 
		{
			arr[i] = im.get(plane,	lx,ly);
			++i;
		}
		
		/* Sort	them with a	variation of bubble	sort (is fastest with
			* low data	amount,	and	we don't need to sort the whole
			* array
			*/
		for	(curSort=0;	curSort<=index;++curSort){
			for	(i=curSort;	i<valCount;	++i) {
				if (arr[curSort] > arr[i]) {
					foo	= arr[curSort];
					arr[curSort] = arr[i];
					arr[i] = foo;
				}
			}
		}
		
		FIL[x][y] =	arr[index];
	}
	
	if (plane == 1)	{
		FREE_IMAGE(im.R);
		im.R =	FIL;
	}
	else
		if (plane == 2)	{
			FREE_IMAGE(im.G);
			im.G =	FIL;
		}
		else {
			FREE_IMAGE(im.B);
			im.B =	FIL;
		}

	delete [] arr;
}

// *************************************************************
// Perform a histogram stretch
// indexPerc: if positiv, then ignore the first indexPerc percent
// data left and right (robust stretch).
// *************************************************************

void histogramStretch (Image &im, int plane, float indexPercLeft,
	float indexPercRight) 
{
	int	minv,maxv;
	Histogram<unsigned int> *H;

	if ((indexPercLeft>0)||(indexPercRight>0))
		H = buildHistogram<unsigned int>(im,false,0);

	if (indexPercLeft<=0)
		minv = imageMin(im, plane);		
	else
		minv = H->getIndexPercLeft (indexPercLeft);				
	
	if (indexPercRight<=0)	
		maxv = imageMax(im, plane);
	else
		maxv = H->getIndexPercRight (indexPercRight);
		
	if ((indexPercLeft>0)||(indexPercRight>0))
		delete H;

	cerr <<	"min: "	<< minv	<< ", max: " <<	maxv <<	"\n";

	if ((minv==0)&&(maxv==255))
		return;

	for	(int y = 0 ; y<	im.ysize; y++) 	
	for	(int x = 0 ; x < im.xsize; x++)	
	{
		im.set (plane, x, y, (byte) scale_value (im.get(plane,x,y), minv,	maxv, 0, 255) /	1) ;
	}
}


// *************************************************************
// Brighten the image
// *************************************************************

void brightenImage (Image &im) 
{
   	int imax = imageMax(im, 1);
   	int mul	= imax > 0 ? (255 / imax)	: 0;
   	if (mul>1)
   		multiply (im, PLANE_RED, mul);
}

