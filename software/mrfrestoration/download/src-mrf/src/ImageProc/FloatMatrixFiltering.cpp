/***************************************************************************
                          FloatMatrixFiltering.cpp  -  description
                             -------------------
    begin                : Fri Mar 14 2003
    copyright            : (C) 2003 by Christian Wolf
    email                : wolf@rfv.insa-lyon.fr
 ***************************************************************************/
 
// From the IMAGE module
#include <FloatMatrix.h>
#include <FloatMatrixColor.h>

// From the this module
#include "FilterMask.h"
#include "ImageProc.h"
 
// *************************************************************
// Filter a FloatMatrix	with a given filter	mask
// Seperate source and destination image!!!
// *************************************************************

static void _filterFloatArr (float **src, float **dst, int xsize, int ysize, FilterMask &fm) {
	int	x, y, ix, iy, bx, fx, fy, by;
	float newVal;
	int fmx = fm.xsize,
		fmy = fm.ysize;
	int fmx2 = fmx/2,
		fmy2 = fmy/2;
	int borderx = xsize-fmx2;
	int bordery = ysize-fmy2;
	int endx=xsize-1;
	int endy=ysize-1;

	// *************************************************************
	// The standard case
	// *************************************************************

	for	(y=fmy2; y<(ysize-fmy2); ++y ) {
		for	(x=fmx2; x<(xsize-fmx2); ++x ) {

			// Calculate the filter	boundaries
			bx = x-fmx2;
			by = y-fmy2;

			// convolve
			newVal = 0;
			for	(iy=by,	fy=0; fy<fmy; ++iy, ++fy) {
				for	(ix=bx,	fx=0; fx<fmx; ++ix, ++fx) {
					newVal += (src[ix][iy]*fm.get(fx,fy));
				}
			}

			dst[x][y] = newVal;
		}
	}

	// *************************************************************
	// Border treatment.
	// Basically the same code as the standard
	// case, except that the borders of the image are checked, and
	// if the access exceeds the borders, the border values are
	// taken.
	// *************************************************************

	// horizontal borders
	// ===================

	for (int y=0; y<ysize; ++y) {

		// The left border
   	    // ---------------

		for (int x=0; x<fmx2&&x<xsize; ++x) {

    		// Calculate the filter	boundaries
    		bx = x-fmx2;
    		by = y-fmy2;

    		// convolve
    		newVal = 0;
    		for	(iy=by,	fy=0; fy<fmy; ++iy, ++fy) {
    			int index_y = (iy<0 ? 0 : (iy >= ysize ? ysize-1 : iy));
    			for	(ix=bx,	fx=0; fx<fmx; ++ix, ++fx) {
    				int index_x = (ix<0 ? 0 : (ix >= xsize ? xsize-1 : ix));
    				newVal += (src[index_x][index_y]*fm.get(fx,fy));
    			}
    		}
    		dst[x][y] = newVal;

    	}

    	// The right border
    	// ---------------

    	for (int x=endx; x>=borderx&&x>=0; --x) {

    		// Calculate the filter	boundaries
    		bx = x-fmx2;
    		by = y-fmy2;

    		// convolve
    		newVal = 0;
    		for	(iy=by,	fy=0; fy<fmy; ++iy, ++fy) {
    			int index_y = (iy<0 ? 0 : (iy >= ysize ? ysize-1 : iy));
    			for	(ix=bx,	fx=0; fx<fmx; ++ix, ++fx) {
    				int index_x = (ix<0 ? 0 : (ix >= xsize ? xsize-1 : ix));
    				newVal += (src[index_x][index_y]*fm.get(fx,fy));
    			}
    		}
    		dst[x][y] = newVal;

    	}
	}

	// vertical borders
	// ===================

	for (int x=0; x<xsize; ++x)	 {

		// The upper border
		// -----------------


		for (int y=0; y<fmx2&&y<xsize; ++y) {

    		// Calculate the filter	boundaries

    		bx = x-fmx2;
    		by = y-fmy2;

    		// convolve
    		newVal = 0;
    		for	(iy=by,	fy=0; fy<fmy; ++iy, ++fy) {
    			int index_y = (iy<0 ? 0 : (iy >= ysize ? ysize-1 : iy));
    			for	(ix=bx,	fx=0; fx<fmx; ++ix, ++fx) {
    				int index_x = (ix<0 ? 0 : (ix >= xsize ? xsize-1 : ix));
    				newVal += (src[index_x][index_y]*fm.get(fx,fy));
    			}
    		}
    		dst[x][y] = newVal;
    	}

    	// The lower border
		// -----------------

		for (int y=endy; y>=bordery&&y>=0; --y) {

    		// Calculate the filter	boundaries
    		bx = x-fmx2;
    		by = y-fmy2;

    		// convolve
    		newVal = 0;
    		for	(iy=by,	fy=0; fy<fmy; ++iy, ++fy) {
    			int index_y = (iy<0 ? 0 : (iy >= ysize ? ysize-1 : iy));
    			for	(ix=bx,	fx=0; fx<fmx; ++ix, ++fx) {
    				int index_x = (ix<0 ? 0 : (ix >= xsize ? xsize-1 : ix));
    				newVal += (src[index_x][index_y]*fm.get(fx,fy));
    			}
    		}
    		dst[x][y] = newVal;
    	}
    }
}

// *************************************************************
// Filter a FloatMatrix	with a given filter	mask
// source and destination image are the same!!!
// *************************************************************

static float ** _filterFloatArr (float **arr, int xsize, int ysize, FilterMask &fm) {
	float **FIL;
	float *p;

	FIL	= FloatMatrix::_AllocPlane(xsize,ysize);

	_filterFloatArr (arr, FIL, xsize, ysize, fm);

 	p = *arr;
	delete [] p;
	delete [] arr;

	return FIL;
}

// *************************************************************
// A static function Filtering an image	 with a median Filter
// the array is given as pointer to the pointer, so the
// array can be replaced
// *************************************************************

static void _filterMedian (float **src, float **dst, int xsize, int ysize, int fsize) {
	int	i, x, y, lx, ly, bx, by, ex, ey, curSort, valCount;
	float foo;
	float *arr;

	arr	= new float	[fsize*fsize];

	valCount = fsize*fsize;

	for	(y=0+fsize/2; y<(ysize-fsize/2); ++y ) {
		for	(x=0+fsize/2; x<(xsize-fsize/2); ++x ) {

			/* Calculate the filter	boundaries */
			bx = x-(fsize/2);
			by = y-(fsize/2);
			ex = x+(fsize/2);
			ey = y+(fsize/2);

			/* Read	the	values into	the	array */
			i=0;
			for	(ly=by;	ly<=ey;	++ly) {
				for	(lx=bx;	lx<=ex;	++lx) {
					arr[i] = src[lx][ly];
					++i;
				}
			}

			/* Sort	them with a	variation of bubble	sort (is fastest with
			 * low data	amount,	and	we need	only to	sort half of the array,
			 * since only the value	in the middle is needed.
			 */

			for	(curSort=0;	curSort<=valCount/2	;++curSort)	{
				for	(i=curSort;	i<valCount;	++i) {
					if (arr[curSort] > arr[i]) {
						foo	= arr[curSort];
						arr[curSort] = arr[i];
						arr[i] = foo;
					}
				}
			}
			dst[x][y] = arr[valCount/2];
		}
	}

	delete arr;
}

// *************************************************************
// A static function Filtering an image	with a median Filter
// *************************************************************
/*
static float ** _filterMedian (float **out_image, int xsize, int ysize, int fsize) {
	float **FIL, *p;

	FIL	= FloatMatrix::_AllocPlane(xsize,ysize);

	_filterMedian (out_image, FIL, xsize, ysize, fsize);

	p = *out_image;
	delete [] p;
	delete [] out_image;

	return FIL;
}
*/

// *************************************************************
// The functions actually called from outside: FILTER MASK
// *************************************************************

void filterFloatMatrix (FloatMatrix &im, FilterMask &fm) 
{
	im.changePlane (_filterFloatArr (im.getPlane(), im.xsize, im.ysize, fm));
}

void filterFloatMatrix (FloatMatrix &src, FloatMatrix &dst, FilterMask &fm) 
{
	_filterFloatArr (src.getPlane(), dst.getPlane(), src.xsize, src.ysize, fm);
}

void filterFloatMatrixColor (FloatMatrixColor &im, FilterMask &fm, int plane) 
{
	im.changePlane(plane, _filterFloatArr (im.getPlane(plane), im.xsize, im.ysize, fm));
}

void filterFloatMatrixColor (FloatMatrixColor &im, FilterMask &fm) 
{
	im.changePlane(1, _filterFloatArr (im.getPlane(1), im.xsize, im.ysize, fm));
	im.changePlane(2, _filterFloatArr (im.getPlane(2), im.xsize, im.ysize, fm));
	im.changePlane(3, _filterFloatArr (im.getPlane(3), im.xsize, im.ysize, fm));
}

// *************************************************************
// The functions actually called from outside: MEDIAN
// *************************************************************

void filterFloatMatrixMedian (FloatMatrix &src, FloatMatrix &dst, int fsize) 
{
	_filterMedian (src.getPlane(), dst.getPlane(), src.xsize, src.ysize, fsize);
}

