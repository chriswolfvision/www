// *************************************************************
// Image_Math.cc
// Mathematical routines for Images
//
// author: Christian Wolf
// *************************************************************

// C
#include <stdlib.h>
#include <stdio.h>
#include <errno.h>
#include <string.h>
#include <math.h>

// From the IMAGE library
#include <Image.h>
#include <color.h>
#include <FloatMatrix.h>
#include <FloatMatrixColor.h>

// From this library
#include "ImageProc.h"


// *************************************************************
// Subtract another image and take the abs value
// was operator -=
// *************************************************************

void absoluteImageDiff (Image &im, const Image & other)	{

	for	(int i = 0 ; i < im.xsize; i++)
	for	(int j = 0 ; j < im.ysize; j++)	{
		im.R[i][j]	= abs (im.R[i][j] - other.R[i][j]);
		if (im.type==3) {
			im.G[i][j]	= abs (im.G[i][j] - other.G[i][j]);
			im.B[i][j]	= abs (im.B[i][j] - other.B[i][j]);
		}
	}
}

// *************************************************************
// Add another image with upper bound 255
// *************************************************************

void add (Image &im, const Image & other)	
{
	int foo;
	for	(int i = 0 ; i < im.xsize; i++)
	for	(int j = 0 ; j < im.ysize; j++)	
	{
		foo	= im.R[i][j] + other.R[i][j];
		im.R[i][j] = foo > 255 ? 255 : foo;

		if (im.type==3) {
			foo	= im.G[i][j] + other.G[i][j];
			im.G[i][j] = foo > 255 ? 255 : foo;
			foo = im.B[i][j] + other.B[i][j];
			im.B[i][j] = foo > 255 ? 255 : foo;
		}
	}
}

// *************************************************************
// Subtract another image with lower bound 0
// *************************************************************

void subtract (Image &im, const Image & other)	
{
	int foo;
	for	(int i = 0 ; i < im.xsize; i++)
	for	(int j = 0 ; j < im.ysize; j++)	{
		foo	= im.R[i][j] - other.R[i][j];
		im.R[i][j] = foo < 0 ? 0 : foo;

		if (im.type==3) {
			foo	= im.G[i][j] - other.G[i][j];
			im.G[i][j] = foo < 0 ? 0 : foo;
			foo = im.B[i][j] - other.B[i][j];
			im.B[i][j] = foo < 0 ? 0 : foo;
		}
	}
}

// *************************************************************
// Add a fixed value
// *************************************************************

void addValue (Image &im, int a)	{

	for	(int i = 0 ; i < im.xsize ; i++)
	for	(int j = 0 ; j < im.ysize ; j++)	{
		if (((int) im.R[i][j])+a >	255)
			im.R[i][j]	= 255;
		else
			im.R[i][j]	+= a;
	}
}

// *************************************************************
// Subtract a fixed value
// *************************************************************

void subValue (Image &im, int a)	{

	for	(int i = 0 ; i < im.xsize ; i++)
	for	(int j = 0 ; j < im.ysize ; j++)	{
		if (((int) im.R[i][j])-a <	0)
			im.R[i][j]	= 0;
		else
			im.R[i][j]	-= a;
	}
}

// *************************************************************
// Transpose the image
// *************************************************************

void transposeImage (Image &im) {
	Image *copy;

	copy = new Image(im);

	im.resizeWithDestroy (copy->ysize,	copy->xsize, copy->type);

	for	(int y=0; y<im.ysize; ++y)
		for	(int x=0; x<im.xsize; ++x)
			im.set	(1,	x, y, copy->get(1, y, x));

	delete copy;
}

// *************************************************************
// Mirror the image	in x direction
// *************************************************************

void mirrorX (Image &im) {
	byte temp;

	for	(int y=0; y<im.ysize; ++y) {
		for	(int x=0; x<im.xsize/2; ++x)	{

			temp = im.get (1, x, y);
			im.set	(1,	x, y, im.get(1, im.xsize-x-1, y));
			im.set	(1,	im.xsize-x-1, y,	temp);

			if (im.type == 3) {
				temp = im.get (2, x, y);
				im.set	(2,	x, y, im.get(2, im.xsize-x-1, y));
				im.set	(2,	im.xsize-x-1, y,	temp);

				temp = im.get (1, x, y);
				im.set	(3,	x, y, im.get(1, im.xsize-x-1, y));
				im.set	(3,	im.xsize-x-1, y,	temp);
			}
		}
	}
}

// *************************************************************
// Mirror the image	in Y direction
// *************************************************************

void mirrorY (Image &im) {
	byte temp;

	for	(int y=0; y<im.ysize/2; ++y)	{
		for	(int x=0; x<im.xsize; ++x) {

			temp = im.get (1, x, y);
			im.set	(1,	x, y, im.get(1, x,	im.ysize-y-1));
			im.set	(1,	x, im.ysize-y-1,	temp);

			if (im.type == 3) {
					temp = im.get (2,x,y);
					im.set	(2,	x, y, im.get(2, x,	im.ysize-y-1));
					im.set	(2,	x, im.ysize-y-1,	temp);

					temp = im.get (3,x,y);
					im.set	(3,	x, y, im.get(3, x,	im.ysize-y-1));
					im.set	(3,	x, im.ysize-y-1,	temp);
			}
		}
	}

}

// *************************************************************
// Multiply	with a fixed value
// *************************************************************

void multiply (Image &im, double m) {
	multiply (im, PLANE_RED, m);
	if (im.type==3) {
		multiply (im, PLANE_BLUE, m);
		multiply (im, PLANE_GREEN, m);
	}
}

void multiply (Image &im, byte plane, double m) {
	byte **I = im.getPlane (plane);
	int v;

	for	(int i = 0 ; i < im.xsize ; i++)
	for	(int j = 0 ; j < im.ysize ; j++)	{
		if (((int) I[i][j])*m >	255)
			I[i][j]	= 255;
		else {
			v = (int) rint (((double) I[i][j]) * m);
			if (v<0) v=0;
			if (v>255) v = 255;
			I[i][j]	= v;
		}
	}
}

// *************************************************************
// Multiply	two images
// *************************************************************

void multiplysqrt (Image &ii, Image &ij) {

	if (ii.type==3) {
     	cerr << "Image::multiply(): color image not supported!";
      	exit (1);
 	}
    for	(int i = 0 ; i < ii.xsize ; i++)
	for	(int j = 0 ; j < ii.ysize ; j++)	{
      	float vi = ii.get(PLANE_RED,i,j);
		float vj = ij.get(PLANE_RED,i,j);
  		float vc = sqrt(vi*vj);
    	int ic = (int) rint(vc);
		if (ic<0) ic=0;
  		if (ic>255) ic=255;
    	ii.set(PLANE_RED,i,j,ic);
	}
}

// *************************************************************
// Reverse Video (revvid function of xv)
// *************************************************************

void reverseVideo (Image &im) 
{
	for (int y=0; y<im.ysize; ++y)
	for (int x=0; x<im.xsize; ++x) 
	{
		im.set(PLANE_RED,x,y, 255-im.get(PLANE_RED,x,y));
		if (im.nbColorPlanes()==3) 
		{
			im.set(PLANE_GREEN,x,y, 255-im.get(PLANE_GREEN,x,y));
			im.set(PLANE_BLUE,x,y, 255-im.get(PLANE_BLUE,x,y));
		}
	}
}

// *************************************************************
// Project the pixels of an image into the blue plane only.
// Used to visualize superimposed data.
// *************************************************************

void projectToBlue (Image &im, float darken) {
	bool hascolor;

	if (im.type==3)
		hascolor = true;
	else {
		im.convertGrayScale2R__();
		hascolor = false;
	}

	for (int y=0; y<im.ysize; ++y)
	for (int x=0; x<im.xsize; ++x) {
		if (hascolor)
			im.B[x][y]	= (int)	((0.299*im.R[x][y] + 0.587*im.G[x][y] + 0.114*im.B[x][y]) / darken);
		else
			im.B[x][y] = (byte) (im.R[x][y] / darken);

		im.R[x][y] = 0;
		im.G[x][y] = 0;
	}
}

// *************************************************************
// Divide by a fixed value
// *************************************************************

void divide (Image &im, int plane, double d) {
	int foo;

	for	(int i = 0 ; i < im.xsize ; i++)
	for	(int j = 0 ; j < im.ysize ; j++)	{
		foo = ( d==0 ? 255 : (int) rint (((double) im.get (plane,i,j)) / d));

		if (foo > 255)
			im.set (plane, i, j, 255);
		else
			im.set (plane, i, j, foo);
	}
}

void divide (Image &im, double d) {
	divide (im, 1,d);
	if (im.type==3) {
		divide (im, 2,d);
		divide (im, 3,d);
	}
}


// *************************************************************
// Calculate the power of the image
// *************************************************************

void power (Image &im, double alpha) {
	double nv;

	for	(int i = 0 ; i < im.xsize ; i++)
	for	(int j = 0 ; j < im.ysize ; j++)	{
		nv = pow (im.R[i][j], alpha);
		if (nv < 0)	nv = 0;
		if (nv > 255) nv = 255;
		im.R[i][j]	= (byte) nv	/ 1;
	}
}

// *************************************************************
// Gamma Correction
// *************************************************************

void gammaCorrection (Image &im, int plane, double gamma) {

	if (gamma==0)
		return;

	for	(int i = 0 ; i < im.xsize ; i++)
	for	(int j = 0 ; j < im.ysize ; j++)	{
		im.set	(plane,	i, j, gammaCorrection4GrayValue	(im.get(plane,i,j),gamma));
	}
}


// *************************************************************
// Take	the	maximum	of two images
// *************************************************************

void maxImage (Image &im, Image &other) {
	bool iscolor=false;

	if ((im.type!=3) && (other.type==3)) {
		im.convertGrayScale2R__();
		im.setPlaneToValue (PLANE_GREEN, 0);
		im.setPlaneToValue (PLANE_BLUE, 0);
	}

	if (im.type==3 && other.type==3)
		iscolor=true;

	for	(int i = 0 ; i < im.xsize ; i++)	{
		for	(int j = 0 ; j < im.ysize ; j++)	{

			im.R[i][j]	= other.R[i][j]	> im.R[i][j] ?
				other.R[i][j] :	im.R[i][j];

			// Color image
			if (iscolor) {
				im.G[i][j]	= other.G[i][j]	> im.G[i][j] ?
					other.G[i][j] :	im.G[i][j];
				im.B[i][j]	= other.B[i][j]	> im.B[i][j] ?
					other.B[i][j] :	im.B[i][j];
			}
		}
	}
}

// *************************************************************
// Take	the	minimum	of two images
// *************************************************************

void minImage (Image &im, Image &other) {
	bool iscolor=false;

	if ((im.type!=3) && (other.type==3)) {
		im.convertGrayScale2R__();
		im.setPlaneToValue (PLANE_GREEN, 0);
		im.setPlaneToValue (PLANE_BLUE, 0);
	}

	if (im.type==3 && other.type==3)
		iscolor=true;

	for	(int i = 0 ; i < im.xsize ; i++)	{
		for	(int j = 0 ; j < im.ysize ; j++)	{

			im.R[i][j]	= other.R[i][j]	< im.R[i][j] ?
				other.R[i][j] :	im.R[i][j];

			// Color image
			if (iscolor) {
				im.G[i][j]	= other.G[i][j]	< im.G[i][j] ?
					other.G[i][j] :	im.G[i][j];
				im.B[i][j]	= other.B[i][j]	< im.B[i][j] ?
					other.B[i][j] :	im.B[i][j];
			}
		}
	}
}

// *************************************************************
// Calculate the minimum of	an image
// *************************************************************

byte imageMin (Image &im, byte plane) {
	byte min=im.get(plane,0,0);

	for	(int i = 0 ; i < im.xsize ; i++)	{
		for	(int j = 0 ; j < im.ysize ; j++)	{
			byte c = im.get (plane, i,	j);
			if (c<min)
				min	= c;
		}
	}
	return min;
}

// *************************************************************
// Calculate the maximum of	an image
// *************************************************************

byte imageMax (Image &im, byte plane) {
	byte max=im.get(plane,0,0);

	for	(int i = 0 ; i < im.xsize ; i++)	{
		for	(int j = 0 ; j < im.ysize ; j++)	{
			byte c = im.get (plane, i,	j);
			if (c>max)
				max	= c;
		}
	}
	return max;
}

// *************************************************************
// Calculate the sum of all the pixels of an image plane
// *************************************************************

double ImageSum	(Image &im, byte plane) {
	double sum=0;
	for	(int i = 0 ; i < im.xsize ; i++)
	for	(int j = 0 ; j < im.ysize ; j++)
		sum += im.get (plane, i,	j);
	return sum;
}

// *************************************************************
// glide a window across the image and
// create two maps: mean and standard deviation.
// *************************************************************

double calcLocalStats (Image &im, int winx, int winy,
	FloatMatrix *&out_m, FloatMatrix *&out_s) {

	FloatMatrix *map_m;
	FloatMatrix *map_s;
	double m,s,max_s, sum, sum_sq, foo;
	int wxh	= winx/2;
	int wyh	= winy/2;
	int x_firstth= wxh;
	int y_lastth = im.ysize-wyh-1;
	int y_firstth= wyh;
	double winarea = winx*winy;

	map_m = new FloatMatrix (im.xsize, im.ysize);
	map_s = new FloatMatrix (im.xsize, im.ysize);

	max_s = 0;
	for	(int j = y_firstth ; j<=y_lastth; j++) {

		// Calculate the initial window at the beginning of the line
		sum = sum_sq = 0;
		for	(int wy=0 ; wy<winy; wy++)
			for	(int wx=0 ; wx<winx; wx++) {
				foo = im.R[wx][j-wyh+wy];
				sum    += foo;
				sum_sq += foo*foo;
			}
		m  = sum / winarea;
		s  = sqrt ((sum_sq - (sum*sum)/winarea)/winarea);
		if (s > max_s)
			max_s = s;
		map_m->set(x_firstth, j, m);
		map_s->set(x_firstth, j, s);

		// Shift the window, add and remove	new/old values to the histogram
		for	(int i=1 ; i <= im.xsize-winx; i++) {

			// Remove the left old column and add the right new column
			for (int wy=0; wy<winy; ++wy) {
				foo = im.R[i-1][j-wyh+wy];
				sum    -= foo;
				sum_sq -= foo*foo;
				foo = im.R[i+winx-1][j-wyh+wy];
				sum    += foo;
				sum_sq += foo*foo;
			}
			m  = sum / winarea;
			s  = sqrt ((sum_sq - (sum*sum)/winarea)/winarea);
			if (s > max_s)
				max_s = s;
			map_m->set(i+wxh, j, m);
			map_s->set(i+wxh, j, s);
		}
	}

	out_m = map_m;
	out_s = map_s;
	return max_s;
}
