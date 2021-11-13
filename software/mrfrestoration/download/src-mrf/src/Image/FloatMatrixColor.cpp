/**********************************************************
 * FloatMatrixColor.cpp
 *
 * Christian Wolf
 * wolf@rf
 **********************************************************/

#include "stdlib.h"
#include "stdio.h"
#include "math.h"
#include <string>
#include <string.h>

using namespace std;

#include "FloatMatrixColor.h"
#include "Image.h"
#include "color.h"

// *********************************************************************
// Constructors
// *********************************************************************

// Plain image
FloatMatrixColor::FloatMatrixColor ()	{
	data1 =	NULL;
	data2 =	NULL;
	data3 =	NULL;
}

FloatMatrixColor::FloatMatrixColor (Image	&other)	{
	_alloc (other.xsize, other.ysize);

	for	(int y=0; y<ysize; ++y)	{
		for	(int x=0; x<xsize; ++x)	{
			set(1,x,y,other.get(1,x,y));

			if (other.type == 3) {
				set(2,x,y,other.get(2,x,y));
				set(3,x,y,other.get(3,x,y));
			}
		}
	}
}

// *********************************************************************
// Copy constructor
// *********************************************************************

FloatMatrixColor::FloatMatrixColor (const FloatMatrixColor &other)	{
    _alloc (other.xsize, other.ysize);

	if ((xsize==0)||(ysize==0))
    	return;

   	memcpy (*data1, *(other.data1), xsize*ysize*sizeof(float));
    memcpy (*data2, *(other.data2), xsize*ysize*sizeof(float));
    memcpy (*data3, *(other.data3), xsize*ysize*sizeof(float));
}

// *********************************************************************
// Help	function for the constructors
// *********************************************************************

void FloatMatrixColor::_alloc (int sx, int sy) {
	xsize	= sx;
	ysize	= sy;
	data1 =	_AllocPlane(sx,sy);
	data2 =	_AllocPlane(sx,sy);
	data3 =	_AllocPlane(sx,sy);
}

// *********************************************************************
// Destructor
// *********************************************************************

FloatMatrixColor::~FloatMatrixColor () {
	if (data2 != NULL) {
		float *p=*data2;
		delete [] p;
		delete [] data2;
	}
	if (data3 != NULL) {
		float *p=*data3;
		delete [] p;
		delete [] data3;
	}
}

// *********************************************************************
// Read	image from a float file
// *********************************************************************

void FloatMatrixColor::read (char *filename,	int	sx,	int	sy)	{
	FILE *fp;

	this->xsize	= sx;
	this->ysize	= sy;

	if ((fp	= fopen(filename, "r"))	== NULL) {
		printf("Cannot open	file %s!\n",filename);
		exit(1);
	}
	data1 =	_AllocPlane(sx,sy);
	data2 =	_AllocPlane(sx,sy);
	data3 =	_AllocPlane(sx,sy);	

	fread (*data1, sx*sy, sizeof(float),	fp);
	if (ferror(fp))	{
		printf ("Could not read	from file %s\n", filename);
		exit (1);
	}
	fread (*data2, sx*sy, sizeof(float),	fp);
	if (ferror(fp))	{
		printf ("Could not read	from file %s\n", filename);
		exit (1);
	}
	fread (*data3, sx*sy, sizeof(float),	fp);

	if (ferror(fp))	{
		printf ("Could not read	from file %s\n", filename);
		exit (1);
	}

	fclose (fp);
}

// *********************************************************************
// Read	image from a float file
// *********************************************************************

void FloatMatrixColor::write	(char *filename) {
	FILE *fp;

	if ((fp	= fopen(filename, "w"))	== NULL) {
		printf("Cannot open	file %s!\n",filename);
		exit(1);
	}

	fwrite (*data1, xsize*ysize,	sizeof(float), fp);
	if (ferror(fp))	{
		printf ("Could not read	from file %s\n", filename);
		exit (1);
	}
	fwrite (*data2, xsize*ysize,	sizeof(float), fp);
	if (ferror(fp))	{
		printf ("Could not read	from file %s\n", filename);
		exit (1);
	}
	fwrite (*data2, xsize*ysize,	sizeof(float), fp);
	if (ferror(fp))	{
		printf ("Could not read	from file %s\n", filename);
		exit (1);
	}

	fclose (fp);
}

// *********************************************************************
// Set all data	to zero
// *********************************************************************

void FloatMatrixColor::setZero () {
	setZero(1);
	setZero(2);
	setZero(3);
}

// *************************************************************
// Combines	the	3 colour planes, weights them and stores them
// into	one	of the 3 planes.
// *************************************************************

void FloatMatrixColor::combine (int plane, double w1, double	w2,	double w3) {
	int	x,y;

	for	(y=0; y<ysize; ++y)
		for	(x=0; x<xsize; ++x)
			set(plane,x,y,w1*get(1,x,y)	+ w2*get(2,x,y)	+ w3*get(3,x,y));
}

// *************************************************************
// Convert into	LUV
// *************************************************************

void FloatMatrixColor::convertRGB2LUV ()	{
	int	x,y;
	double r,g,b,cl,cu,cv;

	init_color ();

	for	(y=0; y<ysize; ++y)	{
		for	(x=0; x<xsize; ++x)	{

			r =	get(1,x,y);
			g =	get(2,x,y);
			b =	get(3,x,y);
			r /= MAX_RED;
			g /= MAX_GREEN;
			b /= MAX_BLUE;

			rgb2luv	(r,g,b,&cl,&cu,&cv);

			set	(1,x,y,cl);
			set	(2,x,y,cu);
			set	(3,x,y,cv);
		}
	}
}

// *************************************************************
// Convert into	LUV
// *************************************************************

void FloatMatrixColor::convertRGB2LAB ()	{
	int	x,y;
	double r,g,b,cl,ca,cb;

	for	(y=0; y<ysize; ++y)	{
		for	(x=0; x<xsize; ++x)	{

			r =	get(1,x,y);
			g =	get(2,x,y);
			b =	get(3,x,y);
			r /= MAX_RED;
			g /= MAX_GREEN;
			b /= MAX_BLUE;

			rgb2lab	(r,g,b,&cl,&ca,&cb);

			set	(1,x,y,cl);
			set	(2,x,y,ca);
			set	(3,x,y,cb);
		}
	}
}


// *************************************************************
// Convert into	RGB
// *************************************************************

void FloatMatrixColor::convertLUV2RGB ()	{
	int	x,y;
	double r,g,b,l,u,v;

	init_color ();

	for	(y=0; y<ysize; ++y)	{
		for	(x=0; x<xsize; ++x)	{

			l =	get(1,x,y);
			u =	get(2,x,y);
			v =	get(3,x,y);

			luv2rgb	(l,u,v,&r,&g,&b);

			r *= MAX_RED;
			g *= MAX_GREEN;
			b *= MAX_BLUE;

			if (r>255) r=255;
			if (g>255) g=255;
			if (b>255) b=255;

			set	(1,x,y,r);
			set	(2,x,y,g);
			set	(3,x,y,b);
		}
	}
}

// ************************************************************************
// Conversion functions for the Gaussian color model
// After JM Geusebroek et al, "Color Invariants", PAMI 23(12)2001 pp1338-1350
// ************************************************************************

void FloatMatrixColor::convertRGB2GCM ()	{
	double r,g,b,e,el,ell;

	for	(int y=0; y<ysize; ++y)	{
		for	(int x=0; x<xsize; ++x)	{

			r =	get(1,x,y);
			g =	get(2,x,y);
			b =	get(3,x,y);

             	e = (0.06*r 	+ 0.63*g + 0.27*b ) / 255.0;
			el = (0.3*r 	+ 0.04*g - 0.35*b ) / 255.0;
			ell = (0.34*r- 0.6*g  + 0.17*b ) / 255.0;

			set	(1,x,y,e);
			set	(2,x,y,el);
			set	(3,x,y,ell);
		}
	}
}

// *************************************************************
// Scale a band;
// *************************************************************

void FloatMatrixColor::scalePlane (int plane, double	omin, double omax,
			double nmin, double	nmax) {

	double color;

	for	(int y=0; y<ysize; ++y)	{
		for	(int x=0; x<xsize; ++x)	{

			color =	get	(plane,	x, y);
			color =	scale_value	(color,	omin, omax,	nmin, nmax);
			set	(plane,	x, y, color);
		}
	}
}

// *************************************************************
// Set plane to	zero
// *************************************************************

void FloatMatrixColor::setZero (int plane) {
	for	(int y=0; y<ysize; ++y)
		for	(int x=0; x<xsize; ++x)
			set	(plane,x,y,0);
}

// *************************************************************
// Add a byte Image
// *************************************************************

void FloatMatrixColor::add(Image &i) {
	for (int p=1; p<=3; ++p)
	for (int y=0; y<ysize; ++y)
	for (int x=0; x<xsize; ++x) 
		set (p,x,y, get(p,x,y) + i.get(p,x,y));

}

// *************************************************************
// Subtract a byte Image
// *************************************************************

void FloatMatrixColor::sub (Image &i) {
	for (int p=1; p<=3; ++p)
	for (int y=0; y<ysize; ++y)
	for (int x=0; x<xsize; ++x)
		set (p,x,y, get(p,x,y) - i.get(p,x,y));

}

// *************************************************************
// Add a FloatMatrix Image -> partly
// *************************************************************

void FloatMatrixColor::add(FloatMatrixColor &other, int xoffset, int yoffset) {
	for (int y=0; y<other.ysize; ++y) {
		for (int x=0; x<other.xsize; ++x) {			
			set (1,x+xoffset,y+yoffset, get(1,x+xoffset,y+yoffset) + other.get(1,x,y));
			set (2,x+xoffset,y+yoffset, get(2,x+xoffset,y+yoffset) + other.get(2,x,y));
			set (3,x+xoffset,y+yoffset, get(3,x+xoffset,y+yoffset) + other.get(3,x,y));
		}
	}		  
}

// *************************************************************
// DEBUG PRINT
// *************************************************************

void FloatMatrixColor::debugPrint ()	{
	for	(int y=0; y<ysize; ++y)
		for	(int x=0; x<xsize; ++x)
		printf ("(%6.2f,%6.2f,%6.2f)", get(1,x,y), get(2,x,y), get(3,x,y));
}


