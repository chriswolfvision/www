/**********************************************************
 * FloatMatrix.cc
 *
 * Christian Wolf, e9226297@stud1.tuwien.ac.at
 * Beginn: 17.6.1999
 **********************************************************/

// C
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

// C++
#include <iostream>

using namespace std;

// Own classes
#include "FloatMatrix.h"
#include "Image.h"
#include "color.h"

/************************************************************************
 * A absolute function for double values
 ************************************************************************/

static double d_abs (double a)	{
	return a < 0 ? -1*a	: a;
}

// *********************************************************************
// Constructors
// *********************************************************************

// Plain image
FloatMatrix::FloatMatrix ()	{
	data1 =	NULL;
}

// Non initialized but malloced	image
FloatMatrix::FloatMatrix (int sx, int sy) {
	this->xsize	= sx;
	this->ysize	= sy;
	data1 =	_AllocPlane (sx,sy);
}

// Malloced	image read from	file
FloatMatrix::FloatMatrix (char *filename) {
	data1 =	NULL;
	this->read (filename);
}

FloatMatrix::FloatMatrix (const Image &other)	{
	this->xsize	= other.xsize;
	this->ysize	= other.ysize;
    data1 =	_AllocPlane(xsize,ysize);

	for	(int y=0; y<ysize; ++y)	{
		for	(int x=0; x<xsize; ++x)	{
			set(x,y,other.get(1,x,y));
		}
	}
}

// *********************************************************************
// Copy	Constructor
// *********************************************************************

FloatMatrix::FloatMatrix (const FloatMatrix &other){
	xsize=other.xsize;
	ysize=other.ysize;
	data1 =	_AllocPlane (xsize,ysize);
	memcpy (*data1, *(other.data1), xsize*ysize*sizeof(float));
}

// *********************************************************************
// Destructor
// *********************************************************************

FloatMatrix::~FloatMatrix () {
	if (data1 != NULL) {
		float *p=*data1;
		delete [] p;
		delete [] data1;
	}
}

/**********************************************************************
 **********************************************************************
 **********************************************************************
 *
 * FILE FORMAT FOR THE FLOAT IMAGES:
 * =================================
 * Header in ASCII FORMAT
 * 1. Line: Magic number: "4081973"
 * 2. Line: Format "F2" for grayscale images, "F3" for color images
 * 3. Line: Dimensions "<xsize> x <ysize>"
 *
 * The rest of the data consists of binary written floats
 *
 **********************************************************************
 **********************************************************************
 **********************************************************************/

// *********************************************************************
// Read	image from a float file
// *********************************************************************

void FloatMatrix::read (const char *filename)	{
	FILE *fp;
	char buf[256];
	int xs, ys;

	if ((fp	= fopen(filename, "r"))	== NULL) {
		printf("Cannot open	file %s!\n",filename);
		exit(1);
	}

	if (fgets(buf,256,fp)!=NULL) {
		buf[strlen(MAGIC_NUMBER_FI)]='\0';
		if (strcmp(buf,MAGIC_NUMBER_FI)!=0) {
			cerr << "Cannot read floatfile " << filename << ": wrong magic number!\n"
				 << "[" << buf << "]" << endl;
			exit (1);
		}
	}

	if (fgets(buf,256,fp)!=NULL) {
		buf[2]='\0';
		if (strcmp(buf,"F2")!=0) {
			cerr << "Wrong format in float file!\n";
			exit (1);
		}
	}

	if (fgets(buf,256,fp)!=NULL) {
		sscanf(buf, "%d %d", &xs, &ys);

		if ((xs<=0) || (ys<=0)) {
			cerr << "Wrong dimensions in float file!\n";
			exit (1);
		}
	}

	xsize = xs;
	ysize = ys;
	data1 =	_AllocPlane(xs,ys);

	fread (*data1, xs*ys, sizeof(float),	fp);
	if (ferror(fp))	{
		printf ("Could not read	from file %s\n", filename);
		exit (1);
	}

	fclose (fp);
}

// *********************************************************************
// Read	image from a float file
// *********************************************************************

void FloatMatrix::write	(const char *filename) {
	FILE *fp;

	if ((fp	= fopen(filename, "w"))	== NULL) {
		printf("Cannot open	file %s!\n",filename);
		exit(1);
	}
	
	if (fprintf (fp,"%s\nF2\n%d %d\n",MAGIC_NUMBER_FI,xsize,ysize)<0) {
		cerr << "Error writing float image to file " << filename << endl;
		exit (1);
	}

	fwrite (*data1, xsize*ysize, sizeof(float), fp);
	if (ferror(fp))	{
		printf ("Could not read	from file %s\n", filename);
		exit (1);
	}

	fclose (fp);
}


// *********************************************************************
// Set all data1 to	zero
// *********************************************************************

void FloatMatrix::setZero () {
	float *p;
	
	if (data1==NULL)
		return;

	p = *data1;
	for	(int i=0; i<xsize*ysize; ++i) {
		p[i] = 0;
	}
}

// *********************************************************************
// Set all data1 to	a given value
// *********************************************************************

void FloatMatrix::setToValue (float val) {
	float *p;

	if (data1==NULL)
		return;

	p = *data1;
	for	(int i=0; i<xsize*ysize; ++i) {
		p[i] = val;
	}
}

// *********************************************************************
// Take the logarithm of all values
// *********************************************************************

void FloatMatrix::log() {
	float *p;
	
	if (data1==NULL)
		return;

	p = *data1;
	for	(int i=0; i<xsize*ysize; ++i) {
		float v = p[i];
		if (v)
			p[i] = ::log(v);
		else
			p[i] = 0;
	}
}

// *********************************************************************
// Take the square root of all values
// *********************************************************************

void FloatMatrix::sqrt() {
	float *p;
	
	if (data1==NULL)
		return;

	p = *data1;
	for	(int i=0; i<xsize*ysize; ++i) {
		float v = p[i];
		if (v>=0)
			p[i] = ::sqrt(v);
		else
			p[i] = 0;
	}
}

// *********************************************************************
// Multiply
// *********************************************************************

void FloatMatrix::multiply(float m) {
	float *p;
	
	if (data1==NULL)
		return;

	p = *data1;
	for	(int i=0; i<xsize*ysize; ++i) {
		p[i] = p[i]*m;
	}
}

// *************************************************************
// Calculate the absolute value
// *************************************************************

void FloatMatrix::abs() {
	float *p;
	if (data1==NULL)
		return;
	p = *data1;
	for	(int i=0; i<xsize*ysize; ++i)
			p[i] = d_abs(p[i]);
}

// *************************************************************
// Calculate the maximum
// *************************************************************

float FloatMatrix::max () {
	float max=get(0,0);
	for	(int i = 0 ; i < xsize ; i++)	{
		for	(int j = 0 ; j < ysize ; j++)	{
			float c = get (i, j);
			if (c>max)
				max	= c;
		}
	}
	return max;
}

// *************************************************************
// Calculate the minimum
// *************************************************************

float FloatMatrix::min () {
	float min=get(0,0);
	for	(int i = 0 ; i < xsize ; i++)	{
		for	(int j = 0 ; j < ysize ; j++)	{
			float c = get (i, j);
			if (c<min)
				min	= c;
		}
	}
	return min;
}

// *************************************************************
// Calculate the minimum of the absolute values
// *************************************************************

float FloatMatrix::minAbs () {
	float min=fabs(get(0,0));
	for	(int i = 0 ; i < xsize ; i++)	{
		for	(int j = 0 ; j < ysize ; j++)	{
			float c = fabs(get (i, j));
			if (c<min)
				min	= c;
		}
	}
	return min;
}

// *************************************************************
// Copy the contents of another one
// *************************************************************

void FloatMatrix::copy (FloatMatrix &other) {
	if (xsize!=other.xsize || ysize!=other.ysize) {
		cerr << "Internal error in FloatMatrix::copy()!!!\n";
		exit (1);
	}
	for	(int y=0; y<ysize; ++y)
	for	(int x=0; x<xsize; ++x)	{
		set(x,y, other.get(x,y));
	}
}

// *************************************************************
// Scale
// *************************************************************

void FloatMatrix::colorScale (double omin, double omax, double nmin, double nmax) {
	double color;

	for	(int y=0; y<ysize; ++y)	{
		for	(int x=0; x<xsize; ++x)	{
			color =	get	(x, y);
			color =	scale_value	(color,	omin, omax,	nmin, nmax);
			set	(x, y, color);
		}
	}
}

// *************************************************************
// Scale
// *************************************************************

void FloatMatrix::reverseVideo () 
{
	double max=this->max();

	for	(int y=0; y<ysize; ++y)
	for	(int x=0; x<xsize; ++x)	
	{
		set	(x, y, max-get(x,y));
	}
}

// *************************************************************
// ReScale the whole range
// *************************************************************

void FloatMatrix::colorReScale (double nmin, double nmax) {
	double omin, omax;
	omin = omax = get(0,0);
	for	(int y=0; y<ysize; ++y)	{
		for	(int x=0; x<xsize; ++x)	{
			float v = get(x,y);
			if (v<omin)
				omin = v;
			else
				if (v>omax)
					omax = v;
		}
	}
	colorScale (omin, omax, nmin, nmax);
}

// *************************************************************
// Help function for the interpolation method
// *************************************************************

inline double weight_func (int x, int x0) {
	double d = ::sqrt(x*x + x0*x0);
	return d == 0 ? 1 : 1/d;
}

// *************************************************************
// Bi-linear interpolation
// *************************************************************

FloatMatrix *FloatMatrix::interpolate (int factor) {
	FloatMatrix *newim;
	int nys=factor*ysize, nxs=factor*xsize;
	int old_x_l, old_x_r,
		old_y_u, old_y_d;
	int new_x_l, new_x_r,
		new_y_u, new_y_d;

#ifdef CHECK_CODE
	if (factor<1) {
		cerr <<	"Internal error in FloatMatrix::interpolate()!\n";
		CORE_DUMP;
	}
#endif

	newim = new FloatMatrix (nxs, nys);
	
	for (int y=0; y<nys; ++y) {
	
		double ry = ((double) y / (double) factor); 	
		
		old_y_u = (int) floor(ry);
		old_y_d = (int) ceil(ry);
		new_y_u = 4*old_y_u;
		new_y_d = 4*old_y_d;
		
		for (int x=0; x<nxs; ++x) {	
			double L=0;
			double w1, w2, w3, w4, weight_tot=0;
					
			double rx = ((double) x / (double) factor);

			// Fast special case if we hit the pixel exactly
			if (rx==rint(rx) && ry==rint(ry)) {
				newim->set(x,y, get((int)rx,(int)ry));
				continue;
			}
							
			old_x_l = (int) floor(rx);
			old_x_r = (int) ceil(rx);
			new_x_l = 4*old_x_l;
			new_x_r = 4*old_x_r;						
			
			if (old_x_r >= xsize) old_x_r = xsize-1;
			if (old_y_d >= ysize) old_y_d = ysize-1;
			
			// The weights
			w1 = weight_func (x-new_x_l, y-new_y_u);
			w2 = weight_func (x-new_x_r, y-new_y_u);
			w3 = weight_func (x-new_x_r, y-new_y_u);
			w4 = weight_func (x-new_x_l, y-new_y_d);
			weight_tot = w1 + w2 + w3 + w4;

			L = (w1 * get (old_x_l,old_y_u) +
				 w2 * get (old_x_r,old_y_u) +
                 w3 * get (old_x_r,old_y_d) +
                 w4 * get (old_x_l,old_y_d)) / weight_tot;   			    			
			newim->set(x,y,L);
	    }	
	}
	
	return newim;
}


// *********************************************************************
// Print the values	to the file, readable!
// *********************************************************************

void FloatMatrix::print	(FILE *fp) {

	fprintf	(fp, "MATRIX %d	x %d:\n", xsize, ysize);

	for	(int y=0; y<ysize; ++y)	{
		for	(int x=0; x<xsize; ++x)	{
			fprintf	(fp,"%5.2f",data1[x][y]);
		}
		fprintf	(fp, "\n");
	}

}

// *************************************************************
// Copy the border pixels
// *************************************************************

void FloatMatrix::copyBorderPixels(int border_x, int border_y) {
	float ul, ur, ll, lr;
	int x;

    // VERTICAL FLANKS
    for (x=border_x; x<xsize-border_x; ++x) {
 		float upper = get(x,border_y);
 		float lower = get(x,ysize-border_y-1);
 		for (int y=0; y<border_y; ++y) {
 			set (x,y, upper);
 			set (x,ysize-y-1, lower);

 		}
 	}

 	// HORIZONTAL FLANKS
 	for (int y=border_y; y<ysize-border_y; ++y) {
 		float left = get(border_x,y);
 		float right = get(xsize-border_x-1,y);
 		for (int x=0; x<border_x; ++x) {
 			set (x, y , left);
 			set (xsize-x-1, y, right);
 		}
 	}

 	ul = get (border_x, border_y);
 	ur = get (xsize-border_x-1, border_y);
 	ll = get (border_x, ysize-border_y-1);
 	lr = get (xsize-border_x-1, ysize-border_y-1);
 	for (x=0; x<border_x; ++x) {
	 	for (int y=0; y<border_y; ++y) {
			set (x,y, ul);
			set (xsize-1-x,y, ur);
			set (x,ysize-1-y, ll);
			set (xsize-1-x,ysize-1-y, lr);
	 	}
	}
}

// *************************************************************
// 
// *************************************************************

FloatMatrix *FloatMatrix::cutSubImage (Rect r, int growX, int growY) 
{
	FloatMatrix *im;
	Rect b;
	int width, height;

	b =	r;
	b.growAndClip (growX, growY, xsize, ysize);

	// Check the size
	width = b.width();

	height = b.height();
	if ((width==0) || (height==0))
		return NULL;
		
#ifdef CHECK_CODE
	if (r.right>=xsize || r.bottom>=ysize) {
		cerr << "Warning in FloatMatrix::cutSubImage(), "
			 << "has been corrected\nrect: " << r
			 << ", image size: " << xsize << "x" << ysize << endl;

		r.right=xsize-1;
		r.bottom=ysize-1;
	}
#endif		
	
	// Create the image		
	im = new FloatMatrix (width, height);
	im->setZero();
	
	// Copy	the	part of	the	big	original image
	for	(int y=0; y<b.height();	++y) 
	for	(int x=0; x<b.width(); ++x)	
		im->set	(x,y, get(b.left+x,b.top+y));
	
	return im;
}		

/****************************************************************
 * Order the coefficents after a FFT transform
 ****************************************************************/

void FloatMatrix::reOrderAfterFFT ()	{
	int	offset_horiz = xsize / 2;
	int	offset_vert	 = ysize/2;
	FloatMatrix copy	(*this);

	for	(int y=0; y<ysize/2; ++y) {
		for	(int x=0; x<xsize/2; ++x) {
			data1[x][y] =							copy.data1[x+offset_horiz][y+offset_vert];
			data1[x+offset_horiz][y+offset_vert] =	copy.data1[x][y];
			data1[x+offset_horiz][y] =				copy.data1[x][y+offset_vert];
			data1[x][y+offset_vert] =				copy.data1[x+offset_horiz][y];
		}
	}
}

// *************************************************************
// Calculate the power of the image
// *************************************************************

void FloatMatrix::power (double alpha) {
	for	(int i = 0 ; i < xsize ; i++)
		for	(int j = 0 ; j < ysize ; j++)
			set	(i,j,::pow(get(i,j), alpha));
}


// *********************************************************************
// Calculate the magnitude from the real and imaginary fft images
// *********************************************************************

FloatMatrix * FloatMatrix::FFTMagnitude (const FloatMatrix &im_image) {
	FloatMatrix *rv = new FloatMatrix (xsize, ysize);
	for	(int y=0; y<ysize; ++y)	{
		for	(int x=0; x<xsize; ++x)	{
			float re = get (x,y);
			float im = im_image.get (x,y);
			rv->set (x,y, ::sqrt (re*re+im*im));
		}
	}
	return rv;
}

// *********************************************************************
// Calculate the phase from the real and imaginary fft images
// *********************************************************************


FloatMatrix * FloatMatrix::FFTPhase (const FloatMatrix &im_image) {
	FloatMatrix *rv = new FloatMatrix (xsize, ysize);
	for	(int y=0; y<ysize; ++y)	{
		for	(int x=0; x<xsize; ++x)	{
			float re = get (x,y);
			float im = im_image.get (x,y);
			if (re==0)
				rv->set(x,y,0);
			else
				rv->set (x,y, atan(im/re));
		}
	}
	return rv;
}

// *********************************************************************
// GLOBAL FUNCTION - complex multiply
// *********************************************************************

void complexMultiplyFloat (
	FloatMatrix &inre1, FloatMatrix &inim1,
	FloatMatrix	&inre2,	FloatMatrix &inim2,
	FloatMatrix	&outre,	FloatMatrix &outim) {

	int	xsize =	inre1.xsize;
	int	ysize =	inre1.ysize;
	for	(int y=0; y	< ysize; y++) {
		for	(int x=0; x	< xsize; x++) {
			outre.data1[x][y] = inre1.data1[x][y]*inre2.data1[x][y] -
								inim1.data1[x][y]*inim2.data1[x][y];
			outim.data1[x][y] = inre1.data1[x][y]*inim2.data1[x][y] +
								inim1.data1[x][y]*inre2.data1[x][y];
		}
	}
}

