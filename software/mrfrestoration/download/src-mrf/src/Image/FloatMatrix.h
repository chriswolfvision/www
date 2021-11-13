/**********************************************************
 * FloatMatrix.h
 *
 * Christian Wolf, chriswolf@gmx.at
 * Beginn: 17.6.1999
 **********************************************************/

#ifndef	_WOLF_FLOATMATRIX_H_
#define	_WOLF_FLOATMATRIX_H_

// C
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

// C++
#include <iostream>

// From the main module
#include <CIL.h>

// From this module
#include "Rect.h"

#define MAGIC_NUMBER_FI "24081973FI"

class Image;
class FilterMask;

// is included _AFTER_ the class definition
class FloatMatrixIterator;

// **************************************************************
// The class
// **************************************************************						  

class FloatMatrix {

	friend
		void complexMultiplyFloat (FloatMatrix &inre1, FloatMatrix &inim1,
			FloatMatrix	&inre2,	FloatMatrix &inim2,
			FloatMatrix	&outre,	FloatMatrix &outim);

	// The iterator	has	full access	to the class
	friend 
		class FloatMatrixIterator;		

	public:
	
		typedef	FloatMatrixIterator iterator;
		typedef float PixelType;

		// Constructors and destructor
		FloatMatrix	();
		FloatMatrix	(int sx, int sy);
		FloatMatrix	(char *filename);
		FloatMatrix (const FloatMatrix &other);
		FloatMatrix (const Image &other);
		template <class T> FloatMatrix (const T &other, int xs, int ys);
		virtual	~FloatMatrix ();

		// Access methods
		inline float get (int x, int y) const ;
		inline void  set (int x, int y, float val);
		inline void inc (int x, int y)			    { ++data1[x][y]; }
		virtual	void setZero ();
		void setToValue (float val);
		
		// **********************************************
		// Return iterators and access functions
		// **********************************************

		iterator iterPos (int x, int y);

		// Access and change the whole data plane
		inline float ** getPlane();
		inline void changePlane (float **data);

		// Data manipulation
		template <class T> void paste (T &other, int xpos, int ypos);
		template <class T> void paste (T &other,int dxpos,int dypos,int sxpos,int sypos,int xlen,int ylen);
		void copy (FloatMatrix &other);
		void copyBorderPixels (int border_x, int border_y);
		FloatMatrix *cutSubImage (Rect r, int growX, int growY);

		// Mathematics
		void  log();
		void  sqrt();
		void  abs();
		float max();
		float min();
		float minAbs();
		void  multiply(float m);
		void power (double alpha);
		
		// Image processing
		void reverseVideo();
		void colorScale  (double omin, double omax, double nmin, double nmax);
		void colorReScale(                          double nmin, double nmax);
  		FloatMatrix *interpolate (int factor);

    	// Operations in the Fourier domain
		FloatMatrix * FFTMagnitude (const FloatMatrix &im);
		FloatMatrix * FFTPhase (const FloatMatrix &im);
		void reOrderAfterFFT ();

		// for debug
		void printPointer(char *s) {
			cerr << "debug: " << s << ": "
			     << (int) data1 << " / "
			     << (int) *data1 << endl;
		}

    	// Input/output
		virtual	void read  (const char *filename);
		virtual	void write (const char *filename);
		virtual	void print (FILE *fp);				        
		
	public :	// CLASS METHODS

    	static float ** _AllocPlane (int ysize, int xsize);
	
	public:		// DATA
		
		int	xsize, ysize;		
	
	protected:	// DATA
		
		float **data1;
};

#include "FloatMatrixIterator.h"

// **************************************************************
// Access and change the whole data plane
// **************************************************************

inline float ** FloatMatrix::getPlane() {
	return data1;
}

inline void FloatMatrix::changePlane (float **d) {
	if (data1!=NULL) {
		float *p=*data1;
		delete [] p;
		delete [] data1;
	}
	data1=d;
}

// **************************************************************
// Global Functions
// **************************************************************

void complexMultiplyFloat (FloatMatrix &inre1, FloatMatrix &inim1,
						  FloatMatrix &inre2, FloatMatrix &inim2,
						  FloatMatrix &outre, FloatMatrix &outim);

// **************************************************************
// Access methods
// **************************************************************

inline float FloatMatrix::get (int x, int y) const {
#ifdef PEDANTIC_CHECK_CODE
	if (x<0 || x>=xsize || y<0 || y>=ysize) {
		cerr << "Coords out of bounds in FloatMatrix::get [" << xsize << "," << ysize << "]: x="
			 << x << " y=" << y << endl;
		CORE_DUMP;
	}
#endif
	return data1[x][y];
}
inline void FloatMatrix::set (int x, int y, float val) {
#ifdef PEDANTIC_CHECK_CODE
	if (x<0 || x>=xsize || y<0 || y>=ysize) {
		cerr << "Coords out of bounds in FloatMatrix::set [" << xsize << "," << ysize << "]: x="
			 << x << " y=" << y << endl;
		CORE_DUMP;
	}
#endif
	data1[x][y] = val;
}
				  
// *************************************************************
// Constructor:
// Create a FloatMatrix from a byte image, but with a different
// size:
//
// Resample the image to fit it to a size for FFT -
// A very simple averaging method.
// *************************************************************

template <class T>
FloatMatrix::FloatMatrix (const T &other, int xs, int ys)
	:FloatMatrix (xs, ys)	{

	FloatMatrix *cnt;

	// Init
	cnt = new FloatMatrix (xs, ys);
	this->setZero();
	 cnt->setZero();

	// Travers the pixels of the source image and add them
	// to their corresponding destination locations.
	for (int y=0; y<ysize; ++y) {
		int ny = (int) rint (y * ys / ysize);
		if (ny>=ysize)
			ny=ysize-1;

		for (int x=0; x<xsize; ++x) {
			int nx = (int) rint (x * xs / xsize);
			if (nx>=xsize)
				nx=xsize-1;

			set(nx,ny, get(nx,ny) + other.get(x,y));
			cnt->inc(nx,ny);
		}
	}

	// Divide by the count to get the average value
	for (int y=0; y<ys; ++y) {
		for (int x=0; x<xs; ++x) {
			float c = cnt->get(x,y);
			set(x,y, (c==0 ? 0 : get(x,y) / c));
		}
	}

	// Clean up
	delete cnt;
}

// *************************************************************
// Paste another Image or FloatMatrix at a specified position
// - The whole image
// *************************************************************

template <class T>
void FloatMatrix::paste (T &other, int xpos, int ypos) {
	paste (other, xpos, ypos, 0, 0, other.xsize, other.ysize);
}

// *************************************************************
// Paste another Image or FloatMatrix at a specified position
// - Parts of the image
// *************************************************************

template <class T>
void FloatMatrix::paste (T &other,int dxpos,int dypos,int sxpos,int sypos,int xlen,int ylen) {
	int rxlen, rylen;

#ifdef CHECK_CODE
	if ((dxpos >= xsize) || (dypos >= ysize)) {
		cerr <<	"Error in template FloatMatrix::paste()!!!\n"
			 << "dxpos: " << dxpos << " xsize: " << xsize << endl
			 << "dypos: " << dypos << " ysize: " << ysize << endl;
		cerr << "Provoking core dump...\n";
		CORE_DUMP;
	}
#endif

	// Check if we have enough space in the image.
	// If not, cut the source image
	rxlen = xsize-dxpos;
	rylen = ysize-dypos;
	rxlen = rxlen < xlen ? rxlen : xlen;
	rylen = rylen < ylen ? rylen : ylen;

	for	(int y=0; y<rylen; ++y) {
		for	(int x=0; x<rxlen; ++x) {
			set(dxpos+x,dypos+y, (float) other.get(sxpos+x,sypos+y));
		}
	}
}

#endif

