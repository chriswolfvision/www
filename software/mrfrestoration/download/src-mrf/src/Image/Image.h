#ifndef	_WOLF_IMAGE_H_
#define	_WOLF_IMAGE_H_

// C
#include <stdlib.h>
#include <stdio.h>

// C++
#include <functional>
#include <iostream>

// From the main module
#include <CIL.h>

// For the PORT to WIndows and VISUAL C++
#include <visualbug.h>

using namespace std;

typedef	unsigned char byte;

class FloatMatrix;
class FloatMatrixColor;
class FilterMask;
class CComponentList;
class Rect;

// is included _AFTER_ the class definition
class ImageIterator;

enum ImagePlane {
PLANE_RED=1,
PLANE_GREEN,
PLANE_BLUE
};

#define OCR_PAGE_WIDTH	2500
#define OCR_PAGE_HEIGHT	3500

// How to stack the sub images when a new image is created from a set
// of images. The method using markers is implemented in OCRImage!!!!
enum StackMethod {
	STACKMETHOD_VERTICAL=0,
	STACKMETHOD_HORIZONTAL,
	STACKMETHOD_VERT_TILED,
	STACKMETHOD_OCR
};

#define TRAVERS_PIXELS for (int y=0; y<ysize; ++y) for (int x=0; x<xsize; ++x)

class Image	{

	// The iterator	has	full access	to the class
	friend class ImageIterator;		

	public:

		typedef	ImageIterator iterator;
		typedef unsigned char PixelType;

		// **********************************************
		// Constructors, Destructor and Assignment
		// **********************************************

		// Constructors
		Image () { xsize=0;	ysize=0; type=0; }
		Image (int xs, int ys, int colorPlanes);
  		Image (int xs, int ys); // Default 1 plane
		Image (FloatMatrix &other);
		Image (FloatMatrixColor &other);
		Image (int xsize, int ysize, byte *one_dimensional_array, int nrPlanes);
		Image (byte	*one_dimensional_arrarys[],	long int xsize[], long int ysize[]);
		Image (const char *filename);
		Image (const char *filename, bool color);

		// A list of files - paste them all together
		Image (const char *filename, StackMethod method, int border, byte bordervalue, bool do_order);

		/*
		template <class T>
		Image (Matrix<T> &m);
		*/

		// Copy	constructor
		Image (const Image &other);

		// Destructor
		~Image ();

		// Assignment
		Image &	operator= (const Image &i);
		Image &	operator= (const FloatMatrixColor &other);

		// **********************************************
		// Iterators
		// **********************************************

		iterator iterPos (int x, int y);
		
		// **********************************************
		// Accessfunctions
		// **********************************************		

		// These acessors are "image" style, i.e.
		// the coordinate order is x,y
		inline byte get (int p, int x,	int	y) const;
		inline byte get (       int x,	int	y) const;
		inline void set (int p, int x,	int	y, byte	v);
		inline void set (       int x,	int	y, byte	v);
		
		// This acessor is "math-matrix" style, i.e.
		// the coordinate order is y,x		
		inline unsigned char & operator() (int y, int x);	
		
		int	nbColorPlanes()	{ return type==3 ? 3 : 1; }

		// **********************************************
		// I/O
		// **********************************************

		void read (const char *filename);
		void readGray (const char *filename);
		void readColor (const char *filename);
		void write (const char *filename);
		void writeGray (const char *filename);
		void writeColor	(const char *filename);

#ifdef UNIX_ARCHITECTURE		
		void writeToTempFile (char *filename);
#endif		

#ifdef HAVE_LIBJPEG		
		bool readJPEG (const char *filename);
		
		// Class method!!!!
		static bool check4JPEG (const char *filename);
#endif		

		// **********************************************
		// Conversion
		// **********************************************

		void convertRGB2GrayScale ();
		void convertR__2GrayScale ();
		void convertGrayScale2RGB ();
		void convertGrayScale2R__ ();
		void convertRGB2LUV	();
		void convertLUV2RGB	();

		// Same function, return the converted image.
		Image *convertRGB2GrayScaleReturn ();

		void copyBorderPixels(int border_x, int border_y);

		Image *cutSubImage (Rect r, int growX,	int	growY);


		// **********************************************
		// Connected Components
		// **********************************************

		// Mark	the	borders	of connected components	in a binary	image
		void componentBorders();

		// **********************************************
		// Data oriented operations
		// **********************************************

		void setZero (int plane);
		inline void setZero ();
		
		void resizeWithDestroy (int	xsize, int ysize, int colorPlanes);		
		
		// Decrease the size of an image by factor 2
		Image 			*reScale2x2 ();
		Image 			*subSample (int factor);

		// interpolate an image, i.e. increase its resolution
		Image           *interpolate (int factor);
		Image 			*interpolateBlocks (int factor);
		FloatMatrixColor *interpolateJolion (FloatMatrixColor *avgImage,
									 	    FloatMatrixColor *varImage,
									        int xoffset, int yoffset, int factor);
		FloatMatrixColor *interpolateBicubicRobust (FloatMatrixColor *avgImage,
									 	   		   FloatMatrixColor *varImage,
									               int xoffset, int yoffset, int factor);

        void crop (Rect &r, int growX, int growY);
        bool toNextPowerOf2 ();
		void copyPlane (int	dst, int src);
		void setPlaneToValue (int plane, byte val);
		void setBorderToValue (int plane, int bordersize, byte value);
		void superImpose (Image	&other,	int	val);
		void paste (Image &other, int xpos, int ypos);
		void paste (Image &other, int dxpos, int dypos,
		int sxpos, int sypos, int sxlen, int sylen);

		void grayvalueChange (byte oldr, byte newr);
		void colorChange (byte oldr, byte oldg,	byte oldb,
						  byte newr, byte newg,	byte newb);

		void separateColumns (int maxheight, int border, byte background);


		// **********************************************
		// Drawing stuff
		// **********************************************

		// Draw	with high contrast
		void mark (int x, int y);
		void mark (int plane, int x, int y, bool docross);
		// void mark (PixelList &pl, int hoffset, int voffset);

		void markDashes (int xb, int yb, int count, int dashlength);

		int drawText (int xb, int yb, const char *string, bool isBig);
		int drawText (int xb, int yb, int value, bool isBig);			

		void print (FILE *);

	public:

		// **********************************************
		// Class functions
		// **********************************************

		static bool isColorImage (const char *filename);

		inline byte** getPlane (int	dst) {
			return (dst==1 ? R : (dst == 2 ? G : B));
		}

	protected:

		// **********************************************
		// Help	functions
		// **********************************************

		void _alloc	(int xs, int ys, int colorPlanes);
		void _free ();
		void _copy (const Image	&other);
		void _copyFloat	(const FloatMatrix &other);
		void _copyFloat	(const FloatMatrixColor &other);

#ifdef UNIX_ARCHITECTURE		
		void _writeColor(int filedescriptor);
		void _writeGray (int filedescriptor);
#endif
		void _writeColor(FILE *fp);
		void _writeGray (FILE *fp);		

	public:

		unsigned char **R, **G,	**B;
	
	public:		// DATA

		int	xsize, ysize;
		int	type;

		
	public:			// CLASS MEMBERS
	
		static void CreateOCRPages (const char *filename, char *name_template, bool order, int border, byte bordervalue, bool doTile);

		static Image *GetDigits();
		static Image *GetDigitsBig();
		static int digitwidths[];
		static int digitwidthsBig[];
		static int digitdistance;		
		static int digitdistanceBig;
		
	private:		// CLASS MEMBERS
	
	
		static char *GetAppDir();
		static Image *_Digits;
		static Image *_DigitsBig;
};

// To sort the images by width if necessary
class image_xless : public binary_function<Image *, Image *, bool> {
	public:
		bool operator()(const Image	*x, const Image *y)	{
			return x->xsize < y->xsize;
		}
};

// To sort the images by height if necessary
class image_ymore : public binary_function<Image *, Image *, bool> {
	public:
		bool operator()(const Image	*x, const Image *y)	{
			return x->ysize > y->ysize;
		}
};

#include "ImageIterator.h"

unsigned char **CREATE_IMAGE (int ysize, int xsize);
void FREE_IMAGE	(unsigned char **im);
void CLEAR_IMAGE (int xsize, int ysize,	byte **im);

// **************************************************************
// ACESSORS
// These acessors are "image" style, i.e.
// the coordinate order is x,y
// **************************************************************

inline byte Image::get (int x, int	y) const 
{
#ifdef PEDANTIC_CHECK_CODE
	if (x<0 || x>=xsize || y<0 || y>=ysize) 
	{	
		CORE_DUMP_ON_EXCEPTION;
		throw EBoundsError (x,y,xsize,ysize);
	}
#endif		
	return R[x][y];
}

inline byte Image::get (int p, int x, int y) const 
{
#ifdef PEDANTIC_CHECK_CODE
	if (x<0 || x>=xsize || y<0 || y>=ysize) 
	{	
		CORE_DUMP_ON_EXCEPTION;
		throw EBoundsError (x,y,xsize,ysize);
	}
#endif		
	switch (p) 
	{
		case 1:	return	R[x][y];
		case 2:	return	G[x][y];
		case 3:	return	B[x][y];
		default: break;
	}
	return 0;
}

inline void Image::set (int p, int x,	int	y, byte	v) 
{
#ifdef PEDANTIC_CHECK_CODE
	if (x<0 || x>=xsize || y<0 || y>=ysize) 
	{	
		CORE_DUMP_ON_EXCEPTION;
		throw EBoundsError (x,y,xsize,ysize);
	}
#endif	
	switch (p) 
	{
		case 1:	R[x][y]	= v; break;
		case 2:	G[x][y]	= v; break;
		case 3:	B[x][y]	= v; break;
		default: break;
	}
}

inline void Image::set (int x, int y, byte v) 
{
#ifdef PEDANTIC_CHECK_CODE
	if (x<0 || x>=xsize || y<0 || y>=ysize) 
	{	
		CORE_DUMP_ON_EXCEPTION;
		throw EBoundsError (x,y,xsize,ysize);
	}
#endif		
	R[x][y]	= v;
}

// **************************************************************
// ACESSORS
// This acessor is "math-matrix" style, i.e.
// the coordinate order is y,x		
// **************************************************************

inline unsigned char & Image::operator() (int y, int x)
{
#ifdef PEDANTIC_CHECK_CODE
	if (x<0 || x>=xsize || y<0 || y>=ysize) 
	{	
		CORE_DUMP_ON_EXCEPTION;
		throw EBoundsError (x,y,xsize,ysize);
	}
#endif		
	return R[x][y];
}

inline void Image::setZero () 
{
	setZero(1);
	if (type==3) 
	{
		setZero(2);
		setZero(3);
	}
}

#endif

