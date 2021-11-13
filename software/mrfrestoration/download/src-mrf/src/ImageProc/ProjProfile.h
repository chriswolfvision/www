#ifndef	_WOLF_PROJPROFILE_H_
#define	_WOLF_PROJPROFILE_H_

// C++
#include <stdlib.h>
#include <stdio.h>
#include <iostream>

// From the IMAGE library
#include <Rect.h>
#include <Image.h>

// From the MATHEMATICS library
#include "DataSerie.h"

class PixelList;

// ----------------------------------------------------------
// Vertical	= vertical projection box,
//			-> horizontal	projection lines
// ----------------------------------------------------------

class ProjProfile : public DataSerie {

	public:

		// Constructors
		ProjProfile	() : DataSerie () {}
		ProjProfile	(int s);
		
		inline ProjProfile (Image &binaryImage, bool isVertical, bool fuzzy);
		inline ProjProfile (Image &binaryImage, Rect & window, bool isVertical, bool fuzzy);
         
		// Destructor
		~ProjProfile () {};
		
		void inc (int bin) { ++values[bin]; }		

		void max (ProjProfile &other);
		
		// bool check4regularCuts (Image *im, int vbeg, int vend);

		// Output
		void print (FILE *fp);
		
		friend ostream & operator << (ostream & os, ProjProfile &p);

	private:
	
		void _alloc (Image &binaryImage, Rect & window, bool isVertical, bool fuzzy);
};

// ***************************************************************
// Constructor
// ***************************************************************

inline ProjProfile::ProjProfile (Image &im, bool isVertical, bool fuzzy) {
	Rect r (0,im.ysize-1,0,im.xsize-1);
	_alloc (im, r, isVertical, fuzzy);
}

// ***************************************************************
// Constructor
// ***************************************************************

inline ProjProfile::ProjProfile (Image	&im, Rect &win, bool isVertical, bool fuzzy) {
	_alloc (im, win, isVertical, fuzzy);
}


#endif

