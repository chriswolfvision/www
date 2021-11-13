/**********************************************************
 * TemplatesHistogramming.h
 * Histogram related routines
 * Templates
 *
 * Christian Wolf
 * Beginn: 30.5.2005
 **********************************************************/
 
// From this module
#include "ImageProc.h"

// *************************************************************
// Build a histogram from an Image
// *************************************************************

template <class T>
Histogram<T> *buildHistogram (const Image & im, 
	bool doIgnore, unsigned char ignoreValue) 
{
	unsigned char v;
	Histogram<T> *h = new Histogram<T> (256, 0, 255, true);

    h->clear();

    for (int y=0; y<im.ysize; ++y)
   	for (int x=0; x<im.xsize; ++x)
   	{
   		v=im.get(x,y);
   		if (!doIgnore || v!=ignoreValue)
   			++h->bins[v];
   	}

	return h;
}

// *************************************************************
// Filter a histogram with a given filter	mask.
// The mask must of course be 1D
// *************************************************************

template <class T>
void filterHistogram (Histogram<T> &h, FilterMask &fm) {
	T *FIL, *todelete;

	int	ix, bx, fx;
	double newVal;

#ifdef CHECK_CODE
	if (fm.ysize!=1) {
		cerr << "Error in Histogram<T>::filter(): Filtermask must be 1D!!!\n";
		exit (1);
	}
#endif

	FIL = new T [h.size];

	// Treat the borders
	for (int i=0; i<fm.xsize/2; ++i) {

		FIL[i] = h.bins[i];
		FIL[h.size-1-i] = h.bins[h.size-1-i];
	}

	// Travers the histogram bins
	for	(int x=0+fm.xsize/2; x<(h.size-fm.xsize/2); ++x ) {

		// the left border
		bx = x-(fm.xsize/2);

    	// convolve
    	newVal = 0;
        for	(ix=bx,	fx=0; fx<fm.xsize; ++ix, ++fx)
			newVal += (((double)h.bins[ix])*fm.get(fx,0));

		FIL[x]=	(int) (newVal / 1);


	}

	// Copy result and clean up
	todelete = h.bins;
	h.bins = FIL;
	delete [] todelete;
}
