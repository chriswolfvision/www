#ifndef	_WOLF_FILTERMASK_H_
#define	_WOLF_FILTERMASK_H_

#include <iostream>
#include <vector>
#include <stdlib.h>

#include "visualbug.h"

using namespace std;

class DataSerie;

enum FilterType 
{	
	FLTTPE_IDENTITY=0,
	FLTTPE_GAUSS_3X3
};

class FilterMask {

	public:

		// Constructors
		FilterMask() : xsize(0), ysize(0), xSizeFrozen(false) {}
		FilterMask(int x, int y);
		FilterMask(double v1, double v2, double	v3);
		FilterMask(const DataSerie &ds);
		
		// Special Filters (Gaussians etc.)
		FilterMask (FilterType);

		// Assignment operator
		void operator=(const FilterMask &o);
		
		void clear() { for (int	i=0; i<xsize*ysize;	++i) M[i]=0; }
		double get(int x, int y) const { return M[y*xsize+x];	}
		void   set(int x, int y, double	v) { M[y*xsize+x]=v; }

		void normalize();

		int	xsize, ysize;
		
		friend        ostream & operator << (ostream & os, FilterMask &h);	
		friend inline FilterMask & operator << (FilterMask &h, double d);	
        friend inline FilterMask & operator << (FilterMask &h, char c);

	private:

		vector<double> M;
		
		bool xSizeFrozen;
};

// ********************************************************************
// Loading operator which makes the construction of filter mask
// in the form like "fm << 1 << 2 << 4 << 2 << 1"
// possible.
// 2D filters can be constructed using newlines
// ********************************************************************

inline FilterMask & operator << (FilterMask &h, double v) {
	h.M.push_back(v);		
	if (h.ysize==0)
		h.ysize=1;
	if (!h.xSizeFrozen)
		++h.xsize;
	return h;
}

inline FilterMask & operator << (FilterMask &h, char c) {
	if (c!='\n') {
		cerr << "Error using FilterMask::operator <<!!!\n";
		exit (1);		
	}	
	++h.ysize;
	h.xSizeFrozen = true;	
	return h;
}

#endif

