// C
#include <math.h>

// From the main module
#include <CIL.h>

// From the MATHEMATICS module
#include <DataSerie.h>

// From this module
#include "FilterMask.h"

/************************************************************************
 * Constructor:	Empty Mask
 ************************************************************************/

FilterMask::FilterMask(int x, int y)
{
	M.resize(x*y);
	xsize =	x;
	ysize =	y;
	xSizeFrozen = true;
}

/************************************************************************
 * Constructor:	Full specified 3x1 Mask
 ************************************************************************/

FilterMask::FilterMask(double v1, double v2, double	v3)
{
	xsize =	3;
	ysize =	1;	
	M.push_back (v1);
	M.push_back (v2);
	M.push_back (v3);
	xSizeFrozen = true;
}

/************************************************************************
 * Constructor:	Out ouf a data serie
 ************************************************************************/

FilterMask::FilterMask(const DataSerie &ds){
	xsize =	ds.size();
	ysize =	1;	
    for (int i=0; i<xsize; ++i)
    	M.push_back(ds[i]);
    xSizeFrozen = true;
}

/************************************************************************
 * Constructor:	special filters
 ************************************************************************/

FilterMask::FilterMask (FilterType type)
{	
	switch (type)
	{
		case FLTTPE_IDENTITY:
			xsize=ysize=3;
			M.resize(9);
			M[0] = 0.; M[1] = 0.; M[2] = 0.;
			M[3] = 0.; M[4] = 1.; M[5] = 0.;
			M[6] = 0.; M[7] = 0.; M[8] = 0.;
			break;
	
		case FLTTPE_GAUSS_3X3:
			xsize=ysize=3;
			M.resize(9);
			M[0] = 0.0625; M[1] = 0.1250; M[2] = 0.0625;
			M[3] = 0.1250; M[4] = 0.2500; M[5] = 0.1250;
			M[6] = 0.0625; M[7] = 0.1250; M[8] = 0.0625;
			break;
		
		default:
			ERR_THROW ("Unknown FilterMask type: " << type);
	}
}


/************************************************************************
 * Assignment Operator
 ************************************************************************/

void FilterMask::operator=(const FilterMask &o){	
	if (this==&o)
		return;
				
	xsize =	o.xsize;
	ysize =	o.ysize;	
	M.clear();
    for (int y=0; y<ysize; ++y)
    for (int x=0; x<xsize; ++x)
    	M.push_back(o.get(x,y));
    xSizeFrozen = true;
}


/************************************************************************
 * Normalize the filter
 ************************************************************************/

void FilterMask::normalize(){
	double sum = 0;

	for	(int y=0; y<ysize; ++y)	{
		for	(int x=0; x<xsize; ++x)	{
			sum	+= get (x,y);
		}
	}

	for	(int y=0; y<ysize; ++y)	{
		for	(int x=0; x<xsize; ++x)	{
			set(x,y,get(x,y)/sum);
		}
	}
}

// *********************************************************************
// OUTPUT
// *********************************************************************

ostream & operator << (ostream & os, FilterMask &h) {
	for (int y=0; y<h.ysize; ++y) {
		os << "fm[" << y << "] = [";
		
		for (int x=0; x<h.xsize; ++x){
			os << " " << h.get(x,y);
		}
		os << " ];\n";
	}
	return os;	
}

// *************************************************************
// Filter a DataSerie with a given filter	mask.
// The mask must of course be 1D
// *************************************************************

void filterDataSerie (DataSerie &d, FilterMask &fm) {
	vector<double> *FIL;
	int	ix, bx, fx;
	double newVal;
	int s=d.size();

#ifdef CHECK_CODE
	if (fm.ysize!=1) {
		cerr << "Error in DataSerie::filter(): Filtermask must be 1D!!!\n";
		exit (1);
	}
#endif

	// Make a copy of the data
	FIL = new vector<double> (d.values);

	// Treat the borders
	for (int i=0; i<fm.xsize/2; ++i) {
		(*FIL)[i] = d.values[i];
		(*FIL)[s-1-i] = d.values[s-1-i];
	}


	// Travers the data
	for	(int x=0+fm.xsize/2; x<(s-fm.xsize/2); ++x ) {

		// the left border
		bx = x-(fm.xsize/2);

    	// convolve
    	newVal = 0;
        for	(ix=bx,	fx=0; fx<fm.xsize; ++ix, ++fx)
			newVal += (((double)d.values[ix])*fm.get(fx,0));

		(*FIL)[x]=	(int) rint (newVal);
	}

	// Copy result and clean up
	d.values = *FIL;
	d.haveStats = d.haveMedian = false;
	delete FIL;
}
