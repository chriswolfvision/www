#include "BicubicPolynom.h"

using namespace std;

double *			BicubicPolynom::values=NULL;
int     			BicubicPolynom::int_factor=0;
BicubicPolynom *	BicubicPolynom::pInstance=NULL;

// ****************************************************************************
// The constructor: Create the table.
// We _always_ only need values in the interval [-2,2]
// The precision of the desired values depends on the interpolation factor.
// ****************************************************************************

BicubicPolynom::BicubicPolynom (int ifactor) {
	cerr << "Precalculating polynom values ..."; cerr.flush();
	int_factor = ifactor;
	values = new double [4*int_factor+1];
	for (int i=0; i<=4*int_factor; ++i) {
		values[i] = calcPolynom (-2.0 + (double) i / (double) int_factor);
	}
	cerr << "Done.\n"; cerr.flush();
}
