#include "ProjProfile.h"

#include "visualbug.h"
using namespace std;

// ***************************************************************
// Constructor
// ***************************************************************

ProjProfile::ProjProfile (int s) {
	values.resize(s);	
}

// ***************************************************************
// Constructor
// Fuzzy: false ... binary image, each pixel counts or not
//        true .... grayscale image, each pixel counts by its
//                  grayvalue
// ***************************************************************

void ProjProfile::_alloc (Image	&im, Rect &win, bool isVertical, bool fuzzy) {
	int	ns = isVertical	? win.height() : win.width();

	values.resize (ns);
	setZero();

	// Travers the pixels of the image
	for	(int y=win.top; y<=win.bottom; ++y) {
		for	(int x=win.left; x<=win.right; ++x) {
			int gv = im.get(PLANE_RED,x,y);
			int contrib = fuzzy ? gv : 1;
			
			if (gv) {
				if (isVertical)	{
					values[y-win.top] += contrib;
				}
				else {
					values[x-win.left] += contrib;
				}
			}
		}
	}
}

// ***************************************************************
// Calculate the maximum of two profiles
// ***************************************************************

void ProjProfile::max (ProjProfile &other) {
	int min_size = size() < other.size() ? size() : other.size();
	
	for (int i=0; i<min_size; ++i)
		if (values[i] < other.values[i])
			values[i] =	other.values[i];
}

// ***************************************************************
// OUTPUT
// ***************************************************************

void ProjProfile::print	(FILE *fp) {
	fprintf	(fp, "Prof (%d): ",	size());
	for	(int i=0; i<size(); ++i)
		fprintf	(fp, "%3.0f,", values[i]);
	fprintf	(fp, "\n");
}

// ***************************************************************
// OUTPUT
// ***************************************************************

ostream & operator << (ostream & os, ProjProfile &p) {
	for (int i=0; i<p.size(); ++i)
		os << p[i] << endl;
	return os;	
}
