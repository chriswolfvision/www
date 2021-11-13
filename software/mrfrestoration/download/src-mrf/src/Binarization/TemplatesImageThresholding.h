/**********************************************************
 * TemplatesImageTresholding.h
 * Routines for Image Thresholding
 * Templates
 *
 * Christian Wolf
 * Beginn: 30.5.2005
 **********************************************************/
 
// From this module
#include "Binarization.h"

inline void thresholdFisherWindowed (Image &i, int kmin, int kmax,  int winx, int winy) 
{
	FloatMatrix *surface = surfaceFisherWindowed (i, kmin, kmax, winx, winy);
	thresholdWithSurface (i, *surface);
	delete surface;
}

inline FloatMatrix * surfaceNiblackImproved (Image &i, NiblackVersion version, int win_x, int win_y, 		double k, double R) 
{
    FloatMatrix *dummy1, *dummy2, *rv;
	rv = surfaceNiblackImproved (i, version, win_x, win_y, k, R, dummy1, dummy2);
	delete dummy1;
	delete dummy2;
	return rv;
}

inline void thresholdNiblackImproved (Image &i, NiblackVersion version, int win_x, int win_y, double k, 	double R) 
{
	FloatMatrix *surface = surfaceNiblackImproved (i, version, win_x, win_y, k, R);
	thresholdWithSurface (i, *surface);
	delete surface;
}

