/**********************************************************
 * Binarization.h
 * 
 *
 * Christian Wolf
 **********************************************************/

#ifndef _WOLF_BINARIZATION_H_
#define _WOLF_BINARIZATION_H_

// From the IMAGE library
#include "Image.h"
#include "FloatMatrix.h"

// From the MATH library
#include "Histogram.h"
#include "Matrix.h"

// Template Routines are included at 
// the end of this file!!!!!
				  
enum NiblackVersion 
{
	NIBLACK_CLASSIC=0,
    NIBLACK_SAUVOLA,
	NIBLACK_WOLF1,
    NIBLACK_WOLF2,
    NIBLACK_WOLF_2007
};

// Fixed thresholds
void thresholdFixed(Image &i, int	t);
void thresholdFixedLower (Image &i, int t);
void thresholdFixedUpper (Image &i, int t);

// Histogram based (Fisher/Otsu) thresholds
int  thresholdFisher(Image &i, int min, int max, int deltaToMax=0);
void thresholdFisherWindowed (Image &i, int kmin, int kmax, int wxsize);

// Local adaptive thresholding methods
void thresholdWithSurface (Image &i, FloatMatrix &surface);
inline void thresholdFisherWindowed (Image &i, int kmin, int kmax,  int winx, int winy);
inline void thresholdNiblackImproved (Image &i, NiblackVersion version, int win_x, int win_y, double k, double R);

// The core functions creating the threshold surfaces
FloatMatrix * surfaceFisherWindowed (Image &i, int kmin, int kmax,  int winx, int winy);
FloatMatrix * surfaceNiblackImproved (Image &i, NiblackVersion version, int win_x, int win_y, double k, double R);
FloatMatrix * surfaceNiblackImproved (Image &i, NiblackVersion version, int win_x, int win_y, double k, double R,
FloatMatrix *& out_map_m, FloatMatrix *& out_map_s);

void postProcessingYanoBruck (Image &i, Image &orig_input, double Tp, double alpha_deriche, int medsize);

#include "TemplatesImageThresholding.h"

#endif
