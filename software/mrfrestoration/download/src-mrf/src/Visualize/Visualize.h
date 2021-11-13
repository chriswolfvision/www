/***************************************************************************
 * Visualize.h
 *
 * Author: Christian Wolf
 * Begin: 17.2.2003
 ***************************************************************************/

#ifndef _WOLF_VISUALIZE_H_
#define _WOLF_VISUALIZE_H_

// From the IMAGE module
#include <Image.h>
#include <FloatMatrix.h>

void visualizeDirections (Image &dircodes, FloatMatrix *intensity, char *filename);

// Visualizes the histogram for each class of a segmented image
void plotClassDistributions (Image *I, Image *segI, ostream &st, 
	unsigned int noErosions=0, int ignoreLabel=-1, float factor=1.0);

// Several functions which crate function plots of various statistical laws
// with given parameters in the given interval	
void plotGaussian (ostream &st, float m, float ss, 
	float min, float max, float step, float factor=1., bool plotXAxis=true);
void plotLogNormal (ostream &st, float m, float ss,
	float min, float max, float step);
#ifdef HAVE_NUM_RECIPES
void plotGeneralizedGaussian (ostream &st, float m, float ss, float p, 
	float min, float max, float step, bool plotXAxis);
#endif



#endif // protection

