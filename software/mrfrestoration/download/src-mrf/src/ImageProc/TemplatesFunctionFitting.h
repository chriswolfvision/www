/***************************************************************************
 * TemplatesFunctionFitting.h
 *
 * Author: Christian Wolf
 * Begin: 2.7.2005
 ***************************************************************************/
 
#ifndef _WOLF_TEMPLATES_FUNCTION_FITTING_H_
#define _WOLF_TEMPLATES_FUNCTION_FITTING_H_

// C
#include <math.h>

// C++
#include <iostream>

// From the main module
#include <CIL.h>

// From the IMAGE module
#include <Image.h>

// From the NUMERICAL RECIPES library
#ifdef HAVE_NUM_RECIPES
#include <NumericalRecipes.h>
#endif

// From this module
#include "MathFuncs.h"
#include "DataSerie.h"

/***************************************************************************
 * Help function: create a data serie from a segmented image
 ***************************************************************************/
 
template <class T> 
inline void dataSerieFromImage (T &im, Image &segI, unsigned char label, DataSerie &outserie)
{
	// Collect the values
	for (int y=0; y<im.ysize; ++y)
	for (int x=0; x<im.xsize; ++x)
		if (segI.get(x,y)==label)
			outserie.add(im.get(x,y));				
	// cerr << "dataserie: " << outserie.size() << " values." << endl;
}

/***************************************************************************
 * Fit a gaussian distribution to the values in an image 
 * specified by a specific label in a label image
 ***************************************************************************/
 
template <class T>
void fitGaussian (T &im, Image &segI, unsigned char label, float &outmean, float &outvar)
{
	DataSerie s;
	dataSerieFromImage (im,segI,label,s);	
							
	// fitGeneralizedGaussian (s, outmean, outvar, outshape);
	outmean=s.getMean();
	outvar=s.getVariance();
}

/***************************************************************************
 * Fit a generalized gaussian distribution to the values in an image 
 * specified by a specific label in a label image
 ***************************************************************************/
 
#ifdef HAVE_NUM_RECIPES  
template <class T>
void fitGeneralizedGaussian (T &im, Image &segI, unsigned char label, 
	float &outmean, float &outvar, float &outshape)
{
	DataSerie s;
	dataSerieFromImage (im,segI,label,s);	
							
	// fitGeneralizedGaussian (s, outmean, outvar, outshape);
	fitGeneralizedGaussian (s, outmean, outvar, outshape);
}
#endif // HAVE_NUM_RECIPES

/***************************************************************************
 * Fit a log normal distribution to the values in an image 
 * specified by a specific label in a label image
 ***************************************************************************/
 
template <class T>
void fitLogNormal (T &im, Image &segI, unsigned char label, float inPosition,
	float &outMean, float &outVar)
{
	DataSerie s;
	dataSerieFromImage (im,segI,label,s);
	
	// cerr << "dataserie:" << endl << s;
						
	// fitGeneralizedGaussian (s, outmean, outvar, outshape);
	fitLogNormal (s, inPosition, outMean, outVar);
}

#endif // protection
