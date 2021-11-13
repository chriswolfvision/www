/**********************************************************
 * ImageProc.h
 * Prototypes for the image processing routines
 *
 * Christian Wolf
 **********************************************************/

#ifndef _WOLF_IMAGEPROC_H_
#define _WOLF_IMAGEPROC_H_

// From the IMAGE module
#include <Image.h>
#include <FloatMatrix.h>

// From the MATH module
#include <Histogram.h>
#include <Matrix.h>
#include <ConfusionMatrix.h>

// From the color module
#include <color.h>

// From this module
#include "FilterMask.h"
#include "Moments.h"
// Template Routines are included at 
// the end of this file!!!!!
			  
/**************************************************************************	
 *	Basic Math stuff		   
 **************************************************************************/

void absoluteImageDiff (Image &i, const Image & other);
void addValue (Image &i, int a);
void subValue (Image &i, int a);
void gammaCorrection (Image &i,int plane, double	gamma);
void multiply (Image &i, byte plane, double m);
void multiply (Image &i, double m);
void multiplysqrt (Image &i, Image &j);
void divide (Image &i,int plane, double d);
void divide (Image &i, double d);
void add (Image &im, const Image & other);
void subtract (Image &i,const Image &other);
void power (Image &i, double alpha);
void projectToBlue (Image &i,float darken);
void reverseVideo (Image &i);
void transposeImage (Image &,Image &);
void mirrorX(Image &i);
void mirrorY(Image &i);
void maxImage (Image &i,Image &other);
void minImage (Image &i,Image &other);
byte imageMin (Image &i,byte plane);
byte imageMax (Image &i,byte plane);
double ImageSum (Image &i,byte plane);

double calcLocalStats (Image &i, int winx, int winy, FloatMatrix *&out_m, FloatMatrix *&out_s);

/**************************************************************************	
 *	Filtering methods
 **************************************************************************/

void filterMedian (Image &i, int plane, int size);
void filterGeneralizedMedian (Image &im, int plane, int fsize, double k, bool horiz, bool vert);
void filterMask	(Image &i, FilterMask	&fm);
void filterDeriche (Image &i, float alpha);
void filterDericheGetDerivatives (FloatMatrixColor &fim, int plane, float alpha,
	FloatMatrix &outDerivX, FloatMatrix &outDerivY);
void filterSobel (Image &i, bool horiz, bool vertical);
Image *filterSobelReturn (Image &im, bool horiz, bool vertical);
void filterPrewittHoriz	(Image &i);

// Filtering for float matrices
void filterFloatMatrix (FloatMatrix &src, FloatMatrix &dst, FilterMask &fm);
void filterFloatMatrix (FloatMatrix &srcdst, FilterMask &fm);
void filterFloatMatrixMedian (FloatMatrix &src, FloatMatrix &dst, int fsize);
void filterFloatMatrixColor (FloatMatrixColor &im, FilterMask &fm);
void filterFloatMatrixColor (FloatMatrixColor &im, FilterMask &fm, int plane);

// Filter other stuff than images
template <class T> void filterHistogram (Histogram<T> &h, FilterMask &fm);
                   void filterDataSerie (DataSerie &d, FilterMask &fm);

/**************************************************************************	
 *	Histogram operations
 **************************************************************************/
 
template <class T> Histogram<T> *buildHistogram (const Image & im,
	bool doIgnore=false, unsigned char ignoreValue=-1);
void histogramStretch(Image &i, int plane, 
	float indexPercLeft, float indexPercRight);
void brightenImage (Image &i);

/**************************************************************************	
 *	Mathematical morphology
 **************************************************************************/

void erode (Image &im, int Nit, bool horizOnly, unsigned char gvObject);
void erode (Image &i, int iterations);
void erodeHoriz	(Image &i, int iterations) ;
void erodeGray (Image &im, int Nit, bool horizOnly);
void dilate (Image &im, int Nit, bool horizOnly, unsigned char gvObject);
void dilate	(Image &i, int iterations) ;
void dilateHoriz (Image &i, int iterations) ;
void levelCurve4 (Image &i);
void levelCurve4Horiz (Image &i);
void ferret	(Image &i, unsigned int minX, unsigned int maxX, unsigned int minY,	unsigned int maxY);

// Special operations for text
void dilateForConnectingCharacters (Image &i, int Nit, byte dilategrayvalue=255);
void erodeForConnectingCharacters (Image &i, int Nit, byte value2bfiltered=255);
void erodeForRemovingBridges(Image &i, int maxheightbridge);

/**************************************************************************	
 * Thresholding routines are in the
 * "Binarization" module
 **************************************************************************/

/**************************************************************************	
 *	Segmentation
 **************************************************************************/

/* If the last argument is NULL, then it must be cast
 * since the argument is a template:
 * (vector<PixelType> *) NULL
 */
template <class T>
void segmentKMeans(T &i, Image &o, unsigned int noClasses, ColorSpaceCode colspc, 
	vector<float> *initValues=NULL,
	Image **out_segmented_visualize=NULL);
	
/**************************************************************************	
 *	Image Handling
 **************************************************************************/
	
Image * ColorSwap (Image *src, char *fname);

/**************************************************************************	
 *	Statistical routines
 **************************************************************************/
 
ConfusionMatrix * CreateConfusionMatrix (Image *seg, Image *gt);

double getProb (Image &i, byte val);
double getProb (Image &i, byte val1, byte val2, int d_x, int d_y);
Image       *uniformNoise (int xsize, int ysize, int planes);
Image       *gaussianNoise (int xsize, int ysize, int planes, float mean, float stddev);
FloatMatrix *gaussianNoise (int xsize, int ysize,             float mean, float stddev);

void getDistributions (Image &i, double khigh, double klow, int winx, int winy,
	FloatMatrix *&outim_t, FloatMatrix *& outim_b, FloatMatrix *& outim_var);
	
/**************************************************************************	
 *	Function and distribution fitting
 **************************************************************************/	

template <class T>
void fitGaussian (T &im, Image &segI, unsigned char label, float &outmean, float &outvar);
#ifdef HAVE_NUM_RECIPES		
template <class T> 
void fitGeneralizedGaussian (T &I, Image &segI, unsigned char label, 
		float &outmean, float &outvar, float &outshape);
#endif
template <class T> 
void fitToto (T &I, Image &segI, unsigned char label, 
		float &outmean, float &outvar);
	
/**************************************************************************	
 *	Astronomical image processing
 **************************************************************************/	
	
// Drawing routines
void drawEllipseFromMoments (FloatMatrix &i, Moments &m, double r_square, double v);
	
// Astronomy
void galaxyShapeSingleBand (FloatMatrix &i, FloatMatrix &o);			

/**************************************************************************	
 *	Include a couple of template routines defined outside
 **************************************************************************/	

#include "TemplatesImageSegmentation.h"
#include "TemplatesHistogramming.h"
#include "TemplatesFunctionFitting.h" 

#endif
